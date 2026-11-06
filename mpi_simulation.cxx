#ifdef USEMPI
#include <mpi.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <new>
#include <numeric>
#include <set>
#include <sstream>
#include <utility>
#include "mpi_simulation.hpp"
#include "mpi_partition.hpp"
#include "mpi_local_fields.hpp"
#include "mpi_transport.hpp"
#include "mpi_halo_kernels.hpp"
#include "bc.hpp"
#include "earthquake_state.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "markerset.hpp"
#include "matprops.hpp"
#include "output.hpp"
#include "phasechanges.hpp"
#include "remeshing.hpp"
#include "rheology.hpp"

namespace {
struct Fields {
    struct Field { double* data; int stride, columns; };
    std::array<Field,64> fields;
    int n=0,columns=0;
    void add(double_vec* v) { if(v) add(v->data(),v->size(),1); }
    template<int N> void add(Array2D<double,N>* v) { if(v) add(v->data(),v->size(),N); }
    void add(double* data,int rows,int cols) {
        if(n==int(fields.size())) die(EXIT_INTERNAL_ASSERT,"too many MPI fields");
        fields[n++]={data,rows,cols};columns+=cols;
    }
};
Fields node_fields(Variables& v) {
    Fields f;
#define F(x) f.add(v.x);
    MPI_NODE_VECTORS(F) MPI_NODE_ARRAYS(F)
#undef F
    return f;
}
Fields element_fields(Variables& v) {
    Fields f;
#define F(x) f.add(v.x);
    MPI_ELEMENT_VECTORS(F) MPI_ELEMENT_ARRAYS(F)
#undef F
    return f;
}
void pack(const Fields& f,const int_vec& rows,double_vec& out) {
    out.resize(std::size_t(rows.size())*f.columns); std::size_t k=0;
    if(rows.empty())return;
    for(int i=0;i<f.n;++i) {
        mpi_pack_columns(f.fields[i].data,f.fields[i].stride,f.fields[i].columns,rows.data(),rows.size(),out.data()+k);
        k+=std::size_t(rows.size())*f.fields[i].columns;
    }
}
void unpack(const Fields& f,const int_vec& rows,const double* in) {
    std::size_t k=0;
    if(rows.empty())return;
    for(int i=0;i<f.n;++i) {
        mpi_unpack_columns(f.fields[i].data,f.fields[i].stride,f.fields[i].columns,rows.data(),rows.size(),in+k);
        k+=std::size_t(rows.size())*f.fields[i].columns;
    }
}
int count(std::size_t n) {
    if(n>INT_MAX) die(EXIT_RUNTIME_RESOURCE,"MPI message exceeds INT_MAX entries");
    return n;
}
void metadata(MpiState& a,Variables& v) {
#define F(x) a.value(v.x)
    F(time); F(last_remesh_time); F(dt); F(steps); F(nremesh); F(max_vbc_val);
    F(max_global_vel_mag); F(global_dt_min); F(compensation_pressure); F(bottom_temperature);
    F(nnode); F(nelem); F(vbc_val_z1_loading_period);
    for(int i=0;i<nbdrytypes;++i) {
        a.value(v.vbc_types[i]);a.value(v.vbc_values[i]);a.value(v.stress_bc_values[i]);
        a.value(v.bnodes[i]);a.value(v.bfacets[i]);
    }
    for(int i=0;i<nbdrytypes_hydro;++i) { a.value(v.hbc_types[i]);a.value(v.stress_bc_types[i]); }
    for(int i=0;i<nbdrytypes*nbdrytypes;++i) a.value(v.edge_slot[i]);
    F(bnormals);F(bcflag);F(edge_vec);F(vbc_vertical_div_x0);F(vbc_vertical_div_x1);
    F(vbc_vertical_ratio_x0);F(vbc_vertical_ratio_x1);a.value(v.surfinfo.max_surf_vel);
#undef F
}
void local_archive(MpiState& a,Variables& v) {
    metadata(a,v);
#define F(x) a.value(v.x);
    MPI_NODE_VECTORS(F) MPI_NODE_ARRAYS(F)
    MPI_ELEMENT_VECTORS(F) MPI_ELEMENT_ARRAYS(F)
    F(connectivity) F(elemmarkers) F(markers_in_elem) F(hydrous_elemmarkers)
    F(log_table) F(tan_table) F(sin_table)
#undef F
    a.value(v.support.arr_data);a.value(v.support.idx_data);a.value(v.support.lidx_data);v.support.rebind();
    int nm=v.markersets.size();a.value(nm);
    if(a.reading()) while(int(v.markersets.size())>nm) {delete v.markersets.back();v.markersets.pop_back();}
    for(int j=0;j<nm;++j) {
        std::string name=a.reading()?"":v.markersets[j]->get_name();a.value(name);
        if(a.reading() && j==int(v.markersets.size())) v.markersets.push_back(new MarkerSet(name));
        v.markersets[j]->mpi_state(a);
    }
    if(a.reading()) {
        if(!v.ntmp) v.ntmp=new double_vec;
        if(!v.etmp) v.etmp=new double_vec;
        if(!v.etmp_int) v.etmp_int=new int_vec;
        if(!v.tmp_result) v.tmp_result=new elem_cache;
        v.ntmp->resize(v.nnode);v.etmp->resize(v.nelem);
        v.etmp_int->resize(v.nelem);v.tmp_result->resize(v.nelem,false);
    }
}
void release(Variables& v) {
    delete v.mat;v.mat=nullptr;
#define F(x) delete v.x;v.x=nullptr;
    MPI_NODE_VECTORS(F) MPI_NODE_ARRAYS(F)
    MPI_ELEMENT_VECTORS(F) MPI_ELEMENT_ARRAYS(F)
    F(connectivity) F(elemmarkers) F(markers_in_elem) F(hydrous_elemmarkers)
    F(log_table) F(tan_table) F(sin_table) F(bnormals) F(bcflag)
    F(ntmp) F(etmp) F(etmp_int) F(tmp_result)
#undef F
    for(int i=0;i<nbdrytypes;++i) {delete v.bnodes[i];delete v.bfacets[i];v.bnodes[i]=nullptr;v.bfacets[i]=nullptr;}
    for(auto* m:v.markersets) delete m;
    v.markersets.clear();
}
template<class T> void select(const std::vector<T>* src,std::vector<T>*& dst,const int_vec& ids) {
    if(!src) {delete dst;dst=nullptr;return;}
    if(!dst) dst=new std::vector<T>;
    dst->resize(ids.size());
    for(int i=0;i<int(ids.size());++i) (*dst)[i]=(*src)[ids[i]];
}
template<class T,int N> void select(const Array2D<T,N>* src,Array2D<T,N>*& dst,const int_vec& ids) {
    if(!src) {delete dst;dst=nullptr;return;}
    if(!dst) dst=new Array2D<T,N>;
    dst->resize(ids.size(),false);
    for(int d=0;d<N;++d) for(int i=0;i<int(ids.size());++i) (*dst)[i][d]=(*src)[ids[i]][d];
}
template<class T> void copy(const T* src,T*& dst) {
    if(!src) {delete dst;dst=nullptr;return;}
    if(!dst) dst=new T;
    *dst=*src;
}
template<class T,int N> void assign(Array2D<T,N>& dst,const Array2D<T,N>& src) {
    dst.resize(src.size(),false); std::copy(src.data(),src.data()+src.num_elements(),dst.data());
}
template<class T,int N> void copy(const Array2D<T,N>* src,Array2D<T,N>*& dst) {
    if(!src) {delete dst;dst=nullptr;return;}
    if(!dst) dst=new Array2D<T,N>;
    assign(*dst,*src);
}

void make_local(const Variables& g,const MpiPartition& p,Variables& v) {
#define F(x) v.x=g.x
    F(time);F(last_remesh_time);F(dt);F(steps);F(nremesh);F(max_vbc_val);
    F(max_global_vel_mag);F(global_dt_min);F(compensation_pressure);F(bottom_temperature);
    F(vbc_val_z1_loading_period);F(edge_vec);F(vbc_vertical_div_x0);F(vbc_vertical_div_x1);
    F(vbc_vertical_ratio_x0);F(vbc_vertical_ratio_x1);F(surfinfo.max_surf_vel);
#undef F
    v.nnode=p.nodes.size();v.nelem=p.elements.size();
#define F(x) select(g.x,v.x,p.nodes);
    MPI_NODE_VECTORS(F) MPI_NODE_ARRAYS(F)
#undef F
#define F(x) select(g.x,v.x,p.elements);
    MPI_ELEMENT_VECTORS(F) MPI_ELEMENT_ARRAYS(F)
    F(connectivity) F(elemmarkers) F(hydrous_elemmarkers)
#undef F
    for(int e=0;e<v.nelem;++e) for(int j=0;j<NODES_PER_ELEM;++j)
        (*v.connectivity)[e][j]=p.node_local.at((*v.connectivity)[e][j]);
    select(g.bcflag,v.bcflag,p.nodes);copy(g.bnormals,v.bnormals);
    for(int i=0;i<nbdrytypes;++i) {
        v.vbc_types[i]=g.vbc_types[i];v.vbc_values[i]=g.vbc_values[i];v.stress_bc_values[i]=g.stress_bc_values[i];
        if(!v.bnodes[i]) v.bnodes[i]=new int_vec;
        v.bnodes[i]->clear();
        for(int n:*g.bnodes[i]) {auto at=p.node_local.find(n);if(at!=p.node_local.end()) v.bnodes[i]->push_back(at->second);}
        if(!v.bfacets[i]) v.bfacets[i]=new std::vector<std::pair<int,int>>;
        v.bfacets[i]->clear();
        for(auto facet:*g.bfacets[i]) {auto at=p.element_local.find(facet.first);
            if(at!=p.element_local.end()) v.bfacets[i]->emplace_back(at->second,facet.second);}
    }
    for(int i=0;i<nbdrytypes_hydro;++i) {v.hbc_types[i]=g.hbc_types[i];v.stress_bc_types[i]=g.stress_bc_types[i];}
    std::copy(g.edge_slot,g.edge_slot+nbdrytypes*nbdrytypes,v.edge_slot);
    // Copy each support in GLOBAL element order, including full stars of owned nodes.
    v.support.idx_data.assign(v.nnode+1,0);v.support.arr_data.clear();v.support.lidx_data.clear();
    for(int n=0;n<v.nnode;++n) {
        int gn=p.nodes[n];const int* patch=g.support.patch(gn);const int* loc=g.support.local(gn);
        for(int i=0;i<g.support.size(gn);++i) {
            auto at=p.element_local.find(patch[i]);
            if(at!=p.element_local.end()) {v.support.arr_data.push_back(at->second);v.support.lidx_data.push_back(loc[i]);}
            else if(n<p.owned_nodes) die(EXIT_INTERNAL_ASSERT,"owned node has an incomplete support");
        }
        v.support.idx_data[n+1]=v.support.arr_data.size();
    }
    v.support.rebind();copy(g.log_table,v.log_table);copy(g.tan_table,v.tan_table);copy(g.sin_table,v.sin_table);
    if(!v.markers_in_elem) v.markers_in_elem=new int_vec2D;
    v.markers_in_elem->assign(v.nelem,int_vec{});
    while(v.markersets.size()>g.markersets.size()) {delete v.markersets.back();v.markersets.pop_back();}
    for(int set=0;set<int(g.markersets.size());++set) {
        const MarkerSet& source=*g.markersets[set];
        if(set==int(v.markersets.size())) v.markersets.push_back(new MarkerSet(source.get_name()));
        int_vec ids,elems;
        for(int m=0;m<source.get_nmarkers();++m) {
            auto at=p.element_local.find(source.get_elem(m));if(at==p.element_local.end())continue;
            if(set==0) (*v.markers_in_elem)[at->second].push_back(ids.size());
            ids.push_back(m);elems.push_back(at->second);
        }
        v.markersets[set]->mpi_subset(source,ids,elems);
    }
}
const char* op_name(MpiOp op) {
    switch(op) {
#define F(x) case MpiOp::x:return #x
        F(State);F(Temperature);F(PorePressure);F(StrainRate);F(Stress);F(Force);F(Velocity);
        F(Coordinate);F(Rotation);F(Volume);F(Dvoldt);F(Edvoldt);F(Nmd);F(Mass);F(TimeStep);
        F(OldMeanStress);F(Friction);F(VelocityBC);F(Residual);F(Surface);F(SwapVolume);
        F(Remesh);F(PhaseChange);F(Hydrous);F(Collect);F(Geometry);F(Boundary);F(Sample);F(Average);F(Metrics);
        F(Report);F(Stop);F(Count);
#undef F
    } return "unknown";
}
bool management(MpiOp op) {
    return op==MpiOp::State || op==MpiOp::Remesh || op==MpiOp::Surface || op==MpiOp::PhaseChange ||
           op==MpiOp::Hydrous || op==MpiOp::Report || op==MpiOp::Collect || op==MpiOp::Geometry || op==MpiOp::Boundary || op==MpiOp::Sample;
}
struct Command {
    int operation,steps;
    double time,dt,max_vbc,max_vel,global_dt,surface_vel,argument;
    bool hydraulic,pt_jump,body_adjustment;
};
}

struct MpiSimulation::Impl {
    explicit Impl(Param& p):param(p) {local.func_time={};packing.func_time={};}
    ~Impl() { release(local);release(packing); }
    Param& param;
    Variables* master=nullptr;
    Variables local{},packing{};
    MpiPartition part;
    std::vector<MpiPartition> parts; // only rank zero retains global partition maps
    struct Buffers {double_vec send,recv;};
    std::vector<Buffers> peer;
    std::vector<MPI_Request> requests;
    std::vector<MpiPeerBuffer> transfers;
    int rank=0,size=1,global_n=0,global_e=0,last_step=-1,generation=0;
    bool active=false,inside=false;
    MpiOp current=MpiOp::State;
    double scalar=0,measurements[2]={};
    std::uint64_t version=0,root_version=0,surface_version=~std::uint64_t(0);
    std::uint64_t calls[int(MpiOp::Count)]={},payload_bytes=0,state_bytes=0,root_bytes=0,all_payload_bytes=0;
    double elapsed[int(MpiOp::Count)]={},communication_s=0;
    double operation_comm[int(MpiOp::Count)]={},wait_s=0,interior_window_s=0;
    double_vec rank_profile;
    bool overlap_enabled=true,pending_active=false;
    int pending=0;
    bool pending_nodes=false;
    Fields pending_fields;
    double posted_at=0;
    std::uint64_t geometry_version=~std::uint64_t(0),overlap_calls=0;
    std::vector<unsigned char> bytes;
    double_vec send,recv,scatter;
    int_vec counts,offsets,queries_n,queries_e;
    array_t average_coord;
    tensor_t average_strain,average_stress;
    double_vec average_plastic;
    double average_time=0;

    void command(const Command& c) {
        local.steps=c.steps;local.time=c.time;local.dt=c.dt;local.max_vbc_val=c.max_vbc;
        local.max_global_vel_mag=c.max_vel;local.global_dt_min=c.global_dt;
        local.surfinfo.max_surf_vel=c.surface_vel;
        param.control.has_hydraulic_diffusion=c.hydraulic;param.control.PT_jump=c.pt_jump;
        param.ic.has_body_force_adjustment=c.body_adjustment;
    }
    void trace() {
        if(rank!=0)return;
        const char* path=std::getenv("DES_MPI_PARTITION_TRACE");
        std::ofstream file;
        if(path) {
            std::ostringstream name;name<<path<<"."<<std::setw(4)<<std::setfill('0')<<generation<<".csv";
            file.open(name.str());if(!file) die(EXIT_IO_OPEN,"cannot create MPI partition trace");
            file<<"entity,id,owner,x,z";
#ifdef THREED
            file<<",y";
#endif
            for(int j=0;j<NODES_PER_ELEM;++j) file<<",n"<<j;
            file<<"\n"<<std::setprecision(17);
        }
        for(int r=0;r<size;++r) {
            const auto& p=parts[r];
            std::cout<<"[mpi][partition] generation="<<generation<<" step="<<master->steps<<" rank="<<r
                     <<" owned_nodes="<<p.owned_nodes<<" owned_elements="<<p.owned_elements
                     <<" owned_markers="<<p.owned_markers<<" ghost_nodes="<<p.nodes.size()-p.owned_nodes
                     <<" ghost_elements="<<p.elements.size()-p.owned_elements<<" neighbors="<<p.halos.size()
                     <<" interior_nodes="<<p.interior_nodes.size()<<" interface_nodes="<<p.interface_nodes.size()<<'\n';
            if(!file)continue;
            for(int i=0;i<p.owned_nodes;++i) {int n=p.nodes[i];file<<"node,"<<n<<','<<r<<','<<(*master->coord)[n][0]<<','<<(*master->coord)[n][NDIMS-1];
#ifdef THREED
                file<<','<<(*master->coord)[n][1];
#endif
                for(int j=0;j<NODES_PER_ELEM;++j) file<<',';
                file<<'\n';}
            for(int i=0;i<p.owned_elements;++i) {int e=p.elements[i];double x[NDIMS]={};
                for(int j=0;j<NODES_PER_ELEM;++j) for(int d=0;d<NDIMS;++d) x[d]+=(*master->coord)[(*master->connectivity)[e][j]][d]/NODES_PER_ELEM;
                file<<"element,"<<e<<','<<r<<','<<x[0]<<','<<x[NDIMS-1];
#ifdef THREED
                file<<','<<x[1];
#endif
                for(int j=0;j<NODES_PER_ELEM;++j) file<<','<<(*master->connectivity)[e][j];
                file<<'\n';}
        }
    }
    void distribute(bool rebalance) {
        const double start=MPI_Wtime();
        if(rank==0) {
            // Serial remeshing leaves this inactive hydraulic scratch array at
            // its old size. Compact selection/collection still visits it.
            if(!param.control.has_hydraulic_diffusion && master->old_mean_stress)
                master->old_mean_stress->resize(master->nelem,0.0);
            global_n=master->nnode;global_e=master->nelem;
            if(rebalance || parts.empty()) {parts=mpi_partition_mesh(*master,size);++generation;trace();}
        }
        MPI_Bcast(&global_n,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&global_e,1,MPI_INT,0,MPI_COMM_WORLD);
        for(int r=0;r<size;++r) {
            std::uint64_t n=0;
            if(rank==0) {
                make_local(*master,parts[r],packing);MpiState a(bytes,true);parts[r].archive(a);local_archive(a,packing);n=bytes.size();
                if(r) {MPI_Send(&n,1,MPI_UINT64_T,r,20,MPI_COMM_WORLD);MPI_Send(bytes.data(),count(n),MPI_BYTE,r,21,MPI_COMM_WORLD);state_bytes+=n;}
            }
            if(rank==r) {
                if(r) {MPI_Recv(&n,1,MPI_UINT64_T,0,20,MPI_COMM_WORLD,MPI_STATUS_IGNORE);bytes.resize(n);MPI_Recv(bytes.data(),count(n),MPI_BYTE,0,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}
                MpiState a(bytes,false);part.archive(a);local_archive(a,local);a.finish();
                if(!local.mat) local.mat=new MatProps(param,local);else local.mat->refresh_elem_cache();
                local.mpi=master?master->mpi:local.mpi;
                if(rebalance && !local.markersets.empty()) local.markersets[0]->check_marker_elem_consistency(local);
            }
        }
        double hydro=rank==0?master->mat->hydro_diff_max:0;
        MPI_Bcast(&hydro,1,MPI_DOUBLE,0,MPI_COMM_WORLD);local.mat->hydro_diff_max=hydro;
        peer.resize(part.halos.size());requests.resize(2*part.halos.size());transfers.resize(part.halos.size());
        last_step=-1;root_version=version;surface_version=geometry_version=~std::uint64_t(0);
        communication_s+=MPI_Wtime()-start;
    }
    void gather_rows(bool nodes,const Fields& from,const Fields& to,bool surface=false,const int_vec* query=nullptr) {
        int_vec rows,global_rows;
        const int_vec* local_rows=nodes?&part.node_rows:&part.element_rows;
        if(surface) local_rows=nodes?&part.surface_owned_nodes:&part.surface_owned_elements;
        if(query) {
            const auto& map=nodes?part.node_local:part.element_local;
            const int owned=nodes?part.owned_nodes:part.owned_elements;
            for(int id:*query) {auto at=map.find(id);if(at!=map.end() && at->second<owned) rows.push_back(at->second);}
            local_rows=&rows;
        }
        pack(from,*local_rows,send);
        if(rank==0) {
            counts.assign(size,0);offsets.assign(size,0);
            for(int r=0;r<size;++r) {
                const auto& p=parts[r];int n=nodes?p.owned_nodes:p.owned_elements;
                if(surface)n=nodes?p.surface_owned_nodes.size():p.surface_owned_elements.size();
                if(query) {n=0;const auto& map=nodes?p.node_local:p.element_local;const int owned=nodes?p.owned_nodes:p.owned_elements;
                    for(int id:*query) {auto at=map.find(id);if(at!=map.end() && at->second<owned) ++n;}}
                counts[r]=count(std::size_t(n)*from.columns);
                if(r)offsets[r]=offsets[r-1]+counts[r-1];
            }
            recv.resize(std::size_t(offsets.back())+counts.back());
        }
        double start=MPI_Wtime();
        MPI_Gatherv(send.data(),count(send.size()),MPI_DOUBLE,recv.data(),counts.data(),offsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        communication_s+=MPI_Wtime()-start;
        if(rank==0) for(int r=0;r<size;++r) {
            const auto& p=parts[r];const auto& ids=nodes?p.nodes:p.elements;
            const int_vec* selected=nodes?&p.node_rows:&p.element_rows;
            if(surface) selected=nodes?&p.surface_owned_nodes:&p.surface_owned_elements;
            global_rows.clear();
            if(query) {const auto& map=nodes?p.node_local:p.element_local;int owned=nodes?p.owned_nodes:p.owned_elements;
                for(int id:*query) {auto at=map.find(id);if(at!=map.end() && at->second<owned) global_rows.push_back(id);}}
            else for(int i:*selected)global_rows.push_back(ids[i]);
            unpack(to,global_rows,recv.data()+offsets[r]);if(r)root_bytes+=std::uint64_t(counts[r])*sizeof(double);
        }
    }
    void scatter_surface(bool nodes,const Fields& from,const Fields& to) {
        const auto& own=nodes?part.surface_nodes:part.surface_elements;
        if(rank==0) {
            counts.assign(size,0);offsets.assign(size,0);scatter.clear();
            for(int r=0;r<size;++r) {
                const auto& p=parts[r];const auto& rows=nodes?p.surface_nodes:p.surface_elements;const auto& ids=nodes?p.nodes:p.elements;
                int_vec global;for(int i:rows)global.push_back(ids[i]);pack(from,global,send);
                offsets[r]=scatter.size();counts[r]=count(send.size());scatter.insert(scatter.end(),send.begin(),send.end());
                if(r)root_bytes+=send.size()*sizeof(double);
            }
        }
        recv.resize(std::size_t(own.size())*to.columns);double start=MPI_Wtime();
        MPI_Scatterv(scatter.data(),counts.data(),offsets.data(),MPI_DOUBLE,recv.data(),count(recv.size()),MPI_DOUBLE,0,MPI_COMM_WORLD);
        communication_s+=MPI_Wtime()-start;unpack(to,own,recv.data());
    }
    void full() {
        if(root_version==version)return;
        Fields gn,ge;if(rank==0) {gn=node_fields(*master);ge=element_fields(*master);}
        gather_rows(true,node_fields(local),gn);gather_rows(false,element_fields(local),ge);
        if(rank==0) {master->mat->refresh_elem_cache();master->mat->hydro_diff_max=local.mat->hydro_diff_max;}root_version=version;
    }
    void geometry() {
        if(geometry_version==version || root_version==version)return;
        Fields gn,ge,ln,le;ln.add(local.coord);le.add(local.volume);
        if(rank==0){gn.add(master->coord);ge.add(master->volume);}
        gather_rows(true,ln,gn);gather_rows(false,le,ge);geometry_version=version;
    }
    Fields surface_n(Variables& v) {Fields f;f.add(v.coord);f.add(v.volume_n);return f;}
    Fields surface_e(Variables& v) {Fields f;f.add(v.volume);f.add(v.plstrain);f.add(v.stress);f.add(v.strain);f.add(v.strain_rate);return f;}
    void boundary() {
        if(surface_version==version)return;
        Fields gn,ge;if(rank==0){gn=surface_n(*master);ge=surface_e(*master);}
        gather_rows(true,surface_n(local),gn,true);gather_rows(false,surface_e(local),ge,true);surface_version=version;
    }
    void query() {
        if(rank==0) {MpiState a(bytes,true);a.value(queries_n);a.value(queries_e);}
        std::uint64_t n=bytes.size();MPI_Bcast(&n,1,MPI_UINT64_T,0,MPI_COMM_WORLD);
        if(rank) bytes.resize(n);
        MPI_Bcast(bytes.data(),count(n),MPI_BYTE,0,MPI_COMM_WORLD);
        if(rank){MpiState a(bytes,false);a.value(queries_n);a.value(queries_e);a.finish();}
        Fields gn,ge;if(rank==0){gn=node_fields(*master);ge=element_fields(*master);}
        gather_rows(true,node_fields(local),gn,false,&queries_n);gather_rows(false,element_fields(local),ge,false,&queries_e);
    }
    void average() {
        const int interval=param.mesh.quality_check_step_interval;
        if(interval<=1)return;
        if(local.steps%interval==1) {
            average_time=local.time;assign(average_coord,*local.coord);assign(average_strain,*local.strain);
            assign(average_stress,*local.stress);average_plastic.assign(local.nelem,0.0);
        } else {
            for(int c=0;c<NSTR;++c) for(int e=0;e<part.owned_elements;++e) average_stress[e][c]+=(*local.stress)[e][c];
            for(int e=0;e<part.owned_elements;++e) average_plastic[e]+=(*local.delta_plstrain)[e];
        }
    }
    void collect_averages() {
        if(param.mesh.quality_check_step_interval<=1)return;
        Fields gn,ge,ln,le;ln.add(&average_coord);le.add(&average_strain);le.add(&average_stress);le.add(&average_plastic);
        if(rank==0) {
            Output& o=*master->output;o.time0=average_time;o.coord0.resize(global_n,false);
            o.strain0.resize(global_e,false);o.stress_avg.resize(global_e,false);o.delta_plstrain_avg.resize(global_e);
            gn.add(&o.coord0);ge.add(&o.strain0);ge.add(&o.stress_avg);ge.add(&o.delta_plstrain_avg);
        }
        gather_rows(true,ln,gn);gather_rows(false,le,ge);
    }
    void execute(MpiOp op,double arg) {
        current=op;inside=true;double started=MPI_Wtime(),comm_before=communication_s;Variables& v=local;
        if(!management(op) && last_step!=v.steps) {v.mat->refresh_elem_cache();last_step=v.steps;}
        switch(op) {
        case MpiOp::State:distribute(true);break;
        case MpiOp::Temperature:update_temperature(param,v,*v.temperature,*v.tmp_result);exchange(MpiEntity::Node,{{v.temperature->data(),1}});break;
        case MpiOp::PorePressure:update_pore_pressure(param,v,*v.ppressure,*v.dppressure,*v.ntmp,*v.tmp_result,*v.stress,*v.old_mean_stress);
            exchange(MpiEntity::Node,{{v.ppressure->data(),1},{v.dppressure->data(),1}});break;
        case MpiOp::StrainRate:update_strain_rate(v,*v.strain_rate);break;
        case MpiOp::Dvoldt:compute_dvoldt(v,*v.ntmp,*v.etmp);exchange(MpiEntity::Node,{{v.ntmp->data(),1}});break;
        case MpiOp::Edvoldt:compute_edvoldt(v,*v.ntmp,*v.edvoldt);break;
        case MpiOp::Stress:update_stress(param,v,*v.stress,*v.stressyy,*v.dpressure,*v.viscosity,*v.strain,*v.plstrain,
            *v.delta_plstrain,*v.strain_rate,*v.ppressure,*v.dppressure,*v.vel,*v.dyn_fric_coeff,*v.state_variable);break;
        case MpiOp::Nmd:NMD_stress(v,*v.stress,*v.ntmp,*v.etmp);break;
        case MpiOp::Force:update_force(param,v,*v.force,*v.force_residual,*v.tmp_result);break;
        case MpiOp::Velocity:update_velocity(v,*v.vel);exchange(MpiEntity::Node,{{v.vel->data(),NDIMS}});break;
        case MpiOp::VelocityBC:apply_vbcs(param,v,*v.vel);exchange(MpiEntity::Node,{{v.vel->data(),NDIMS}});break;
        case MpiOp::Coordinate:update_coordinate(v,*v.coord);exchange(MpiEntity::Node,{{v.coord->data(),NDIMS}});break;
        case MpiOp::Rotation:rotate_stress(v,*v.stress,*v.strain,arg);break;
        case MpiOp::Volume:compute_volume(v,*v.volume);exchange(MpiEntity::Element,{{v.volume->data(),1}});break;
        case MpiOp::Mass:compute_mass(param,v,arg,*v.volume_n,*v.mass,*v.tmass,*v.hmass,*v.ymass,*v.tmp_result);
            exchange(MpiEntity::Node,{{v.volume_n->data(),1},{v.mass->data(),1},{v.tmass->data(),1},{v.hmass->data(),1},{v.ymass->data(),1}});break;
        case MpiOp::OldMeanStress:update_old_mean_stress(param,v,*v.stress,*v.old_mean_stress);break;
        case MpiOp::Friction:refresh_rsf_friction(param,v,*v.dyn_fric_coeff,*v.state_variable);break;
        case MpiOp::TimeStep:scalar=compute_dt(param,v,arg!=0);if(rank==0){master->max_global_vel_mag=v.max_global_vel_mag;master->global_dt_min=v.global_dt_min;}break;
        case MpiOp::Residual:scalar=calculate_residual_force(v,*v.force_residual);if(rank==0)master->l2_residual=scalar;break;
        case MpiOp::SwapVolume:std::swap(v.volume,v.volume_old);if(rank==0)std::swap(master->volume,master->volume_old);break;
        case MpiOp::Surface: {
            bool whole=v.steps%param.mesh.quality_check_step_interval==0 || (param.control.surface_process_option!=0 && param.control.surface_process_option!=1);
            if(whole)full();else boundary();
            if(rank==0) surface_processes(param,*master,*master->coord,*master->stress,*master->strain,*master->strain_rate,
                *master->plstrain,*master->volume,*master->volume_n,master->surfinfo,master->markersets,*master->elemmarkers,*master->markers_in_elem);
            if(whole)distribute(false);
            else {Fields gn,ge;if(rank==0){gn=surface_n(*master);ge=surface_e(*master);}
                scatter_surface(true,gn,surface_n(v));scatter_surface(false,ge,surface_e(v));}
            break;
        }
        case MpiOp::Remesh:full();if(rank==0)remesh(param,*master,int(arg));distribute(true);break;
        case MpiOp::PhaseChange:full();if(rank==0){phase_changes(param,*master);master->mat->refresh_elem_cache();}distribute(false);break;
        case MpiOp::Hydrous:full();if(rank==0)advect_hydrous_markers(param,*master,arg,*master->markersets[master->hydrous_marker_index],*master->hydrous_elemmarkers);distribute(false);break;
        case MpiOp::Collect:full();if(arg!=0)collect_averages();break;
        case MpiOp::Geometry:geometry();break;
        case MpiOp::Boundary:boundary();break;
        case MpiOp::Sample:query();break;
        case MpiOp::Average:average();break;
        case MpiOp::Metrics:compute_earthquake_metrics(v,measurements[0],measurements[1]);break;
        case MpiOp::Report: {
            MPI_Reduce(&payload_bytes,&all_payload_bytes,1,MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);
            double sample[3*int(MpiOp::Count)+3];
            for(int i=0;i<int(MpiOp::Count);++i) {sample[3*i]=double(calls[i]);sample[3*i+1]=elapsed[i];sample[3*i+2]=operation_comm[i];}
            sample[3*int(MpiOp::Count)]=wait_s;sample[3*int(MpiOp::Count)+1]=interior_window_s;
            sample[3*int(MpiOp::Count)+2]=double(overlap_calls);
            if(rank==0)rank_profile.resize(size*(3*int(MpiOp::Count)+3));
            MPI_Gather(sample,3*int(MpiOp::Count)+3,MPI_DOUBLE,rank_profile.data(),3*int(MpiOp::Count)+3,MPI_DOUBLE,0,MPI_COMM_WORLD);
            break;
        }
        default:die(EXIT_INTERNAL_ASSERT,"unknown MPI operation");
        }
        if(op!=MpiOp::Collect && op!=MpiOp::Geometry && op!=MpiOp::Boundary && op!=MpiOp::Sample && op!=MpiOp::Average && op!=MpiOp::Metrics) ++version;
        ++calls[int(op)];elapsed[int(op)]+=MPI_Wtime()-started;
        operation_comm[int(op)]+=communication_s-comm_before;inside=false;
    }
    void begin_exchange(MpiEntity entity,std::initializer_list<std::pair<double*,int>> arrays) {
        if(pending_active)die(EXIT_INTERNAL_ASSERT,"overlapping halo buffers reused before completion");
        const bool nodes=entity==MpiEntity::Node;Fields f;
        for(auto a:arrays) f.add(a.first,nodes?local.nnode:local.nelem,a.second);
        double start=MPI_Wtime();
        for(int i=0;i<int(part.halos.size());++i) {
            const auto& h=part.halos[i];const auto& receiving=nodes?h.recv_nodes:h.recv_elements;
            peer[i].recv.resize(std::size_t(receiving.size())*f.columns);
            pack(f,nodes?h.send_nodes:h.send_elements,peer[i].send);
            transfers[i]={h.rank,peer[i].send.data(),peer[i].recv.data(),count(peer[i].send.size()),count(peer[i].recv.size())};
            payload_bytes+=peer[i].send.size()*sizeof(double);
        }
        pending=mpi_post_buffers(transfers.data(),transfers.size(),nodes?30:31,requests.data());
        pending_nodes=nodes;pending_fields=f;pending_active=true;posted_at=MPI_Wtime();
        communication_s+=MPI_Wtime()-start;
    }
    void finish_exchange() {
        if(!pending_active)die(EXIT_INTERNAL_ASSERT,"halo completion without a posted exchange");
        const double started=MPI_Wtime();interior_window_s+=started-posted_at;
        mpi_finish_buffers(pending,requests.data());wait_s+=MPI_Wtime()-started;
        for(int i=0;i<int(part.halos.size());++i)
            unpack(pending_fields,pending_nodes?part.halos[i].recv_nodes:part.halos[i].recv_elements,peer[i].recv.data());
        communication_s+=MPI_Wtime()-started;pending_active=false;
    }
    void exchange(MpiEntity entity,std::initializer_list<std::pair<double*,int>> arrays) {
        begin_exchange(entity,arrays);finish_exchange();
    }
};

MpiSimulation::MpiSimulation(Param& p):impl_(new Impl(p)) {
    int provided;MPI_Init_thread(nullptr,nullptr,MPI_THREAD_FUNNELED,&provided);
    if(provided<MPI_THREAD_FUNNELED)die(EXIT_RUNTIME_RESOURCE,"MPI_THREAD_FUNNELED required");
    MPI_Comm_rank(MPI_COMM_WORLD,&impl_->rank);MPI_Comm_size(MPI_COMM_WORLD,&impl_->size);
    if(is_root()) {
        const char* option=std::getenv("DES_MPI_OVERLAP");
        if(option && std::string(option)!="0" && std::string(option)!="1")die(EXIT_CONFIG_VALUE,"DES_MPI_OVERLAP must be 0 or 1");
        impl_->overlap_enabled=!option || std::string(option)=="1";
    }
    MPI_Bcast(&impl_->overlap_enabled,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
    impl_->local.mpi=this;
    if(is_root())std::cout<<"[mpi] compact spatial partitions; ranks="<<impl_->size<<"; neighbor halos; centralized mesh management\n";
}
MpiSimulation::~MpiSimulation() {
    if(is_root()){Command c{};c.operation=int(MpiOp::Stop);MPI_Bcast(&c,sizeof(c),MPI_BYTE,0,MPI_COMM_WORLD);}
    MPI_Finalize();
}
bool MpiSimulation::is_root()const{return impl_->rank==0;}
bool MpiSimulation::distributed()const{return impl_->size>1;}
bool MpiSimulation::in_kernel()const{return impl_->inside && !management(impl_->current);}
int MpiSimulation::owned_nodes()const{return impl_->part.owned_nodes;}
int MpiSimulation::owned_elements()const{return impl_->part.owned_elements;}
int MpiSimulation::global_nodes()const{return impl_->global_n;}
double MpiSimulation::result()const{return impl_->scalar;}
void MpiSimulation::attach(Variables& v){impl_->master=&v;impl_->active=true;synchronize(v);}
void MpiSimulation::synchronize(Variables& v){dispatch(MpiOp::State,v);}
bool MpiSimulation::dispatch(MpiOp op,const Variables& v,double arg) {
    auto& m=*impl_;if(!distributed() || !m.active || m.inside || m.master!=&v)return false;
    Command c{};c.operation=int(op);c.steps=v.steps;c.time=v.time;c.dt=v.dt;c.max_vbc=v.max_vbc_val;
    c.max_vel=v.max_global_vel_mag;c.global_dt=v.global_dt_min;c.surface_vel=v.surfinfo.max_surf_vel;c.argument=arg;
    c.hydraulic=m.param.control.has_hydraulic_diffusion;c.pt_jump=m.param.control.PT_jump;c.body_adjustment=m.param.ic.has_body_force_adjustment;
    MPI_Bcast(&c,sizeof(c),MPI_BYTE,0,MPI_COMM_WORLD);m.command(c);m.execute(op,arg);return true;
}
int MpiSimulation::worker_loop()try {
    for(;;){Command c{};MPI_Bcast(&c,sizeof(c),MPI_BYTE,0,MPI_COMM_WORLD);if(MpiOp(c.operation)==MpiOp::Stop)return 0;
        impl_->command(c);impl_->execute(MpiOp(c.operation),c.argument);}
}catch(const std::bad_alloc&){die(EXIT_RUNTIME_RESOURCE,"MPI worker allocation failed");}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';die(EXIT_INTERNAL_UNREACHABLE,"MPI worker exception");}
void MpiSimulation::exchange(MpiEntity e,std::initializer_list<std::pair<double*,int>> f){impl_->exchange(e,f);}
void MpiSimulation::begin_exchange(MpiEntity e,std::initializer_list<std::pair<double*,int>> f){
    impl_->begin_exchange(e,f);
    if(!impl_->overlap_enabled)impl_->finish_exchange();else ++impl_->overlap_calls;
}
void MpiSimulation::finish_exchange(){impl_->finish_exchange();}
bool MpiSimulation::overlap()const{return impl_->overlap_enabled;}
MpiNodeRows MpiSimulation::gather_rows(int pass)const {
    const auto& rows=pass==0?impl_->part.interior_nodes:impl_->part.interface_nodes;
    return {rows.data(),int(rows.size())};
}
void MpiSimulation::collect(const Variables& v,bool a){dispatch(MpiOp::Collect,v,a?1:0);}
void MpiSimulation::geometry(const Variables& v){dispatch(MpiOp::Geometry,v);}
void MpiSimulation::boundary(const Variables& v){dispatch(MpiOp::Boundary,v);}
void MpiSimulation::sample(const Variables& v,const int_vec& nodes,const int_vec& elements){impl_->queries_n=nodes;impl_->queries_e=elements;dispatch(MpiOp::Sample,v);}
bool MpiSimulation::metrics(const Variables& v,double& s,double& p){if(!dispatch(MpiOp::Metrics,v))return false;s=impl_->measurements[0];p=impl_->measurements[1];return true;}
void MpiSimulation::minimum(double* data,int n){double t=MPI_Wtime();MPI_Allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);impl_->communication_s+=MPI_Wtime()-t;}
double MpiSimulation::maximum(double v){double r,t=MPI_Wtime();MPI_Allreduce(&v,&r,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);impl_->communication_s+=MPI_Wtime()-t;return r;}
double MpiSimulation::sum(double v){double r,t=MPI_Wtime();MPI_Allreduce(&v,&r,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);impl_->communication_s+=MPI_Wtime()-t;return r;}
int MpiSimulation::sum(int v){int r;double t=MPI_Wtime();MPI_Allreduce(&v,&r,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);impl_->communication_s+=MPI_Wtime()-t;return r;}
void MpiSimulation::maximum(int* v,int n){double t=MPI_Wtime();MPI_Allreduce(MPI_IN_PLACE,v,n,MPI_INT,MPI_MAX,MPI_COMM_WORLD);impl_->communication_s+=MPI_Wtime()-t;}
void MpiSimulation::report() {
    if(impl_->master) dispatch(MpiOp::Report,*impl_->master);
    const auto& m=*impl_;std::cout<<std::fixed<<std::setprecision(9);
    for(int i=0;i<int(MpiOp::Count);++i)if(m.calls[i])std::cout<<"[mpi][kernel] name="<<op_name(MpiOp(i))<<" calls="<<m.calls[i]<<" seconds="<<m.elapsed[i]<<'\n';
    std::cout<<"[mpi][communication] seconds="<<m.communication_s<<" halo_sent_bytes="<<m.all_payload_bytes
             <<" root_collection_bytes="<<m.root_bytes<<" state_bytes="<<m.state_bytes<<'\n';
    const int stride=3*int(MpiOp::Count)+3;
    for(int r=0;r<m.size && !m.rank_profile.empty();++r) {
        const double* p=m.rank_profile.data()+r*stride;
        for(int i=0;i<int(MpiOp::Count);++i)if(p[3*i])
            std::cout<<"[mpi][rank-profile] rank="<<r<<" operation="<<op_name(MpiOp(i))<<" calls="<<p[3*i]
                     <<" elapsed_s="<<p[3*i+1]<<" communication_s="<<p[3*i+2]
                     <<" outside_communication_s="<<std::max(0.0,p[3*i+1]-p[3*i+2])<<'\n';
        std::cout<<"[mpi][overlap] rank="<<r<<" calls="<<p[stride-1]<<" wait_s="<<p[stride-3]
                 <<" work_window_s="<<p[stride-2]<<'\n';
    }
}
#endif
