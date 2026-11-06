#ifdef USEMPI
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include "mpi_partition.hpp"
#include "utils.hpp"
#ifdef USEMETIS
#include <metis.h>
#endif

void MpiPartition::index() {
    node_local.clear(); element_local.clear();
    for (int i=0;i<int(nodes.size());++i) node_local.emplace(nodes[i],i);
    for (int i=0;i<int(elements.size());++i) element_local.emplace(elements[i],i);
    node_rows.resize(owned_nodes); element_rows.resize(owned_elements);
    std::iota(node_rows.begin(),node_rows.end(),0);
    std::iota(element_rows.begin(),element_rows.end(),0);
}
void MpiPartition::archive(MpiState& a) {
    a.value(owned_nodes); a.value(owned_elements); a.value(owned_markers);
    a.value(nodes); a.value(elements); a.value(surface_nodes); a.value(surface_elements);
    a.value(interior_nodes); a.value(interface_nodes);
    a.value(surface_owned_nodes); a.value(surface_owned_elements);
    int n=halos.size(); a.value(n); if(a.reading()) halos.resize(n);
    for (auto& h:halos) h.archive(a);
    if(a.reading()) index();
}

namespace {
#ifdef USEMETIS
void metis_partition(const Variables& v,int ranks,int_vec& owner) {
    // A face-only dual graph omits dependencies in DES node-star gathers.
    // Connect elements sharing any vertex. The volume objective is a proxy;
    // actual node/element halo bytes still depend on the node-owner assignment.
    std::vector<idx_t> offsets(v.nelem+1,0),adjacency,part(v.nelem);
    int_vec neighbors;
    for(int e=0;e<v.nelem;++e) {
        neighbors.clear();
        for(int j=0;j<NODES_PER_ELEM;++j) {
            const int n=(*v.connectivity)[e][j];
            const int* patch=v.support.patch(n);
            neighbors.insert(neighbors.end(),patch,patch+v.support.size(n));
        }
        std::sort(neighbors.begin(),neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(),neighbors.end()),neighbors.end());
        for(int other:neighbors) if(other!=e) adjacency.push_back(other);
        if(adjacency.size()>std::size_t(IDX_MAX))
            die(EXIT_RUNTIME_RESOURCE,"METIS adjacency exceeds idx_t; use a 64-bit METIS build");
        offsets[e+1]=adjacency.size();
    }
    idx_t vertices=v.nelem,constraints=1,parts=ranks,objective=0,options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING]=0;
    options[METIS_OPTION_SEED]=0;
    options[METIS_OPTION_OBJTYPE]=METIS_OBJTYPE_VOL;
    options[METIS_OPTION_UFACTOR]=30; // allow 3% element-count imbalance
    const int status=METIS_PartGraphKway(&vertices,&constraints,offsets.data(),adjacency.data(),
        nullptr,nullptr,nullptr,&parts,nullptr,nullptr,options,&objective,part.data());
    if(status!=METIS_OK) die(EXIT_MESH_QUALITY,"METIS_PartGraphKway failed");
    for(int e=0;e<v.nelem;++e) {
        if(part[e]<0 || part[e]>=ranks) die(EXIT_INTERNAL_ASSERT,"invalid METIS owner");
        owner[e]=part[e];
    }
    std::cout<<"[mpi][metis] objective=volume value="<<objective<<" imbalance_tolerance=0.03\n";
}
#endif

void bisect(int_vec& ids, int b,int e,int rb,int re,
            const std::vector<std::array<double,NDIMS>>& centers,int_vec& owner) {
    if (re-rb==1) { for(int i=b;i<e;++i) owner[ids[i]]=rb; return; }
    std::array<double,NDIMS> lo{},hi{};
    if(b<e) {lo=hi=centers[ids[b]];}
    for(int i=b;i<e;++i) for(int d=0;d<NDIMS;++d) {
        lo[d]=std::min(lo[d],centers[ids[i]][d]); hi[d]=std::max(hi[d],centers[ids[i]][d]);
    }
    int axis=0;
    for(int d=1;d<NDIMS;++d) if(hi[d]-lo[d]>hi[axis]-lo[axis]) axis=d;
    std::sort(ids.begin()+b,ids.begin()+e,[&](int a,int c){
        return centers[a][axis]<centers[c][axis] ||
              (centers[a][axis]==centers[c][axis] && a<c);
    });
    const int rm=(rb+re)/2, middle=b+std::int64_t(e-b)*(rm-rb)/(re-rb);
    bisect(ids,b,middle,rb,rm,centers,owner); bisect(ids,middle,e,rm,re,centers,owner);
}
}

std::vector<MpiPartition> mpi_partition_mesh(const Variables& v,int nr) {
    const auto started=std::chrono::steady_clock::now();
    const char* requested=std::getenv("DES_MPI_PARTITIONER");
    const std::string method=requested?requested:"rcb";
    if(method!="rcb" && method!="metis")
        die(EXIT_CONFIG_VALUE,"DES_MPI_PARTITIONER must be rcb or metis");
#ifndef USEMETIS
    if(method=="metis") die(EXIT_UNSUPPORTED_LIB,"METIS requested; rebuild with usempi=1 usemetis=1");
#endif
    int_vec elem_owner(v.nelem),node_owner(v.nnode,nr);
    // METIS coarsening is inappropriate for tiny meshes with empty partitions.
    const bool geometric=method=="rcb" || v.nelem<=nr;
    if(geometric) {
        int_vec ids(v.nelem);std::iota(ids.begin(),ids.end(),0);
        std::vector<std::array<double,NDIMS>> centers(v.nelem);
        for(int e=0;e<v.nelem;++e) for(int d=0;d<NDIMS;++d) {
            double sum=0;for(int j=0;j<NODES_PER_ELEM;++j) sum+=(*v.coord)[(*v.connectivity)[e][j]][d];
            centers[e][d]=sum/NODES_PER_ELEM;
        }
        bisect(ids,0,v.nelem,0,nr,centers,elem_owner);
    }
#ifdef USEMETIS
    else metis_partition(v,nr,elem_owner);
#endif
    int_vec node_count(nr,0);
    std::vector<std::pair<int,int_vec>> shared;
    for(int n=0;n<v.nnode;++n) {
        std::set<int> candidates;
        const int* patch=v.support.patch(n);
        for(int i=0;i<v.support.size(n);++i) candidates.insert(elem_owner[patch[i]]);
        if(candidates.empty()) candidates.insert(0);
        if(candidates.size()==1) {node_owner[n]=*candidates.begin();++node_count[node_owner[n]];}
        else shared.emplace_back(n,int_vec(candidates.begin(),candidates.end()));
    }
    for(const auto& item:shared) {
        int r=item.second.front();
        for(int candidate:item.second) if(node_count[candidate]<node_count[r]) r=candidate;
        node_owner[item.first]=r;++node_count[r];
    }
    std::vector<char> surf_e(v.nelem,0),surf_n(v.nnode,0);
    for(int n:*v.surfinfo.top_nodes) {
        const int* patch=v.support.patch(n);
        for(int j=0;j<v.support.size(n);++j) surf_e[patch[j]]=1;
    }
    for(int e=0;e<v.nelem;++e) if(surf_e[e])
        for(int j=0;j<NODES_PER_ELEM;++j) surf_n[(*v.connectivity)[e][j]]=1;
    for(int bd=0;bd<nbdrytypes;++bd) for(int n:*v.bnodes[bd]) surf_n[n]=1;
    std::vector<MpiPartition> parts(nr);
    for(int r=0;r<nr;++r) {
        auto& p=parts[r]; std::vector<char> needed_e(v.nelem,0),needed_n(v.nnode,0);
        for(int e=0;e<v.nelem;++e) if(elem_owner[e]==r) {
            p.elements.push_back(e);needed_e[e]=1;
            if(v.markers_in_elem) p.owned_markers+=(*v.markers_in_elem)[e].size();
        }
        p.owned_elements=p.elements.size();
        for(int n=0;n<v.nnode;++n) if(node_owner[n]==r) {
            p.nodes.push_back(n);needed_n[n]=1;
            const int* patch=v.support.patch(n);
            for(int j=0;j<v.support.size(n);++j) needed_e[patch[j]]=1;
        }
        p.owned_nodes=p.nodes.size();
        for(int e=0;e<v.nelem;++e) if(needed_e[e]) {
            if(elem_owner[e]!=r) p.elements.push_back(e);
            for(int j=0;j<NODES_PER_ELEM;++j) needed_n[(*v.connectivity)[e][j]]=1;
        }
        for(int n=0;n<v.nnode;++n) if(needed_n[n] && node_owner[n]!=r) p.nodes.push_back(n);
        p.index();
        for(int i=0;i<p.owned_nodes;++i) {
            const int n=p.nodes[i];const int* patch=v.support.patch(n);bool interior=true;
            for(int j=0;j<v.support.size(n);++j) if(elem_owner[patch[j]]!=r) {interior=false;break;}
            (interior?p.interior_nodes:p.interface_nodes).push_back(i);
        }
        for(int i=0;i<int(p.nodes.size());++i) if(surf_n[p.nodes[i]]) {
            p.surface_nodes.push_back(i); if(i<p.owned_nodes) p.surface_owned_nodes.push_back(i);
        }
        for(int i=0;i<int(p.elements.size());++i) if(surf_e[p.elements[i]]) {
            p.surface_elements.push_back(i); if(i<p.owned_elements) p.surface_owned_elements.push_back(i);
        }
    }
    std::vector<std::map<int,MpiHalo>> links(nr);
    for(int r=0;r<nr;++r) {
        auto& p=parts[r];
        for(int i=p.owned_nodes;i<int(p.nodes.size());++i) {
            int s=node_owner[p.nodes[i]]; links[r][s].rank=s; links[s][r].rank=r;
            links[r][s].recv_nodes.push_back(i);
            links[s][r].send_nodes.push_back(parts[s].node_local.at(p.nodes[i]));
        }
        for(int i=p.owned_elements;i<int(p.elements.size());++i) {
            int s=elem_owner[p.elements[i]];links[r][s].rank=s;links[s][r].rank=r;
            links[r][s].recv_elements.push_back(i);
            links[s][r].send_elements.push_back(parts[s].element_local.at(p.elements[i]));
        }
    }
    for(int r=0;r<nr;++r) for(auto& h:links[r]) parts[r].halos.push_back(std::move(h.second));
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::cout<<"[mpi][partitioner] requested="<<method<<" actual="<<(geometric?"rcb":"metis")
             <<" seconds="<<elapsed<<'\n';
    return parts;
}
#endif
