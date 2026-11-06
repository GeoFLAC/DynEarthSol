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
    a.value(surface_owned_nodes); a.value(surface_owned_elements);
    int n=halos.size(); a.value(n); if(a.reading()) halos.resize(n);
    for (auto& h:halos) h.archive(a);
    if(a.reading()) index();
}

namespace {

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
    int_vec elem_owner(v.nelem),node_owner(v.nnode,nr);
    {
        int_vec ids(v.nelem);std::iota(ids.begin(),ids.end(),0);
        std::vector<std::array<double,NDIMS>> centers(v.nelem);
        for(int e=0;e<v.nelem;++e) for(int d=0;d<NDIMS;++d) {
            double sum=0;for(int j=0;j<NODES_PER_ELEM;++j) sum+=(*v.coord)[(*v.connectivity)[e][j]][d];
            centers[e][d]=sum/NODES_PER_ELEM;
        }
        bisect(ids,0,v.nelem,0,nr,centers,elem_owner);
    }
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
    std::cout<<"[mpi][partitioner] requested=rcb actual=rcb seconds="<<elapsed<<'\n';
    return parts;
}
#endif
