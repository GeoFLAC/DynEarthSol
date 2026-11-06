// Production SoA halo indexing + OpenACC device allocation + the real MPI
// transport. --host explicitly selects CPU fallback for testing without a GPU.
#define DES_MPI_DEVICE_KERNELS
#include <openacc.h>
#include <mpi.h>
#ifdef OMPI_MAJOR_VERSION
#include <mpi-ext.h>
#endif
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include "mpi_halo_kernels.hpp"
#include "mpi_transport.hpp"
#include "constants.hpp"
#include "exit_codes.hpp"

[[noreturn]] void die(ExitCode code,const char* message) {
    std::cerr<<"[accelerator halo exit "<<int(code)<<"] "<<message<<'\n';
    int initialized=0;MPI_Initialized(&initialized);
    if(initialized)MPI_Abort(MPI_COMM_WORLD,int(code));
    std::abort();
}
[[noreturn]] void die(ExitCode code) {die(code,"accelerator halo failure");}

template<class T> T* allocate(int count) {
    if(!count)return nullptr;
    auto* p=static_cast<T*>(acc_malloc(sizeof(T)*count));
    if(!p)die(EXIT_RUNTIME_RESOURCE,"OpenACC allocation failed");return p;
}
struct Peer {
    int rank,nsend,nrecv,stride,columns;
    std::vector<int> send_ids,recv_ids;
    double *field,*target,*send,*receive;
    int *sending,*receiving;
    Peer(int own,int other,int width):rank(other),stride(67),columns(width) {
        nsend=(own+2*other)%5==0?0:5+(own+other)%9;
        nrecv=(other+2*own)%5==0?0:5+(own+other)%9;
        if(own==other)nsend=nrecv=7;
        for(int i=0;i<nsend;++i)send_ids.push_back((7*i+3)%stride);
        for(int i=0;i<nrecv;++i)recv_ids.push_back((9*i+1)%stride);
        field=allocate<double>(stride*columns);target=allocate<double>(stride*columns);
        send=allocate<double>(nsend*columns);receive=allocate<double>(nrecv*columns);
        sending=allocate<int>(nsend);receiving=allocate<int>(nrecv);
        if(nsend)acc_memcpy_to_device(sending,send_ids.data(),nsend*sizeof(int));
        if(nrecv)acc_memcpy_to_device(receiving,recv_ids.data(),nrecv*sizeof(int));
    }
    ~Peer(){for(auto* p:{field,target,send,receive})if(p)acc_free(p);if(sending)acc_free(sending);if(receiving)acc_free(receiving);}
};
int main(int argc,char** argv) {
    const bool host=argc==2 && std::strcmp(argv[1],"--host")==0;
    if(!host) {
        const int devices=acc_get_num_devices(acc_device_nvidia);
        if(!devices){std::cerr<<"BLOCKED: no NVIDIA device exposed; use --host only for CPU validation\n";return EXIT_UNSUPPORTED_LIB;}
        const char* local=std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");int index=local?std::atoi(local):0;
        if(index>=devices){std::cerr<<"Need all local GPUs visible and one GPU per local rank\n";return EXIT_CONFIG_VALUE;}
        acc_set_device_num(index,acc_device_nvidia);
    } else acc_set_device_type(acc_device_host);
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(!host) {
#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
        if(!MPIX_Query_cuda_support())die(EXIT_UNSUPPORTED_LIB,"MPI does not support CUDA device pointers");
#else
        die(EXIT_UNSUPPORTED_LIB,"GPU test requires Open MPI's CUDA capability query");
#endif
    }
    for(int columns:{1,NDIMS,NSTR,NODES_PER_ELEM*NDIMS}) {
        std::vector<Peer*> peers;std::vector<MpiPeerBuffer> buffers;std::vector<MPI_Request> requests(2*size);
        // One rank still exercises kernels, using a self transfer.
        for(int other=0;other<size;++other)if(size==1 || other!=rank)peers.push_back(new Peer(rank,other,columns));
        for(int iteration=0;iteration<3;++iteration) {
            buffers.clear();
            for(auto* p:peers) {
                std::vector<double> input(p->stride*columns),empty(input.size(),-1);
                for(int c=0;c<columns;++c)for(int i=0;i<p->stride;++i)
                    input[c*p->stride+i]=iteration*10000000+rank*100000+p->rank*1000+c*p->stride+i;
                acc_memcpy_to_device(p->field,input.data(),input.size()*sizeof(double));
                acc_memcpy_to_device(p->target,empty.data(),empty.size()*sizeof(double));
                if(p->nsend)mpi_pack_columns(p->field,p->stride,columns,p->sending,p->nsend,p->send);
                buffers.push_back({p->rank,p->send,p->receive,p->nsend*columns,p->nrecv*columns});
            }
            // Kernels are synchronous: GPU packing is complete before MPI reads.
            int pending=mpi_post_buffers(buffers.data(),buffers.size(),95,requests.data());
            mpi_finish_buffers(pending,requests.data());
            for(auto* p:peers) {
                if(p->nrecv)mpi_unpack_columns(p->target,p->stride,columns,p->receiving,p->nrecv,p->receive);
                std::vector<double> output(p->stride*columns),expected(output.size(),-1);
                acc_memcpy_from_device(output.data(),p->target,output.size()*sizeof(double));
                for(int c=0;c<columns;++c)for(int i=0;i<p->nrecv;++i)
                    expected[c*p->stride+p->recv_ids[i]]=iteration*10000000+p->rank*100000+rank*1000+c*p->stride+(7*i+3)%p->stride;
                if(output!=expected)die(EXIT_INTERNAL_ASSERT,"device pack/MPI/unpack differs from independent CPU reference");
            }
        }
        for(auto* p:peers)delete p;
    }
    if(rank==0)std::cout<<"PASS: "<<size<<" ranks; "<<(host?"OpenACC HOST (not GPU verification)":"OpenACC NVIDIA device buffers")
                        <<"; irregular SoA rows, asymmetric empty messages, scalar/vector/tensor widths\n";
    MPI_Finalize();return EXIT_OK;
}
