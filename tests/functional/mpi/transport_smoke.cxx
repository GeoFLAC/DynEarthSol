// Compile normally for the CPU test, or with -DDES_TEST_CUDA and CUDA runtime
// headers/libs to exercise the identical MPI transport with device buffers.
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include "mpi_transport.hpp"
#include "utils.hpp"
#ifdef DES_TEST_CUDA
#include <cuda_runtime.h>
#ifdef OMPI_MAJOR_VERSION
#include <mpi-ext.h>
#endif
#endif

// This standalone executable does not link the solver's IO/error dependencies.
[[noreturn]] void die(ExitCode code,const char* message) {
    std::cerr<<"[transport test exit "<<int(code)<<"] "<<message<<'\n';
    int initialized=0;MPI_Initialized(&initialized);
    if(initialized)MPI_Abort(MPI_COMM_WORLD,int(code));
    std::abort();
}
[[noreturn]] void die(ExitCode code) {die(code,"transport test failure");}
#ifdef DES_TEST_CUDA
void check(cudaError_t status) {if(status!=cudaSuccess)die(EXIT_RUNTIME_RESOURCE,cudaGetErrorString(status));}
#endif
struct Buffer {
    int peer;
    double_vec send,receive;
    double *source=nullptr,*destination=nullptr;
    Buffer(int rank,int other):peer(other) {
        send.resize((rank+2*other)%5==0?0:1024+17*rank+other);
        receive.resize((other+2*rank)%5==0?0:1024+17*other+rank);
#ifdef DES_TEST_CUDA
        if(!send.empty())check(cudaMalloc(reinterpret_cast<void**>(&source),send.size()*sizeof(double)));
        if(!receive.empty())check(cudaMalloc(reinterpret_cast<void**>(&destination),receive.size()*sizeof(double)));
#else
        source=send.data();destination=receive.data();
#endif
    }
    ~Buffer() {
#ifdef DES_TEST_CUDA
        if(source)check(cudaFree(source));if(destination)check(cudaFree(destination));
#endif
    }
};
int main(int argc,char** argv) {
#ifdef DES_TEST_CUDA
    int devices=0;check(cudaGetDeviceCount(&devices));
    if(!devices)die(EXIT_RUNTIME_RESOURCE,"no CUDA devices");
    const char* env=std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");int local_rank=env?std::atoi(env):0;
    int device=devices==1?0:local_rank;
    if(device>=devices)die(EXIT_CONFIG_VALUE,"need one available GPU per local rank");
    check(cudaSetDevice(device));
#endif
    MPI_Init(&argc,&argv);
    int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
#ifdef DES_TEST_CUDA
#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
    if(!MPIX_Query_cuda_support())die(EXIT_UNSUPPORTED_LIB,"MPI runtime is not CUDA-aware");
#else
    die(EXIT_UNSUPPORTED_LIB,"this test requires Open MPI's CUDA capability query");
#endif
    MPI_Comm local;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
    int local_size;MPI_Comm_size(local,&local_size);cudaDeviceProp properties{};check(cudaGetDeviceProperties(&properties,device));
    std::vector<char> ids(16*local_size);MPI_Allgather(properties.uuid.bytes,16,MPI_BYTE,ids.data(),16,MPI_BYTE,local);
    for(int i=0;i<local_size;++i)for(int j=i+1;j<local_size;++j)
        if(std::memcmp(ids.data()+16*i,ids.data()+16*j,16)==0)die(EXIT_CONFIG_VALUE,"multiple ranks selected the same GPU");
    MPI_Comm_free(&local);
#endif
    {
        std::vector<Buffer*> storage;std::vector<MpiPeerBuffer> peers;std::vector<MPI_Request> requests(2*size);
        for(int other=0;other<size;++other)if(other!=rank)storage.push_back(new Buffer(rank,other));
        for(int iteration=0;iteration<3;++iteration) {
            peers.clear();
            for(auto* b:storage) {
                for(std::size_t i=0;i<b->send.size();++i)b->send[i]=iteration*10000000+rank*100000+b->peer*1000+i;
#ifdef DES_TEST_CUDA
                if(!b->send.empty())check(cudaMemcpy(b->source,b->send.data(),b->send.size()*sizeof(double),cudaMemcpyHostToDevice));
#endif
                peers.push_back({b->peer,b->source,b->destination,int(b->send.size()),int(b->receive.size())});
            }
            mpi_exchange_buffers(peers.data(),peers.size(),77,requests.data());
            for(auto* b:storage) {
#ifdef DES_TEST_CUDA
                if(!b->receive.empty())check(cudaMemcpy(b->receive.data(),b->destination,b->receive.size()*sizeof(double),cudaMemcpyDeviceToHost));
#endif
                for(std::size_t i=0;i<b->receive.size();++i)
                    if(b->receive[i]!=iteration*10000000+b->peer*100000+rank*1000+i)
                        die(EXIT_INTERNAL_ASSERT,"received data differ from the sending rank");
            }
        }
        for(auto* b:storage)delete b;
    }
    if(rank==0)std::cout<<"PASS: "<<size<<" ranks; unequal and zero-length peer messages; "
#ifdef DES_TEST_CUDA
                      <<"CUDA device buffers\n";
#else
                      <<"CPU buffers\n";
#endif
    MPI_Finalize();return EXIT_OK;
}
