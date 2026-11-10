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
#ifdef DES_TEST_CUDA
#include <cuda_runtime.h>
#endif
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
    Peer(int own,int other,int width,int extent):rank(other),stride(extent),columns(width) {
        const int count=(extent==67?5:4096)+(own+other)%9;
        nsend=(own+2*other+width)%5==0?0:count;
        nrecv=(other+2*own+width)%5==0?0:count;
        if(own==other)nsend=nrecv=extent==67?7:8192;
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
    // CUDA-aware MPI may initialize its CUDA context inside MPI_Init_thread.
    // Select the rank's device first, then verify the placement with MPI itself.
    int selected=-1;
    if(!host) {
        const char* local=std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");
        if(!local)die(EXIT_CONFIG_VALUE,"GPU mode requires Open MPI local-rank information");
        selected=std::atoi(local);
        const int devices=acc_get_num_devices(acc_device_nvidia);
        if(!devices)die(EXIT_UNSUPPORTED_LIB,"no NVIDIA device exposed; --host only validates CPUs");
        if(selected>=devices)die(EXIT_CONFIG_VALUE,"need one visible GPU per local rank");
        acc_set_device_num(selected,acc_device_nvidia);
        acc_init(acc_device_nvidia);
    } else acc_set_device_type(acc_device_host);
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    if(provided<MPI_THREAD_FUNNELED)die(EXIT_UNSUPPORTED_LIB,"MPI_THREAD_FUNNELED unavailable");
    int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm local;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
    int local_rank,local_size;MPI_Comm_rank(local,&local_rank);MPI_Comm_size(local,&local_size);
    if(!host) {
        const int devices=acc_get_num_devices(acc_device_nvidia);
        if(local_size>devices || selected!=local_rank)
            die(EXIT_CONFIG_VALUE,"Open MPI device selection differs from shared-communicator placement");
    }
    char uuid[33]="host",bus[32]="host";
    if(!host) {
#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
        if(!MPIX_Query_cuda_support())die(EXIT_UNSUPPORTED_LIB,"MPI does not support CUDA device pointers");
#else
        die(EXIT_UNSUPPORTED_LIB,"GPU test requires Open MPI's CUDA capability query");
#endif
#ifdef DES_TEST_CUDA
        int device=-1;
        if(cudaGetDevice(&device)!=cudaSuccess || device!=local_rank)
            die(EXIT_INTERNAL_ASSERT,"CUDA and OpenACC selected different devices");
        cudaDeviceProp properties{};
        if(cudaGetDeviceProperties(&properties,device)!=cudaSuccess ||
           cudaDeviceGetPCIBusId(bus,sizeof(bus),device)!=cudaSuccess)
            die(EXIT_RUNTIME_RESOURCE,"cannot identify selected CUDA device");
        const char* hex="0123456789abcdef";
        for(int i=0;i<16;++i) {
            unsigned char byte=static_cast<unsigned char>(properties.uuid.bytes[i]);
            uuid[2*i]=hex[byte>>4];uuid[2*i+1]=hex[byte&15];
        }
        std::vector<char> ids(33*local_size);
        MPI_Allgather(uuid,33,MPI_CHAR,ids.data(),33,MPI_CHAR,local);
        for(int i=0;i<local_size;++i)for(int j=i+1;j<local_size;++j)
            if(std::memcmp(ids.data()+33*i,ids.data()+33*j,33)==0)
                die(EXIT_CONFIG_VALUE,"multiple local ranks selected the same physical GPU");
#else
        die(EXIT_UNSUPPORTED_LIB,"GPU verification requires DES_TEST_CUDA and CUDA runtime linkage");
#endif
    }
    char hostname[MPI_MAX_PROCESSOR_NAME];int host_length;
    MPI_Get_processor_name(hostname,&host_length);
    // Serialize small identity records so launcher output can be audited reliably.
    for(int current=0;current<size;++current) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==current)std::cout<<"[accelerator][rank] rank="<<rank<<" host="<<hostname
            <<" local_rank="<<local_rank<<" local_size="<<local_size
            <<" uuid="<<uuid<<" pci="<<bus<<std::endl;
    }
    MPI_Comm_free(&local);
    for(int extent:{67,32771})for(int columns:{1,NDIMS,NSTR,NODES_PER_ELEM*NDIMS}) {
        std::vector<Peer*> peers;std::vector<MpiPeerBuffer> buffers;std::vector<MPI_Request> requests(2*size);
        // One rank still exercises kernels, using a self transfer.
        for(int other=0;other<size;++other)if(size==1 || other!=rank)peers.push_back(new Peer(rank,other,columns,extent));
        for(int iteration=0;iteration<3;++iteration) {
            buffers.clear();
            for(auto* p:peers) {
                std::vector<double> input(p->stride*columns),empty(input.size(),-1);
                for(int c=0;c<columns;++c)for(int i=0;i<p->stride;++i)
                    input[c*p->stride+i]=iteration*10000000+rank*100000+p->rank*1000+c*p->stride+i;
                acc_memcpy_to_device(p->field,input.data(),input.size()*sizeof(double));
                acc_memcpy_to_device(p->target,empty.data(),empty.size()*sizeof(double));
                if(p->nsend)mpi_pack_device_columns(p->field,p->stride,columns,p->sending,p->nsend,p->send);
                buffers.push_back({p->rank,p->send,p->receive,p->nsend*columns,p->nrecv*columns});
            }
            // Kernels are synchronous: GPU packing is complete before MPI reads.
            int pending=mpi_post_buffers(buffers.data(),buffers.size(),95,requests.data());
            mpi_finish_buffers(pending,requests.data());
            for(auto* p:peers) {
                if(p->nrecv)mpi_unpack_device_columns(p->target,p->stride,columns,p->receiving,p->nrecv,p->receive);
                std::vector<double> output(p->stride*columns),expected(output.size(),-1);
                acc_memcpy_from_device(output.data(),p->target,output.size()*sizeof(double));
                for(int c=0;c<columns;++c)for(int i=0;i<p->nrecv;++i)
                    expected[c*p->stride+p->recv_ids[i]]=iteration*10000000+p->rank*100000+rank*1000+c*p->stride+(7*i+3)%p->stride;
                if(output!=expected) {
                    for(std::size_t i=0;i<output.size();++i)if(output[i]!=expected[i]) {
                        std::cerr<<"halo mismatch rank="<<rank<<" peer="<<p->rank<<" extent="<<extent
                            <<" columns="<<columns<<" iteration="<<iteration<<" index="<<i
                            <<" actual="<<output[i]<<" expected="<<expected[i]<<std::endl;
                        break;
                    }
                    die(EXIT_INTERNAL_ASSERT,"device pack/MPI/unpack differs from independent CPU reference");
                }
            }
        }
        for(auto* p:peers)delete p;
    }
    if(rank==0)std::cout<<"PASS: "<<size<<" ranks; "<<(host?"OpenACC HOST (not GPU verification)":"OpenACC NVIDIA device buffers")
                        <<"; small/large irregular SoA rows, asymmetric empty messages, scalar/vector/tensor widths\n";
    MPI_Finalize();return EXIT_OK;
}
