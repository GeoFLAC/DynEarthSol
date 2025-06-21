#ifndef DYNEARTHSOL3D_KNN_CU
#define DYNEARTHSOL3D_KNN_CU

#ifdef ACC

class CudaKNN
{
public:
    CudaKNN(const Param& param, const double3_vec& points) : 
        points(points), resolution(param.mesh.resolution) {};
    ~CudaKNN() {};

    void search_grid(const double3_vec& queries, neighbor_vec& neighbors, 
            int k, double resoTimes = 3);
private:
    const double3_vec& points;
    const int resolution;
};

#endif

#endif // DYNEARTHSOL3D_KNN_CU