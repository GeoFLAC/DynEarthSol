#ifndef DYNEARTHSOL3D_REMESHING_HPP
#define DYNEARTHSOL3D_REMESHING_HPP

int bad_mesh_quality(const Param&, const Variables&, int&);
void remesh(const Param&, Variables&, int);
void initialize_elem_size_n(const Variables&, double_vec &);

#endif
