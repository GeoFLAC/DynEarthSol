#ifndef DES_MPI_LOCAL_FIELDS_HPP
#define DES_MPI_LOCAL_FIELDS_HPP
#ifdef USEMPI
// Fields stored on compact compute meshes and reconstructed for output/management.
#ifdef USEMMG
// Size hints are populated only in MMG builds. Without MMG their root vectors
// stay empty and must not be indexed by compact node IDs or packed into halos.
#define MPI_MMG_NODE_VECTORS(F) F(init_elem_size_n) F(refined_elem_size_n)
#else
#define MPI_MMG_NODE_VECTORS(F)
#endif
#define MPI_NODE_VECTORS(F) F(volume_n) F(mass) F(tmass) F(hmass) F(ymass) \
 F(temperature) F(ppressure) F(dppressure) F(dppressure_zero) F(fluid_source) F(fmass) \
 MPI_MMG_NODE_VECTORS(F)
#define MPI_NODE_ARRAYS(F) F(coord) F(coord0) F(vel) F(force) F(force_residual)
#define MPI_ELEMENT_VECTORS(F) F(volume) F(volume_old) F(edvoldt) F(plstrain) \
 F(delta_plstrain) F(plstrain_remesh) F(stressyy) F(dpressure) F(viscosity) \
 F(old_mean_stress) F(radiogenic_source) F(dyn_fric_coeff) F(state_variable)
#define MPI_ELEMENT_ARRAYS(F) F(stress) F(strain) F(strain_rate)
#endif
#endif
