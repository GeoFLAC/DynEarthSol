#ifndef DYNEARTHSOL3D_NONLOCAL_HPP
#define DYNEARTHSOL3D_NONLOCAL_HPP

// Nonlocal (integral-type) regularization of plastic-strain softening (mat.nonlocal_option 1).
//
// Each step's plastic-strain increment is averaged over the neighbourhood and accumulated into
// the nonlocal plastic strain kbar:
//   kbar_e += sum_j a_ej * dpls_j,   a_ej = w(r_ej) V_j / sum_i w(r_ei) V_i,
// with the Gaussian weight w(r) = exp(-r^2 / l^2) cut off at r = 3 l (l: mat.nonlocal_length_in_m).
// The normalization renormalizes the weights at boundaries. Weakening reads
//   khat = m * kbar + (1 - m) * k
// (m: mat.nonlocal_m; m > 1 is the over-nonlocal form of Vermeer & Brinkgreve 1994), where k is
// the local plastic strain. kbar starts from the local k after the initial weak zone, so the
// prescribed seed is not smeared, and is a state variable: remapped at remesh, corrected with k
// at the surface, checkpointed.
//
// kbar is lagged by one step (update_stress reads the value accumulated after the previous
// step), so within a step only the local term moves the yield surface: update_stress scales the
// softening modulus by (1 - m).
//
// Only softening is regularized: where khat >= pls1 the softening modulus is 0, the averaging
// has nothing to act on, and flow there localizes as without this option.

// After the initial weak zone: kbar = k.
void init_nonlocal_plstrain(const Param& param, Variables& var);
// After update_stress: add the averaged increment of this step to kbar.
void accumulate_nonlocal_plstrain(const Param& param, Variables& var);

#endif
