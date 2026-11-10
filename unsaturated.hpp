#ifndef DYNEARTHSOL3D_UNSATURATED_HPP
#define DYNEARTHSOL3D_UNSATURATED_HPP

// Retention laws, Bishop pressure and the pointwise nodal content update.
// Pore pressure p is compression positive, suction is -p; alpha is in 1/Pa.
// The nonlinear step is one bracketed scalar root per node.

#include <algorithm>
#include <cfloat>
#include <cmath>

// Model 0 uses linear saturated storage.
enum RetentionModel { retention_saturated = 0, retention_van_genuchten = 1,
                      retention_brooks_corey = 2, retention_gardner = 3 };

// n is the van Genuchten n (> 1) or the Brooks-Corey pore-size index (> 0).
// sr is the residual saturation.
struct Retention
{
    int model;
    double alpha, n, sr;
};

// Pressure at and above which the material is saturated.
#pragma acc routine seq
inline double retention_p_sat(const Retention &r)
{
    return r.model == retention_brooks_corey ? -1 / r.alpha : 0.0;
}

#pragma acc routine seq
inline bool retention_wet(const Retention &r, double p)
{
    return r.model == retention_saturated || p >= retention_p_sat(r);
}

// Effective saturation Se and, if requested, dSe/dp >= 0, from one evaluation.
// On the unsaturated branch x = alpha*suction > 0 (> 1 for Brooks-Corey).
#pragma acc routine seq
inline double retention_se(const Retention &r, double p, double *dse = nullptr)
{
    if (retention_wet(r, p)) {
        if (dse) *dse = 0;
        return 1.0;
    }
    const double x = -r.alpha * p;
    double se, rate;  // rate = -dSe/dx
    switch (r.model) {
    case retention_van_genuchten: {
        const double y = std::pow(x, r.n);
        se = std::pow(1 + y, 1 / r.n - 1);
        rate = (r.n - 1) * se * y / (x * (1 + y));
        break;
    }
    case retention_brooks_corey:
        se = std::pow(x, -r.n);
        rate = r.n * se / x;
        break;
    default:
        se = std::exp(-x);
        rate = se;
    }
    if (dse) *dse = r.alpha * rate;
    return se;
}

#pragma acc routine seq
inline double retention_sat(const Retention &r, double p)
{
    return r.sr + (1 - r.sr) * retention_se(r, p);
}

// Mualem (van Genuchten), Burdine (Brooks-Corey) or exponential (Gardner).
#pragma acc routine seq
inline double retention_kr(const Retention &r, double p)
{
    if (retention_wet(r, p)) return 1.0;
    if (r.model == retention_van_genuchten) {
        // Se^(1/m) = 1/(1+y), so 1 - (1 - Se^(1/m))^m = -expm1(-m*log1p(1/y)),
        // accurate both near saturation and in the dry tail.
        const double y = std::pow(-r.alpha * p, r.n), m = 1 - 1 / r.n;
        const double inner = -std::expm1(-m * std::log1p(1 / y));
        return std::pow(1 + y, -0.5 * m) * inner * inner;
    }
    const double se = retention_se(r, p);
    return r.model == retention_brooks_corey ? std::pow(se, 3 + 2 / r.n) : se;
}

// Unsaturated-branch inverse of Se, for 0 < se < 1.
#pragma acc routine seq
inline double retention_p_of_se(const Retention &r, double se)
{
    switch (r.model) {
    case retention_van_genuchten:
        return -std::pow(std::expm1(-std::log(se) / (1 - 1 / r.n)), 1 / r.n) / r.alpha;
    case retention_brooks_corey:
        return -std::pow(se, -1 / r.n) / r.alpha;
    default:
        return std::log(se) / r.alpha;
    }
}

// Bishop equivalent pore pressure pi = Se*p (chi = Se), so that
// sigma' = sigma + alpha_B*pi. Equal to p when saturated.
#pragma acc routine seq
inline double bishop_pressure(const Retention &r, double p)
{
    return retention_se(r, p) * p;
}

// A node's lumped content is M(p) = hmass*p + sum_k w_k*S_k(p). Entry k is
// one (patch element, material) pair: w_k = f*V*phi/NODES_PER_ELEM with f the
// element's marker fraction of that material, and S_k is that material's law.
// Patch provides size(), law(k) and weight(k); entries with w_k = 0 (absent
// materials) are skipped, so they never set the saturation threshold.

// Pore content sum_k w_k*(1 - Sr_k)*Se_k(x) and, optionally, its x-derivative.
#pragma acc routine seq
template <class Patch>
inline double pore_content(const Patch &patch, double x, double *slope = nullptr)
{
    double c = 0, dc = 0;
    for (int k = 0; k < patch.size(); ++k) {
        const Retention &r = patch.law(k);
        if (r.model == retention_saturated || patch.weight(k) == 0) continue;
        const double w = patch.weight(k) * (1 - r.sr);
        double dse;
        c += w * retention_se(r, x, slope ? &dse : nullptr);
        if (slope) dc += w * dse;
    }
    if (slope) *slope = dc;
    return c;
}

// Highest saturation threshold among the patch's present unsaturated laws, or
// -HUGE_VAL when there is none: p and x at or above it are saturated for every
// entry, where the update is linear.
#pragma acc routine seq
template <class Patch>
inline double patch_p_sat(const Patch &patch)
{
    double p_sat = -HUGE_VAL;
    for (int k = 0; k < patch.size(); ++k) {
        const Retention &r = patch.law(k);
        if (r.model != retention_saturated && patch.weight(k) != 0)
            p_sat = std::max(p_sat, retention_p_sat(r));
    }
    return p_sat;
}

// Content change M(x) - M(p), for the fluid balance.
#pragma acc routine seq
template <class Patch>
inline double content_change(const Patch &patch, double hmass, double p, double x)
{
    return hmass * (x - p) + (pore_content(patch, x) - pore_content(patch, p));
}

// Pressure x with M(x) - M(p) = dm (hmass > 0). If p and trial are saturated,
// return trial; otherwise safeguarded Newton brackets p and p + dm/hmass,
// with bisection as fallback. Sets (never clears) invalid for nonfinite input
// or a residual that does not reach round-off.
#pragma acc routine seq
template <class Patch>
inline double advance_content(const Patch &patch, double hmass, double p, double dm, int &invalid)
{
    const double trial = p + dm / hmass;
    const double p_sat = patch_p_sat(patch);
    if (p >= p_sat && trial >= p_sat) return trial;
    double pores = 0, widest = 0;
    int guide = 0;
    for (int k = 0; k < patch.size(); ++k) {
        const Retention &r = patch.law(k);
        if (r.model == retention_saturated || patch.weight(k) == 0) continue;
        const double w = patch.weight(k) * (1 - r.sr);
        pores += w;
        if (w > widest) { widest = w; guide = k; }
    }

    const double c0 = pore_content(patch, p);
    // g is a difference of contents of size ~pores, so its round-off floor is
    // ~eps*pores; stopping there is what the fluid balance can resolve.
    const double tol = 4 * DBL_EPSILON * (std::fabs(dm) + pores);
    double lo = std::min(p, trial), hi = std::max(p, trial);
    // First guess: the dominant element's inverse with storage neglected.
    double x = 0.5 * (lo + hi);
    {
        const Retention &r = patch.law(guide);
        const double se = retention_se(r, p) + dm / pores;
        if (se > 0 && se < 1) x = std::min(std::max(retention_p_of_se(r, se), lo), hi);
    }
    for (int it = 0; it < 100; ++it) {
        double slope;
        const double g = hmass * (x - p) + (pore_content(patch, x, &slope) - c0) - dm;
        double next = x - g / (hmass + slope);
        // At the floor, the free Newton correction removes the one-sided
        // residual that would otherwise accumulate in the fluid balance.
        if (std::fabs(g) <= tol) return std::min(std::max(next, lo), hi);
        if (g < 0) lo = x; else hi = x;
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        if (next == x) return x;
        x = next;
    }
    invalid = 1;
    return x;
}

#endif
