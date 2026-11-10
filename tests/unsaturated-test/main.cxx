// Backend check for unsaturated.hpp: the same routines run serially on the
// host and inside a parallel loop (OpenACC on the device, or OpenMP), and the
// results are compared. Also reports throughput for saturated, unsaturated
// and mixed-material nodes. Exit status is nonzero on any failure.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "unsaturated.hpp"

namespace {

const int max_patch = 4;

// A node patch stored in flat arrays: entries [c*max_patch, c*max_patch + count).
struct FlatPatch
{
    const Retention *laws;
    const int *id;
    const double *w;
    int first, count;
    #pragma acc routine seq
    int size() const { return count; }
    #pragma acc routine seq
    const Retention &law(int k) const { return laws[id[first + k]]; }
    #pragma acc routine seq
    double weight(int k) const { return w[first + k]; }
};

const Retention law_table[] = {
    {retention_van_genuchten, 1 / 5e3, 1.6, 0.1},
    {retention_van_genuchten, 1 / 2e4, 3.0, 0.0},
    {retention_brooks_corey, 1 / 3e3, 0.8, 0.05},
    {retention_gardner, 2 / 9810.0, 0.0, 0.0},
    {retention_saturated, 1e-4, 2.0, 0.0},
};
const int nlaw = sizeof(law_table) / sizeof(law_table[0]);

struct Cases
{
    std::vector<int> id, count;
    std::vector<double> w, hmass, p, dm;
    int size() const { return static_cast<int>(p.size()); }
};

// kind 0: saturated nodes; 1: one unsaturated law; 2: mixed patches.
Cases make_cases(int kind, int n, unsigned seed)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> uni(0, 1);
    std::normal_distribution<double> nrm(0, 1);
    Cases c;
    c.id.assign(n * max_patch, 0);
    c.w.assign(n * max_patch, 0);
    for (int i = 0; i < n; ++i) {
        const double h = 1e-10 * (0.5 + uni(gen));
        int count = 1;
        double p, dm;
        if (kind == 0) {
            count = 2;
            c.id[i * max_patch] = 0;
            c.id[i * max_patch + 1] = 2;
            p = 1e6 * (1 + uni(gen));
            dm = 1e-5 * (2 * uni(gen) - 1);  // |dm/h| < 2e5 < p: stays saturated
        } else if (kind == 1) {
            c.id[i * max_patch] = i % 4;
            p = -std::pow(10.0, 6 * uni(gen));
            dm = 2.5e-4 * nrm(gen);
        } else {
            count = 4;
            for (int k = 0; k < 4; ++k) c.id[i * max_patch + k] = (i + k) % nlaw;
            p = (uni(gen) < 0.8 ? -std::pow(10.0, 5 * uni(gen)) : 3e3 * uni(gen));
            dm = 5e-4 * nrm(gen);
        }
        for (int k = 0; k < count; ++k) c.w[i * max_patch + k] = 0.01 + 0.04 * uni(gen);
        c.count.push_back(count);
        c.hmass.push_back(h);
        c.p.push_back(p);
        c.dm.push_back(dm);
    }
    return c;
}

// Updates every case into x, change and invalid; parallel unless serial.
void advance_all(const Cases &c, bool serial, std::vector<double> &x,
                 std::vector<double> &change, std::vector<int> &invalid)
{
    const int n = c.size();
    x.assign(n, 0);
    change.assign(n, 0);
    invalid.assign(n, 0);
    const Retention *laws = law_table;
    const int *id = c.id.data(), *count = c.count.data();
    const double *w = c.w.data(), *hm = c.hmass.data(), *p = c.p.data(), *dm = c.dm.data();
    double *xo = x.data(), *co = change.data();
    int *bad = invalid.data();
    if (serial) {
        for (int i = 0; i < n; ++i) {
            const FlatPatch patch = {laws, id, w, i * max_patch, count[i]};
            xo[i] = advance_content(patch, hm[i], p[i], dm[i], bad[i]);
            co[i] = content_change(patch, hm[i], p[i], xo[i]);
        }
        return;
    }
#ifdef ACC
    #pragma acc parallel loop gang vector copyin(laws[0:nlaw], id[0:n*max_patch], count[0:n], \
        w[0:n*max_patch], hm[0:n], p[0:n], dm[0:n]) copy(xo[0:n], co[0:n], bad[0:n])
#else
    #pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i) {
        const FlatPatch patch = {laws, id, w, i * max_patch, count[i]};
        xo[i] = advance_content(patch, hm[i], p[i], dm[i], bad[i]);
        co[i] = content_change(patch, hm[i], p[i], xo[i]);
    }
}

// se, dse, sat, kr and Bishop pressure for each law at each pressure.
void evaluate_all(const std::vector<double> &pv, bool serial, std::vector<double> &out)
{
    const int n = static_cast<int>(pv.size()), m = nlaw * n * 5;
    out.assign(m, 0);
    const Retention *laws = law_table;
    const double *p = pv.data();
    double *o = out.data();
    if (serial) {
        for (int j = 0; j < nlaw * n; ++j) {
            const Retention &r = laws[j / n];
            const double x = p[j % n];
            double dse;
            o[5*j] = retention_se(r, x, &dse);
            o[5*j+1] = dse;
            o[5*j+2] = retention_sat(r, x);
            o[5*j+3] = retention_kr(r, x);
            o[5*j+4] = bishop_pressure(r, x);
        }
        return;
    }
#ifdef ACC
    #pragma acc parallel loop gang vector copyin(laws[0:nlaw], p[0:n]) copy(o[0:m])
#else
    #pragma omp parallel for
#endif
    for (int j = 0; j < nlaw * n; ++j) {
        const Retention &r = laws[j / n];
        const double x = p[j % n];
        double dse;
        o[5*j] = retention_se(r, x, &dse);
        o[5*j+1] = dse;
        o[5*j+2] = retention_sat(r, x);
        o[5*j+3] = retention_kr(r, x);
        o[5*j+4] = bishop_pressure(r, x);
    }
}

int failures = 0;

void check(bool ok, const char *what, double value)
{
    std::printf("%-58s %-4s %.3g\n", what, ok ? "PASS" : "FAIL", value);
    if (!ok) ++failures;
}

double rel(double a, double b)
{
    if (a == b) return 0;
    return std::fabs(a - b) / std::max(std::fabs(a), std::fabs(b));
}

} // namespace

int main()
{
#ifdef ACC
    std::printf("backend: OpenACC\n");
#elif defined(_OPENMP)
    std::printf("backend: OpenMP\n");
#else
    std::printf("backend: serial (parallel loop compiled serially)\n");
#endif

    // Constitutive values across 11 decades of suction and the saturated side.
    std::vector<double> pv;
    for (int i = 0; i <= 400; ++i) pv.push_back(-std::pow(10.0, -3 + 11.0 * i / 400));
    for (double x : {0.0, 1e3, 1e6}) pv.push_back(x);
    std::vector<double> host, dev;
    evaluate_all(pv, true, host);
    evaluate_all(pv, false, dev);
    double worst = 0;
    bool finite = true;
    for (std::size_t j = 0; j < host.size(); ++j) {
        worst = std::max(worst, rel(host[j], dev[j]));
        finite = finite && std::isfinite(dev[j]);
    }
    check(finite && worst <= 1e-12, "retention values: max relative host/parallel difference", worst);

    const char *names[] = {"saturated", "unsaturated", "mixed"};
    for (int kind = 0; kind < 3; ++kind) {
        const Cases c = make_cases(kind, 200000, 17 + kind);
        std::vector<double> xh, ch, xd, cd;
        std::vector<int> bh, bd;
        advance_all(c, true, xh, ch, bh);
        advance_all(c, false, xd, cd, bd);
        double resid = 0, diff = 0;
        int bad = 0, inexact = 0;
        for (int i = 0; i < c.size(); ++i) {
            resid = std::max(resid, std::fabs(cd[i] - c.dm[i]));
            diff = std::max(diff, std::fabs(xd[i] - xh[i]) / (1 + std::fabs(xh[i])));
            bad += bd[i];
            if (kind == 0 && xd[i] != c.p[i] + c.dm[i] / c.hmass[i]) ++inexact;
        }
        char what[96];
        std::snprintf(what, sizeof what, "%s: invalid flags", names[kind]);
        check(bad == 0, what, bad);
        std::snprintf(what, sizeof what, "%s: max content residual (parallel)", names[kind]);
        check(resid <= 1e-15, what, resid);
        std::snprintf(what, sizeof what, "%s: max pressure host/parallel difference", names[kind]);
        check(diff <= 1e-9, what, diff);
        if (kind == 0) {
            check(inexact == 0, "saturated: cases differing from p + dm/hmass", inexact);
        }

        // Throughput: repeat the parallel update, including transfers for ACC.
        const int repeat = 5;
        const auto t0 = std::chrono::steady_clock::now();
        for (int r = 0; r < repeat; ++r) advance_all(c, false, xd, cd, bd);
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("%-58s %.1f ns/node\n", (std::string(names[kind]) + ": parallel throughput").c_str(),
                    1e9 * t / (repeat * c.size()));
    }

    // Non-finite input must be flagged, never silently accepted.
    {
        Cases c = make_cases(1, 4, 99);
        c.p[0] = NAN;
        c.dm[1] = INFINITY;
        std::vector<double> x, ch;
        std::vector<int> b;
        advance_all(c, false, x, ch, b);
        check(b[0] == 1 && (b[1] == 1 || !std::isfinite(x[1])), "non-finite input flagged", b[0]);
    }

    std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "OK", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
