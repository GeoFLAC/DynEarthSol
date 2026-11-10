# Plastic tunnel: Salençon comparison

**The associated-flow tunnel agrees with Salençon's solution and is a
worked example in the manual.** With 16 unloading steps and mixed stress,
the 0.05 m mesh gives hoop-stress RMS error 1.57% of $P$, displacement RMS
error 0.21% of the wall convergence and a plastic radius within 0.34%; all
errors fall as the mesh is refined. The radial strain error, 3.88% of
$P/G$, misses a 3% target and lies almost entirely within half a radius
of the wall. Two numerical limits remain and are stated on the page: the
result depends on the number of unloading steps, and without mixed stress
the coarse mesh diverges. Observations and a proposed mechanism are described
in [Cause of the two numerical limits](#cause-of-the-two-numerical-limits).

![Setup](../../benchmarks/poromechanics/results/runs/application-qualification-20261108/fig/salencon/setup.svg)

A quarter annulus has opening radius $a=1$ m and outer radius $R=20a$.
Initially all three normal stresses are $-P$, with $P=500$ kPa. Inner
support pressure decreases from $P$ to zero while outer pressure stays
fixed. Straight cuts have zero normal displacement. Gravity, fluid flow,
heat and mesh motion are off. Strain and displacement are measured from
the initially compressed state, rather than from an unloaded solid.

The material is perfectly plastic Mohr–Coulomb: $c=100$ kPa,
$\phi=30^\circ$, $G=10$ GPa and $K=65/3$ GPa, giving $\nu=0.3$.
Associated flow ($\psi=30^\circ$) and zero dilation ($\psi=0$) are
separate cases. There is no initial weak zone or strain softening.

## Reference and derivation

The infinite-domain stress and displacement reference is the cylindrical
solution of Salençon (1969), reproduced in the
[Itasca verification example](https://docs.itascacg.com/itasca940/3dec/docproject/source/examples/CylindricalHoleInAnIMCM.html).
The finite-annulus extension below removes outer-boundary truncation
error from the numerical comparison. Reference stress in this paragraph
is **compression positive**; DES output has the opposite sign.

Define the friction and dilation factors and compressive strength:

$$
N=\frac{1+\sin\phi}{1-\sin\phi},\qquad
M=\frac{1+\sin\psi}{1-\sin\psi},\qquad q=2c\sqrt N.
$$

Radial equilibrium, $s_r'+(s_r-s_\theta)/r=0$, and yield,
$s_\theta=Ns_r+q$, give, inside the plastic radius $b$,

$$
C=\frac{P_i+q/(N-1)}{a^{N-1}},\qquad
s_r=Cr^{N-1}-\frac{q}{N-1},\quad
s_\theta=NCr^{N-1}-\frac{q}{N-1}.
$$

In the elastic exterior, write $s_r=A-B/r^2$ and $s_\theta=A+B/r^2$.
Matching radial stress $s_b$ at $b$ and pressure $P$ at $R$ gives

$$
A=\frac{PR^2-s_b b^2}{R^2-b^2},\qquad
B=\frac{b^2R^2(P-s_b)}{R^2-b^2}.
$$

The interface is also at yield, so solve the scalar equation

$$
Cb^{N-1}-\frac{q}{N-1}
=\frac{(2P-q)R^2+qb^2}{(N+1)R^2-(N-1)b^2}.
$$

Here $b=1.39575056$ m. The infinite-domain value is $1.39405010$ m,
obtained by taking $R\to\infty$ in the same equation.

For displacement, return to tension-positive stress $\sigma=-s$ and
write $\lambda=K-2G/3$. Elastic displacement relative to the initial
compression is

$$
u_e(r)=\frac{(P-A)r}{2(\lambda+G)}-\frac{B}{2Gr}.
$$

On the intermediate out-of-plane principal-stress branch,
$\epsilon_r^p=-M\epsilon_\theta^p$ and $\epsilon_y^p=0$.
Plane-strain Hooke's law gives the elastic increments

$$
e_r^e=\frac{(1-\nu)(\sigma_r+P)-\nu(\sigma_\theta+P)}{2G},
\qquad
e_\theta^e=\frac{(1-\nu)(\sigma_\theta+P)-\nu(\sigma_r+P)}{2G}.
$$

Thus $u'+Mu/r=e_r^e+Me_\theta^e$. Define
$a_1=1-\nu-M\nu$ and $a_2=M(1-\nu)-\nu$, followed by

$$
D_0=\frac{(1+M)(1-2\nu)[P+q/(N-1)]}{2G},
\qquad D_1=\frac{C(a_1+Na_2)}{2G}.
$$

Integrating and matching $u(b)=u_e(b)$ gives

$$
J(r)=\frac{D_0r^{M+1}}{M+1}-\frac{D_1r^{M+N}}{M+N},
\qquad
u_p(r)=\frac{b^Mu_e(b)+J(r)-J(b)}{r^M}.
$$

Finally $\epsilon_r=u'$, $\epsilon_\theta=u/r$ and
$\sigma_{yy}=-P+\nu(\sigma_r+\sigma_\theta+2P)$ on this branch.
The reference checks boundary traction, yield, equilibrium, displacement
derivatives, matching and the infinite-domain displacement for three
dilation angles. It rejects parameters that leave the assumed principal
branch. The comparison uses the finite solution, rather than fitting a
material parameter or plastic radius.

## Current spatial evidence

![Finest associated-flow result and load-increment sensitivity](../../benchmarks/poromechanics/results/runs/salencon-audit-20261108/convergence/fig/salencon/result.svg)

An additional diagnostic profile uses the finest executed unstructured mesh,
0.025 m (39,999 nodes and 79,207 triangles), with 32 unloading increments.
It completes loading; radial-strain RMS error is 4.44%, against the 3%
target. The guided manual instead selects the 0.05 m, 16-step case described
above, keeping the three-mesh sequence at one increment count. The
64-increment and fixed-load controls below retain the limits of agreement;
neither choice establishes load-increment convergence.

The radial-strain error is concentrated near the opening: 99.998% of its
squared error lies at $r<1.5a$. For $2a\le r\le5a$, radial-strain RMS
error is 0.0162% of $P/G$. This localization explains the close overall
profiles; the unchanged acceptance band remains $1.05a\le r\le5a$.

<details>
<summary>Earlier spatial study at 16 unloading increments</summary>

For associated flow with 16 unloading increments, mixed stress enabled
and graded refinement in $r<3a$, the frozen solver after both material
corrections gives the following errors. All three meshes finish loading.
Stress and strain RMS errors use every element centroid in
$1.05a\le r\le5a$; displacement uses all nodes in $a\le r\le5a$.

| Inner spacing, m | Elements | Hoop stress RMS / $P$ | Out-of-plane stress RMS / $P$ | Displacement RMS / wall contraction | Radial strain RMS / $(P/G)$ |
| --- | --- | --- | --- | --- | --- |
| 0.2 | 1,321 | 8.83% | 3.28% | 2.66% | 18.02% |
| 0.1 | 5,106 | 4.88% | 1.85% | 0.808% | 9.80% |
| 0.05 | 19,898 | 1.57% | 0.613% | 0.213% | 3.88% |

The finest plastic radius is 1.39106 m, versus 1.39575 m analytically
(0.336% error). Its accepted PT force residual is
$1.34\times10^{-13}PR$; force assembled independently from the saved,
objectively rotated stress is $6.47\times10^{-9}PR$. These sample different
points in the physical step. The remaining failed field criterion is
radial strain, whose target remains 3%; no field tolerance was loosened.
These are fixed-increment spatial results, not evidence of load-increment
convergence. Zero-dilation qualification is separate.

</details>

### Strain quantity and load-increment check

The stored DES strain is a tensor total strain relative to the compressed
initial state. Its shear entry is $\epsilon_{xz}$, not engineering shear.
`update_strain_rate` first forms the symmetric P1 velocity gradient;
`compute_dvoldt`, `compute_edvoldt` and `update_stress` replace its element
divergence by a nodally averaged value. For an element increment
$\mathbf E_e=\operatorname{sym}\nabla\Delta\mathbf u_e$, define

$$
\theta_n=\frac{\sum_{e\ni n} A_e\operatorname{tr}\mathbf E_e}
                 {\sum_{e\ni n} A_e},\qquad
\bar\theta_e=\frac13\sum_{n\in e}\theta_n,\qquad
\bar{\mathbf E}_e=\mathbf E_e+
 \frac{\bar\theta_e-\operatorname{tr}\mathbf E_e}{2}\mathbf I_2.
$$

The physical step then rotates accumulated strain using the accepted
velocity's spin. This strain correction is also present when the optional
mixed-stress projection is disabled. An independent reconstruction from
every saved displacement increment reproduces stored strain to
$6.78\times10^{-15}$ in units of $P/G$ across the three original meshes.
An independent central difference of analytical displacement reproduces
$\epsilon_r=du/dr$ to $3.15\times10^{-10}$ in the same units. The hoop
comparison uses $\epsilon_\theta=u/r$; $P/G=5\times10^{-5}$. The reported
radial-strain error is therefore not a component or normalization mistake.

Every associated-flow sequence below reaches force equilibrium, but at
fixed spacing the field errors increase with more increments. The material is rate independent,
so these are refinements of the support-pressure increment $P/n$; there
is no viscous or thermal timescale to change by extending the run.

| Inner spacing, m | Unloading increments | Hoop stress RMS / $P$ | Displacement RMS / wall contraction | Radial strain RMS / $(P/G)$ |
| --- | --- | --- | --- | --- |
| 0.05 | 16 | 1.57% | 0.213% | 3.88% |
| 0.05 | 32 | 3.09% | 0.372% | 8.91% |
| 0.05 | 64 | 6.71% | 1.559% | 21.06% |
| 0.025 | 32 | 1.24% | 0.113% | 4.44% |
| 0.025 | 64 | 4.71% | 1.197% | 16.74% |

The numerical displacement difference between successive sequences grows
from 0.272% to 1.386% of the reference wall contraction. Thus this mesh
does not demonstrate load-increment convergence, despite small force
residuals. Refining the mesh must be assessed at matching increment counts;
selecting only the closest-fitting sequence would conceal this sensitivity.

At the same 64 increments, halving the spacing to 0.025 m reduces hoop
error from 6.71% to 4.71% and radial-strain error from 21.06% to 16.74%.
It improves agreement but does not meet the unchanged field criteria.

Zero dilation is a separate regime. With mixed stress enabled at 0.05 m,
the 32-increment sequence stops at the final equilibrium solve; the
64-increment sequence stops at increment 56. Both report exit 54
(nonconvergence). An ordinary-explicit restart of the failed final
32-increment step reaches a stationary state under the same tractions.
After 20,000 relaxation steps the maximum speed is already small,
$6.37\times10^{-8}$ m/s, but independently assembled force still reaches
58.2 N/m. Extending the relaxation lowers that force below 0.005 N/m;
a subsequent 4,000-step hold finishes normally, stays below 0.005 N/m
and changes stress by at most $1.8\times10^{-7}P$. The force acceptance
limit is $10^{-7}PR=1$ N/m. Thus low speed alone is insufficient, but
the settled control establishes an equilibrium that this PT solve misses.

At that state hoop and radial stress RMS errors are 3.53% and 1.41% of
$P$, and radial strain RMS error is 3.30% of $P/G$. This is a single-mesh
zero-dilation control, rather than a spatial-convergence qualification.
It uses inertial scaling $10^5$, characteristic speed $10^{-5}$ m/s,
automatic timesteps, mixed stress and zero final support pressure.
Its relaxation clock is numerical; it does not predict excavation time.
The associated-flow profiles above do not qualify every zero-dilation
loading path. A geometry-aligned radial-mesh control also completes
associated loading but has larger errors; alignment alone does not
resolve the sensitivity. All controls remain in the audit summary.

<details>
<summary>Fixed-load holding and zero-strain return diagnostics</summary>

### Constant-load control

![Increment sensitivity and fixed-load holding](../../benchmarks/poromechanics/results/runs/salencon-audit-20261108/convergence/fig/salencon/controls.svg)

Restart the original 0.05 m, 16-increment final state and advance four
steps with both wall tractions unchanged. Coordinates, stress, strain and
plastic strain match the saved source exactly at restart; all physical
configuration sections are identical. The four hold steps run without
further restarts. There is no viscosity, fluid flow, heating, weakening or
changed load to prescribe additional deformation.

The maximum additional wall contraction is 2.31%, 4.02%, 5.41% and 6.61%
of the analytical wall contraction after the four steps. Its wall-node
mean reaches 5.05%. All four independently assembled force residuals
remain below $4.35\times10^{-9}PR$. Thus the small force residual does
not certify a stationary history for this discrete state.

A driver-only diagnostic separates the existing material and projection
kernels on the completed coarse HDF5 control (834 elements). At exactly
zero strain, the local associated-flow return agrees with an independent
elastic-energy projection onto all six Mohr–Coulomb halfspaces and three
tensile caps to $4.66\times10^{-16}P$; all returned stresses satisfy yield.
The reference enumerates up to three independent active constraints and
checks nonnegative multipliers and every inequality, rather than copying
the production branch selection.

The input to that return is the previously saved averaged stress. Returning
it changes stress by up to $0.257P$ and creates a force imbalance of
$1.80\times10^{-3}PR$ without any strain increment. The following nodal
projection, also reconstructed independently, leaves
$1.30\times10^{-3}PR$. A subsequent equilibrium solve therefore has a
nonzero correction to make even at unchanged tractions. This isolates a
readjustment mechanism in the interaction between projected stress and
constitutive history in this configuration. It does not identify a local
return error, establish a defect in every mixed-stress use, or invalidate
earlier verification. Uniform affine material tests cannot expose this
heterogeneous history interaction.

</details>

### What the earlier failure meant

The earlier report combined different configurations and revisions.
Without mixed stress, the current coarse PT case becomes non-finite at
stage 14 of 16; ordinary explicit relaxation of that same load increment
also diverges, including with a tenfold smaller timestep. Changing that
increment to elastic response converges in both methods. This narrows the
unstable configuration but does not establish its cause or invalidate
DES's existing verification.

With mixed stress enabled, ordinary explicit relaxation of the same
plastic load step converges: the independently assembled final residual
is $2.59\times10^{-9}PR$. Parent and current solvers preserve all 2,355
physical payloads in 81 files exactly for that control. PT also completes
the full 16-stage unloading. The coarse ungraded mesh has large field
errors; the graded sequence above shows their reduction. Solver
nonconvergence, completed runs with excessive mesh error, and unqualified
loading regimes must be reported separately.

Mixed stress is an existing optional formulation. Its hydrostatic
projection changes the saved stress after the constitutive return.
Evaluating local Mohr–Coulomb slack on that averaged stress is therefore
a diagnostic, not a direct test that the local return mapping failed.
The checker retains the measured slack but applies its local yield gate
only without mixed stress. Independent affine tests below check the
constitutive return directly. The earlier blanket attribution to a
mixed-stress bug was not supported by the comparison.

<details>
<summary>Material-point checks and earlier failure controls</summary>

An independent four-node affine material test identifies one increment
size issue without a free mechanical degree of freedom. With a single
increment, the parent solver clips an elastic predictor at the tensile cutoff before the
shear correction even though the reference final stress is compressive.
The largest stress errors are 56.5% / 96.8% of $P$ for zero/associated
dilation. Splitting the same strain into 16 or 64 increments recovers the
reference within $5.4\times10^{-15}P$. This qualifies those refined
material paths in the parent only. The corrected solver now checks the
uncut shear return first and recovers the reference even in one increment:
the largest errors are $3.50\times10^{-16}P$ / $5.83\times10^{-16}P$
for zero/associated dilation. All six 1/16/64-increment cases pass in
binary and HDF5, with maximum error $5.36\times10^{-15}P$.
The driver now requires the large-increment comparison to pass as well
as the refined paths. This alone does not establish spatial accuracy.

## Two-face shear-edge return

A second four-node test imposes equal in-plane extensions
$\epsilon_{xx}=\epsilon_{zz}=1.2\times10^{-5}$, with
$\epsilon_{yy}=0$, on the same initially compressed material. Every
velocity is prescribed, mesh motion and mixed-stress smoothing are off,
and the motion has no spin. There is no free force residual to converge.
The configured tensile limit is 1 GPa; the actual cap is also limited to
$c/\tan\phi=173.205$ kPa. The trial maximum is 100 kPa and does not reach
that cap. Symmetry requires
$\sigma_{xx}=\sigma_{zz}$; these are the two equal maximum principals.

The independent reference enforces both active yield faces. In the order
$(xx,zz,yy)$, let $\mathbf C=\lambda\mathbf 1\mathbf 1^T+2G\mathbf I$,
and form the elastic trial stress from the initial $-P\mathbf1$.
The two rows of $\mathbf A$ are $(N,0,-1)$ and $(0,N,-1)$; the rows of
$\mathbf B$ replace $N$ by the dilation factor $M$. Solve

$$
(\mathbf A\mathbf C\mathbf B^T)\boldsymbol\eta
  =\mathbf A\boldsymbol\sigma^{\mathrm{tr}}-q\mathbf1,
\qquad
\boldsymbol\sigma=\boldsymbol\sigma^{\mathrm{tr}}
  -\mathbf C\mathbf B^T\boldsymbol\eta.
$$

Both multipliers are positive, and the resulting $yy$ stress remains the
minimum principal, so the assumed two-face branch is consistent. This
small two-by-two solve is a reference, separate from the production return.

For zero dilation the exact $(xx,zz,yy)$ stress is
$(81.282,81.282,-102.564)$ kPa. The one-step DES result is
$(100.000,76.603,-116.603)$ kPa in the parent solver. Its actual sorted principals give
$f=\sigma_{\min}-N\sigma_{\max}+q=-70.192$ kPa, whereas admissibility
requires $f\ge0$. The parent associated-flow case also fails. Using 16 or
64 increments reduces the error but leaves both cases outside yield.
Thus the single-face correction can reorder the principals without
returning to the adjacent face/edge. This material defect is independent
of stress smoothing; it does not explain every spatial error by itself.

The corrected solver retains the single-face arithmetic when its returned
principals remain ordered. Otherwise one shared routine in `rheology.cxx`
tries the two adjacent active shear faces, checks nonnegative multipliers,
ordering and the tensile cap, and restores the repeated eigenvalue without
an arbitrary preferred axis. Plastic increments use the combined flow
directions and each existing branch's strain measure. The flow convention
is the [Mohr–Coulomb nonassociated potential](https://docs.itascasoftware.com/itasca980/common/models/mohr/doc/modelmohr.html);
the simultaneous two-face solve is the tested extension here.

All six right-edge cases (both dilation angles, 1/16/64 increments) now
match the independent reference within $1.17\times10^{-16}P$ in binary and
HDF5, without a yield violation at the test tolerance. Direct production
tests also cover both repeated-principal edges and rotated directions in
2D plane strain and 3D. The current serial C++20 and OpenMP C++17 suites
pass 7,579 / 11,171
checks. The two dry core comparisons remain bit-exact, and 340 existing
pressure-strength files retain 11,091 physical payloads exactly.

This qualifies the constant-strength two-face component. The inherited
within-step hardening tangent is retained; evolving-strength accuracy,
general apex and mixed tensile/shear returns remain unqualified. The
additional tensile-predictor correction below addresses a separate
branch-selection error without qualifying those remaining branches.

## A tensile predictor can return by shear

The elastic trial's maximum alone does not determine whether plastic
tensile flow is needed. Before clipping, both dimensional implementations
now test the uncut shear return. It is accepted only with a nonnegative
plastic multiplier, ordered principals, satisfied shear yield and a final
maximum below the tensile cap. Otherwise the remaining composite/tensile
branches retain their existing behavior. The same shear routine supplies
stress and plastic history in both paths.

An associated trial with sorted principals $(169,185,200)$ kPa has the
same $c=100$ kPa, $\phi=\psi=30^\circ$, $K$ and $G$ as above.
Its elastic maximum exceeds $T=173.205$ kPa, yet the independent two-face
solution is $(165.347657,170.585939,170.585939)$ kPa. Both shear
multipliers are positive and the tensile constraint is inactive. The
parent plane-strain clipping and 3D domain selector miss this state;
the corrected production test checks all stresses and the plastic
increment in both dimensions. The ordinary pressure-strength cases keep
all 11,091 physical payloads unchanged.

The two-face-only spatial controls remain unsuccessful. With no mixed-stress smoothing,
both parent and corrected solvers become non-finite at stage 14 of 16 on
the 0.2 m mesh. Reducing `PT_CFL` from 0.25 to 0.1 does not change that
failure stage. The corrected 128-increment run stops at stage 109 with
nonconvergence. These are numerical failures, not evidence of physical
tunnel collapse; the analytical equilibrium still exists.

Before the tensile-predictor correction, the mixed-stress run on the
graded 0.2 m mesh stopped at stage 15. After both corrections it completes,
as do the two finer meshes reported above. The earlier failure remains in
`mc-return-20261108`; it is not the current mixed-stress result.

After the tensile-predictor correction, the unsmoothed 16-stage coarse
case still becomes non-finite at stage 14 with PT option 1. PT option 0
stops at stage 12 with nonconvergence. The corrected local material paths
therefore establish two constitutive improvements, not a spatial
qualification or a claim that changing the PT option resolves the tunnel.

The checker still reports saved-stress yield slack in
$1.05a\le r\le5a$. For mixed stress, that quantity describes projected
stress and is not used as a local return gate. At the finest current mesh,
374 of 8,496 projected stresses have negative slack at $10^{-8}P$;
that observation alone does not establish a constitutive defect.
For slope applications, the pressure sampled by mechanics and the stress
projection must likewise be identified before interpreting a yield audit
or attributing nonconvergence to physical collapse.

</details>

## Cause of the two numerical limits

Without mixed stress, stage 14 of the 0.2 m mesh diverges in PT and in
ordinary explicit relaxation alike. Saving the explicit run every 50 steps
shows how. The pressure starts smooth; within a few hundred steps the
triangles with an edge on the wall fall toward tension while their
neighbours with a vertex on the wall rise in compression. By step 1013
the edge elements sit at the Mohr–Coulomb apex (+173 kPa, $c/\tan\phi$,
no shear strength) beside neighbours at −11 MPa, about twenty times the
in-situ stress. The alternating element pressure and small nodal force
residual are consistent with a spurious pressure mode in the nodally averaged
volumetric formulation. Pressure-dependent strength offers a growth
mechanism: elements with less compression are weaker and yield more. The
out-of-plane stress remains the intermediate principal until the apex, so
the two-face returns above are not involved.

Mixed stress smooths the element pressure after the plastic return and so
suppresses this mode. The smoothed stress, however, is no longer on the
yield surface, and the next step returns it again. This is a plausible
mechanism for the fixed-load drift and increasing error with more steps
(hoop 1.57%, 3.09% and 6.71% for 16, 32 and 64 steps at 0.05 m).
Post-return projection is an existing optional formulation; departure from
local yield after that projection does not demonstrate a failed local
return. A causal correction has not been qualified. Trial-pressure
projection or a separate stress for forces are possible experiments,
requiring the fixed-load hold, 16/32/64-step sequence and unsmoothed stage
14 controls, as well as DES regressions. These observations do not yet
establish one common cause or a remedy for both limits.

## Reproduce

From the repository, use fresh output directories. The tunnel command is
expected to return nonzero while qualification fails.

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
tunnel_solver=benchmarks/poromechanics/results/runs/shear-cap-20261108/candidate/dynearthsol2d
tunnel_manifest=benchmarks/poromechanics/results/runs/shear-cap-20261108/qualification.json
tunnel_root=benchmarks/poromechanics/results/runs/salencon-audit-20261108
python3 benchmarks/poromechanics/salencon_reference.py
python3 tests/functional/hydromechanics/check_plastic_path.py \
  "$tunnel_solver" "$tunnel_root/replay-plastic-path" \
  --binary-manifest "$tunnel_manifest"
python3 tests/functional/hydromechanics/check_salencon.py \
  "$tunnel_solver" "$tunnel_root/replay-salencon" \
  --spacings .2 .1 .05 --dilation 30 --mixed --refined --load-steps 16 \
  --binary-manifest "$tunnel_manifest"
```

The edge diagnostic is small and 2D. It now passes in strict mode with
the corrected frozen solver:

```sh
python3 tests/functional/hydromechanics/check_plastic_path.py \
  "$tunnel_solver" "$tunnel_root/replay-mc-edges" --edges \
  --binary-manifest "$tunnel_manifest"
```

The parent failure cases are under
`application-qualification-20261108/mc-edge-final/{binary,hdf}/`; the existing
spatial data's separate admissibility audit is `salencon-yield-audit.json`.
Corrected binaries, cases, source patch and hash-bound qualification are
under `mc-return-20261108/` (`edges-2d`, `edges-hdf-2d`, `qualification.json`).
The subsequent predictor correction is in `shear-cap-20261108/`:
`verified-{binary,hdf}-{paths,edges}` retain the strict material checks;
`qualification.json` records the frozen binaries and unsuccessful tunnel
controls. It supersedes the parent large-increment failure as the current
material result, while retaining that failure for comparison.

The initial audit is `salencon-audit-20261108/audit.json`. It binds the
earlier control scripts, solver hashes and raw outputs. `replay.py` independently
recomputes the three graded-mesh reports and the qualification decision
without rerunning the solver. The original reports remain unchanged;
`refinement-replay/qualification.json` records the corrected interpretation
of mixed-stress slack. The subsequent `convergence/` directory contains
the paired load/mesh studies, strain reconstruction, fixed-load hold and
driver-only return diagnostic. Its `audit.json` binds the commands,
source/binary hashes and raw artifacts; `summary.json` lists every
completed or failed configuration. The selected plots can be redrawn:

```sh
python3 "$tunnel_root/convergence/summarize.py"
python3 "$tunnel_root/convergence/plot.py"
```

For a fresh 64-increment comparison, use a new output directory:

```sh
python3 tests/functional/hydromechanics/check_salencon.py \
  "$tunnel_solver" "$tunnel_root/replay-steps64" \
  --spacings .05 .025 --dilation 30 --mixed --refined --load-steps 64 \
  --binary-manifest "$tunnel_manifest"
```

Remaining acceptance work is near-wall strain accuracy, load-increment
convergence, a settled fixed-load history and zero-dilation qualification.
The alternating pressure mode narrows the instability diagnosis;
its proposed mechanism and a correction remain unqualified. The elastic [Kirsch component](kirsch.md) and the
[sealed poroelastic borehole study](borehole.md) have separate evidence;
broader injection/loading modes remain open.
