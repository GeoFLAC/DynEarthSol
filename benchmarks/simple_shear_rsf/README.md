# Local EP/RSF Benchmark

This directory reproduces the constitutive verification described in the DES
v2.0 paper revision. The main benchmark has nine two-element simple-shear
cases:

- static-friction elastoplasticity at 10, 20, and 30 degrees
- steady-state RSF at \(V_0=10^{-6},10^{-5},10^{-4}\) m/s
- aging-law RSF at \(D_c=10^{-3},3\times10^{-3},10^{-2}\) m

Each case uses a fixed 0.01 s step for 200,000 steps. The top boundary moves at
\(10^{-5}\) m/s, giving a final time of 2000 s and an engineering shear strain
of 0.02. The RSF cases explicitly select the total-strain-invariant rate
measure. In this deforming cell,

\[
V(t)=\frac{v_x}{\sqrt{(1-v_x t/H)^2+1}}.
\]

The runner can also reproduce the paper's separate zero-velocity healing test:
1546 fixed steps of three days with \(D_c=10^{-2}\) m and
\(V_0=4\times10^{-9}\) m/s.

## Files

- simple_shear_base.cfg: shared CFG template
- run_simple_shear_benchmark.py: case generation and execution
- benchmark_reference.py: paper reference solutions and monitor readers
- plot_simple_shear_benchmark.py: paper-layout figure and metrics CSV
- check_simple_shear_benchmark.py: non-plotting regression check

## Full Paper Run

~~~bash
cd benchmarks/simple_shear_rsf
python3 run_simple_shear_benchmark.py --clean --aging-transient --healing
python3 plot_simple_shear_benchmark.py -o simple_shear_benchmark.pdf
~~~

This creates nine main case directories, three higher-frequency aging-law
monitor directories for the state inset, and one healing directory under
runs/. The plotting command creates the requested figure and
simple_shear_benchmark_metrics.csv.

The stress comparison is the two-element mean absolute shear stress. At every
sample after \(t=0\), the reported error is

\[
\frac{|\overline{\sigma}_{xy}^{DES}
      -\overline{\sigma}_{xy}^{ref}|}
     {\overline{\sigma}_{xy}^{ref}}.
\]

## Regression Check

The default check runs one steady-state case, one aging-law case, and the
zero-velocity healing test:

~~~bash
cd benchmarks/simple_shear_rsf
python3 check_simple_shear_benchmark.py
~~~

Run all nine simple-shear cases plus healing with:

~~~bash
python3 check_simple_shear_benchmark.py --all
~~~

The checker verifies pointwise stress error, recovers the rate used by the RSF
law from the monitored friction and state outputs, and iterates the discrete
zero-velocity healing update independently. Temporary results are removed
unless --keep-output is supplied.

The default maximum pointwise stress-error tolerance is \(2\times10^{-5}\) as
a fraction. The paper reports a largest pointwise error of
\(1.16\times10^{-3}\) percent.

## Check Existing Results

~~~bash
python3 plot_simple_shear_benchmark.py \
  --check --no-save --max-relative-error 2e-5
~~~

NumPy and Matplotlib are required only for plotting; the runner and regression
checker use the Python standard library.
