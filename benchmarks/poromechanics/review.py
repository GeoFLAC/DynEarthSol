"""Setup drawings and reference derivations for the portable benchmark review.

Each admitted result needs a setup, reference equations and a stated validity
limit. Math is rendered locally so the HTML remains readable without a CDN.
The drawings are explanatory diagrams, never meshes or numerical results.
"""

import html
import json
import re
import shutil

import prepare_figures as pf
import des_options
from manual_pages import MANUAL
from setup_drawings import DRAW, schematic


SPECS = {
    'terzaghi': dict(
        title='Terzaghi consolidation',
        kind='vertical', top='Instantaneous load; drained top',
        body='Saturated\nelastic column', base='Sealed base; vertical displacement fixed',
        geometry='H = 10 m; width = 1 m\n2D plane strain; drawing not to scale',
        labels=[('Load and drainage', 'At $t=0^+$ a compression $P=100$ kPa is applied; the top is drained, $p=0$.'),
                ('Constraints', 'Sides: zero normal displacement and zero normal water flux. Base: vertical displacement fixed, water flux zero.'),
                ('Initial state and comparison', 'Initially unloaded. The reference starts just after the instantaneous undrained load. Compare pressure and downward top settlement.')],
        equations=[
            r'S_0=\phi/K_f+(\alpha-\phi)/K_s,\quad H_c=K+4G/3,\quad S_c=S_0+\alpha^2/H_c',
            r'c_v=\kappa/(\mu S_c),\quad k_m=(2m+1)\pi/(2H),\quad p_0=\alpha P/(H_cS_c)',
            r'p(d,t)=p_0\sum_{m=0}^{\infty}\frac{4\sin(k_md)}{(2m+1)\pi}\exp(-c_vk_m^2t)',
            r'\delta(t)=\frac{PH}{H_c}-\frac{\alpha Hp_0}{H_c}\sum_{m=0}^{\infty}\frac{8\exp(-c_vk_m^2t)}{[(2m+1)\pi]^2}'
        ],
        definitions=r"Depth $d$ runs downward from the drained top ($0$) to the base ($H$). $P$ is positive compression and $\delta$ positive downward settlement; $\kappa$ is intrinsic permeability, $\mu$ viscosity, $\phi$ porosity, and $\alpha$ the Biot coefficient. $K$, $G$, $K_f$ and $K_s$ are the drained bulk, shear, fluid and grain moduli.",
        limit='Small-strain, homogeneous, saturated, laterally confined elasticity; constant properties, no gravity or heating. Frame 0 is written before the load is applied and is not compared. This example uses the PT solver; a sudden load with GVS is not tested here.',
        derivation=[
            (r'Equilibrium makes the total vertical stress uniform: $\sigma_{zz}=-P$. With tension-positive strain, $H_c\varepsilon_{zz}-\alpha p=-P$. Substitute this strain into the fluid content $\zeta=\alpha\varepsilon_{zz}+S_0p$.',
             r'\zeta=-\alpha P/H_c+S_cp,\qquad S_cp_t=(\kappa/\mu)p_{dd}'),
            (r'Zero fluid-content change during the instantaneous load gives $p_0$. Thereafter $p(0,t)=0$ and $\partial_dp(H,t)=0$. The mixed-boundary eigenfunctions are $\sin(k_md)$; expanding the uniform initial pressure gives the coefficients above. Integrating compressive strain gives the settlement, including its instantaneous undrained part.',
             r'\delta(t)=-\int_0^H\epsilon_{zz}\,dd=PH/H_c-(\alpha/H_c)\int_0^Hp(d,t)\,dd')],
        reference='The series and settlement follow directly from the stated one-dimensional Biot initial/boundary problem; evaluated in terzaghi_reference.py.'),
    'gardner': dict(
        title='Gardner transient infiltration', kind='vertical', top='Downward infiltration I = 0.9 Ks',
        body='Partly saturated\nGardner soil', base='Water table: p = 0',
        geometry='H = 1 m; width = 0.05 m\n3D column; drawing not to scale',
        labels=[('Water boundaries', 'Constant downward inflow starts at the top at t = 0. Base pressure is zero; all side faces are sealed.'),
                ('Initial profile and gravity', r'Hydrostatic initial pressure $p=-\rho_fgz$, with elevation $z$ above the base. Gravity points downward.'),
                ('Mechanical conditions', 'Side-normal displacement and base vertical displacement fixed; traction-free top. Compare K/Ks derived from pressure. Bishop coupling is not tested here (see the Bishop pages).')],
        equations=[
            r'w=S_e=k_r=\exp(\alpha_p p),\quad A=\alpha_p\rho_fg,\quad Z=Az,\quad L=AH,\quad T=AK_st/\Delta\theta',
            r'w_b(Z)=b+(1-b)e^{-Z},\quad a=I_a/K_s,\quad b=I_b/K_s,\quad \Delta\theta=\phi(1-S_r)',
            r'w=w_b-4(b-a)e^{(L-Z)/2-T/4}\sum_{n=1}^{\infty}\frac{\sin(\lambda_nZ)\sin(\lambda_nL)e^{-\lambda_n^2T}}{1+L/2+2\lambda_n^2L}',
            r'\sin(\lambda_nL)+2\lambda_n\cos(\lambda_nL)=0,\qquad p=\log(w)/\alpha_p'
        ],
        definitions=r"$K_s=\kappa\rho_f g/\mu$ is saturated hydraulic conductivity. $\alpha_p$ is the pressure retention coefficient in $\mathrm{Pa^{-1}}$ and $A$ its head-based value in $\mathrm{m^{-1}}$. Here $A=10\ \mathrm{m^{-1}}$, $I_a=0$, $I_b=2.5\times10^{-6}\ \mathrm{m/s}$, $K_s=2.77778\times10^{-6}\ \mathrm{m/s}$, $\phi=0.4$ and $S_r=0$. $\lambda_n$ is the positive root in $((n-1/2)\pi/L,n\pi/L)$.",
        limit='This exact series solves ideal Gardner Richards flow with retention storage only. DES also retains elastic/compressible storage, so it is a reference-limit comparison, not an exact solution of every active DES term. The displayed error includes that model difference. It does not validate ponding or a seepage-face switch.',
        derivation=[
            (r'For unsaturated $p\le0$, $\theta=\theta_r+\Delta\theta\,w$ and the upward water flux is $q_z=-K_s(w_z/A+w)$. Substitution into $\theta_t+\partial_zq_z=0$ and nondimensionalization give a linear drift-diffusion equation.',
             r'w_T=w_{ZZ}+w_Z,\quad w(0,T)=1,\quad w_Z(L,T)+w(L,T)=b'),
            (r'Before the flux change the column is at $w_a=a+(1-a)e^{-Z}$. Subtract the new steady solution wb and remove the drift. Separation of variables gives a sine series and the Robin eigenvalue equation; projecting the initial difference supplies the coefficient in the reference above.',
             r'w-w_b=e^{-Z/2-T/4}f,\quad f_T=f_{ZZ},\quad f(0,T)=0,\quad f_Z(L,T)+f(L,T)/2=0')],
        reference='Srivastava and Yeh (1991), homogeneous-column solution; evaluated in unsaturated_reference.py. The reference implementation is independently checked against a matrix exponential.',
        url='https://tian.hwr.arizona.edu/publications/data/a9ab4ba8091091d98c03d35241dcbf0c.pdf'),
    'van_genuchten': dict(
        title='van Genuchten steady infiltration', kind='vertical', top='Downward infiltration I = 0.5 Ks',
        body='van Genuchten\nretention + Mualem kr', base='Water table: p = 0',
        geometry='H = 1 m; width = 0.05 m\n2D column; drawing not to scale',
        labels=[('Water boundaries', r'Constant downward inflow $I=1.38889\times10^{-6}$ m/s at the top. Base pressure zero; sides sealed.'),
                ('Retention law', r'Pressure-based $\alpha_p=5\times10^{-4}$ Pa$^{-1}$, $n=2$, residual saturation $S_r=0.1$. Gravity points downward.'),
                ('Mechanical conditions and comparison', 'Side-normal displacement and base vertical displacement fixed; traction-free top. Start hydrostatic and compare the final steady pressure profile.')],
        equations=[
            r's=\max(-p,0),\quad m=1-1/n,\quad S_e=[1+(\alpha_p s)^n]^{-m},\quad S=S_r+(1-S_r)S_e',
            r'k_r=S_e^{1/2}\left[1-(1-S_e^{1/m})^m\right]^2',
            r'J=I/K_s,\quad \gamma=\rho_fg,\quad \frac{dp}{dz}=\gamma\left(\frac{J}{k_r(p)}-1\right),\quad p(0)=0',
            r'z=\frac{1}{\gamma}\int_{p(z)}^0\frac{k_r(\xi)}{k_r(\xi)-J}\,d\xi'
        ],
        definitions=r"Elevation $z$ runs upward from the water-table base. $I$ is positive downward infiltration, $K_s=\kappa\rho_f g/\mu$ saturated hydraulic conductivity, $J=0.5$ and $\gamma=9810\ \mathrm{Pa/m}$. The integral is an implicit analytical relation, evaluated by high-accuracy integration of its equivalent first-order ODE.",
        limit='Homogeneous, steady one-dimensional flow, constant gravity and prescribed flux. This is not a closed-form pressure profile. Storage drops out of the steady equation; the deliberately reduced fluid modulus shortens the run (a larger stable flow step) and does not test transient storage.',
        reference='The implicit relation is derived from steady Darcy flow; evaluated by DOP853 (rtol 10⁻¹², atol 10⁻⁹ Pa) in unsaturated_reference.py.'),
    'brooks_corey': dict(
        title='Brooks–Corey steady infiltration', kind='vertical', top='Downward infiltration I = 0.5 Ks',
        body='Brooks–Corey\nretention\nBurdine kr', base='Water table: p = 0',
        geometry='H = 1 m; width = 0.05 m\n2D column; drawing not to scale',
        labels=[('Water boundaries', r'Constant downward inflow $I=1.38889\times10^{-6}$ m/s at the top. Base pressure zero; sides sealed.'),
                ('Air entry', r'Pressure-based $\alpha_p=10^{-3}$ Pa$^{-1}$, $n=1.5$ and $S_r=0.05$. The material stays saturated down to $p=-1/\alpha_p=-1000$ Pa.'),
                ('Mechanical conditions and comparison', 'Side-normal displacement and base vertical displacement fixed; traction-free top. Start hydrostatic and compare steady pressure, including the air-entry transition.')],
        equations=[
            r's=\max(-p,0),\quad S_e=\min[1,\max(\alpha_p s,1)^{-n}],\quad S=S_r+(1-S_r)S_e',
            r'k_r=S_e^{3+2/n},\qquad p_{\rm entry}=-1/\alpha_p',
            r'J=I/K_s,\quad \gamma=\rho_fg,\quad \frac{dp}{dz}=\gamma\left(\frac{J}{k_r(p)}-1\right),\quad p(0)=0',
            r'z=\frac{1}{\gamma}\int_{p(z)}^0\frac{k_r(\xi)}{k_r(\xi)-J}\,d\xi'
        ],
        definitions=r"Elevation $z$ runs upward from the water-table base. $I$ is positive downward infiltration, $K_s=\kappa\rho_f g/\mu$ saturated hydraulic conductivity and $J=0.5$. In the saturated zone $p=-\gamma(1-J)z$; air entry occurs at $z_{\rm entry}=1/[\alpha_p\gamma(1-J)]\simeq0.204\ \mathrm{m}$.",
        limit='The steady relation is implicit and its unsaturated part is evaluated numerically, rather than treated as a closed-form profile. This test checks the conductivity and air-entry transition. The fluid modulus is deliberately reduced to shorten the run; it does not alter the steady reference and does not qualify transient storage.',
        reference='Derived from steady Darcy flow, with the saturated segment treated as part of the same pressure-dependent law; evaluated in unsaturated_reference.py.'),
    'undrained_wave': dict(
        title='Physical-inertia undrained wave', kind='wave', top='Sealed; vertical displacement fixed',
        body='Saturated\nelastic column', base='Sealed; vertical displacement fixed',
        geometry='H = 1 m; width = 0.025 m\n3D column; drawing not to scale',
        labels=[('Water and displacement boundaries', 'All water boundaries sealed; κ = 0 throughout. End vertical displacement and side-normal displacement fixed.'),
                ('Initial disturbance', r'$v_z(z,0)=V\sin[\pi(z+H)/H]$, $V=10^{-6}$ m/s. Initially zero pressure, stress and strain.'),
                ('Physical regime and comparison', 'Physical mixture density; no gravity, damping, heating, PT or GVS. Compare velocity and pressure over one period.')],
        equations=[
            r'S_0=\phi/K_f+(\alpha-\phi)/K_s,\quad H_u=K+4G/3+\alpha^2/S_0,\quad \rho_b=(1-\phi)\rho_s+\phi\rho_f',
            r'c_u=\sqrt{H_u/\rho_b},\quad k=\pi/H,\quad \omega=kc_u,\quad T=2H/c_u',
            r'v_z=V\sin[k(z+H)]\cos(\omega t),\quad \epsilon_{zz}=(V/c_u)\cos[k(z+H)]\sin(\omega t)',
            r'p=-\frac{\alpha V}{S_0c_u}\cos[k(z+H)]\sin(\omega t),\quad \sigma_{zz}=H_u\epsilon_{zz}'
        ],
        definitions=r"$z$ is the DES coordinate, $-H$ at the base and $0$ at the top. $V$ is the initial velocity amplitude, $\rho_b$ physical mixture density, $H_u$ the undrained constrained modulus and $T$ the wave period. The selected material gives $T\simeq1.162\ \mathrm{ms}$.",
        limit='Linear small-strain, homogeneous, completely sealed, zero-permeability elasticity. The input file sets inertial_scaling, but with is_quasi_static = false DES uses the physical (unscaled) mass, so the value has no effect. This wave is a dynamic coupling check and does not establish earthquake-cycle or GVS equilibrium accuracy.',
        derivation=[
            ('Sealed local fluid content remains zero because it starts at zero. Eliminate pressure from the uniaxial effective-stress relation, giving an undrained modulus in the momentum equation.',
             r'\zeta=\alpha\epsilon_{zz}+S_0p=0,\quad \sigma_{zz}=H_u\epsilon_{zz},\quad \rho_b u_{tt}=H_u u_{zz}'),
            ('The fundamental sine mode vanishes at both fixed ends. Choose its amplitude to match the initial velocity V and differentiate displacement for strain; the zero-content relation then gives pressure.',
             r'u_z=(V/\omega)\sin[k(z+H)]\sin(\omega t)')],
        reference='Direct eigenmode solution of the stated sealed one-dimensional Biot wave problem; evaluated in undrained_wave.py.'),
    'ogata_banks': dict(
        title='Groundwater heat transport: Ogata–Banks', kind='horizontal',
        labels=[('Steady groundwater flow', r'$p(x)=2-x$ Pa, $\kappa/\mu=1$ m²/(Pa s), $q_x=1$ m/s; gravity zero. Pressure is prescribed at both ends.'),
                ('Temperature and mechanical conditions', r'Initially $T_0=273.15$ K. At $t=0$ the left end is set to $T_{in}=274.15$ K. Sides and right end have zero conductive heat flux; heat can advect out. Normal motion is fixed on every face.'),
                ('Reference regime', 'Constant properties and steady Darcy flux; fixed mesh, no mechanical pressure feedback. Deliberately scaled coefficients isolate transport; this is not an aquifer application.')],
        equations=[
            r'C=(1-\phi)\rho_sc_s+\phi\rho_fc_f,\quad v_T=\rho_fc_fq_x/C,\quad D=\lambda/C',
            r'\theta=(T-T_0)/(T_{\rm in}-T_0),\qquad \theta_t+v_T\theta_x=D\theta_{xx}',
            r'\theta(x,t)=\frac{1}{2}\left[\operatorname{erfc}\left(\frac{x-v_Tt}{2\sqrt{Dt}}\right)+e^{v_Tx/D}\operatorname{erfc}\left(\frac{x+v_Tt}{2\sqrt{Dt}}\right)\right]'
        ],
        definitions=r"$C$ is bulk volumetric heat capacity, $c_f$ and $c_s$ fluid and solid specific heat capacities, $\lambda$ thermal conductivity and $q_x$ Darcy flux. Here $\phi=0.3$, $C=4.7\times10^6\ \mathrm{J/(m^3\,K)}$, $v_T=0.851064\ \mathrm{m/s}$ and $D=0.005\ \mathrm{m^2/s}$; the final front position is $v_Tt=0.4\ \mathrm{m}$.",
        limit=r'The analytical solution is semi-infinite: $\theta(0,t)=1$, $\theta(x,0)=0$ for $x>0$ and $\theta(\infty,t)=0$. The numerical domain ends at $x=1$ m with zero conductive flux, so the comparison uses early times before the front reaches that boundary. Constant-capacity, divergence-free flow is tested; general evolving-storage enthalpy conservation is not established. The intentionally tiny fluid modulus makes the hydraulic storage large without changing the steady flux.',
        derivation=[
            ('For constant C and divergence-free steady flux, the bulk sensible-heat balance reduces to the equation above. Dividing by C shows that the thermal front travels at vT, rather than at qx or qx/φ.',
             r'C T_t+\rho_fc_fq_xT_x=\lambda T_{xx}'),
            ('Laplace-transform the normalized equation in time with its zero initial condition. Reject the growing spatial root and impose the unit step at x = 0. Inverting this transform gives the two complementary-error-function terms above.',
             r'\widehat\theta(x,s)=s^{-1}\exp\left[\frac{v_T-\sqrt{v_T^2+4Ds}}{2D}x\right]')],
        reference='Ogata and Banks (1961), USGS Professional Paper 411-A; the transport solution is applied here to normalized temperature. Evaluated stably using erfcx in check_heat_transport.py.',
        url='https://pubs.usgs.gov/publication/pp411A')
}

SPECS['point_sources'] = dict(
    title='Closed-box fluid and heat sources', kind='sources',
    labels=[('Closed domain', 'Unit 3D box; all water and conductive-heat boundaries sealed. Normal motion fixed on every face. Initial pressure 0 Pa, temperature 300 K.'),
            ('Two signed sources', 'Point 1: (0.37, 0.43, −0.53) m, fluid +8×10⁻⁸ m³/s and heat +4000 W, active 0.25–2.75 s. Point 2: (0.61, 0.57, −0.41) m, fluid −3×10⁻⁸ m³/s and heat −1500 W, active 0.8–2.1 s.'),
            ('Isolated balance', 'κ = λ = 0; mechanical pressure feedback and thermal expansion disabled. This case compares integrated volume and heat; it does not claim a spatial point-source solution.')],
    equations=[
        r'\ell_i(t)=\max[0,\min(t,b_i)-a_i],\quad \Delta V_f(t)=\sum_i R_i[\ell_i(t)-\ell_i(0)],\quad \Delta E(t)=\sum_i Q_i[\ell_i(t)-\ell_i(0)]',
        r'\sum_n m_{h,n}[p_n(t)-p_n(0)]=\Delta V_f(t),\quad \sum_n m_{T,n}[T_n(t)-T_n(0)]=\Delta E(t)',
        r'm_{h,n}=\sum_{e\ni n}V_e S_B/N_e,\quad m_{T,n}=\sum_{e\ni n}V_e C/N_e,\quad S_B=S_0+\alpha^2/K_*'
    ],
    definitions=r'$R_i$ and $Q_i$ are signed fluid-volume and heat rates, positive for supply; $a_i$ and $b_i$ are start and end times in seconds. The plotted run has $C=(1-\phi)\rho_sc_s+\phi\rho_fc_f=2.95\times10^6$ J/(m³ K), $K_*=K$ in 3D. In 2D the in-plane mean-stress split uses $K_*=K+G/3$ and rates are per unit out-of-plane thickness. At 4 s the exact net supply is 1.61×10⁻⁷ m³ and 8050 J.',
    limit='This is an exact integral balance for an isolated component test with fixed capacity and mechanical pressure feedback disabled. It does not qualify the spatial singular solution, nonlinear storage, the enthalpy of injected fluid (see the injection-temperature page) or conservative remesh transfer. Fluid volume and heat sources are independent inputs. Separate forced-remesh runs check source relocation and continuation.',
    derivation=[
        ('A point rate is shared among the nodes of the element that contains it by the element shape functions. They sum to one, so the nodal loads add up to the requested signed rate. A source window that starts or ends inside a time step is integrated over its exact overlap with the step, so every source acts exactly once.',
         r'f_n^{(i)}=R_iN_n(\mathbf{x}_i),\quad \sum_nN_n(\mathbf{x}_i)=1,\quad \int_0^t r_i(\tau)\,d\tau=R_i[\ell_i(t)-\ell_i(0)]'),
        ('Summing the lumped nodal update cancels internal transport contributions. With sealed boundaries, fixed capacities and the mechanical source disabled, only the prescribed source integrals remain. Apply the same argument to thermal mass and heat rate. A total-amount schedule uses rate equal to amount divided by (bi−ai).', None)],
    reference='Direct time integral of the supplied schedules plus the assembled closed-domain storage/heat balances; independently evaluated from simplex volumes and saved fields in check_point_sources.py.')

REFERENCE_ONLY = ()
SPECS['theis'] = dict(
    title='Theis confined-aquifer pumping reference', kind='radial',
    labels=[('Geometry and initial state', 'Infinite homogeneous confined aquifer of thickness B. A fully penetrating line well at r = 0 starts extracting total rate Q at t = 0. Initial head drawdown is zero.'),
            ('Boundaries', 'Top and bottom sealed. Drawdown tends to zero as r tends to infinity. The dotted line marks an observation radius r, not an outer boundary. Exclude the singular well axis.'),
            ('Scope', 'Constant transmissivity and storativity; no wellbore storage or explicit mechanical field. Positive Q extracts water; negative Q injects. This appendix supplies an analytical reference, not a completed DES radial-flow qualification.')],
    equations=[
        r'\mathcal{T}_p=\kappa B/\mu,\quad S_p=B S_B,\quad S_B=S_0+\alpha^2/K',
        r'C_Q=Q/(4\pi\mathcal{T}_p),\quad u=r^2S_p/(4\mathcal{T}_p t),\quad s_p(r,t)=Q E_1(u)/(4\pi\mathcal{T}_p)',
        r'E_1(u)=\int_u^\infty e^{-v}/v\,dv,\quad s_p(r,0)=0\quad(r>0)',
        r's_p(r,t)=\sum_{t>t_k}\Delta Q_k\,\mathcal{G}(r,t-t_k)',
        r'\mathcal{G}(r,\tau)=\frac{E_1[r^2S_p/(4\mathcal{T}_p\tau)]}{4\pi\mathcal{T}_p}\quad(\tau>0)'
    ],
    definitions=r'$s_p=p_{\rm initial}-p$ is pressure drawdown in Pa and $\mathcal{G}=E_1(u)/(4\pi\mathcal{T}_p)$ the drawdown of a unit-rate well. $\mathcal{T}_p$ has units $\mathrm{m^3/(Pa\,s)}$, $S_p$ has units $\mathrm{m/Pa}$, and $Q$ is total discharge in $\mathrm{m^3/s}$. Here $S_0=\phi/K_f+(\alpha-\phi)/K_s$. Disabling mechanical pressure feedback retains the selected fixed-stress storage $S_B$. Terms are zero before each rate change $\Delta Q_k$ at time $t_k$. This pressure formulation remains valid at zero gravity; it uses no conversion from pressure to head.',
    limit='Infinite, homogeneous, isotropic, confined aquifer with constant coefficients and a fully penetrating line well; no wellbore storage, skin or finite-radius boundary. The finite DES box is compared before its sealed far sides affect the observation window.',
    derivation=[
        (r'Away from the well, integrate pressure diffusion through thickness $B$ and use zero initial drawdown. The inward discharge at $r\to0$ supplies the line-sink boundary. Integrating the radial diffusion Green function in time gives $E_1$.',
         r'S_p\,\partial_t s_p=\mathcal{T}_p\,r^{-1}\partial_r(r\partial_r s_p),\quad -2\pi r\mathcal{T}_p\,\partial_r s_p\longrightarrow Q'),
        ('Differentiate the solution to check the well limit. Storage inside an observation circle supplies the remaining discharge; linearity permits pumping/shut-in and multiwell superposition.',
         r'-2\pi r\mathcal{T}_p\,\partial_r s_p=Qe^{-u}')],
    reference='Theis (1935), evaluated with the exponential integral in theis_reference.py; the same module superposes changes of rate.',
    url='https://water.usgs.gov/ogw/pubs/Theis-1935.pdf')
SPECS['averaged_aquifer'] = dict(
    title='Generalized Bear–Corapcioglu aquifer reference', kind='aquifer',
    labels=[('Confined annulus', 'Aquifer thickness B, well radius rw and outer radius R. The well is fully penetrating; top and bottom are hydraulically sealed. Initially pressure, strain and displacement increments are zero.'),
            ('Illustrated boundaries', 'Extraction Q at the well face; fixed radial displacement Ur = 0 at rw and R; p = 0 at R. The reference also separately supports a sealed outer face and traction-free radial boundaries.'),
            ('Mechanical approximation', 'Incremental vertical and vertical-shear total stresses vanish; horizontal displacement is uniform through the thickness. This is a thickness-averaged model. A thin-slab DES comparison exists; full spatial/time and thin-layer qualification remain open.')],
    equations=[
        r'S_a=S_0+\frac{\alpha^2}{\lambda+G},\quad c=\frac{2G\alpha}{\lambda+2G},\quad D=\frac{\kappa}{\mu S_a},\quad \lambda=K-2G/3',
        r'S_a p_t+c h_t=\frac{\kappa}{\mu}\left(p_{rr}+\frac{p_r}{r}\right),\quad e_h=\frac{\alpha p}{2(\lambda+G)}+h(t)',
        r'\epsilon_{zz}=\frac{\alpha p-\lambda e_h}{\lambda+2G},\quad \delta=-B\epsilon_{zz},\quad U_r=\frac{C_w+\alpha J/[2(\lambda+G)]+h(r^2-r_w^2)/2}{r}',
        r'J=\int_{r_w}^r\xi p(\xi,t)\,d\xi,\quad \widehat p=A I_0(qr)+B_0 K_0(qr)-(c/S_a)\widehat h,\quad q^2=s\mu S_a/\kappa'
    ],
    definitions='Sa is pressure-based storage (1/Pa), S0 = φ/Kf + (α−φ)/Ks, λ the drained Lamé coefficient, eh areal strain, h(t) the uniform radial-equilibrium integration function and δ positive downward subsidence. Cw follows the inner mechanical boundary; for the illustrated fixed well face Cw = 0. Iν and Kν are modified Bessel functions; hats denote the Laplace transform in time, with variable s. B0 is a transform coefficient, distinct from layer thickness B.',
    limit='This generalizes the averaged confined model to finite well/domain radii and compressible fluid/grains. The finite annulus is solved analytically in Laplace space and inverted numerically by a 32-node Talbot method; it is not a closed-form time history. The model assumes a thin layer and cannot reproduce the full near-well 3D transient or phreatic response. Finite boundaries generally make h nonzero, so a finite-domain pressure is not automatically Theis.',
    derivation=[
        (r'Zero incremental vertical stress gives $\varepsilon_{zz}$. Eliminating it from the horizontal stresses and integrating radial equilibrium gives $e_h=\alpha p/[2(\lambda+G)]+h(t)$. Substitution into $\zeta=\alpha(e_h+\varepsilon_{zz})+S_0p$ gives the storage and $h$ coupling in the equations below.',
         r'\epsilon_v=\alpha p/(\lambda+G)+2Gh/(\lambda+2G),\quad \zeta=S_ap+ch'),
        (r'Two hydraulic conditions and the outer radial mechanical condition determine the Laplace-space constants $A$, $B_0$ and $\hat h$ (defined in the theory note) after the inner mechanical boundary eliminates $C_w$. The illustrated fixed inner face sets $C_w=0$; a traction-free inner face instead sets $C_w=r_w^2(\lambda_*+2G)h/(2G)$. At $R$ choose either $U_r=0$ or $\sigma_{rr}=0$.',
         r'2\pi r_wB(\kappa/\mu)\widehat p_r(r_w)=Q/s,\quad \widehat p(R)=0,\quad \widehat U_r(R)=0'),
        (r'For the optional traction boundaries, $\lambda_*$ is the plane-stress Lamé coefficient. The integral $J$ of the theory note is evaluated analytically from the Bessel basis. This specifies the small linear system before numerical inverse transformation.',
         r'\lambda^*=2G\lambda/(\lambda+2G),\quad \sigma_{rr}=(\lambda^*+2G)h-2GU_r/r'),
        (r'In the infinite-domain line-well limit $h=0$. With $C_Q=Q\mu/(4\pi\kappa B)$ and $u=r^2/(4Dt)$, pressure and displacement reduce to closed forms. This retains fluid/grain compressibility; the published incompressible-grain formulas of Bear and Corapcioglu are a further limit of these.',
         r'p=-C_QE_1(u),\quad e_h=\epsilon_{zz}=\frac{\alpha p}{2(\lambda+G)},\quad U_r=-\frac{\alpha C_QDt}{(\lambda+G)r}[uE_1(u)-e^{-u}+1]')],
    reference='Rederived thickness-averaged Biot extension of Bear and Corapcioglu (1981); averaged_aquifer_reference.py. Independent checks cover radial finite-volume comparison, sealed balance, boundary conditions and limiting closed forms. The full theory note is bundled separately.')

# The DES slab is compared with the infinite line-well limit. The optional
# finite-annulus system remains in the derivation and the full theory appendix.
SPECS['averaged_aquifer']['equations'] = [
    r'\lambda=K-2G/3,\quad S_a=S_0+\alpha^2/(\lambda+G),\quad D=\kappa/(\mu S_a)',
    r'C_Q=Q\mu/(4\pi\kappa B),\quad u=r^2/(4Dt),\quad p(r,t)=-C_QE_1(u)',
    r'e_h=\epsilon_{zz}=\frac{\alpha p}{2(\lambda+G)},\quad\delta=-B\epsilon_{zz}',
    r'U_r=-\frac{\alpha C_QDt}{(\lambda+G)r}[uE_1(u)-e^{-u}+1]']
SPECS['averaged_aquifer']['definitions'] = r'$p$ is the pressure change (negative during extraction), $S_a$ is pressure storage in $\mathrm{Pa^{-1}}$, $S_0=\phi/K_f+(\alpha-\phi)/K_s$, $e_h=\epsilon_{xx}+\epsilon_{yy}$ is areal strain, and $\delta$ is positive downward subsidence. The plotted reference is the infinite-domain line-well limit. The full finite-annulus Laplace solution, which includes an additional spatially uniform strain function $h(t)$, is derived below and in the theory appendix.'

for key in ('van_genuchten', 'brooks_corey'):
    SPECS[key]['derivation'] = [
        (r'At steady state $q_z=-I$ is constant. Darcy flow with $z$ upward includes both pressure and gravity. Rearranging and separating variables gives the integral shown above, starting from $p=0$ at $z=0$.',
         r'-I=-\frac{\kappa k_r(p)}{\mu}(p_z+\gamma),\quad dz=-\frac{k_r(p)}{\gamma[k_r(p)-J]}\,dp'),
        (r'For $0<J<1$ the pressure decreases with elevation toward the suction at which $k_r=J$. The integral diverges at that limiting pressure. Storage and elastic compressibility do not enter this steady relation. The analytical relation is evaluated with an independent high-accuracy ODE solver; no DES output or fitting parameter is used.', None)]


SPECS['convection'] = dict(
    title='Onset of hydrothermal convection',
    kind='convection',
    labels=[('Heating and water boundaries', r'A square porous box is heated from below: $T=T_0+\Delta T$ at the base and $T=T_0$ at the top. All four sides are sealed, $q_n=0$; the sides are insulated, $\partial T/\partial n=0$.'),
            ('Initial state', r'Steady conduction $T_c$ with the hydrostatic pressure of the temperature-dependent water density, plus a 1 mK Gaussian seed centred at $x=0.25$ m, mid-depth.'),
            ('What is compared', r'The amplitude of the temperature mode $\cos(\pi x/W)\sin(\pi z/H)$. Its growth or decay rate $\sigma$ is fitted over the second half of each run and compared with linear stability theory. A seed-free run must stay at rest.')],
    equations=[
        r'\rho_f=\rho_{f0}[1-\alpha_f(T-T_0)],\quad \mathbf{q}=-(\kappa/\mu)(\nabla p+\rho_fg\mathbf{e}_z)',
        r'C\,\partial_tT+\rho_{f0}c_f\,\mathbf{q}\cdot\nabla T=\lambda\nabla^2T,\quad \nabla\cdot\mathbf{q}=0',
        r'Ra=\rho_{f0}^2c_fg\alpha_f\Delta T\kappa H/(\mu\lambda),\quad a=\pi H/W',
        r'\theta\propto\cos(\pi x/W)\sin(\pi z/H)e^{\sigma t}',
        r'\sigma=\frac{\lambda}{CH^2}\left(\frac{Ra\,a^2}{a^2+\pi^2}-a^2-\pi^2\right),\quad Ra_c=\frac{(a^2+\pi^2)^2}{a^2}=4\pi^2'],
    definitions=r'Here $z$ is height above the base, and $T_0$ is the temperature at which the water density is $\rho_{f0}$. $\kappa$ is intrinsic permeability, $\mu$ viscosity, $\alpha_f$ the thermal expansivity of water, $c_f$ its specific heat, $\lambda$ the bulk conductivity and $C$ the bulk heat capacity. The runs use $W=H=1$ m, $\Delta T=10$ K, $C=10^6\ \mathrm{J\,m^{-3}K^{-1}}$ and $\lambda=1\ \mathrm{W\,m^{-1}K^{-1}}$, so the conductive time $CH^2/\lambda$ is $10^6$ s; $\kappa$ sets $Ra$.',
    limit=r'Linear stability of the incompressible Boussinesq limit: small disturbances, constant properties, and a rigid matrix without pore-pressure feedback. DES keeps a small water storage; pressure diffuses about $10^3$ times faster than heat, so the incompressible reference applies. 3D cells, temperature-dependent viscosity and the Elder transient-plume problem are not yet verified.',
    derivation=[
        (r'Conduction with no flow is the base state, $T_c=T_0+\Delta T(1-z/H)$, with the hydrostatic pressure of $\rho_f(T_c)$. Perturb it by $\theta$ and $\mathbf{q}$. Taking the curl of Darcy’s law twice removes the pressure: the vertical Darcy flux $w$ is driven only by horizontal buoyancy gradients. The heat balance couples $\theta$ back to $w$ through the base gradient $-\Delta T/H$.',
         r'\nabla^2w=\frac{\kappa\rho_{f0}g\alpha_f}{\mu}\nabla_h^2\theta\quad C\theta_t-\rho_{f0}c_f\frac{\Delta T}{H}w=\lambda\nabla^2\theta'),
        (r'Impermeable isothermal top and base, and sealed insulated sides, admit $\theta,w\propto\cos(ax/H)\sin(\pi z/H)\,e^{\sigma t}$. With $K^2=a^2+\pi^2$, substitution gives the growth rate above. The square box ($a=\pi$) has the smallest single-cell onset, $Ra_c=4\pi^2\approx39.48$.',
         r'\hat{w}=\frac{\kappa\rho_{f0}g\alpha_f}{\mu}\frac{a^2}{K^2}\hat\theta\quad C\sigma=\frac{\lambda}{H^2}\left(Ra\frac{a^2}{K^2}-K^2\right)')],
    reference='Horton and Rogers (1945); Lapwood (1948); Nield and Bejan, Convection in Porous Media, section 6.2. Evaluated in convection_reference.py; no DES output or fitted parameter enters the reference.')

SPECS['convection']['url'] = 'https://empslocal.ex.ac.uk/people/staff/tej202/docs/teaching/Lapwood48.pdf'

SPECS['mandel'] = dict(
    title='Mandel consolidation under a rigid plate',
    kind='mandel',
    labels=[('Load and plates', r'Rigid, frictionless, impermeable plates squeeze the strip $-a<x<a$ with a total force $2F$ per unit length, applied at $t=0^+$. DES models the top plate with a rigid platen: every top node shares one vertical velocity and slides freely.'),
            ('Drainage and symmetry', r'The faces $x=\pm a$ are drained ($p=0$) and traction-free. One quarter is meshed: symmetry at $x=0$ and at the mid-plane, both sealed.'),
            ('What is compared', r'Pressure at every node and the uniform plate settlement, against the series solution. At the centre the pressure first rises above its undrained value $p_0$: the Mandel–Cryer effect.')],
    equations=[
        r'\tan\alpha_n=\frac{1-\nu}{\nu_u-\nu}\alpha_n,\quad A_n=\frac{\sin\alpha_n}{\alpha_n-\sin\alpha_n\cos\alpha_n}',
        r'p_0=\frac{FB(1+\nu_u)}{3a},\quad E_n(t)=e^{-\alpha_n^2ct/a^2}',
        r'p(x,t)=2p_0\sum_{n=1}^\infty A_n[\cos(\alpha_nx/a)-\cos\alpha_n]E_n(t)',
        r'\frac{u_z(b,t)}{b}=-\frac{F(1-\nu)}{2Ga}+\frac{F(1-\nu_u)}{Ga}\sum_{n=1}^\infty A_n\cos\alpha_n E_n(t)',
        r'c=\frac{\kappa}{\mu[S_0+\alpha^2/(K+4G/3)]},\quad S_0=\phi/K_f+(\alpha-\phi)/K_s',
        r'B=\alpha M/K_u,\quad K_u=K+\alpha^2M,\quad M=1/S_0',
        r'\nu=\frac{3K-2G}{2(3K+G)},\quad\nu_u=\frac{3K_u-2G}{2(3K_u+G)}'],
    definitions=r'Here $a$ and $b$ are the half width and half height, $\nu$ and $\nu_u$ the drained and undrained Poisson ratios, $G$ the shear modulus, $\kappa$ the intrinsic permeability and $\mu$ the viscosity. The undrained pressure is $p_0=FB(1+\nu_u)/(3a)$. The runs use $a=b=1$ m, $F=Pa$ with $P=10^5$ Pa, $K=G=10^7$ Pa, $\phi=0.3$, $\alpha=1$, $K_f=2.17$ GPa, $K_s=37$ GPa and $\kappa=1.03\times10^{-12}\ \mathrm{m^2}$, so $c=0.0239\ \mathrm{m^2/s}$ and $a^2/c=41.8$ s.',
    limit=r'Plane strain, small strain, linear isotropic poroelasticity with constant properties; ideal rigid frictionless plates. DES couples pressure and mechanics sequentially; the early-time error is the one-step lag of that split and falls with the time step. No gravity.',
    derivation=[
        (r'A rigid frictionless plate makes $\varepsilon_{zz}(t)$ uniform; $\sigma_{xx}=0$ at the drained faces and equilibrium make $\sigma_{xx}=0$ everywhere, so the problem is one-dimensional in $x$. The fluid balance then contains the unknown plate strain rate, and the plate force fixes it through the integral of $\sigma_{zz}$ over the width.',
         r'\int_0^a\sigma_{zz}\,dx=-F'),
        (r'Separation of variables with the drained condition $p(a,t)=0$ and this integral constraint gives the roots $\alpha_n$ and the series above; the initial undrained state fixes the coefficients. The same series gives the plate settlement.', None)],
    reference='Mandel (1953); Abousleiman et al. (1996), Mandel’s problem revisited. Evaluated independently in mandel_reference.py with the same material values.',
    url='https://doi.org/10.1680/geot.1996.46.2.187')

SPECS['cryer'] = dict(
    title='Cryer consolidation of a sphere',
    kind='cryer',
    labels=[('Load and drainage', r'A poroelastic sphere of radius $R$ receives a uniform compressive normal stress $P$ at $t=0^+$ on its drained surface ($p=0$). DES applies it as a normal stress on the curved facets.'),
            ('Mesh and symmetry', r'One octant is meshed: the three flat faces are symmetry planes (normal motion fixed, sealed), the curved surface is subdivided on the sphere.'),
            ('What is compared', r'Pressure at every node, and especially at the centre, against the series solution. The centre pressure first rises above its undrained value $p_0$ (Mandel–Cryer effect).')],
    equations=[
        r'S=\frac{1}{M}+\frac{\alpha^2}{K+4G/3},\quad c=\frac{\kappa}{\mu S},\quad \beta=\frac{4G\alpha^2}{3KS(K+4G/3)}',
        r'p_0=\frac{\alpha P}{KS(1+\beta)}',
        r'3\beta(\sin\xi_n-\xi_n\cos\xi_n)=(1+\beta)\,\xi_n^2\sin\xi_n',
        r'\frac{p(r,t)}{p_0}=\sum_n w_n\left(1-\frac{R\sin(\xi_nr/R)}{r\sin\xi_n}\right)e^{-\xi_n^2ct/R^2}',
        r'G(\xi)=3(1-\xi\cot\xi)/\xi^2,\quad w_n=-\frac{2(1+\beta)}{\beta\xi_nG^{\prime}(\xi_n)}'],
    definitions=r'Here $K$, $G$ are drained moduli, $\alpha$ the Biot coefficient, $1/M=\phi/K_f+(\alpha-\phi)/K_s$, $\kappa$ the intrinsic permeability and $\mu$ the viscosity. The runs use $R=1$ m, $P=10^5$ Pa and the material of the Mandel page, so $c=0.0239\ \mathrm{m^2/s}$, $R^2/c=41.8$ s and $\beta=1.33$.',
    limit=r'Small strain, linear isotropic poroelasticity with constant properties; the curved surface is a faceted approximation of the sphere. DES couples pressure and mechanics sequentially, so the earliest frames carry the one-step lag of that split. No gravity.',
    derivation=[
        (r'Spherical symmetry makes the displacement irrotational, so $(K+4G/3)\,\varepsilon=\alpha p+f(t)$. The surface traction fixes $f$ through the volume-average pressure $\langle p\rangle$, and the fluid balance becomes a diffusion equation forced by $d\langle p\rangle/dt$.',
         r'p_t-c\nabla^2p=-\beta\,\frac{d\langle p\rangle}{dt}'),
        (r'With $q=\sqrt{s/c}$, $x=qR$ and $g(x)=3(x\coth x-1)/x^2$, the Laplace transform is closed. Its poles $s_n=-c\xi_n^2/R^2$ give the roots and their residues the weights $w_n$. The series is checked against a numerical (Stehfest) inversion of the transform.',
         r'\frac{\hat p}{p_0}=\frac{(1+\beta)\left(1-\frac{R\sinh qr}{r\sinh x}\right)}{s\,(1+\beta(1-g(x)))}')],
    reference='Cryer (1963); Detournay and Cheng (1993). Derived and evaluated in cryer_reference.py with the same material values; the series agrees with a Stehfest inversion to 1.3e-4.',
    url='https://academic.oup.com/qjmam/article-abstract/16/4/401/1871972')

SPECS['theis'].update(
    labels=[('Geometry and well', r'A quarter box, $0\le x,y\le 60$ m, of thickness $B=5$ m. A fully penetrating well on the corner edge pumps $Q/4$, a quarter of the total rate $Q=10^{-4}\ \mathrm{m^3/s}$, for 80 s and then shuts in for 80 s.'),
            ('Boundaries and mechanics', r'All faces are sealed; the symmetry planes $x=0$, $y=0$ make the corner well radial. The far sides stay beyond the drawdown: its reference value there is below $10^{-3}$ of the pressure scale. The skeleton is stiff and mechanical pressure feedback is off, so the storage is the fixed-stress value $S_B$.'),
            ('What is compared', r'The pressure drawdown $s_p=p_{\rm initial}-p$ at every node with $10\le r\le 30$ m, at eight times through pumping and recovery, against the superposed Theis solution. Errors are fractions of $C_Q=Q/(4\pi\mathcal{T}_p)$.')],
    limit=r'Confined, homogeneous, isotropic aquifer with constant properties; a line well without wellbore storage; the singular axis $r<10$ m is excluded. Pressure only: the stiff skeleton removes poroelastic coupling, which the averaged-aquifer page addresses.')
SPECS['averaged_aquifer'].update(
    labels=[('DES slab', r'A thin quarter slab, $400\times400$ m and $B=10$ m thick, with a fully penetrating well on its corner edge pumping $Q/4$ ($Q=0.01\ \mathrm{m^3/s}$). Symmetry rollers at $x=0$, $y=0$, a roller base and a traction-free top; the far sides hold the initial pressure and zero normal displacement.'),
            ('Reference', r'The thickness-averaged (generalised Bear–Corapcioglu) model: vanishing incremental vertical and vertical-shear total stress, horizontal displacement uniform through the thickness. It applies while the slab is thin compared with the radius, up to $O((B/r)^2)$.'),
            ('What is compared', r'Nodal pressure and element areal and vertical strains for $15\le r\le120$ m at four times, while the pressure front stays inside the slab, on structured meshes.')],
    limit=r'Linear poroelasticity with constant properties; the thin-slab geometry approaches the averaged model only up to $O((B/r)^2)$, so a small model difference remains at the finest mesh. Two-way coupling with PT; no gravity.')

SPECS['rudnicki'] = dict(
    title='Fluid injection into a poroelastic solid',
    kind='rudnicki',
    labels=[('Source', r'Water is injected at a constant volume rate $Q$ from $t=0$ at a point in an unbounded poroelastic solid. DES uses one octant with the scheduled point source at its corner node, injecting $Q/8$.'),
            ('Boundaries', r'The planes $x=0$, $y=0$ and $z=0$ are symmetry planes (normal displacement zero, sealed). The far faces of the 60 m box hold $p=0$ and zero normal displacement; the mesh is refined within 10 m of the source.'),
            ('What is compared', r'Pressure at the nodes and the radial and hoop strains at the elements for $2.5\le r\le8$ m, at all saved times while the diffusion front stays well inside the box.')],
    equations=[
        r'S=\frac{1}{M}+\frac{\alpha^2}{K+4G/3},\quad c=\frac{\kappa}{\mu S},\quad a=2\sqrt{ct}',
        r'p(r,t)=\frac{Q\mu}{4\pi\kappa r}\,\operatorname{erfc}\frac{r}{a}',
        r'u_r=\frac{A}{r^2}\left[\frac{r^2}{2}\operatorname{erfc}\frac{r}{a}+\frac{a^2}{4}\operatorname{erf}\frac{r}{a}-\frac{ar}{2\sqrt{\pi}}e^{-r^2/a^2}\right],\quad A=\frac{\alpha Q\mu}{4\pi\kappa(K+4G/3)}',
        r'\varepsilon_{\theta\theta}=\frac{u_r}{r},\quad \varepsilon_{rr}=\frac{\alpha p}{K+4G/3}-\frac{2u_r}{r}',
        r'\sigma_{rr}=-4G\,\frac{u_r}{r},\quad \sigma_{\theta\theta}=2G\left(\frac{u_r}{r}-\frac{\alpha p}{K+4G/3}\right)'],
    definitions=r'Here $r$ is the distance from the source, $K$ and $G$ are drained moduli, $\alpha$ the Biot coefficient, $1/M=\phi/K_f+(\alpha-\phi)/K_s$, $\kappa$ the intrinsic permeability and $\mu$ the viscosity; stresses are total changes, tension positive. The runs use $Q=10^{-4}\ \mathrm{m^3/s}$, $K=2\times10^8$ Pa, $G=10^8$ Pa, $\alpha=0.8$, $\phi=0.25$, $K_f=2.2$ GPa, $K_s=1$ GPa and $\kappa=10^{-13}\ \mathrm{m^2}$, so $c=0.0387\ \mathrm{m^2/s}$.',
    limit=r'Linear isotropic poroelasticity with constant properties and small strain; an unbounded solid approximated by a 60 m box, whose fixed far faces perturb the strain by about $(r/L)^3$, with $L=60$ m (under 0.3% at $r=8$ m). The singular region $r<2.5$ m is excluded. No fault or friction is modelled: the stresses are those a fault would feel.',
    derivation=[
        (r'For a point source in an unbounded solid the displacement is irrotational and vanishes at infinity, so the volumetric strain follows the pressure, $(K+4G/3)\,\varepsilon=\alpha p$. The fluid content $\alpha\varepsilon+p/M$ is then $Sp$, and the pressure obeys a diffusion equation with the coupled coefficient $c$.',
         r'S\,p_t=\frac{\kappa}{\mu}\nabla^2p+Q\,\delta(\mathbf{x})'),
        (r'Its continuous point-source solution is the erfc profile above. The radial displacement is $u_r=r^{-2}\int_0^r\varepsilon s^2ds$, which integrates in closed form; the strains and total stresses follow from $u_r$ and $p$.',
         r'\sigma_{ij}=2G\,(u_{i,j}-\delta_{ij}\,\varepsilon),\ u_i=\phi_{,i}')],
    reference='Rudnicki (1986), Fluid mass sources and point forces in linear elastic diffusive solids, Mechanics of Materials 5, 383-393. Evaluated in rudnicki_reference.py, whose closed-form displacement is checked against direct integration.')

SPECS['bishop'] = dict(
    title='Bishop effective stress: lowering the water table',
    kind='bishop',
    labels=[('Column and loading', r'A laterally confined elastic column, $0\le z\le H$ with $H=10$ m, starts drained and at rest with a hydrostatic water table at $z_0=9$ m. The base pressure is then lowered to that of a table at $z_1=5$ m; the sealed top and sides let the column drain to the new hydrostatic state.'),
            ('Mechanics', r'Side-normal and base displacements are fixed; the top is traction-free. Gravity acts on the partly saturated bulk weight. The PT solver keeps mechanical equilibrium, with pressure feedback on.'),
            ('What is compared', r'The change in vertical displacement between the two drained states, integrated from the element strains, against the exact solution. Three retention laws are run, and a column with van Genuchten below $z=7$ m and Gardner above, so interface nodes are shared by two laws.')],
    equations=[
        r"\sigma'=\sigma+\alpha\,S_e(p)\,p,\qquad \chi=S_e",
        r'\rho=\rho_s(1-\phi)+\rho_f\,\phi\,S,\qquad S=S_r+(1-S_r)\,S_e',
        r'p_k(z)=-\rho_fg\,(z-z_k),\qquad \Delta\pi=S_e(p_1)p_1-S_e(p_0)p_0',
        r'\Delta\sigma_{zz}(z)=-\rho_fg\,\phi\,(1-S_r)\int_z^H\Delta S_e\,dz^\prime',
        r'u_z(z)=\frac{1}{K+4G/3}\left[\alpha\int_0^z\Delta\pi\,dz^\prime-\rho_fg\,\phi\,(1-S_r)\int_0^H\min(z^\prime,z)\,\Delta S_e\,dz^\prime\right]'],
    definitions=r'Pressure $p$ is compression positive and suction is $-p$; stress is tension positive. $S_e$ is the effective saturation of the van Genuchten, Brooks–Corey or Gardner law and $S$ the saturation. The runs use $\phi=0.4$, $\alpha=0.9$, $K=2\times10^7$ Pa, $G=10^7$ Pa, $\rho_f=1000\ \mathrm{kg/m^3}$, and retention $\alpha_p=1/(2\times10^4)$ Pa$^{-1}$, $n=1.6$, $S_r=0.1$ (van Genuchten); $1/10^4$ Pa$^{-1}$, $\lambda=0.8$, $S_r=0.05$ (Brooks–Corey); $0.5/(\rho_fg)$ Pa$^{-1}$ (Gardner).',
    limit=r'Both states are drained and at rest, so storage and permeability do not enter; the pore fluid is made compressible only to allow a larger stable flow step, which shortens the drainage run. This checks the Bishop stress and the partly saturated weight, not the transient; the transient coupling is checked separately by a sealed cell against its exact ordinary differential equation. Small strain, linear elasticity, no hysteresis.',
    derivation=[
        (r'Laterally confined uniaxial strain gives $\varepsilon_{zz}=(\Delta\sigma_{zz}+\alpha\Delta\pi)/(K+4G/3)$. Vertical equilibrium with the changed bulk weight and a traction-free top gives the total stress change as the weight lost above $z$.',
         r'\frac{d\Delta\sigma_{zz}}{dz}=\rho_fg\,\phi\,\Delta S,\qquad \Delta\sigma_{zz}(H)=0'),
        (r'Integrating the strain from the fixed base and exchanging the order of integration in the weight term gives $u_z(z)$ above. The saturated limit is the classical $u_z=\alpha\Delta p\,z/(K+4G/3)$.', None)],
    reference='Bishop (1959); Lu, Godt and Wu (2010), $\chi=S_e$. Evaluated in unsaturated_reference.py, checked against closed forms, an independent double integral and reversibility; the two-law column is integrated numerically in check_bishop_column.py.',
    url='https://doi.org/10.1029/2009WR008646')

SPECS['hantush'] = dict(
    title='Hantush–Jacob pumping from a leaky aquifer',
    kind='hantush',
    labels=[('Geometry and well', r'A quarter box, $0\le x,y\le75$ m. An aquifer $B=2.5$ m thick lies under an aquitard one element thick ($b^\prime=h$) whose top is held at $p=0$. A well on the corner edge, screened through the aquifer only, pumps $Q/4$ ($Q=10^{-4}\ \mathrm{m^3/s}$) for 200 s and then shuts in.'),
            ('Layers', r'The aquitard is a second material whose vertical leakance $k_z^\prime/(\mu b^\prime)$ is the same on every mesh, so $\lambda=15$ m throughout; it stores a thousandth of the aquifer water per volume. Vertical permeability is ten times horizontal in both layers, so the aquifer has negligible vertical resistance, as the theory assumes.'),
            ('What is compared', r'The pressure drawdown at the aquifer nodes with $10\le r\le30$ m, at eight times through pumping and recovery, against the superposed Hantush–Jacob solution. The skeleton is stiff and pressure feedback is off.')],
    equations=[
        r'\mathcal{T}_p=\kappa B/\mu,\quad S_p=BS_B,\quad \lambda^2=\frac{\kappa B\,b^\prime}{\kappa_z^\prime}',
        r'S_p\,\partial_ts_p=\mathcal{T}_p\left(\partial_r^2s_p+\frac{1}{r}\partial_rs_p\right)-\frac{\mathcal{T}_p}{\lambda^2}\,s_p',
        r's_p(r,t)=\frac{Q}{4\pi\mathcal{T}_p}\,W\!\left(\frac{r^2S_p}{4\mathcal{T}_pt},\frac{r}{\lambda}\right)',
        r'W(u,\beta)=\int_u^\infty\frac{1}{y}\exp\!\left(-y-\frac{\beta^2}{4y}\right)dy',
        r'W(u,0)=E_1(u),\qquad W(0,\beta)=2K_0(\beta)'],
    definitions=r'$s_p=p_{\rm initial}-p$ is the pressure drawdown, $\kappa$ and $\kappa_z^\prime$ the horizontal aquifer and vertical aquitard permeabilities, and $S_B$ the fixed-stress storage of the Theis page. Shut-in is superposed as in Theis. Errors are fractions of $C_Q=Q/(4\pi\mathcal{T}_p)$.',
    limit=r'Hantush and Jacob neglect aquitard storage, horizontal aquitard flow and vertical resistance in the aquifer. Their dimensionless measures here are $S^\prime_B b^\prime/S_p\le0.12\%$, $(b^\prime/\lambda)^2/10\le0.28\%$ and $B^2/(20\lambda^2)=0.14\%$, the first two falling with $h$. $S^\prime_B$ is aquitard storage per volume and $S_p=BS_B$ aquifer storage per area. Line well, no wellbore storage; the axis $r<10$ m is excluded.',
    derivation=[
        (r'With no aquitard storage the aquitard pressure is linear in $z$, so it leaks a flux $(\kappa_z^\prime/\mu b^\prime)\,s_p$ into the aquifer. Integrating over the aquifer thickness adds this sink to the Theis equation.',
         r'q_{\rm leak}=\frac{\kappa_z^\prime}{\mu\,b^\prime}\,s_p=\frac{\mathcal{T}_p}{\lambda^2}\,s_p'),
        (r'The Laplace transform in time and the line-sink limit at $r\to0$ give the leaky well function $W$, which tends to Theis at early time and to the steady cone $2K_0(r/\lambda)$ at late time. hantush_reference.py checks both limits, the identity $W(u,\beta)+W(\beta^2/4u,\beta)=2K_0(\beta)$ and the equation itself.', None)],
    reference='Hantush and Jacob (1955), Non-steady radial flow in an infinite leaky aquifer, Trans. AGU 36, 95-100. Evaluated in hantush_reference.py by quadrature.',
    url='https://doi.org/10.1029/TR036i001p00095')

SPECS['mctigue'] = dict(
    title='Thermal pressurization: McTigue surface heating',
    kind='mctigue',
    labels=[('Column and heating', r'A laterally confined column $0\le z\le1$ m (depth $z$) of fluid-saturated rock starts at $T_0$. At $t=0^+$ its top is raised to $T_0+\Delta T$, $\Delta T=10$ K, and held there; the top is drained ($p=0$) and traction-free.'),
            ('Other boundaries', r'Sides and base are sealed, insulated and fixed normally. Heat conducts without advection and gravity is off. Both diffusion lengths stay well inside the column, so it behaves as a half-space.'),
            ('What is compared', r'Pressure and temperature at the nodes, vertical strain in each element layer and the surface heave, at eight times to 40 ks. DES applies a boundary temperature present at initialization as part of the initial state, so the run restarts from a checkpoint at $T_0$ with the hot surface.')],
    equations=[
        r'\varepsilon=\frac{\alpha p+K\alpha_sT}{H},\qquad H=K+\frac{4}{3}G',
        r'\zeta=\alpha\varepsilon+S_0p-\beta T,\qquad \beta=\phi\alpha_f+(\alpha-\phi)\alpha_s',
        r'S\,p_t-\frac{k}{\mu}p_{zz}=S\gamma\,T_t,\quad S=S_0+\frac{\alpha^2}{H},\quad \gamma=\frac{\beta-\alpha K\alpha_s/H}{S},\quad c=\frac{k}{\mu S}',
        r'T=\Delta T\,\operatorname{erfc}\frac{z}{2\sqrt{\kappa t}}',
        r'p=\frac{\gamma\Delta T}{1-c/\kappa}\left[\operatorname{erfc}\frac{z}{2\sqrt{\kappa t}}-\operatorname{erfc}\frac{z}{2\sqrt{ct}}\right]',
        r'w(t)=\frac{2}{\sqrt{\pi}}\left[\frac{K\alpha_s\Delta T}{H}\sqrt{\kappa t}+\frac{\alpha\gamma\Delta T}{H(1-c/\kappa)}\left(\sqrt{\kappa t}-\sqrt{ct}\right)\right]'],
    definitions=r'Stress is tension positive and $p$ compression positive. $\alpha_s$ and $\alpha_f$ are volumetric thermal expansions, $\kappa$ the thermal diffusivity, $S_0=\phi/K_f+(\alpha-\phi)/K_s$ and $w$ the upward surface heave. The runs use $K=G=1$ GPa, $\alpha=0.8$, $\phi=0.2$, $K_f=2.2$ GPa, $K_s=5$ GPa, $k=10^{-19}\ \mathrm{m^2}$, $\mu=10^{-3}$ Pa s, $\alpha_s=3\times10^{-5}$ K$^{-1}$, $\alpha_f=2.07\times10^{-4}$ K$^{-1}$ and $\kappa=10^{-6}\ \mathrm{m^2/s}$, so $\gamma=101$ kPa/K and $c/\kappa=0.21$.',
    limit=r'Linear thermo-poroelasticity with constant properties, small strain and conduction only (no Darcy heat advection, no temperature-dependent fluid properties). A finite column stands for the half-space while both diffusion lengths stay well inside it. DES couples pressure and mechanics sequentially, so the earliest frames carry the one-step lag.',
    derivation=[
        (r'Uniaxial strain with a traction-free surface keeps the vertical total stress zero, which fixes the strain from $p$ and $T$. Substituting it into the fluid content and Darcy flow gives a pressure diffusion equation forced by the heating rate.',
         r'\sigma_{zz}=H\varepsilon-\alpha p-K\alpha_sT=0'),
        (r'Each erfc term solves its own diffusion equation, so a difference of the thermal and hydraulic erfc solves the forced equation and vanishes at the drained surface and at $t=0$. Integrating the strain over depth gives the heave. mctigue_reference.py checks the equation, the conditions, the heave integral and an independent Crank–Nicolson solution.', None)],
    reference='McTigue (1986), Thermoelastic response of fluid-saturated porous rock, J. Geophys. Res. 91(B9), 9533-9542. Evaluated in mctigue_reference.py.',
    url='https://doi.org/10.1029/JB091iB09p09533')

SPECS['bishop_cell'] = dict(
    title='Bishop coupling in a sealed cell under compression',
    kind='bishop_cell',
    labels=[('Cell and loading', r'A homogeneous block, 1 m high, of partly saturated Gardner soil at uniform suction is compressed uniaxially at a constant rate to a vertical strain of $-10^{-3}$ over 1000 s. All faces are sealed and the permeability is negligible, so no water moves.'),
            ('Two states', r'Moderate suction, $p_0=-5$ kPa with $\alpha_p=10^{-4}$ Pa$^{-1}$ ($\alpha_ps=0.5$), where $\chi\Pi^\prime>0$; and strong suction, $p_0=-2$ kPa with $\alpha_p=10^{-3}$ Pa$^{-1}$ ($\alpha_ps=2$, strong relative to the retention scale), where the Bishop pressure $\Pi=S_ep$ decreases with suction ($\Pi^\prime<0$). The skeleton is soft, $K=G=10$ kPa, so the stabilising storage is large.'),
            ('What is compared', r'The mean pore pressure along each run’s own strain path, against the exact solution of the rate equation; PT at 10 to 320 steps and the explicit strain form, in 2D and 3D, and the same scheme without the correction.')],
    equations=[
        r'\left(S_0+\phi(1-S_r)\frac{dS_e}{dp}\right)\dot p=-\alpha\,\chi\,\dot\varepsilon_v,\qquad \chi=S_e',
        r'\Pi=S_ep,\qquad \Pi^\prime=S_e+p\,\frac{dS_e}{dp},\qquad L=\alpha^2/K_b',
        r'D=\overline{p-p_{\rm old}},\qquad E=\overline{\Pi(p)-\Pi(p_{\rm old})}',
        r'q_{\rm FS}=\frac{\alpha\chi\,\Delta\sigma_m}{K_b\,\Delta t}+\frac{L(\chi E-D)}{\Delta t}',
        r'S_0\Delta p+\phi\Delta S=-\alpha\chi\Delta\varepsilon_v+L(D-\Delta p)'],
    definitions=r'$S_0=\phi/K_f+(\alpha-\phi)/K_s$, $K_b$ is the mean-stress modulus of the fixed-stress split, and the overbars are element means of the increments of the last physical step, already applied by the mechanical solve. $q_{\rm FS}$ is the fixed-stress source per unit volume; the mechanical stress change $\Delta\sigma_m$ already contains $-\alpha E$.',
    limit=r'Sealed and homogeneous, so flow and spatial accuracy are not tested; the column pages cover them. The remaining PT error is the one-step lag $L(D-\Delta p)$, first order in the step and weighted by the soft skeleton’s large $L$. Small strain, linear elasticity, no hysteresis.',
    derivation=[
        (r'Without flow, fluid content changes only by storage, desaturation and the Bishop-weighted skeleton strain, which gives the rate equation above. The fixed-stress source reads the stress change of the last mechanical step, $K_b\Delta\varepsilon-\alpha E$, while the scheme stores $S_0+L$.',
         r'(S_0+L)\Delta p+\phi\Delta S=-\alpha\chi\Delta\varepsilon_v+\chi LE\ \mathrm{(uncorrected)}'),
        (r'Adding $L(\chi E-D)$ cancels the Bishop increment inside the stress change and leaves $L(D-\Delta p)$, which vanishes as consecutive increments approach each other. Clipping $\chi\Pi^\prime$ at zero is not used: where $\Pi^\prime<0$ it would converge to a different continuum.', None)],
    reference='Bishop (1959); Lu, Godt and Wu (2010). The rate equation is integrated with a high-order ODE solver in check_bishop_cell.py; the correction follows the fixed-stress split of Kim, Tchelepi and Juanes (2011) with the Bishop increment.',
    url='https://doi.org/10.1016/j.cma.2010.12.022')

SPECS['nusselt'] = dict(
    numerical_reference=True,
    title='Steady hydrothermal convection: Nusselt number',
    kind='nusselt',
    labels=[('Cell', r'The Horton–Rogers–Lapwood square cell of the onset page, $0\le x,z\le H$: impermeable walls, base at $T_0+\Delta T$, top at $T_0$, insulated sides, Boussinesq Darcy flow with heat conduction and advection.'),
            ('Run', r'Each run starts from conduction with a seeded roll at $\mathrm{Ra}=100$ or $200$ and continues until the heat transport stops changing. Storage does not enter the steady state, so a compressible pore fluid shortens the transient.'),
            ('What is compared', r'The steady Nusselt number, from the volume form below and from the temperature gradient at the top, against an independent steady solution on three fine grids, extrapolated.')],
    equations=[
        r'\nabla^2\psi=-\mathrm{Ra}\,\theta_x,\qquad u\,\theta_x+w\,\theta_z=\nabla^2\theta,\qquad u=\psi_z,\ w=-\psi_x',
        r'\mathrm{Ra}=\frac{\rho_f^2c_f\,g\,\alpha_f\,\Delta T\,k\,H}{\mu\,\lambda}',
        r'\mathrm{Nu}=-\frac{1}{H}\int_0^H\left.\frac{\partial\theta}{\partial z}\right|_{z=H}dx=1+\frac{\rho_fc_fH}{\lambda\,\Delta T}\langle q_zT\rangle'],
    definitions=r'$\theta=(T-T_0)/\Delta T$, lengths are scaled by $H$ and $\psi$ is the dimensionless streamfunction; $q_z$ is the vertical Darcy flux and $\langle\cdot\rangle$ the volume mean. The second form of Nu follows from the steady heat balance at every height, averaged over the cell; it is the one evaluated on DES fields, with $q$ rebuilt from the saved pressure and temperature.',
    limit=r'Steady state only: the transient depends on storage, which is chosen for speed. Single-roll solution, as seeded; Boussinesq, constant properties, no mechanical coupling.',
    derivation=[
        (r'Taking the curl of Darcy’s law with Boussinesq buoyancy removes the pressure and gives the streamfunction equation; the energy balance is steady advection–diffusion.', None),
        (r'convection_steady.py solves both with second-order central differences on $32^2$, $64^2$ and $128^2$ grids by relaxed Picard iteration and extrapolates Nu (observed order 1.99). Below onset it returns the conductive $\mathrm{Nu}=1$; at Ra = 100 it gives 2.6459, matching published values near 2.65.', None)],
    reference='Horton and Rogers (1945); Lapwood (1948); steady values as in Caltagirone (1975). Computed independently in convection_steady.py.',
    url='https://doi.org/10.1017/S0022112075000936')

SPECS['regional'] = dict(
    title='Application: seasonal pumping and land subsidence',
    kind='regional',
    labels=[('Basin and pumping', r'A confined aquifer of soft basin sediment, $B=50$ m thick, is modelled as a quarter of a 6 km square slab with a well field on its corner edge. It pumps $Q=0.02\ \mathrm{m^3/s}$ for six months, rests for six months, and repeats once.'),
            ('Boundaries', r'Symmetry rollers at $x=0$ and $y=0$, a roller base and a traction-free top. The outer sides hold $p=0$ and zero normal displacement, so the basin is closed. Gravity is off: all fields are changes from the initial state.'),
            ('What is compared', r'Pressure at the nodes and subsidence $B\,|\varepsilon_{zz}|$ of the elements for $250\le r\le1500$ m at 17 times, and the subsidence at three observation wells, against the averaged-aquifer model of the closed basin.')],
    equations=[
        r'S_a\,p_t-\frac{\kappa}{\mu}\nabla^2p=-\alpha\,\dot h(t)\cdot\frac{2G}{\lambda+2G},\qquad S_a=S_0+\frac{\alpha^2}{\lambda+G}',
        r'e_h=\frac{\alpha\,p}{2(\lambda+G)}+h(t),\qquad \varepsilon_{zz}=\frac{\alpha p-\lambda e_h}{\lambda+2G},\qquad \delta=-B\,\varepsilon_{zz}',
        r'p(R,t)=0,\qquad u_r(R,t)=0\ \Rightarrow\ h(t)\ne0',
        r'Q(t)=Q\ \mathrm{for}\ 0\leq t<T\ \mathrm{and}\ 2T\leq t<3T,\ \mathrm{else}\ 0,\qquad T=1/2\ \mathrm{yr}'],
    definitions=r'$e_h$ is the areal strain, $h(t)$ the uniform areal strain of the closed basin and $\delta$ the subsidence (positive down). The square is represented by the circle of the same area, $R=2L/\sqrt\pi=3.39$ km (half side $L=3$ km). The runs use $K=5$ MPa, $G=3$ MPa, $\alpha=0.9$, $\phi=0.3$, $K_f=2.2$ GPa, $K_s=37$ GPa and $\kappa=10^{-12}\ \mathrm{m^2}$.',
    limit=r'Linear elastic, thin-slab averaged model: horizontal displacement uniform through the thickness, up to $O((B/r)^2)$. No gravity, no overlying aquitard and no inelastic compaction of clays, which dominate most real subsidence. The square basin is approximated by a circle of equal area in the reference.',
    derivation=[
        (r'Averaging Biot poroelasticity over the thickness with zero vertical and shear stress gives the generalised Bear–Corapcioglu model; the areal strain follows the local pressure up to a uniform function $h(t)$ fixed by the outer boundary.', None),
        (r'For a closed basin the fixed outer displacement makes $h(t)$ nonzero, and through the coupling it deepens the drawdown everywhere; averaged_aquifer_reference.py solves this finite annulus exactly in Laplace space (Talbot inversion), and the response to the schedule is superposed over the rate changes.', None)],
    reference='Bear and Corapcioglu (1981a, b), generalised to Biot coefficients and a finite basin in averaged_aquifer_reference.py.',
    url='https://doi.org/10.1029/WR017i004p00947')

SPECS['theis_wells'] = dict(
    title='Several wells: pumping, injection and recovery', kind='theis_wells',
    labels=[('Geometry and wells', 'A 160 by 160 by 5 m confined slab contains full-height wells at (60,80), (100,80) and (80,100) m. Their rates are Q, 0.6Q and -0.25Q, with $Q=10^{-4}$ m³/s.'),
            ('Water boundaries and initial state', 'All faces are sealed, initial pressure is zero, and every well stops at 40 s. The comparison ends at 80 s, before the pressure front reaches the far walls.'),
            ('Mechanical and material assumptions', 'Constant permeability and pressure capacity; mechanical feedback, gravity and heat are disabled. Observations remain 10–30 m from the nearest well axis.')],
    equations=[r'\mathcal{T}_p=\kappa B/\mu,\qquad S_p=BS_B,\qquad D=\mathcal{T}_p/S_p',
               r'W(r,t)=E_1[r^2/(4Dt)]\ (t>0),\qquad W(r,t)=0\ (t\leq0)',
               r's_p(\mathbf{x},t)=\sum_i\frac{Q_i}{4\pi\mathcal{T}_p}[W(r_i,t)-W(r_i,t-40)]',
               r'S_B\int_V[p(t)-p(0)]\,dV=-\sum_i Q_i\min(t,40)'],
    definitions=r'$s_p=p(0)-p(t)$ is drawdown in Pa; extraction has positive $Q_i$. $E_1(u)=\int_u^\infty e^{-s}ds/s$. The fixed comparison scale is $C_\Sigma=\sum_i|Q_i|/(4\pi\mathcal{T}_p)$.',
    limit='Infinite-aquifer reference, compared inside a sealed finite slab before its walls matter. Linear storage and constant properties; this does not establish superposition for nonlinear or evolving-property flow.',
    derivation=[('The pressure diffusion equation is linear in both pressure and well rate. A well switched off at 40 s is the sum of a positive step at zero and its negative at 40 s.', None),
                ('Summing those responses over the three well locations gives the pressure reference. Integrating the pressure equation over sealed boundaries gives the signed volume budget independently of the pointwise solution.', None)],
    reference='Theis superposition, evaluated in theis_reference.py and check_theis_wells.py; the saved combined run is also compared with three independent single-well DES runs.',
    url='https://doi.org/10.1029/TR016i002p00519')

SPECS['open_flow'] = dict(
    title='Open boundaries: a draining pressure pulse', kind='open_flow',
    labels=[('Geometry', 'A 10 m column has a 1 m² cross-section. Both end faces hold zero pressure; the lateral faces are sealed.'),
            ('Initial condition', 'Pressure is $P\sin(\pi x/L)$, with $P=1000$ Pa and $L=10$ m. The run restarts from the shipped seed checkpoint, whose frame 0 holds this field.'),
            ('Material and mechanics', 'Diffusivity is 1 m²/s and pressure capacity is constant. Mechanical pressure feedback, gravity and heat are disabled.')],
    equations=[r'p(x,t)=P\sin(\pi x/L)e^{-\lambda t},\qquad \lambda=D\pi^2/L^2',
               r'M(t)=S_B A\int_0^L p\,dx=M_0e^{-\lambda t}',
               r'M_0=2S_BPAL/\pi,\qquad V_{\rm out}(t)=M_0(1-e^{-\lambda t})',
               r'M(t)+V_{\rm out}(t)=M_0'],
    definitions=r'$M$ and $V_{\rm out}$ are fluid volumes, in m³. $V_{\rm out}$ is positive for outward drainage. Stored fluid is computed from the element volumes and nodal pressures, independently of DES’s own fluid balance.',
    limit='Constant-capacity saturated flow on a fixed mesh. Variable storage, unsaturated coupling, deforming volume and remapping need separate open-system fluid budgets.',
    derivation=[('The sine is the first eigenfunction of the diffusion equation with zero pressure at each end. Separation of variables supplies its exponential decay.', r'p_t=Dp_{xx},\qquad p(0,t)=p(L,t)=0'),
                ('Integrating the sine over the column gives the stored volume. Darcy outflow at the two ends supplies its equal and opposite loss.', r'dM/dt=-\lambda M,\qquad dV_{\rm out}/dt=\lambda M')],
    reference='Separation-of-variables solution of the stated diffusion problem; pressure, storage and the boundary budget of every time step are checked in check_open_flow.py.')

SPECS['kirsch'] = dict(
    title='Elastic circular opening: Kirsch limit', kind='kirsch',
    labels=[('Geometry and load', 'A quarter annulus has opening radius a = 1 m and outer radius R = 20a. The opening is free; the outer curved face carries isotropic compressive traction P = 100 kPa.'),
            ('Constraints', 'The two straight symmetry cuts have zero normal displacement. This is a dry plane-strain problem with gravity and transport disabled.'),
            ('Material and comparison', 'K = G = 1 GPa. Stress and strain are compared over $1.1a\le r\le5a$ with the finite-annulus reference; the analytical curve remains visible and DES markers are thinned for display.')],
    equations=[r'A_0=-PR^2/(R^2-a^2),\qquad B_0=-Pa^2R^2/(R^2-a^2)',
               r'\sigma_{rr}=A_0-B_0/r^2,\qquad \sigma_{\theta\theta}=A_0+B_0/r^2',
               r'u_r=\frac{A_0r}{2(\lambda+G)}+\frac{B_0}{2Gr},\qquad \lambda=K-2G/3',
               r'\epsilon_{rr}=u_r^\prime,\qquad \epsilon_{\theta\theta}=u_r/r',
               r'\sigma_{yy}=\nu(\sigma_{rr}+\sigma_{\theta\theta})'],
    definitions=r'Stress and displacement are tension/outward positive; $\sigma_{yy}$ is the out-of-plane stress. $\nu=\lambda/[2(\lambda+G)]$. Displacement is measured from an unloaded annulus, rather than from a precompressed excavation state.',
    limit='Isotropic, small-strain, dry elasticity. The finite outer wall changes the opening hoop stress by 0.251% relative to infinite-domain Kirsch. Anisotropic loading, plasticity, excavation history and poroelastic injection remain separate tests.',
    derivation=[(r'Radial equilibrium and the isotropic plane-strain Hooke law give the displacement $Cr+D/r$. The boundary stresses determine the two integration constants.', r'd\sigma_{rr}/dr+(\sigma_{rr}-\sigma_{\theta\theta})/r=0'),
                ('Substitute zero radial stress at the opening and -P at the outer wall. As the wall moves to infinity, the finite solution approaches the isotropically loaded Kirsch limit.', r'\sigma_{rr}(a)=0,\qquad \sigma_{rr}(R)=-P')],
    reference='Finite-annulus solution derived from radial equilibrium; evaluated and checked in kirsch_reference.py and check_kirsch.py.')

SPECS['strength'] = dict(
    title='Pore pressure and material strength', kind='strength',
    labels=[('Material cell', r'Every normal velocity component is prescribed on the unit cell. Normal strains are $e$ in $x$ and $-e$ in $z$, with zero $y$ strain in 3D; 2D uses plane strain.'),
            ('Pressure and initial stress', r'Initial total stress is $-C\mathbf{I}$. Pressure is uniform and constant; permeability and mechanical pressure feedback are zero. Gravity, mesh motion and heat are disabled.'),
            ('Strength law', r'$K=G=1$ GPa, cohesion $c=100$ kPa, friction angle $30^\circ$ and dilation angle zero. EP-RSF is checked in the constant-friction limit $a=b=0$.')],
    equations=[r'\boldsymbol{\sigma}^\prime=\boldsymbol{\sigma}+q\mathbf{I},\qquad q=\alpha\langle S_e(p)p\rangle_e',
               r't_y=c\cos\phi+(C-q)\sin\phi',
               r't(e)=\min(2Ge,t_y),\qquad \epsilon_p=\max[0,e-t_y/(2G)]',
               r'N_\phi=(1+\sin\phi)/(1-\sin\phi)',
               r'A_{\rm total}=2c\sqrt{N_\phi}-(N_\phi-1)q,\qquad T_{\rm total}=T-q'],
    definitions=r'$q=\alpha p$ when saturated. $t$ is half the extreme principal-stress difference and $T$ the tensile cap of the skeleton. DES stores total stress, so it applies the effective-stress criterion by shifting the shear intercept and the tensile cap by $q$, as in the last two equations.',
    limit='A uniform-strain, fixed-pressure, zero-dilation material check on the compressive shear branch. Evolving rate-and-state friction is checked on the rate-and-state page. Tensile failure, mixed-material suction, flow, dilation and slope or injection applications are not tested here.',
    derivation=[('The effective principal stresses are -(C-q)-t and -(C-q)+t. Substituting them in Mohr-Coulomb gives the shear yield load.', r'\sigma_{\min}^\prime-N_\phi\sigma_{\max}^\prime+2c\sqrt{N_\phi}=0'),
                ('Zero dilation keeps mean stress fixed on this path. Stress first follows 2Ge and then stays at the perfect-plastic yield load; the remaining strain is plastic.', None)],
    reference='Closed-form Mohr-Coulomb response derived for the stated imposed path; stress and plastic strain are checked at every saved frame in check_pore_strength.py.')

SPECS['rainfall'] = dict(
    title='Application: rainfall triggers a shallow slide', kind='rainfall',
    labels=[('Slope', r'The slope of the pressure-driven page (2 m of soil on 20° bedrock, 40 m long, strong ends, gravity applied as a load), with a weaker soil: $c=1$ kPa, $\phi=30^\circ$, no dilation.'),
            ('Water', r'Gardner retention ($\alpha=2.4\times10^{-4}$ Pa$^{-1}$, $S_r=0.05$), Bishop $\chi=S_e$, initial suction 5 kPa. Rain $q=10^{-5}$ m/s (36 mm/h) enters through the surface; the bedrock is sealed; the toe drains at $p=0$.'),
            ('What is compared', r'The time and the mid-slope bedrock pore pressure at the onset of basal yield and when the mean plastic strain of the basal 0.6 m reaches 0.001, against the same soil loaded by bedrock pressure.')],
    equations=[r"\theta_t+\nabla\cdot\mathbf{q}=0,\qquad \mathbf{q}=-\frac{k\,k_r(S_e)}{\mu}\left(\nabla p+\rho_fg\,\hat{\mathbf{z}}\right)",
               r"\boldsymbol{\sigma}'=\boldsymbol{\sigma}+\alpha S_e\,p\,\mathbf{I}",
               r"\frac{\sigma'_1-\sigma'_3}{2}+\frac{\sigma'_1+\sigma'_3}{2}\sin\phi=c\cos\phi"],
    definitions=r"Richards flow with Gardner retention and Mualem-type relative permeability; Bishop effective stress with $\chi=S_e$ (tension positive); Mohr–Coulomb yield on the effective stresses.",
    limit=r'No exact solution: an application. DES has no runoff or ponding, so the rain rate is chosen to keep the surface unsaturated, and the surface pore pressure is checked to stay negative throughout. Failure is read at a fixed basal slip: the slip curve converges with mesh, while the pressure at the last time step for which PT still finds equilibrium does not. 2D only.',
    derivation=[(r'Each ingredient is verified separately: unsaturated flow against Srivastava–Yeh, Bishop coupling in the drained column, the sealed cell and the draining column, and effective-stress plasticity against the confined layer. Here they act together on a slope.', None),
                (r'A perched water table leaves suction, hence Bishop strength, above it and drives seepage parallel to the slope. Pressure from the bedrock instead drives an upward seepage force through the whole layer. The rainfall slope therefore needs a higher bedrock pressure to fail.', None)],
    reference='Infinite-slope stability under rainfall, as in Iverson (2000), Water Resources Research 36(7), 1897–1910; here computed fully coupled rather than from an infinite-slope factor of safety.')

SPECS['embankment'] = dict(
    title='Factor of safety of an embankment by strength reduction', kind='embankment',
    labels=[('Slope', r'A homogeneous slope $H=10$ m high with a 2:1 face, on a foundation $H$ deep; toe platform and crest 20 m long; plane strain. The base is fixed and the sides are rollers.'),
            ('Material', r'Mohr–Coulomb, $c=10$ kPa, $\phi=20^\circ$, no dilation; $E=100$ MPa, $\nu=0.3$, $\gamma=20\ \mathrm{kN/m^3}$; dry.'),
            ('Strength reduction', r'Gravity is applied elastically (from an unstressed, gravity-free checkpoint). Each trial restarts from that state with $c/F$ and $\tan\phi/F$ and runs explicitly (no PT) with DES’s default damping for 2000 s of mass-scaled time on the fixed original mesh (small strain, as the factor of safety is defined). The factor of safety lies between the largest $F$ that decelerates toward rest and the smallest that does not.'),
            ('What is compared', r'The DES factor of safety, bracketed on each mesh, and the displacement as $F$ approaches it, against Bishop’s simplified method of slices on the same slope.')],
    equations=[r"c_F=\frac{c}{F},\qquad \tan\phi_F=\frac{\tan\phi}{F}",
               r"FS=\frac{\sum\left[c\,b+(W-u\,b)\tan\phi\right]/m_\alpha}{\sum W\sin\alpha},\qquad m_\alpha=\cos\alpha\left(1+\frac{\tan\alpha\tan\phi}{FS}\right)"],
    definitions=r"$b$, $W$ and $\alpha$ are the width, weight and base inclination of a slice of the sliding mass above a circular slip surface, and $u$ the pore pressure on its base (zero here). Bishop’s FS is the minimum over circles.",
    limit=r'Bishop’s method satisfies moment equilibrium of a circular mechanism, so it is an engineering reference, not an exact solution; strength reduction finds the mechanism itself. Non-associated flow ($\psi=0$) and a finite mesh both shift the DES value. 2D, dry, small strain; on a moving mesh a failing slope slides, flattens and stops, which is a run-out question, not a factor of safety.',
    derivation=[(r'Reducing $c$ and $\tan\phi$ by the same $F$ until the slope cannot stand gives the same factor of safety as limit equilibrium defines it: the ratio of available to mobilised strength on the critical surface. Below it the explicit model dissipates the released energy and comes to rest; above it the sliding mass keeps moving.', None),
                (r'Bishop’s simplified method takes moments about the centre of a trial circle and neglects interslice shear. For a shallow circle in cohesionless sand it reproduces the infinite-slope value $\tan\phi/\tan\beta$, and the critical value converges as the number of slices grows.', None)],
    reference='Bishop (1955), Géotechnique 5(1), 7–17; strength reduction as in Griffiths and Lane (1999), Géotechnique 49(3), 387–403.',
    url='https://doi.org/10.1680/geot.1999.49.3.387')

SPECS['salencon'] = dict(
    title='Plastic tunnel: Salençon solution', kind='salencon',
    labels=[('Geometry and load', r'A quarter annulus with opening radius $a=1$ m and outer radius $R=20a$, in plane strain. The rock starts under uniform compression $P=500$ kPa in all three directions. The support pressure on the opening wall is lowered from $P$ to zero in 16 equal steps, while the outer boundary keeps the pressure $P$.'),
            ('Constraints', 'The two straight symmetry cuts have zero normal displacement. Gravity, fluid flow and heat are off, and the mesh does not move.'),
            ('Material', r'Perfectly plastic Mohr–Coulomb: $c=100$ kPa, $\phi=30^\circ$, associated flow ($\psi=30^\circ$), $G=10$ GPa, $K=21.7$ GPa ($\nu=0.3$). No softening and no weak zone.'),
            ('What is compared', r'Radial, hoop and out-of-plane stress, radial and hoop strain, and radial displacement for $1.05a\le r\le5a$ after the support is removed, on meshes refined near the opening to 0.2, 0.1 and 0.05 m, against the finite-annulus Salençon solution.')],
    equations=[r'N=\frac{1+\sin\phi}{1-\sin\phi},\qquad q=2c\sqrt{N}',
               r's_r=Cr^{N-1}-\frac{q}{N-1},\quad s_\theta=Ns_r+q\quad(a\leq r\leq b)',
               r'C=\frac{P_i+q/(N-1)}{a^{N-1}}',
               r's_r=A-B/r^2,\quad s_\theta=A+B/r^2\quad(b\leq r\leq R)',
               r'Cb^{N-1}-\frac{q}{N-1}=\frac{(2P-q)R^2+qb^2}{(N+1)R^2-(N-1)b^2}'],
    definitions=r'$s_r$ and $s_\theta$ are radial and hoop stress, compression positive (DES output is tension positive). $P_i$ is the support pressure on the wall, here zero at the end, and $b$ the plastic radius: 1.39575 m for $R=20a$ (1.39405 m for an infinite domain). $A$ and $B$ follow from continuity of $s_r$ at $b$ and $s_r=P$ at $R$. The displacement integrates the plastic flow rule inside $b$ and the elastic solution outside it; the derivation gives it.',
    limit='Dry, isotropic loading, small strain, perfectly plastic. The radial strain misses a 3% target (3.88% of $P/G$ on the finest mesh), almost all of it within half a radius of the wall. With mixed stress (needed here) the stored stress readjusts slightly at every step, so the result depends on the number of unloading steps: on the 0.05 m mesh, 16, 32 and 64 steps give hoop errors of 1.57%, 3.09% and 6.71%. Use the 16 steps shipped here, and do not hold the final load for many extra steps. Without mixed stress the coarse mesh diverges near the wall. Zero dilation is not covered.',
    derivation=[(r'Inside the plastic zone, radial equilibrium $s_r^\prime+(s_r-s_\theta)/r=0$ with the yield condition $s_\theta=Ns_r+q$ gives a power law in $r$; the wall pressure $P_i$ fixes $C$. Outside it the Lamé solution holds. Requiring the elastic stress at $b$ to reach yield gives the scalar equation for $b$, solved numerically.', r's_r^\prime+(s_r-s_\theta)/r=0'),
                (r'For the displacement, the plastic strains follow the flow rule $\epsilon_r^p=-M\epsilon_\theta^p$ with $M=(1+\sin\psi)/(1-\sin\psi)$, and the out-of-plane plastic strain is zero. Adding the elastic strains of plane-strain Hooke’s law gives $u^\prime+Mu/r$ as a known function of $r$; integrating from $b$, where $u$ matches the elastic solution, gives $u$ in the plastic zone. Strain and displacement are measured from the initially compressed state.', r'u^\prime+\frac{M}{r}u=e_r^e+Me_\theta^e')],
    reference='Salençon (1969), Contraction quasi-statique d’une cavité à symétrie sphérique ou cylindrique dans un milieu élastoplastique, Annales des Ponts et Chaussées 4, 231–236; the solution as given in the Itasca verification example for a cylindrical hole in a Mohr–Coulomb material.',
    url='https://docs.itascacg.com/itasca940/3dec/docproject/source/examples/CylindricalHoleInAnIMCM.html')

SPECS['drawdown'] = dict(
    title='Application: rapid drawdown of a reservoir against an embankment', kind='drawdown',
    labels=[('Slope', r'The embankment of the strength-reduction page ($H=10$ m, 2:1, $c=10$ kPa, $\phi=20^\circ$, $\psi=0$), saturated: porosity 0.3, saturated unit weight $20\ \mathrm{kN/m^3}$, $k=10^{-13}\ \mathrm{m^2}$, Gardner retention ($\alpha=2.4\times10^{-4}$ Pa$^{-1}$), Bishop $\chi=S_e$.'),
            ('Reservoir', r'Water stands at the crest. DES loads the face and the toe platform with its pressure, and the slope fills through the crest and the base until its pore pressure is hydrostatic.'),
            ('Drawdown', r'At $t=0$ the reservoir drops to the toe platform: the water load is removed, the face and the toe platform drain at $p=0$, and the crest and base are sealed. DES computes the undrained response and the drainage that follows.'),
            ('What is compared', r'At several times, the DES factor of safety by strength reduction (pore pressure held during each trial) against Bishop’s simplified method on the pore pressure and unit weight DES has then.')],
    equations=[r"FS=\frac{\sum\left[c\,b+(W-\chi p\,b)\tan\phi\right]/m_\alpha}{\sum W\sin\alpha},\qquad m_\alpha=\cos\alpha\left(1+\frac{\tan\alpha\tan\phi}{FS}\right)",
               r"W=b\int_{z_b}^{z_s}\gamma\,dz,\qquad \gamma=g\left[\rho_s(1-n)+\rho_f\,n\,S\right]"],
    definitions=r"$p$ and $S$ are DES’s pore pressure and saturation at that time, $\chi=S_e$ (Bishop), $n$ the porosity; the other symbols are those of the strength-reduction page.",
    limit=r'No exact solution: Bishop’s method needs the pore pressure as an input, here taken from DES, so the comparison checks the mechanics given the hydraulics. 2D, small strain; drawdown is instantaneous; the earliest comparison is a few hundred seconds after drawdown, once the mass-scaled mechanics has settled.',
    derivation=[(r'Before drawdown the reservoir’s pressure on the face balances the pore pressure inside: only buoyant weight drives sliding. Removing the water removes that support at once, and the undrained response lowers the pore pressure by about the removed load; the remaining pressure then drains through the face. The factor of safety is lowest just after drawdown and recovers as the slope drains.', None),
                (r'A strength-reduction trial holds the pore pressure: no flow, and a compressible pore fluid so that deformation does not change it. With a stiff fluid the trial would be undrained and excess pore pressure of order 10 kPa would lower the factor of safety below limit equilibrium.', None)],
    reference='Morgenstern (1963), Géotechnique 13(2), 121–131 (stability charts for rapid drawdown); Bishop (1955), Géotechnique 5(1), 7–17.')

SPECS['slope'] = dict(
    title='Application: pore pressure triggers a translational slide', kind='slope',
    labels=[('Slope', r'A soil layer of vertical thickness $D=2$ m on bedrock inclined at $\beta=20^\circ$, 40 m long, plane strain. The base is fixed and supplies water; the surface is free and drained; the vertical ends are rollers. The 4 m at each end are strong (cohesion 10 MPa) so the corners stay elastic.'),
            ('Material', r'Mohr–Coulomb, $c=2$ kPa, $\phi=35^\circ$, no dilation; $K=50$ MPa, $G=5.17$ MPa ($\nu=0.45$); saturated, $\rho_s=2650\ \mathrm{kg/m^3}$, $\phi_n=0.3$, $k=10^{-11}\ \mathrm{m^2}$.'),
            ('Loading', r'Gravity is applied as a load: the slope restarts with gravity from an unstressed, gravity-free checkpoint. The bedrock pressure then rises linearly over 2000 s; the seepage stays quasi-steady.'),
            ('What is compared', r'The bedrock pore pressure at the onset of basal yield, at fixed mean plastic strains of the basal 0.6 m in the middle third of the slope, and at the last time step for which PT still finds equilibrium, against three estimates.')],
    equations=[r"\sigma_n=\gamma D\cos^2\beta,\qquad \tau=\gamma D\sin\beta\cos\beta",
               r"FS=\frac{c+(\sigma_n-u)\tan\phi}{\tau}=1\quad\mathrm{(classical)}",
               r"\tan\phi^*=\sin\phi,\qquad c^*=c\cos\phi\quad\mathrm{(non\text{-}associated,\ }\psi=0\mathrm{)}"],
    definitions=r"$\gamma$ is the saturated unit weight and $u$ the pore pressure on the bedrock plane. The pointwise estimate applies Mohr–Coulomb to the elastic infinite slope: slope-parallel effective stress $K_0(\sigma_n-u)$ with $K_0=\nu/(1-\nu)$, slope-normal $\sigma_n-u$, shear $\tau$.",
    limit=r'No exact solution: an application compared with engineering estimates. The slip curve converges with mesh and step; the pressure at the last time step for which PT still finds equilibrium does not, so it is not read as a collapse pressure. There is no dynamic run-out. 2D, saturated, zero dilation. For the classic embankment by strength reduction see the embankment page.',
    derivation=[(r'For a long slope, statics fix the normal and shear stress on planes parallel to the surface. Gravity applied to an unstressed slope with zero strain along it gives slope-parallel stress $K_0\sigma_n$; a pore-pressure rise reduces the slope-normal effective stress by $u$ and the parallel one by $K_0u$. The first $u$ at which that state reaches Mohr–Coulomb is the pointwise onset, 14.1 kPa.', None),
                (r'A zero-dilation shear band along the bedrock fails with the reduced parameters above, giving 16.5 kPa; the classical factor of safety, which assumes associated flow, gives 20.8 kPa. A non-associated (coaxial) material fails below the associated limit load.', None)],
    reference='Drescher and Detournay (1993), Géotechnique 43(3), 443–456 (limit load of translational mechanisms for non-associative materials, with Davis’s reduced parameters).',
    url='https://experts.umn.edu/en/publications/limit-load-in-translational-failure-mechanisms-for-associative-an/')

SPECS['rsf'] = dict(
    title='Rate-and-state friction under a velocity step', kind='rsf',
    labels=[('Cell', r'The material-strength unit cell: every node prescribed, normal extension rates $\pm\dot e$, plane strain, initial stress $-C\mathbf{I}$ with $C=1$ MPa and a fixed pore pressure of 0.8 MPa ($\alpha=1$), cohesion 100 kPa, no dilation.'),
            ('Friction', r'EP-RSF with the aging law: $a=0.01$, $b=0.015$, $D_c=1$ mm, $\mu_0=\tan30^\circ$, $V_0=V_1$, so the cell starts and stays at steady state; slip rate $V=2w\dot\epsilon_{II}$ (total-strain-rate projection).'),
            ('Rate step', r'After 64 s at $\dot e_1=7.8\times10^{-6}\ \mathrm{s^{-1}}$ the cell restarts from its last checkpoint at $10\dot e_1$ for another 64 s.'),
            ('What is compared', r'Friction, state and principal shear stress at every step against the exact discrete recursion and the continuum solution, for steps of 1, 0.5, 0.25 and 0.125 s.')],
    equations=[r"\mu=\mu_0+a\ln\frac{V}{V_0}+b\ln\frac{V_0\theta}{D_c},\qquad \dot\theta=1-\frac{V\theta}{D_c}",
               r"\theta(t)=\frac{D_c}{V_2}+\left(\frac{D_c}{V_1}-\frac{D_c}{V_2}\right)e^{-V_2(t-t_1)/D_c},\quad t>t_1",
               r"\tau=\min\left(2G\,e,\ c\cos\phi+(C-\alpha p)\sin\phi\right),\quad \tan\phi=\mu"],
    definitions=r"$\theta$ is the state variable, $V$ the slip rate, $w$ the element's minimum altitude and $\tau$ half the principal-stress difference. DES advances $\theta$ by forward Euler, once per physical step.",
    limit=r'Aging law and the total-strain-rate slip projection only; constant pore pressure; no dilation, no elastic interaction (every node prescribed), so stick–slip instability is not tested.',
    derivation=[(r'With $V$ held constant the aging law is linear in $\theta$ and integrates to the exponential above; at $t_1$ friction jumps by $a\ln10$ through the direct effect and relaxes to the new steady state, a net change of $(a-b)\ln10$.', None),
                (r'With every node prescribed and no dilation, the stress is elastic until yield and then equals the Mohr–Coulomb load at the current friction. The discrete reference repeats DES’s forward-Euler update; the continuum one is exact.', None)],
    reference='Dieterich (1979), Journal of Geophysical Research 84(B5), 2161–2168; Ruina (1983), Journal of Geophysical Research 88(B12), 10359–10370.',
    url='https://doi.org/10.1029/JB084iB05p02161')

SPECS['borehole'] = dict(
    title='Fluid injection into a poroelastic borehole', kind='borehole',
    labels=[('Geometry', r'A quarter annulus $a=1$ m, $R=20a$, plane strain, with symmetry rollers on both cuts; meshes refined out to $5a$.'),
            ('Loading', r'At $t=0^+$ the wall takes fluid at $p_w=1$ MPa: pore pressure $p_w$ (a 1 ms ramp) and normal traction $-p_w$. The outer arc is sealed and traction-free. The rock starts unstressed at $p=0$; gravity is off.'),
            ('Material', r'$K=200$ MPa, $G=100$ MPa, $\alpha=1$, $1/M=\phi/K_f=10^{-9}\ \mathrm{Pa^{-1}}$, $k/\mu=10^{-9}\ \mathrm{m^2/(Pa\,s)}$; diffusivity $c=0.25\ \mathrm{m^2/s}$.'),
            ('What is compared', r'Nodal pressure and element total radial and hoop stress for $1.05a\leq r\leq4a$ at 9, 18, 27 and 36 s, relative to $p_w$.')],
    equations=[r"S\,p_t+2\alpha\,\dot C_1=\frac{k}{\mu}\frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial p}{\partial r}\right),\qquad S=\frac{1}{M}+\alpha m,\quad m=\frac{\alpha}{\lambda+2G}",
               r"\sigma_{rr}=-2Gm\frac{J}{r^2}+2(\lambda+G)C_1-2G\frac{C_2}{r^2},\qquad \sigma_{\theta\theta}=2Gm\left(\frac{J}{r^2}-p\right)+2(\lambda+G)C_1+2G\frac{C_2}{r^2}",
               r"J(r)=\int_a^r p\,r'\,dr',\qquad \sigma_{rr}(a)=-p_w,\quad \sigma_{rr}(R)=0"],
    definitions=r"Tension positive. The radial displacement is $u=mJ/r+C_1r+C_2/r$, so the volumetric strain is $mp+2C_1$; $C_1(t)$ and $C_2(t)$ follow from the two traction conditions and are linear in $p_w$ and $J(R)$.",
    limit=r'Axisymmetric loading only (no deviatoric in-situ stress), linear poroelasticity, small strain. The finite annulus is part of the reference; the infinite-medium Detournay–Cheng solution is its $R\to\infty$ limit.',
    derivation=[(r'Radial symmetry makes the displacement irrotational, so equilibrium integrates once: the volumetric strain is $mp$ plus a uniform term $2C_1$. Fluid mass, $\alpha\dot\varepsilon_v+\dot p/M=\nabla\cdot(k\nabla p/\mu)$, then gives the diffusion equation above.', None),
                (r'On a finite annulus $C_1\neq0$ and depends on $J(R)$, which couples the pressure to its own integral; it moves the pressure by about 1% of $p_w$ here. In Laplace space the pressure is a Bessel ($I_0$, $K_0$, $I_1$, $K_1$) solution plus a uniform part $-2\alpha\bar C_1/S$, with $\bar C_1$ from one scalar equation; Stehfest inversion agrees with an independent Crank–Nicolson solution to $5\times10^{-6}$ of $p_w$.', None)],
    reference='Detournay and Cheng (1988), International Journal of Rock Mechanics and Mining Sciences 25(3), 171–182, axisymmetric loading modes; the finite annulus is solved in borehole_reference.py.',
    url='https://experts.umn.edu/en/publications/poroelastic-response-of-a-borehole-in-a-non-hydrostatic-stress-fi/')

SPECS['injection_temperature'] = dict(
    title='Hot-water injection: Avdonin solution and heat ledger', kind='injection_temperature',
    labels=[('Option', r'injection.temperature gives the temperature $T_{in}$ of the fluid of each point source; it requires fluid thermal advection. Unset, the injected fluid carries no heat, and any heat source must be given separately.'),
            ('(a) Avdonin', r'A line source injects $q=2\times10^{-5}\ \mathrm{m^2/s}$ at $T_{in}=333$ K into a layer at $T_0=283$ K. DES models a quarter plane, 8 m square, with the source at its corner, sealed and insulated symmetry sides, drained far sides and a stiff skeleton without pressure feedback.'),
            ('(b) Doublet', r'An 8 m sealed, insulated box with an injector at $T_{in}$ and a producer 4 m away at the same rate $q=10^{-4}\ \mathrm{m^2/s}$. The hot front reaches the producer, so heat leaves at a changing temperature.'),
            ('What is compared', r'(a) Nodal temperature against the similarity solution at four times; stored against injected heat. (b) Stored heat against the independent ledger of injected minus extracted heat at every step.')],
    equations=[r"C\,T_t+\rho_fc_f\,\mathbf{q}\cdot\nabla T=\nabla\cdot(\lambda\nabla T)+\rho_fc_f\,Q\,(T_{in}-T)",
               r"\frac{T-T_0}{T_{in}-T_0}=\frac{\Gamma(\nu,u)}{\Gamma(\nu)},\qquad u=\frac{r^2C}{4\lambda t},\quad \nu=\frac{\rho_fc_f\,q}{4\pi\lambda}",
               r"E(t)-E(0)=\rho_fc_f\,q\int_0^t\left(T_{in}-T_{prod}\right)ds"],
    definitions=r"$C$ is the bulk heat capacity, $Q$ the injected volume rate per unit volume, $\Gamma(\nu,u)$ the upper incomplete gamma function and $E=\sum_nC A_nT_n$ the stored sensible heat. $T_{prod}$ is the temperature at the producer, interpolated with the element's shape functions.",
    limit=r'Fixed heat capacity: heat carried by fluid that goes into storage is not tracked, so a sealed cell storing all its injection is outside the scope. No pressure work, dissipation or thermal expansion. Thermal advection uses the Darcy flux at the start of each step, which loses part of the first step of heat when a source starts from rest and that step needs several thermal substeps.',
    derivation=[(r'The conservative balance $\partial_t(CT)+\nabla\cdot(\rho_fc_f\mathbf qT)=\nabla\cdot(\lambda\nabla T)+\rho_fc_fQT_{in}$ gives the advective form above once fluid mass conservation, $\nabla\cdot\mathbf q=Q$, removes $\rho_fc_fT\,\nabla\cdot\mathbf q$. Extraction removes fluid at the local temperature, which the advective term already carries, so extracting points add nothing.', None),
                (r'For radial flow $\mathbf q=q\,\hat r/(2\pi r)$ the equation has a similarity solution in $u$; the injected heat $\rho_fc_fq(T_{in}-T_0)t$ is all stored, since none leaves. avdonin_reference.py checks both.', None)],
    reference='Avdonin’s radial thermal-injection problem (1964). The similarity solution is derived here and checked against the equation and the heat balance in avdonin_reference.py; the ledger follows from energy conservation alone.')

SPECS['drainage'] = dict(
    numerical_reference=True,
    title='Drainage of a saturated column into the unsaturated state', kind='drainage',
    labels=[('Column', r'A laterally confined sand column $H=1$ m high, saturated and hydrostatic with the water table at its top. Top sealed and traction-free; base fixed.'),
            ('Drainage', r'At $t=0^+$ the base is opened to $p=0$ (a 1 ms ramp of the scheduled boundary pressure). Water drains out, the top desaturates and suction develops.'),
            ('Material', r'After Liakopoulos: $E=1.3$ MPa, $\nu=0.4$, $\alpha=1$, $\phi=0.2975$, $k=4.5\times10^{-13}\ \mathrm{m^2}$; van Genuchten $\alpha_{vG}=1/3000\ \mathrm{Pa^{-1}}$, $n=1.8$, $S_r=0.1$, Mualem permeability; Bishop $\chi=S_e$. Gravity and pressure feedback (two-way coupling) are on.'),
            ('What is compared', r'Nodal pressure and the settlement of the top at 5, 20, 60 and 120 min against an independent solution of the same continuum.')],
    equations=[r"c(p)\,p_t+\alpha S_e\,\dot\varepsilon=-\partial_z q,\qquad q=-\frac{k\,k_r(S_e)}{\mu}\left(\partial_z p+\rho_f g\right)",
               r"\dot\varepsilon=\frac{\dot\sigma+\alpha\,\Pi'(p)\,p_t}{K+\frac{4}{3}G},\qquad \Pi=S_e\,p",
               r"\dot\sigma(z)=-g\rho_f\phi(1-S_r)\int_z^H\frac{dS_e}{dp}\,p_t\,dz'"],
    definitions=r"$c(p)$ is the fluid storage including $\phi\,dS/dp$; $\varepsilon$ the vertical strain, $\sigma$ the total vertical stress, which changes only by the weight of the water that leaves the column above. $k_r$ is the Mualem relative permeability.",
    limit=r'No closed form: the reference is numerical. Uniaxial strain, small strain, constant porosity, no hysteresis. The base opens instantly, so the first minutes are limited by the step.',
    derivation=[(r'Uniaxial strain makes the vertical effective-stress change equal to the constrained modulus times the strain. The Bishop stress adds $\alpha\Pi$, and the total stress follows the weight, so drying above a point unloads it: the storage becomes nonlocal.', None),
                (r'drainage_reference.py integrates this system with vertex-centred finite volumes, the dense coupled storage matrix and an implicit BDF method. Grid differences fall at second order (5.2 to 1.3 Pa at 300 s on 100 and 200 cells, of a 9.8 kPa scale); the comparison uses 200 cells.', None)],
    reference='Liakopoulos (1965) experiment as benchmarked by Schrefler and Zhan (1993), Water Resources Research 29(1), 155–167; reference computed independently in drainage_reference.py.',
    url='https://doi.org/10.1029/92WR01737')

SPECS['pressure_yield'] = dict(
    title='Pressure-driven yielding of a confined layer', kind='pressure_yield',
    labels=[('Layer', r'A saturated layer $H=2$ m deep under a drained surcharge $q=20$ kPa, water table at the top, an isotropic lithostatic initial stress with hydrostatic pressure. Sides and base are fixed normally and sealed.'),
            ('Loading', r'The top pressure rises linearly over 2000 s to 95% of the apex pressure; permeability is high, so the whole layer follows. Gravity acts; pressure feedback is on.'),
            ('Material', r'Elastic–perfectly plastic Mohr–Coulomb: $K=100$ MPa, $G=60$ MPa, $c=10$ kPa, $\phi=20^\circ$, no dilation, $\alpha=1$, $\phi_n=0.3$.'),
            ('What is compared', r'For each element, the pressure rise at the first step with plastic strain against the exact onset $\Delta p^*$ at its centroid depth.')],
    equations=[r"\sigma'_v=-(\gamma'd+q)+\alpha\Delta p,\qquad \sigma'_h=-(\gamma'd+\kappa q)+\kappa\alpha\Delta p",
               r"\frac{\sigma'_1-\sigma'_3}{2}+\frac{\sigma'_1+\sigma'_3}{2}\sin\phi=c\cos\phi\ \mathrm{at}\ \Delta p=\Delta p^*(d)",
               r"\Delta p_{\mathrm{apex}}=\frac{c\cot\phi+q}{\alpha}"],
    definitions=r"Tension positive; $d$ is depth, $\kappa=\nu/(1-\nu)$, $\gamma'=(\rho_b-\alpha\rho_f)g$ and $\sigma'_1\geq\sigma'_3$ are the principal effective stresses. $\Delta p^*(d)$ is the first root of the yield function of the elastic trial state.",
    limit=r'Uniaxial strain, a uniform pressure rise, no dilation, fully saturated. The test stops short of the apex, beyond which no equilibrium exists. 3D is checked only to onset; an inclined (infinite-slope) version needs an initial-stress option DES does not have.',
    derivation=[(r'Under uniaxial strain the horizontal stresses stay equal and change by $\kappa$ times the vertical change; the total vertical stress is fixed by the weight and the surcharge, so a uniform pressure rise moves only the effective stresses, as above.', None),
                (r'$\Delta p^*$ increases with depth, so a yield front moves down. Each depth is independent under uniaxial strain, so the front stays exact after the layer above yields; the yielded state follows the cone until its apex.', None),
                (r'The surcharge keeps the top in compression: without it the yielded top reaches the apex within a few kPa of the onset. pressure_yield_reference.py checks the onset on the yield surface, the elastic state just below it, the surcharge-free closed form and the ordering of onset and apex.', None)],
    reference='Terzaghi effective stress with Mohr–Coulomb yield; the uniaxial-strain stress path as in Jaeger, Cook and Zimmerman, Fundamentals of Rock Mechanics (2007). Computed in pressure_yield_reference.py.')

SPECS['anisotropic_time'] = dict(
    title='Anisotropic flow with changing permeability', subtitle='Two directional decay rates across a prescribed tensor change',
    labels=[('Domain and drainage', r'A 1 m square, $0<x<1$, $-1<z<0$; all four sides hold $p=0$. Gravity, heat and mesh motion are off.'),
            ('Initial pressure', r'Two sine modes have relative amplitudes 1 and 0.35. Their different horizontal wave numbers distinguish the two directional diffusivities.'),
            ('Coefficient change', r'At $t_s=0.5$ s, scalar permeability doubles and the axis-aligned tensor changes. Pressure and storage are continuous. The second phase is a restart at $t_s$ with the new coefficients; a further restart part-way through it reproduces the uninterrupted run.')],
    equations=[r'S_0\,\partial_t p=\mu^{-1}\nabla\!\cdot[\mathbf{k}(t)\nabla p]',
               r'S_0=\phi/K_f+(\alpha-\phi)/K_s',
               r'p=P_0\sum_{n=1}^{2}a_n(t)\sin(n\pi x)\sin[\pi(z+1)]',
               r'\lambda_{n,j}=\frac{\pi^2}{\mu S_0}(n^2 k_{x,j}+k_{z,j})',
               r'\tau_1=\min(t,t_s),\qquad\tau_2=\max(0,t-t_s)',
               r'a_n(t)=a_n(0)\exp[-\lambda_{n,1}\tau_1-\lambda_{n,2}\tau_2]',
               r'\mathbf{k}_1=k_0\,\mathrm{diag}(4,1),\qquad\mathbf{k}_2=k_0\,\mathrm{diag}(2,18)'],
    definitions=r'Lengths in the modes are divided by the 1 m side. $P_0=1000$ Pa, $a_1(0)=1$, $a_2(0)=0.35$, $k_0=10^{-15}$ m², $\mu=1.002\times10^{-3}$ Pa s, $\phi=0.3$, $\alpha=1$, $K_f=2.17$ GPa and $K_s=37$ GPa. ',
    limit='Prescribed, homogeneous, axis-aligned tensors; fixed storage and geometry. Mechanical pressure feedback is off, so deformation does not change the pressure; pressure still loads the solid, so the time step must stay within the elastic stability limit. This tests the flow equation, not coupled deformation or a general fluid-content budget.',
    derivation=[('Each sine mode vanishes on the drained faces. Substitution into the diffusion equation gives its scalar decay rate; the different horizontal wave numbers make both tensor components observable.', None),
                ('Integrate the decay rate in each time interval and carry the amplitude continuously through the coefficient jump. This gives the exponential formula without fitted parameters.', None),
                ('The reference stored fluid volume (pressure storage only) is $S_0\\int p\\,dV=4S_0P_0a_1/\\pi^2$: the second horizontal mode integrates to zero. This excludes deformation storage because its feedback is disabled.', None)],
    reference='Separation of variables for the stated linear diffusion equation; the independent reference is implemented in check_anisotropic_flow.py.')

SPECS['content'] = dict(
    title='Water balance in a deforming column', subtitle='Pressure, deformation and open-boundary flow measured independently',
    labels=[('Geometry', 'A 1 m confined column; 0.1 m wide on the selected 2D mesh. Fixed geometry, zero gravity, elastic skeleton.'),
            ('Mechanical boundaries', r'Base vertical motion and side-normal motion are fixed. The top has zero total normal traction: $\sigma_{zz}=0$.'),
            ('Water boundaries', r'Initially $p=0$. The top stays drained; the base supplies $q_i=10^{-6}$ m/s. Sides are sealed.')],
    equations=[r'S=\phi/K_f+(\alpha-\phi)/K_s,\qquad H_m=K+4G/3',
               r'C_o=S+\alpha^2/H_m,\qquad m=k/\mu,\qquad D=m/C_o',
               r'k_n=(2n+1)\pi/(2H)',
               r'p(d,t)=\frac{q_i}{m}\left[d-\frac{2}{H}\sum_{n=0}^{\infty}\frac{(-1)^n\sin(k_nd)}{k_n^2}e^{-Dk_n^2t}\right]',
               r'\bar p=\frac{q_i}{m}\left[\frac{H}{2}-\frac{2}{H^2}\sum_{n=0}^{\infty}\frac{(-1)^n}{k_n^3}e^{-Dk_n^2t}\right]',
               r'V_s=AHC_o\bar p,\qquad V_o=Aq_it-V_s',
               r'Z=\int_V(Sp+\alpha\varepsilon_v)\,dV',
               r'Z_{\rm split}=Z+\sum_eV_eL\overline{(p_n-p_{n-1})}_e',
               r'L=\alpha^2/B,\qquad B_{2D}=K+G/3,\qquad B_{3D}=K'],
    definitions=r'$d$ is depth below the drain, $A$ cross-section and $H=1$ m. Stress and volume strain are tension-positive. $K=20$ MPa, $G=10$ MPa, $\alpha=0.8$, $\phi=0.25$, $K_f=100$ MPa, $K_s=1$ GPa, $k=10^{-11}$ m² and $\mu=10^{-3}$ Pa s. The steady basal pressure is 100 Pa.',
    derivation=[('Vertical equilibrium of the confined, traction-free column gives strain proportional to pressure. Insert it in physical Biot content to obtain the constrained capacity.', r'\varepsilon_v=\alpha p/H_m,\qquad \zeta=C_op'),
                ('Mass balance with Darcy mobility gives diffusion and its mixed boundary conditions.', r'C_op_t=mp_{dd},\qquad p(0,t)=0,\qquad mp_d(H,t)=q_i'),
                ('Subtract the steady profile. The remaining field has zero value at the top and zero derivative at the base, giving the odd half-wave sine eigenfunctions. Projection of the initial negative steady profile yields the pressure series above.', None),
                ('Integrate the pressure series for stored volume; differentiate it at the top for independent outward Darcy flux.', r'q_o=q_i\left[1-\frac{2}{H}\sum_{n=0}^{\infty}\frac{(-1)^n}{k_n}e^{-Dk_n^2t}\right]'),
                ('With PT, the mechanics of each step uses the pressure increment of the previous step. The increment not yet passed to the mechanics is a small pending term, and adding it to the stored water closes the budget exactly. It is a numerical lag, not physical water, and it halves when the time step is halved.', None)],
    reference='Separation of variables and mesh-independent quadrature for the stated mixed-boundary Biot diffusion problem; fluid_content.py also checks the PDE and endpoint derivatives.',
    limit='Fixed geometry, saturated flow and constant material properties. This example does not check storage that changes with deformation, large-strain mass, mass conservation through remeshing, or heat. A closed budget does not by itself mean an accurate equilibrium; the two are checked separately.')

SPECS['thermal_cell'] = dict(
    title='Thermal pressurization in a sealed cell', subtitle='An exact heat pulse with separate PT and open-energy controls',
    labels=[('Geometry', 'Rigid unit material cell, zero gravity, fixed mesh and constant porosity. Every water boundary is sealed in the selected result.'),
            ('Heating', 'Temperature starts at 273 K. A uniform supplied power raises it at 40 K/s for 0.25 s, followed by a hold through 1 s.'),
            ('Coefficients', 'No conduction, Darcy permeability or solid thermal expansion. The heat is supplied by one point source per node, each proportional to the volume that node represents.')],
    equations=[r'\theta=T-T_0,\qquad S=\phi/K_f+(\alpha-\phi)/K_s',
               r'\beta=\phi\alpha_f,\qquad C_0=(1-\phi)\rho_sc_s+\phi\rho_fc_f',
               r'\theta=40\min(t,0.25)\ {\rm K}',
               r'p=\beta\theta/S,\qquad \Delta\sigma_{ii}=-\alpha p',
               r'E-E_0=Q_{\rm supplied}=VC_0\theta',
               r'(S+L)p_n=\beta\theta_n+Lp_{n-1}',
               r'a=\rho_fc_f\beta,\qquad R_E=-Va\theta^2/2'],
    definitions=r'$K=20$ MPa, $G=10$ MPa, $K_s=1$ GPa, $K_f=100$ MPa, $\alpha=0.98$, $\phi=0.25$ and $\alpha_f=2\times10^{-4}$ K$^{-1}$. $C_0=3.0325$ MJ/(m³ K). $L=\alpha^2/B$ is the PT stabilizer, with $B=K+G/3$ in 2D and $B=K$ in 3D. $R_E$ is the separate drained control’s retained-fluid enthalpy residual.',
    derivation=[('For a sealed rigid cell, total volumetric strain and fluid-content change vanish. This gives the physical pressure directly.', r'\zeta=Sp+\alpha\varepsilon_v-\beta\theta=0'),
                ('The fixed-capacity heat equation integrates the supplied pulse exactly. The selected explicit pressure follows the physical content equation. PT instead retains its pending increment; the displayed recurrence is exact for that split, and its physical-pressure error is first order for small timesteps.', None),
                ('For the additional drained rigid cell, all vertices hold zero pressure. Fluid expansion must leave through the boundary reaction.', r'V_o=V\beta\theta'),
                ('The retained fluid mass changes the heat capacity. Outflow carries enthalpy at the changing cell temperature.', r'C(\theta)=C_0-a\theta\qquad H_o=Va(T_0\theta+\theta^2/2)'),
                ('Use retained mass in the energy change, rather than integrating the original fixed capacity.', r'\Delta E_{\rm ret}=V[(C_0-a\theta)(T_0+\theta)-C_0T_0]'),
                ('Adding outflow enthalpy and subtracting the supplied heat produces the displayed negative half-quadratic deficit. At 10 K it is −10,450 J, or 0.03446% of supplied heat, independent of timestep. There is no strain work or pressure-boundary work in this rigid, zero-pressure control.', None)],
    reference='Linear Biot thermal expansion and integration of the stated heat pulse; the retained-fluid deficit is derived for the implemented linear fluid-content law.',
    limit='This sealed cell checks thermal pressurization and the fixed-capacity heat equation. With PT the pressure can lag noticeably during rapid heating even when the force balance converges. In a drained cell a small, persistent energy deficit remains, so open systems such as thermal plumes are not covered; shear heating is not included.')

SPECS['plastic_content'] = dict(
    title='Water balance through plastic yielding', subtitle='Total volume change remains observable after an elastic stress relation ceases to apply',
    labels=[('Geometry and motion', 'Unit column with a fixed base and fixed side-normal motion. The top moves down at 0.001 m/s for 1 s. Plane strain in 2D; zero gravity and a fixed mesh.'),
            ('Water boundaries', 'All boundaries are sealed in the selected case. Initial pressure is zero; porosity and pressure capacity remain fixed.'),
            ('Yield law', 'Mohr–Coulomb: cohesion 100 Pa, friction 20° and selected dilation 10°. Other controls use zero and negative dilation.')],
    equations=[r'S=\phi_n/K_f+(\alpha-\phi_n)/K_s',
               r'Z=\int_V(Sp+\alpha\varepsilon_v)\,dV',
               r'\bar p(t)=\alpha vt/(HS)',
               r'L=\alpha^2/B,\qquad B_{2D}=K+G/3,\qquad B_{3D}=K',
               r'Z_{\rm split}=Z+VL(\bar p_n-\bar p_{n-1})',
               r'\bar p_n=\frac{\alpha vt_n/H+L\bar p_{n-1}}{S+L}',
               r'(S+L)\Delta p+\alpha\Delta\varepsilon_v-LD=0'],
    definitions=r'$K=20$ MPa, $G=10$ MPa, $\alpha=1$, $\phi_n=0.25$, $K_f=1$ MPa and $K_s=10^{20}$ Pa. $v=10^{-3}$ m/s and $H=1$ m. $D$ is the previous raw pressure increase consumed by mechanics. The physical pressure reaches 4 kPa at 1 s; the 0.01 s PT step leaves about 0.17% pending lag.',
    derivation=[('Integrate total volumetric strain over the actual mesh. The imposed boundary velocity fixes the volume loss even after the stress return becomes plastic.', r'\int_V\varepsilon_v\,dV=-Avt'),
                ('Sealed linear Biot conservation gives the physical mean pressure. Adding the separately recorded pending-pressure contribution gives the exact PT recurrence shown above.', None),
                ('Converting the stress change back to volume with the elastic moduli fails after a plastic return, including the out-of-plane plastic strain in plane strain. DES’s fixed-stress source therefore uses the accepted volumetric strain of the step and the previous raw pressure increment; pressure substeps reuse the same accepted mechanical history.', None),
                ('With thermal expansion the direct strain source uses the full thermal fluid coefficient, without an additional elastic stress conversion.', r'\beta=\phi_n\alpha_f+(\alpha-\phi_n)\alpha_s'),
                ('The partly saturated homogeneous Gardner control has an independent scalar content inversion at each accepted step. It checks the implemented retention-rate law; it does not establish a finite-strain porosity law.', r'(S+L)\Delta p+\phi_n\Delta S_e+\alpha S_{e,\rm old}\Delta\varepsilon_v-LD=0'),
                ('Additional EP-RSF controls check the independent aging update and effective-stress return while pressure evolves. Open drains include measured boundary exchange; they check material content, rather than a separate spatial Darcy-flux solution.', None)],
    reference='Integral Biot conservation with prescribed boundary motion; independent scalar pressure recurrence and local effective-stress/thermal controls in check_plastic_content.py.',
    limit='Fixed porosity and geometry. The pending PT term is a numerical coupling lag, not physical water. Porosity change from plastic strain (including enhancement and healing), remeshing, the energy balance of open systems and shear heating are not covered.')

# Index order: themes from simple to coupled; builders list records in this order.
THEMES = (('groundwater', 'Groundwater flow'), ('unsaturated', 'Unsaturated flow'),
          ('consolidation', 'Poroelastic consolidation'), ('injection', 'Injection-induced stress'),
          ('dynamics', 'Poroelastic dynamics'), ('mechanics', 'Material strength'), ('tunnels', 'Circular openings'),
          ('heat', 'Heat and groundwater'))
SPECS['plastic_heat'] = dict(
    title='Plastic shear heating and sealed pressurization', kind='prescribed pure-shear cell',
    labels=[('Geometry and loading', 'Unit square, −1 ≤ z ≤ 0; plane strain. Prescribed normal velocities give ε̇xx=0.5/s, ε̇zz=−0.5/s, zero spin and zero volume strain. Four nodes and two triangles; mesh motion disabled.'),
            ('Initial state', 'T=300 K, p=0, isotropic compression 30 MPa; K=20 GPa, G=1 GPa, cohesion 10 MPa. Aging EP-RSF: a=0.01, b=0.015, Dc=0.01 m, V0=0.001 m/s, θ0=Dc/V0.'),
            ('Water and heat', 'Every face sealed and insulated. Porosity 0.2, Biot coefficient 0.8, zero permeability/conductivity, grain/fluid moduli 100/2.2 GPa. Fixed mixture capacity C=2,996,000 J/(m³ K); β=4×10⁻⁵/K.'),
            ('Numerical method', 'PT option 1; mixed stress disabled; 1000 physical steps of 0.0002 s. Heat is deposited on the following physical interval; pending energy and pressure are separate numerical reservoirs.')],
    equations=[
        r'\gamma=2\dot\epsilon t,\quad \tau_{\rm EP}=\min(G\gamma,c),\quad Q_{\rm EP}=c(\gamma-c/G)_+',
        r'\theta(t)=D_c/V+(\theta_0-D_c/V)e^{-Vt/D_c}',
        r'\mu(t)=\tan30^\circ+a\ln(V/V_0)+b\ln(V_0\theta/D_c)',
        r'\tau_y(t)=\frac{c+P\mu(t)}{\sqrt{1+\mu(t)^2}},\quad G\gamma(t_y)=\tau_y(t_y)',
        r'Q(t)=\dot\gamma\int_{t_y}^t\tau_y(s)\,ds-\frac{\tau_y(t)^2-\tau_y(t_y)^2}{2G}\quad(t>t_y)',
        r'\Delta T=Q/C,\quad p=\beta\Delta T/S_0,\quad S_0=\phi/K_f+(\alpha-\phi)/K_s',
        r'E_T+E_{\rm pending}=\sum_n V_e\Delta W^p_e,\quad \Delta W^p=\bar\sigma^{\prime}:\Delta\epsilon^p'],
    definitions=r'$\tau=(\sigma_{xx}-\sigma_{zz})/2$ is the principal half-span, conjugate to $\gamma$, not a prescribed fault traction. $P=30$ MPa is effective confinement, $Q$ heat per volume, and $V=2w\dot\epsilon$ with $w=1/\sqrt2$ m. The integral reference uses tight numerical quadrature of the analytical aging solution; $Q=0$ before first yield.',
    derivation=[('Solve the aging ODE at constant slip rate. The decreasing yield strength intersects the elastic predictor once; after this point the cell stays at yield.', r'\dot\theta=1-V\theta/D_c,\quad \tau=\min(G\gamma,\tau_y)'),
                ('Boundary work is τ dγ. Subtract the change in elastic shear energy τ²/(2G) to obtain plastic heat; unit heat fraction is assumed.', r'dQ=\tau\,d\gamma-d(\tau^2/2G)'),
                ('Fixed volume and sealed fluid give zero content. With constant capacity and no heat loss, temperature follows Q/C. DES’s pending-pressure reservoir vanishes under time refinement.', r'\zeta=S_0p-\beta\Delta T=0')],
    limit='A material-level earthquake-heating qualification with constant effective confinement. It does not demonstrate thermal weakening at fixed total normal stress, a propagating rupture, general enthalpy conservation with changing storage, or a complete earthquake cycle.',
    reference='Direct boundary-work/elastic-storage identity and analytical aging law, evaluated in plastic_heat.py. Integration of the nonelastic work proposed in Eh Tan’s PR #88; thermal fault weakening motivates further fault-scale models (Rice, 2006, doi:10.1029/2005JB004006).')

SPECS['fault_heat'] = dict(
    title='Fault-zone heating at fixed normal traction', kind='sealed direct-shear patch',
    labels=[('Geometry and motion', '10 mm square, four nodes and two triangles; plane strain and fixed reference mesh. The base is fixed; both sides impose a linear horizontal velocity from zero at the base to 1 m/s at the top. The top normal platen is free under 50 MPa compression.'),
            ('Initial material', 'Normal stresses −50 MPa, shear stress 25.866 MPa at yield; T=300 K, p=0. Mohr–Coulomb: K=20 GPa, G=10 GPa, cohesion 1 MPa, friction 30°, dilation zero. Physical mixture density 2360 kg/m³; damping and PT off.'),
            ('Water and heat', 'All faces sealed and insulated, gravity zero. Porosity 0.2, Biot coefficient 0.8, grain/fluid moduli 100/2.2 GPa. Fixed mixture heat capacity C=2,996,000 J/(m³ K), β=4×10⁻⁵/K; permeability, conductivity and solid expansion zero.')],
    equations=[r'\gamma=Vt/H,\quad \sigma_{zz}=-N,\quad \ell=K-2G/3',
               r'd=(\sigma_{xx}+N)/2,\quad m=(\sigma_{xx}-N)/2,\quad q=\sqrt{d^2+\tau^2}',
               r'q+\sin\phi(m+\alpha p)-c\cos\phi=0',
               r'Sp\prime+\alpha e=\beta T\prime,\quad CT\prime=rq',
               r'(\ell+2G+\alpha^2/S)e+(Ga-\alpha\beta q/(SC))r=\tau',
               r'[\sin\phi(\ell+G)-Ga]e-Gr=-Gb',
               r'\sigma_{xx}\prime=\ell e-Gra+\tau-\alpha p\prime',
               r'\tau\prime=G-Grb-d,\quad W^p{}\prime=rq'],
    definitions=r'Prime means differentiation with respect to γ, a=d/q, b=τ/q, e=εv′, and r is the plastic multiplier with engineering plastic shear rb and plastic work rq. S=ϕ/Kf+(α−ϕ)/Ks. The reference imposes normal equilibrium and checks that σyy=−N+ℓεv−αp stays the intermediate principal.',
    derivation=[('The independent continuous reference solves normal equilibrium and yield consistency for the volume strain rate and plastic multiplier. It integrates total stress, fluid content, temperature and work; it does not call or duplicate the DES return algorithm.', None),
                ('Zero dilation gives opposite plastic normal strains. The stress equations include the Jaumann spin from simple shear; the sealed fluid law retains normal compliance.', r'\epsilon^p_{xx}{}\prime=ra/2,\quad \epsilon^p_{zz}{}\prime=-ra/2,\quad \gamma^p{}\prime=rb'),
                ('A closed-form stiff-skeleton limit takes q≈τ, r≈1 and εv≈αp/(ℓ+G). It explains exponential thermal weakening; the plotted reference retains the finite-modulus equations above.', r'\Lambda=\beta/[S+\alpha^2/(\ell+G)],\quad \tau\simeq\tau_0e^{-\alpha\sin\phi\Lambda\gamma/C}'),
                ('Pending heat is included in accepted work. Exact power integration is checked only for cases saving every physical interval, including the first 20 intervals after switching slow GVS to physical mass.', r'E_T+E_{\rm pending}=\sum_n V\Delta t_n q_n')],
    limit='A uniform two-triangle local patch with one normal mechanical mode, constant friction and fixed capacity. Saved-frame stress errors do not bound every wave peak. Fast GVS controls all reach the speed cap; they do not prove arbitrary scaling independence. Manual slow-to-fast restart does not establish automatic events, long tectonic clock handling, a spatial rupture or a full earthquake cycle.',
    reference='Independent differential Mohr–Coulomb/Jaumann and sealed fluid/energy solution in fault_heat.py; a tighter integrator tolerance changes normalized stress by less than 1.6×10⁻¹⁰. Thermal-pressurization motivation: Rice (2006), doi:10.1029/2005JB004006. Full derivation and controls are in the supporting note.')

SPECS['convection_3d'] = dict(
    title='Three-dimensional porous convection cells', kind='linear 3D convection',
    labels=[('Geometry', 'Unit cube, rigid saturated matrix and fixed tetrahedral mesh. All faces are impermeable; vertical faces are thermally insulated.'),
            ('Temperature', 'The base is at 283 K and top at 273 K. A 0.1 mK Gaussian perturbs the conductive profile; its projection onto the (1,1,1) mode measures growth or decay.'),
            ('Pressure storage', 'The fluid bulk modulus is deliberately scaled to 10 Pa for an economical explicit test. The analytical solution retains this finite storage; this is not a realistic water-property application.'),
            ('Comparison', 'Three meshes, three Rayleigh numbers, an unseeded rest case and a separate halved physical timestep. No material or seed fitting.')],
    equations=[r'\tau=\kappa_Tt/H^2,\quad \kappa_T=\lambda/C,\quad D=D_p/\kappa_T',
               r'a^2=\pi^2(m^2+n^2),\quad K_m^2=a^2+\ell^2\pi^2',
               r'Ra_c=K_m^4/a^2=9\pi^2/2\quad(m=n=\ell=1)',
               r'\delta T=A(\tau)\cos(\pi x/H)\cos(\pi y/H)\sin[\pi(z/H+1)]',
               r'\dot P=-DK_m^2P-D\pi A',
               r'\dot A=Ra\pi P+(Ra-K_m^2)A',
               r'b_1=(D+1)K_m^2-Ra,\quad b_0=D(K_m^4-Ra\,a^2)',
               r'\lambda_+=\frac{-b_1+\sqrt{b_1^2-4b_0}}{2}',
               r'\sigma=\kappa_T\lambda_+/H^2'],
    definitions=r'Dots mean derivatives with respect to $\tau$. $P$ is pressure-mode amplitude divided by $\rho_fg\alpha_fH$, with temperature units. $D_p=k/(\mu S_B)$, $S_B=\phi/K_f+(\alpha-\phi)/K_s+\alpha^2/K$. $Ra=\rho_f^2c_fg\alpha_f\Delta T kH/(\mu\lambda)$. The cube uses $\kappa_T=10^{-6}$ m²/s; $D$ varies with the prescribed permeability.',
    derivation=[('Linearize Darcy flow about the conductive hydrostatic state. The pressure disturbance has the same horizontal cosines and a vertical cosine. Its diffusion equation couples to the vertical derivative of the temperature mode, giving the first amplitude equation.', None),
                ('The background vertical temperature gradient converts the perturbed vertical Darcy flux into heat advection. Adding thermal diffusion gives the second amplitude equation.', None),
                ('Set the determinant of the two-mode system to zero for its eigenvalues. Its positive or least negative root gives the late growth rate. At zero growth the storage cancels, so onset is the same as in incompressible Darcy convection.', r'\lambda^2+b_1\lambda+b_0=0'),
                ('As pressure diffusion becomes infinitely fast, eliminate P to recover the incompressible result. The implementation checks this limit, zero-Rayleigh diffusion and the exact neutral point independently.', r'\lambda_+\longrightarrow Ra\,a^2/K_m^2-K_m^2\quad(D\to\infty)')],
    reference='Separation of the linear Darcy–Boussinesq pressure and heat equations, evaluated in convection_reference.py. The (1,1,1) mode extends the Horton–Rogers–Lapwood layer reference to a cube with finite pressure storage.',
    limit='Linear 3D cell growth and decay with deliberately scaled storage, fixed properties and a rigid matrix. The reference retains finite storage; it must not be replaced by the incompressible expression at D around 10. This does not establish nonlinear steady cells, realistic compressibility or mechanical coupling.')

SPECS['mass'] = dict(
    title='Finite-density water in a sealed cell',
    labels=[('Motion', 'Unit cell with fixed side-normal motion and base. The top shortens at 0.01 m/s for 0.1 s; the mesh stays fixed and the inventory uses the implied bulk volume.'),
            ('Water and temperature', 'Every face is sealed, initial pressure is zero, and temperature stays at 273 K. Gravity is zero.'),
            ('Material', 'Saturated elasticity with evolving Eulerian Biot porosity and a constant-compressibility fluid EOS; no plasticity or recovery in this selected control.')],
    equations=[r'\theta=T-T_0,\qquad V=V_0e^{\varepsilon_v}',
               r'\rho=\rho_0e^{p/K_f-\alpha_f\theta}',
               r'\phi=\alpha-(\alpha-\phi_0)e^{-\varepsilon_v-p/K_s+\alpha_s\theta+h}',
               r'M=V\rho\phi,\qquad M(p_*,T,\varepsilon_v)=M_0',
               r'L_n=\rho_nV_n\alpha^2/B,\qquad R_n=L_n(p_{n+1}-p_n)',
               r'M_{n+1}+R_n=M_n+R_{n-1}',
               r'B_{2D}=K+G/3,\qquad B_{3D}=K'],
    definitions=r'$\varepsilon_v$ is the tension-positive bulk strain interpreted as log volume in the inventory law; $h$ is recovered pore dilation (zero here). $\alpha=0.98$, $\phi_0=0.25$, $K_f=100$ MPa, $K_s=1$ GPa, $K=20$ MPa, $G=10$ MPa and $M_0=250$ kg per metre thickness. $R$ is an extensive numerical buffer, never physical mass.',
    derivation=[('Integrating constant fluid compressibility and expansivity gives the density law. Integrating the Eulerian porosity rate gives the displayed porosity. The sealed continuum reference is the monotone scalar root M=M0, evaluated independently in mass_reference.py.', None),
                ('When fluid and grain moduli are equal, the implicit pressure has a closed form. The unequal-modulus reference shown here uses a bracketed root instead.', r'p_*=K_f\left\{\ln\frac{\phi_0e^{\alpha_f\theta}+(\alpha-\phi_0)e^{\alpha_s\theta+h}}{\alpha}-\varepsilon_v\right\}'),
                ('The fixed-stress pressure stabilizer leaves a separately recorded numerical mass increment. With sealed walls, physical M differs from M0 by minus R. Refining the physical timestep reduces this defect at first order; an exact M+R ledger alone does not establish exact physical mass conservation.', None)],
    reference='Direct EOS and porosity integration; independent scalar root and step-refinement controls in mass_reference.py and check_mass.py.',
    limit='Selected homogeneous isothermal prescribed-motion cell. The broader suite checks remapping and accepted-path spatial transport, but this figure does not establish finite-strain constitutive accuracy, unsaturated finite-density flow, buoyant hydrostatic balance or a full thermal energy equation.')

SPECS['drained_heat'] = dict(
    title='Heating a drained cell with changing water mass',
    labels=[('Geometry and constraints', 'Rigid unit cell, every face at p=0, no gravity or solid thermal expansion. All thermal faces are insulated.'),
            ('Heat source', 'Uniform supplied power 30 MW per metre thickness for 5 s; no conduction or fluid point injection.'),
            ('Exchange', 'Expanded water leaves the pressure-controlled boundary at the current cell temperature. The original solid inventory and its heat capacity stay fixed.')],
    equations=[r'\theta=T-T_0,\qquad M=M_0e^{-\alpha_f\theta}',
               r'C_s=(1-\phi_0)\rho_sc_sV_0',
               r'E=(C_s+c_fM)T,\qquad dM_o=-dM',
               r'dQ=dE+c_fT\,dM_o=(C_s+c_fM)\,dT',
               r'Q=C_s\theta+\frac{c_fM_0}{\alpha_f}(1-e^{-\alpha_f\theta})',
               r'H_o=Q-[E(T)-E(T_0)]'],
    definitions=r'$M_o$ and $H_o$ are cumulative exported water mass and sensible heat. $T_0=273$ K, $M_0=250$ kg, $C_s=1.9875$ MJ/K, $c_f=4180$ J/(kg K) and $\alpha_f=2\times10^{-4}$ K$^{-1}$. Values refer to unit volume per metre thickness.',
    derivation=[('At fixed zero pressure, zero bulk strain and zero solid expansion, porosity stays fixed. Exponential thermal expansion gives the retained fluid mass directly.', None),
                ('Differentiate retained sensible energy. Its c_f T dM term cancels the outflow enthalpy, leaving a temperature-dependent heat capacity. Integrating from the initial temperature gives Q(T); invert Q=power times time for the temperature history.', None),
                ('The independent reference accounts for water leaving at the changing cell temperature. Subtracting retained energy gain from supplied heat gives its cumulative export. No pressure or strain boundary work occurs in this control.', None)],
    reference='Direct integration of the declared exponential-density sensible-heat model in mass_reference.py.',
    limit='A rigid zero-pressure open-system control. Sensible heat is measured from the stated temperature origin. Pressure work, thermoelastic caloric energy, phase changes and recovery energy are outside this model; spatial thermal-plume accuracy is not established here.')

SPECS['mixing'] = dict(
    title='Warm-water mixing in a drained cell',
    labels=[('Cell', 'Rigid unit cell at p=0, initially 273 K, insulated and without gravity. Uniform distributed sources keep the temperature spatially uniform.'),
            ('Injection', '0.01 m³/s enters at 323 K for 5 s. Its mass rate uses density at the inlet temperature.'),
            ('Outflow', 'Natural boundary reaction carries water at the changing cell temperature; it balances injection and the loss of retained mass from thermal expansion.')],
    equations=[r'\dot M_i=\rho(0,T_i)Q_i,\qquad M(T)=M_0e^{-\alpha_f(T-T_0)}',
               r'\dot M_o=\dot M_i-\dot M',
               r'E=(C_s+c_fM)T',
               r'(C_s+c_fM)\dot T=c_f\dot M_i(T_i-T)',
               r't=\int_{T_0}^{T}\frac{C_s+c_fM(u)}{c_f\dot M_i(T_i-u)}\,du',
               r'E-E_0=\int c_fT_i\,dM_i-\int c_fT\,dM_o'],
    definitions=r'$Q_i$ is inlet volume rate, $\dot M_i$ inlet mass rate and $\dot M_o$ the outward boundary mass rate. $C_s=1.9875$ MJ/K, $c_f=4180$ J/(kg K), $M_0=250$ kg and $\alpha_f=2\times10^{-4}$ K$^{-1}$. The inlet density is lower than the initial cell density.',
    derivation=[('Water balance sets the outward mass rate. In the sensible-energy balance, the retained c_f T dM term cancels the matching part of boundary export, giving the displayed temperature equation.', None),
                ('Separate variables to obtain the independent integral solution. Its integrand is evaluated after a logarithmic distance-to-inlet-temperature substitution, avoiding an endpoint pole. With zero expansivity the solution reduces to an exponential.', r'T=T_i-(T_i-T_0)e^{-c_f\dot M_it/(C_s+c_fM_0)}'),
                ('The reported budget uses actual inlet density, actual retained fluid mass and boundary enthalpy at cell temperature. It includes no numerical fixed-stress buffer in the heat capacity.', None)],
    reference='Independent integral solution and constant-density exponential limit in mass_reference.py.',
    limit='Uniform well-mixed sensible heat at zero pressure and zero strain. This control does not resolve a local well singularity, a thermal front or a buoyant plume.')

SPECS['elder'] = dict(
    numerical_reference=True,
    title="Elder's transient heater problem",
    labels=[('Geometry', 'Symmetric half of a 4H-wide porous layer: 0 ≤ x ≤ 2H, 0 ≤ z ≤ H, H = 1 m. The matrix is rigid; all walls are impermeable.'),
            ('Temperature', 'Initially T = T0 = 273 K. At t = 0 the base segment 0 ≤ x ≤ H is heated to T0 + 10 K; the top stays at T0. The rest of the base and both sidewalls are insulated.'),
            ('Numerics', 'Ra = 400, fixed 50 s physical steps, 0 ≤ τ ≤ 0.1. The displayed DES mesh has h = H/64. An independent incompressible reference uses 256 × 512 cells; DES retains finite pressure storage with fluid bulk modulus 100 Pa.')],
    equations=[r'\theta=(T-T_0)/\Delta T,\quad \tau=\kappa t/(\sigma H^2)',
               r'\nabla^2\psi=-Ra\,\partial_x\theta,\quad u=\partial_z\psi,\quad w=-\partial_x\psi',
               r'\partial_\tau\theta+u\partial_x\theta+w\partial_z\theta=\nabla^2\theta',
               r'\mathcal{H}(\tau)=\int_0^1\int_0^2\theta(x,z,\tau)\,dx\,dz'],
    definitions=r'Coordinates, velocity and stream function in the reference are dimensionless. $\kappa=10^{-6}$ m²/s, $\sigma=1$ and $Ra=400$. $\psi=0$ on the walls imposes no normal flow. $\mathcal{H}$ is scaled heat content, not the complete dimensional energy ledger.',
    derivation=[('Take the curl of incompressible Boussinesq Darcy flow to eliminate pressure. Introduce a stream function to satisfy continuity exactly, giving the Poisson equation above.', None),
                ('Scale heat transport by H, ΔT and σH²/κ. Conduction, advection and the heater boundary condition then give the displayed temperature equation. There is no closed-form solution for this transient.', None),
                ('elder_reference.py solves the stream-function Poisson equation by sparse LU and advances heat with third-order SSP Runge–Kutta and independent stable steps. Central spatial differences are checked on four grids. The 128 and 256 grids approach the same plume branch; coarser grids do not.', None)],
    reference='Elder (1967), Journal of Fluid Mechanics 27, 609–623; Voss and Souza (1987), Water Resources Research 23, 1851–1866; van Reeuwijk et al. (2009), Water Resources Research 45, W04416. Independent reference: elder_reference.py.',
    limit='Symmetric thermal half box through τ = 0.1. Early plume errors reach about 14% in heat content and 0.56 in scaled pointwise temperature. DES uses deliberately scaled finite storage; the reference is incompressible. This comparison does not establish grid-independent nonlinear DES plumes, an exact energy budget, saline transport or later asymmetric branches.')

ORDER = ('theis', 'hantush', 'theis_wells', 'open_flow', 'anisotropic_time', 'point_sources', 'van_genuchten', 'brooks_corey', 'gardner', 'bishop', 'bishop_cell', 'drainage', 'content', 'mass', 'terzaghi', 'mandel',
         'cryer', 'averaged_aquifer', 'regional', 'embankment', 'drawdown', 'slope', 'rainfall', 'rudnicki', 'undrained_wave', 'strength', 'rsf', 'plastic_heat', 'fault_heat', 'plastic_content', 'pressure_yield', 'kirsch', 'salencon', 'borehole', 'ogata_banks', 'thermal_cell', 'drained_heat', 'mixing', 'mctigue', 'convection', 'convection_3d', 'nusselt', 'elder', 'injection_temperature')

SLUGS = dict(elder='elder', mass='mass', drained_heat='drained-heat', mixing='mixing', fault_heat='fault-heat', plastic_heat='plastic-heat', content='content', thermal_cell='thermal-cell', plastic_content='plastic-content', terzaghi='terzaghi', gardner='gardner', van_genuchten='vg',
             theis_wells='theis-wells', open_flow='open-flow', anisotropic_time='anisotropic-time', strength='strength', kirsch='kirsch',
             brooks_corey='bc', undrained_wave='wave', ogata_banks='heat',
             point_sources='sources', theis='theis', averaged_aquifer='aquifer', convection='convection', convection_3d='convection-3d', mandel='mandel', cryer='cryer', rudnicki='rudnicki', bishop='bishop', hantush='hantush', mctigue='mctigue', bishop_cell='bishop-cell', nusselt='nusselt', regional='regional', pressure_yield='pressure-yield', drainage='drainage', injection_temperature='injection-temperature', borehole='borehole', rsf='rsf', embankment='embankment', salencon='salencon', drawdown='drawdown', slope='slope', rainfall='rainfall')

STYLE = """
:root{color-scheme:light;--ink:#202b33;--muted:#52616b;--blue:#22658a}
*{box-sizing:border-box}body{font:16px/1.65 system-ui,sans-serif;color:var(--ink);
max-width:1000px;margin:auto;padding:32px 34px 64px;background:white}
a{color:var(--blue);text-underline-offset:3px}h1{font-size:30px;line-height:1.25;
font-weight:650;margin:20px 0}h2{font-size:21px;margin:0 0 16px;font-weight:600}
h3{font-size:16px}p{margin:12px 0}nav{display:flex;gap:18px;flex-wrap:wrap;
font-size:14px;border-bottom:1px solid #d8dfe3;padding-bottom:15px}
.kicker{font-size:12px;letter-spacing:.1em;color:var(--muted);text-transform:uppercase}
section{margin:36px 0;padding-top:22px;border-top:1px solid #d8dfe3}
figure{margin:18px 0}img{max-width:100%;height:auto}img.setup{width:100%;max-width:500px}
img.result{width:100%}.setup-grid{display:grid;grid-template-columns:1fr 1fr;
gap:30px;align-items:start}figcaption,small{font-size:13px;color:var(--muted)}
dl{margin:0}dt{font-weight:600;margin-top:15px}dd{margin:4px 0 16px}
.equation{display:flex;align-items:center;gap:16px;margin:14px 0;overflow-x:auto}
.equation img{display:block;flex-shrink:0;width:auto;max-width:none}
.eqno{font-size:12px;color:var(--muted);margin-left:auto}.definitions{font-size:14px;
color:var(--muted)}.note{border-left:3px solid #6e99af;padding:4px 0 4px 15px}
table{width:100%;border-collapse:collapse}th,td{padding:14px 12px;border-bottom:1px solid
#d8dfe3;text-align:left;vertical-align:top}th{font-size:13px;color:var(--muted)}
td:first-child{width:34%}details{margin:18px 0}summary{cursor:pointer;color:var(--blue)}
pre{font-size:12px;overflow:auto;padding:14px;background:#f4f6f7;white-space:pre-wrap}
td code{word-break:break-word}h3{margin-top:24px}
footer{margin-top:40px;padding-top:16px;border-top:1px solid #d8dfe3;font-size:13px;color:var(--muted)}
@media(max-width:700px){body{padding:20px 18px}.setup-grid{grid-template-columns:1fr}
h1{font-size:26px}th,td{padding:10px 6px}}
@media print{nav,details{display:none}body{max-width:none;padding:0}section{break-inside:auto}
h2,h3,dt{break-after:avoid}figure,.equation{break-inside:avoid}
a{color:inherit;text-decoration:none}}
"""


def page_head(title, css='style.css'):
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{html.escape(title)} · DES</title><link rel="stylesheet" href="{css}">'
            '</head><body>')


INTRO = ('Worked examples of coupled fluid flow, heat and deformation in DES. Each page sets up one problem, gives the '
         'input file and the commands to run it, says which output to look at, and compares the result with an analytical '
         'or independent reference solution, with the limits of that comparison.')
HOW = ('Every example ships its input file and any mesh or restart files it needs, so it runs as is with DES '
       'built from this source.')

GLOSSARY = ('Terms: DES is DynEarthSol. PT is its pseudo-transient solver, which iterates each time step to mechanical '
            'equilibrium by dynamic relaxation; a physical step is one model time step, as opposed to a PT iteration. '
            'The explicit form is ordinary explicit mechanics with inertial (mass) scaling, and GVS is that form with '
            'global velocity scaling. Pressure feedback means that deformation changes the pore pressure, so flow and '
            'mechanics are coupled both ways. EP is elasto-plastic Mohr–Coulomb plasticity, and EP-RSF the same with '
            'rate-and-state friction.')


def reference_lines(key):
    """Separate definitions without shrinking the equations to fit a wide plate."""
    equations = list(SPECS[key]['equations'])
    if key == 'gardner':
        equations[2:3] = [
            r'F_n=\frac{\sin(\lambda_nZ)\sin(\lambda_nL)e^{-\lambda_n^2T}}{1+L/2+2\lambda_n^2L}',
            r'w=w_b-4(b-a)e^{(L-Z)/2-T/4}\sum_{n=1}^{\infty}F_n']
    if key == 'ogata_banks':
        equations[2:3] = [
            r'a=(x-v_Tt)/(2\sqrt{Dt}),\quad b=(x+v_Tt)/(2\sqrt{Dt})',
            r'\theta=\frac{1}{2}[\operatorname{erfc}(a)+e^{v_Tx/D}\operatorname{erfc}(b)]']
    if key == 'theis':
        equations[-1:] = [
            r's_p(\mathbf{x},t)=\sum_j\sum_k\frac{\Delta Q_{jk}}{4\pi\mathcal{T}_p}E_1(u_{jk})',
            r'u_{jk}=|\mathbf{x}-\mathbf{x}_j|^2S_p/[4\mathcal{T}_p(t-t_{jk})],\quad t>t_{jk}']
    out = []
    for equation in equations:
        # Split separate equations (\qquad) onto their own lines; a \quad keeps a condition with its equation.
        out.extend(part.strip(' ,;') for part in re.split(r'\\qquad', equation) if part.strip())
    return out


def equation_figure(equation, folder, name):
    fig = pf.plt.figure(figsize=(7, .65), facecolor='white')
    item = fig.text(.01, .5, '$'+equation+'$', fontsize=12, va='center', color='#202b33')
    fig.canvas.draw()
    width = item.get_window_extent().width / fig.dpi
    assert width <= 7.1, (name, width, equation)
    fig.savefig(folder/(name+'.svg'), bbox_inches='tight', pad_inches=.04)
    pf.plt.close(fig)


YEAR = 365.2422*86400                     # DES's year, as in Dynearthsol.py
_OPTIONS = {}


def load_catalogues(folder):
    """Use qualified remote help without executing a local solver."""
    for dimension in (2, 3):
        text = (folder/f'des-help{dimension}d.txt').read_text()
        options = des_options.parse_help(text)
        assert 'control.hydraulic_content_model' in options
        _OPTIONS[dimension] = options


def des_catalogue(dimension):
    """DES's options with defaults and help, from the executable built in this checkout."""
    if dimension not in _OPTIONS:
        _OPTIONS[dimension] = des_options.catalogue(pf.ROOT/f'dynearthsol{dimension}d')
    return _OPTIONS[dimension]


def _shown(key, value):
    """A setting's value, with times in years also given in seconds."""
    if key.endswith('_in_yr'):
        numbers = des_options._numbers(value)
        if numbers and all(0 <= x < 1e8 for x in numbers) and any(numbers):
            seconds = [f'{x*YEAR:.6g}' for x in numbers]
            return f'{value} yr = ' + (f'{seconds[0]} s' if len(seconds) == 1 else '[' + ', '.join(seconds) + '] s')
    return value


def _cell(text):
    return html.escape(text, quote=False).replace('|', '&#124;')


def _record_text(key, rec, field):
    """Record text with this page's manual wording applied (exact replacements, so computed numbers stay live)."""
    text = rec[field]
    for old, new in MANUAL.get(key, {}).get('replace', {}).get(field, ()):
        assert old in text, (key, field, old)
        text = text.replace(old, new)
    return text


def manual_sections(key, slug, rec, output):
    """HTML and Markdown for 'Input and how to run' and 'Output'."""
    m = MANUAL[key]
    folder = output/'data'/slug
    commands = ' '.join(command for command, _ in m.get('run', ()))
    dimension = 3 if 'dynearthsol3d' in commands else 2 if 'dynearthsol2d' in commands else int(rec.get('dimension') or 2)
    exe = f'dynearthsol{dimension}d'
    settings = des_options.read(folder/'input.cfg')
    changed, unknown = des_options.changed(folder/'input.cfg', des_catalogue(dimension))
    assert not unknown, (key, unknown)
    catalogue = des_catalogue(dimension)
    model = settings.get('sim.modelname', 'result')
    run = m.get('run') or [(f'{exe} input.cfg', 'Runs the example.')]
    script = '\n'.join(f'# {note}\n{command}' if note else command for command, note in run)
    keys = []
    for item in m['keys']:
        name, note = item if isinstance(item, (tuple, list)) else (item, catalogue[item][2])
        assert name in settings, (key, name)
        keys.append((name, _shown(name, settings[name]), note))
    html_parts = [
        '<section id="input"><h2>2. Input and how to run</h2>',
        f'<p>Every file this example needs is in <a href="data/{slug}/">data/{slug}/</a>. Run from that folder:</p>',
        '<pre>'+html.escape(script)+'</pre>',
        '<h3>Settings that define this example</h3><table><thead><tr><th>Setting</th><th>Value</th><th>What it does here</th></tr></thead><tbody>']
    html_parts += [f'<tr><td><code>{html.escape(n)}</code></td><td><code>{html.escape(v)}</code></td><td>{html.escape(t)}</td></tr>' for n, v, t in keys]
    html_parts += ['</tbody></table>',
        f'<details><summary>All {len(changed)} settings of input.cfg that differ from DES’s defaults</summary>',
        '<table><thead><tr><th>Setting</th><th>Value (default)</th><th>DES help</th></tr></thead><tbody>']
    html_parts += [f'<tr><td><code>{html.escape(n)}</code></td><td><code>{html.escape(_shown(n, v))}</code> ({html.escape(str(d))})</td>'
                   f'<td>{html.escape(h)}</td></tr>' for n, v, d, h in changed]
    html_parts += ['</tbody></table></details>',
        f'<p><a href="data/{slug}/input.cfg">input.cfg</a> · <code>{exe} -h</code> lists every option with its default.</p></section>']
    md = ['## Input and how to run', f'Every file this example needs is in [data/{slug}/](data/{slug}/). Run from that folder:',
          '```bash\n'+script+'\n```', '### Settings that define this example',
          '| Setting | Value | What it does here |\n| --- | --- | --- |\n'
          + '\n'.join(f'| `{n}` | `{_cell(v)}` | {_cell(t)} |' for n, v, t in keys),
          f'[input.cfg](data/{slug}/input.cfg) holds all {len(changed)} settings that differ from DES’s defaults; '
          f'`{exe} -h` lists every option with its default.']

    fields = m['fields']
    first = fields[0][0]
    snippet = (f"from Dynearthsol import Dynearthsol      # in the DES source folder\n"
               f"model = Dynearthsol('{model}')\n"
               f"frame = model.frames[-1]                 # model.time[i] is the time of frame i, in seconds\n"
               f"x = model.read_field(frame, 'coordinate')\n"
               f"value = model.read_field(frame, '{first}')")
    restarted = settings.get('sim.is_restarting', 'no').lower() in des_options.TRUE
    vtk = f'python3 "$DES_SOURCE/2vtk.py" {"-g " if restarted else ""}{model}'
    files = (f'DES writes <code>{model}.info</code> (one line per saved frame: frame, step, time in seconds), '
             f'<code>{model}.save.NNNNNN</code> (the fields of each saved frame) and <code>{model}.chkpt.NNNNNN</code> '
             '(checkpoints to restart from). <code>sim.output_step_interval</code> and '
             '<code>sim.output_time_interval_in_yr</code> set how often a frame is saved.')
    html_parts += ['<section id="output"><h2>3. Output</h2>', '<p>'+files+'</p>',
        '<p>Set <code>DES_SOURCE</code> and <code>PYTHONPATH</code> as described in '
        '<a href="REPRODUCE.md">Running the examples</a> before using the Python reader or converter.</p>',
        '<table><thead><tr><th>Field</th><th>What to look at</th></tr></thead><tbody>']
    html_parts += [f'<tr><td><code>{html.escape(f)}</code></td><td>{html.escape(t)}</td></tr>' for f, t in fields]
    html_parts += ['</tbody></table><p>Read a frame in Python:</p><pre>'+html.escape(snippet)+'</pre>',
        f'<p>View it in ParaView after <code>{html.escape(vtk)}</code>'
        + (' (<code>-g</code> also converts the frames of the run it restarted from).' if restarted else '.') + '</p></section>']
    md += ['## Output', re.sub(r'</?code>', '`', files),
           'Set `DES_SOURCE` and `PYTHONPATH` as described in [Running the examples](REPRODUCE.md) before using the Python reader or converter.',
           '| Field | What to look at |\n| --- | --- |\n' + '\n'.join(f'| `{f}` | {_cell(t)} |' for f, t in fields),
           'Read a frame in Python:', '```python\n'+snippet+'\n```', f'View it in ParaView after `{vtk}`.']
    return html_parts, md


def write_review(records, output):
    """One problem per page, from one shared source for HTML and Markdown."""
    (output/'style.css').write_text(STYLE)
    all_records=[*records,*[dict(key=key,reference_only=True) for key in REFERENCE_ONLY]]
    setup_owners = {}
    for rec in all_records:
        key=rec['key'];s=SPECS[key];slug=SLUGS[key]
        setup_slug = setup_owners.setdefault(DRAW[key], slug)
        setup = f'fig/{setup_slug}/setup'
        assets=output/'fig'/slug
        if assets.exists():shutil.rmtree(assets)
        assets.mkdir(parents=True)
        numerical=not rec.get('reference_only',False)
        if numerical:
            for ext in ('svg','png','pdf'):
                shutil.move(str(output/f'figures/{key}.{ext}'),assets/f'result.{ext}')
            source=output/'data'/key;target=output/'data'/slug
            if source!=target:
                if target.exists():shutil.rmtree(target)
                source.rename(target)
            rec.update(page=slug+'.html',schematic=setup+'.svg',
                       analytical_solution=f'fig/{slug}/e1.svg',derivation=slug+'.html#derivation')
        if setup_slug == slug:
            schematic(key,assets)
        equations=reference_lines(key)
        for i,equation in enumerate(equations,1):equation_figure(equation,assets,f'e{i}')
        rec['equation_files']=[f'fig/{slug}/e{i}.svg' for i in range(1,len(equations)+1)]
        title=rec.get('title',s['title'])
        manual=numerical and key in MANUAL
        parts=[page_head(title),'<nav><a href="index.html">All problems</a><a href="#setup">Setup</a>',
               '<a href="#input">Input</a><a href="#output">Output</a>' if manual else '',
               '<a href="#result">Result</a>' if numerical else '',
               '<a href="#solution">Reference</a><a href="#scope">Limits</a></nav>',
               '<p class="kicker">DES · '+('Application example' if rec.get('status') == 'preliminary' else 'Worked example' if numerical else 'Analytical reference')+'</p>',
               '<h1>'+html.escape(title)+'</h1>']
        md=['---\ntitle: '+json.dumps(title, ensure_ascii=False)+'\n---','[All problems](README.md)','# '+title]
        if numerical:
            problem,result=_record_text(key,rec,'problem'),_record_text(key,rec,'result')
            parts+=['<p>'+html.escape(problem)+'</p>','<p class="note">'+html.escape(result)+'</p>']
            md+=[problem,'**'+result+'**']
        else:
            note='Analytical reference. Full DES radial-flow or subsidence qualification remains unfinished.'
            parts+=['<p class="note">'+note+'</p>'];md+=[note]
        parts+=['<section id="setup"><h2>1. Problem setup</h2><div class="setup-grid">',
            f'<figure><img class="setup" src="{setup}.svg" alt="{html.escape(title)}: geometry and boundary conditions">',
            f'<figcaption>Geometry and imposed conditions. <a href="{setup}.pdf">Vector PDF</a></figcaption></figure><dl>']
        md+=['## Problem setup',f'![Setup]({setup}.png)']
        for label,description in s['labels']:
            parts+=['<dt>'+html.escape(label)+'</dt><dd>'+html.escape(description)+'</dd>']
            md+=['**'+label+'.** '+description]
        parts+=['</dl></div></section>']
        if manual:
            more_html,more_md=manual_sections(key,slug,rec,output);parts+=more_html;md+=more_md
        numerical_reference=s.get('numerical_reference',False)
        number=4 if manual else 2
        if numerical:
            results='Result: DES and the reference' if numerical_reference else 'Result: DES and the analytical solution'
            caption=MANUAL.get(key,{}).get('caption') or rec['scope']
            meaning=_record_text(key,rec,'meaning')
            parts+=[f'<section id="result"><h2>{number}. {results}</h2>',
                f'<figure><img class="result" src="fig/{slug}/result.svg" alt="DES compared with the independent reference">',
                '<figcaption>'+html.escape(caption)+'</figcaption></figure>',
                '<p>'+html.escape(meaning)+'</p>',
                f'<p><a href="fig/{slug}/result.pdf">PDF</a> · <a href="fig/{slug}/result.png">PNG</a></p>']
            md+=['## '+results,f'![Result](fig/{slug}/result.png)',caption,meaning]
            if rec.get('supplements'):
                parts += ['<p>'+' · '.join(f'<a href="appendix/{name}">{html.escape(label)}</a>'
                                         for label, name in rec['supplements'])+'</p>']
                md += [' · '.join(f'[{label}](appendix/{name})' for label, name in rec['supplements'])]
            if rec.get('tables'):
                links = [(label, f'data/{slug}/{name}') for label, name in rec['tables']]
                parts += ['<p>'+' · '.join(f'<a href="{path}">{html.escape(label)}</a>' for label,path in links)+'</p>']
                md += [' · '.join(f'[{label}]({path})' for label,path in links)]
            (output/'data'/slug/'record.json').write_text(json.dumps(dict(problem=rec['problem'],result=rec['result'],scope=rec['scope'],
                meaning=rec['meaning'],metric=rec.get('metric'),metric_label=rec.get('metric_label'),spacing_m=rec.get('spacing_m'),
                evidence=json.loads(rec['evidence']) if isinstance(rec.get('evidence'),str) and rec['evidence'].lstrip()[:1] in '{[' else rec.get('evidence')),
                indent=1,ensure_ascii=False)+'\n')
            study=MANUAL.get(key,{}).get('study')
            if study:
                parts+=['<details><summary>Reproduce the whole study behind this figure</summary><p>From the DES source folder, with the '
                        'executable built there (each runner writes a new output folder and refuses an existing one):</p>',
                        '<pre>'+html.escape('\n'.join(study))+'</pre>',
                        f'<p>The study’s full numbers: <a href="data/{slug}/record.json">record.json</a>.</p></details>']
                md+=['To reproduce the whole study behind this figure, from the DES source folder:','```bash\n'+'\n'.join(study)+'\n```',
                     f'The study’s full numbers: [record.json](data/{slug}/record.json).']
            parts+=['</section>'];number+=1
        heading='Reference solution (numerical)' if numerical_reference else 'Analytical reference'
        parts+=[f'<section id="solution"><h2>{number}. {heading}</h2>']
        md+=['## '+heading]
        for i,equation in enumerate(equations,1):
            parts+=[f'<div class="equation"><img src="fig/{slug}/e{i}.svg" alt="{html.escape(equation)}"><span class="eqno">({i})</span></div>']
            md+=['$$\n'+equation+'\n$$']
        parts+=['<p class="definitions">'+html.escape(s['definitions'])+'</p>']
        md+=[s['definitions']]
        parts+=['<details id="derivation"><summary>Derivation</summary>']
        md+=['### Derivation']
        for i,(description,equation) in enumerate(s['derivation'],1):
            parts+=['<p>'+html.escape(description)+'</p>'];md+=[description]
            if equation:
                for j,line in enumerate(re.split(r'\\qquad',equation),1):
                    line=line.strip(' ,;');name=f'd{i}-{j}';equation_figure(line,assets,name)
                    parts+=[f'<div class="equation"><img src="fig/{slug}/{name}.svg" alt="{html.escape(line)}"></div>']
                    md+=['$$\n'+line+'\n$$']
        parts+=['</details>','<p>'+html.escape(s['reference'])+'</p>'];md+=[s['reference']]
        if 'url' in s:
            parts+=[f'<p><a href="{html.escape(s["url"])}">Original reference</a></p>'];md+=[f'[Original reference]({s["url"]}).']
        if key=='averaged_aquifer':
            parts+=['<p><a href="appendix/averaged-aquifer.md">Full finite-annulus theory note</a></p>'];md+=['[Full finite-annulus theory note](appendix/averaged-aquifer.md).']
        parts+=['</section>','<section id="scope"><h2>Validity and limits</h2><p>'+html.escape(s['limit'])+'</p></section>']
        md+=['## Validity and limits',s['limit']]
        parts+=['<footer><a href="index.html">All problems</a> · <a href="'+slug+'.md">Portable Markdown</a></footer></body></html>']
        (output/(slug+'.html')).write_text('\n'.join(parts)+'\n')
        (output/(slug+'.md')).write_text('\n\n'.join(md)+'\n')
    if not any((output/'figures').iterdir()):(output/'figures').rmdir()
    title='Hydromechanics verification'
    parts=[page_head(title),'<p class="kicker">DES · worked examples</p>','<h1>'+title+'</h1>',
           '<p>'+html.escape(INTRO)+'</p>','<p>'+html.escape(HOW)+' <a href="REPRODUCE.md">Running the examples</a>.</p>',
           '<p class="note">'+html.escape(GLOSSARY)+'</p>',
           '<h2>Examples</h2>']
    md=['# '+title,INTRO,HOW+' See [Running the examples](REPRODUCE.md).',
        GLOSSARY,'Pressure is compression-positive; stress and volumetric strain are tension/expansion-positive.','## Examples']
    for theme, heading in THEMES:
        group=[rec for rec in records if rec.get('theme')==theme]
        if not group: continue
        parts+=[f'<h3>{html.escape(heading)}</h3><table><thead><tr><th>Problem</th><th>Result</th></tr></thead><tbody>']
        md+=[f'### {heading}']
        for rec in group:
            slug=SLUGS[rec['key']]
            parts+=[f'<tr id="{rec["key"]}"><td><a href="{slug}.html">{html.escape(rec["title"])}</a></td><td>{html.escape(rec["result"])}</td></tr>']
            md+=[f'- [{rec["title"]}]({slug}.md).']
        parts+=['</tbody></table>']
    assert {rec.get('theme') for rec in records} <= {t for t, _ in THEMES}, 'every record needs a theme'
    if REFERENCE_ONLY:
        parts+=['<section><h2>Analytical references</h2><p>These references have standalone checks but no DES comparison yet.</p><ul>']
        md+=['## Analytical references','Standalone references without a DES comparison yet.']
        for key in REFERENCE_ONLY:
            slug=SLUGS[key];label=SPECS[key]['title']
            parts+=[f'<li><a href="{slug}.html">{html.escape(label)}</a></li>'];md+=[f'- [{label}]({slug}.md).']
        parts+=['</ul></section>']
    parts+=['<section><h2>Background</h2><p><a href="appendix/physical-inertia.md">Physical-inertia notes</a> · <a href="evidence.html">Verification records</a> · <a href="provenance.json">Source and data versions</a></p></section>',
            '<footer>Pressure is compression-positive; stress and volumetric strain are tension/expansion-positive. <a href="README.md">Markdown index</a></footer></body></html>']
    md+=['## Background','[Physical-inertia notes](appendix/physical-inertia.md). [Verification records](evidence.html). [Source and data versions](provenance.json).']
    (output/'index.html').write_text('\n'.join(parts)+'\n');(output/'README.md').write_text('\n\n'.join(md)+'\n')
    appendix=[page_head('Reference solutions','../style.css'),'<h1>Reference solutions</h1><p>Each problem has its own page.</p>']
    appendix_md=['# Reference solutions']
    for rec in all_records:
        key=rec['key'];slug=SLUGS[key]
        appendix+=[f'<p id="{key}"><a href="../{slug}.html#solution">{html.escape(SPECS[key]["title"])}</a></p>']
        appendix_md+=[f'<a id="{key.replace("_","-")}"></a>',f'[{SPECS[key]["title"]}](../{slug}.md).']
    (output/'appendix/analytical-solutions.html').write_text('\n'.join(appendix)+'\n</body></html>\n')
    (output/'appendix/analytical-solutions.md').write_text('\n\n'.join(appendix_md)+'\n')
