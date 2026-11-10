"""Per-example manual content for the review pages: how to run, the settings that define it, the output to look at.

MANUAL[key] = dict(
    run=[(command, note)],        # shell commands, in order, from inside data/<page>/
    keys=[(setting, note)],       # settings of input.cfg that define the example; values are read from the file
    fields=[(field, note)],       # DES output fields the figure is made from
    caption=...,                  # what each figure panel shows
    study=[command],              # repository commands that reproduce the whole study behind the figure
    replace={field: [(old, new)]})  # plain wording for a record's problem/result/meaning text
"""

MANUAL = {'theis': {'run': [('dynearthsol3d input.cfg',
                    'Runs the 3D quarter box: 80 steps of 2 s, the well pumping for the first 80 s and shut in for '
                    'the last 80 s. DES writes theis.save.* frames every 10 steps (20 s); frame 0 is the initial '
                    'state.')],
           'keys': [('control.has_hydraulic_diffusion', 'turns on Darcy flow of pore pressure; wells need it'),
                    ('control.has_poroelastic_pressure_feedback',
                     'off: deformation does not change the pressure, so the storage stays at the fixed-stress value '
                     'S_B and pressure obeys the Theis diffusion equation'),
                    ('control.gravity',
                     'no gravity, so there is no hydrostatic gradient and the drawdown is compared directly in '
                     'pressure'),
                    ('control.fixed_dt', 'physical step of 2 s; DES splits it into stable flow substeps by itself'),
                    ('sim.max_steps', '80 steps of 2 s end the run at 160 s: 80 s of pumping and 80 s of recovery'),
                    ('sim.output_step_interval',
                     'a frame every 10 steps (20 s), so eight frames are compared with Theis'),
                    ('bc.num_wells', 'one rate-controlled vertical well'),
                    ('bc.well_x',
                     'with bc.well_y, puts the well on the corner edge x = y = 0; the sealed symmetry planes make it '
                     'one quarter of a full well'),
                    ('bc.well_rate',
                     'extraction rate of the quarter well, Q/4 with Q = 1e-4 m³/s, spread uniformly over the 5 m '
                     'screen'),
                    ('bc.well_schedule_time_in_yr',
                     'schedule times 0, 80 s, 80 s (in years); the repeated time makes a step change at 80 s'),
                    ('bc.well_schedule_ratio', 'rate multipliers 1, 1, 0: full rate until 80 s, then shut in'),
                    ('mat.bulk_modulus',
                     'very stiff skeleton (1e12 Pa), so the mechanical storage α²/K is negligible and the solid '
                     'barely moves'),
                    ('mat.porosity', 'with the fluid bulk modulus, sets the pressure storage φ/K_f of the aquifer'),
                    ('mat.hydraulic_perm', 'chosen so that the pressure diffusivity κ/(μS_B) is 1 m²/s')],
           'fields': [('pore pressure',
                       'drawdown s_p = p(frame 0) − p at every node, compared with Theis at nodes 10–30 m from the '
                       'well'),
                      ('coordinate', 'gives the radius r = √(x²+y²) of each node from the well axis'),
                      ('connectivity',
                       'with the coordinates, gives element volumes for the closed-box water budget against the '
                       'pumped volume')],
           'caption': '(a) Radial drawdown profiles s_p/C_Q on the y = 0 symmetry plane at four saved times, two '
                      'during pumping and two during recovery: DES nodes (circles) and the Theis solution with the '
                      'shut-in superposed (lines). The shaded band is the compared range 10 ≤ r ≤ 30 m. (b) Drawdown '
                      'history at r ≈ 10 m and r ≈ 20 m through pumping and recovery; the dotted line marks shut-in. '
                      '(c) Maximum and RMS drawdown error over the compared nodes and times against mesh spacing, '
                      'with a slope-2 guide.',
           'study': ['python3 tests/functional/hydromechanics/check_theis.py dynearthsol3d OUT --spacings 5 2.5 1.25 '
                     '0.625']},
 'hantush': {'run': [('dynearthsol3d input.cfg',
                      'Runs the 3D quarter box with the one-element aquitard on top: 160 steps of 2 s, pumping for '
                      '200 s, then recovery to 320 s. DES writes hantush.save.* frames every 20 steps (40 s). On one '
                      'thread this mesh takes about 45 minutes.')],
             'keys': [('control.has_hydraulic_diffusion', 'turns on Darcy flow of pore pressure; wells need it'),
                      ('control.has_poroelastic_pressure_feedback',
                       'off: deformation does not change the pressure, so only flow is tested'),
                      ('control.gravity', 'no gravity, so the drawdown is compared directly in pressure'),
                      ('sim.max_steps',
                       '160 steps of 2 s end the run at 320 s: 200 s of pumping and 120 s of recovery'),
                      ('bc.hbc_z1',
                       'type 2: the top of the aquitard is held at the pressure given in bc.hbc_pressure, here 0, so '
                       'water leaks down through the aquitard'),
                      ('bc.well_rate', 'extraction rate of the quarter well, Q/4 with Q = 1e-4 m³/s'),
                      ('bc.well_z_top',
                       'the screen starts one element below the top, so the well draws from the aquifer only, not '
                       'the aquitard'),
                      ('bc.well_schedule_ratio',
                       'rate multipliers 1, 1, 0 at 0, 200 s, 200 s: full rate, then shut in'),
                      ('ic.layer_mattypes', 'material 1 (the aquitard) on top of material 0 (the aquifer)'),
                      ('ic.mattype_layer_depths',
                       'the aquitard is the top 0.2 of the 3.125 m height, i.e. one 0.625 m element'),
                      ('mat.hydraulic_perm',
                       'aquifer and aquitard permeabilities; the aquitard value is set so that the leakage factor λ '
                       'is 15 m'),
                      ('mat.permeability_tensor',
                       'vertical permeability ten times horizontal in both layers, so the aquifer has negligible '
                       'vertical resistance'),
                      ('mat.porosity',
                       'aquifer 0.25 and aquitard 1e-4 (with Biot coefficient 0.3), so the aquitard stores almost no '
                       'water, as Hantush–Jacob assume'),
                      ('mat.bulk_modulus',
                       'very stiff skeleton (1e12 Pa) in both layers, so the solid barely moves')],
             'fields': [('pore pressure',
                         'drawdown s_p = p(frame 0) − p at the aquifer nodes 10–30 m from the well, compared with '
                         'Hantush–Jacob'),
                        ('coordinate', 'gives the radius from the well and selects the nodes on the aquifer base')],
             'caption': '(a) Drawdown s_p/C_Q along the aquifer base on the y = 0 plane at four times during '
                        'pumping: DES nodes (circles) and the Hantush–Jacob solution (lines). The dashed grey curve '
                        'is the Theis solution without leakage at the last of these times, and the band marks the '
                        'compared range 10–30 m. (b) Drawdown history at r ≈ 10 m and r ≈ 20 m: Hantush–Jacob '
                        '(solid), Theis without leakage (dashed) and DES (circles); the dotted line marks shut-in. '
                        '(c) Maximum and RMS drawdown error against mesh spacing, with a slope-2 guide.',
             'study': ['python3 benchmarks/poromechanics/hantush_reference.py',
                       'python3 tests/functional/hydromechanics/check_hantush.py dynearthsol3d OUT --spacings 2.5 '
                       '1.25 0.625']},
 'theis_wells': {'run': [('dynearthsol3d input.cfg',
                          'Runs the 3D slab with three wells: 640 steps of 0.125 s to 80 s, all wells stopping at 40 '
                          's. DES writes multiwell.save.* frames every 80 steps (10 s).')],
                 'keys': [('control.has_hydraulic_diffusion', 'turns on Darcy flow of pore pressure; wells need it'),
                          ('control.has_poroelastic_pressure_feedback',
                           'off: deformation does not change the pressure, so the flow equation is linear and '
                           'responses superpose'),
                          ('control.gravity', 'no gravity, so drawdown is compared directly in pressure'),
                          ('control.fixed_dt',
                           'outer step of 0.125 s; DES splits it into stable flow substeps by itself'),
                          ('sim.max_steps',
                           '640 steps of 0.125 s end the run at 80 s, before the pressure front reaches the far '
                           'walls'),
                          ('sim.output_step_interval', 'a frame every 80 steps (10 s)'),
                          ('bc.num_wells',
                           'three vertical rate-controlled wells, each screened through the full 5 m'),
                          ('bc.well_x', 'x positions 60, 100 and 80 m'),
                          ('bc.well_y', 'y positions 80, 80 and 100 m'),
                          ('bc.well_rate',
                           'rates Q, 0.6Q and −0.25Q (Q = 1e-4 m³/s): two wells extract and the negative rate '
                           'injects'),
                          ('bc.well_schedule_time_in_yr',
                           'schedule times 0, 40 s, 40 s (in years); the repeated time makes a step at 40 s'),
                          ('bc.well_schedule_ratio',
                           'multipliers 1, 1, 0: every well runs at its full rate until 40 s and then stops'),
                          ('mat.bulk_modulus', 'very stiff skeleton (1e12 Pa), so the solid barely moves'),
                          ('mat.hydraulic_perm', 'chosen so that the pressure diffusivity is 1 m²/s')],
                 'fields': [('pore pressure',
                             'drawdown p(frame 0) − p, compared with superposed Theis at nodes 10–30 m from the '
                             'nearest well, and summed with element volumes for the fluid budget'),
                            ('coordinate', 'node positions, for the distance to each well and the y = 80 m profile'),
                            ('connectivity',
                             'with the coordinates, gives element volumes for the stored-fluid budget')],
                 'caption': '(a) Drawdown in kPa along the line y = 80 m at the end of pumping (40 s): DES nodes '
                            '(circles) and the superposed Theis solution of the three wells (line). The shaded bands '
                            'within 10 m of the two wells on this line are excluded from the comparison. (b) Maximum '
                            'pressure error relative to C_Σ against mesh spacing, for an outer step of 1 s and, on '
                            'the finer mesh, 0.125 s. (c) Closed fluid budget: the fluid deficit in litres, computed '
                            'from the DES pressure and element volumes (circles), against the time integral of the '
                            'net well rate (line).',
                 'study': ['python3 tests/functional/hydromechanics/check_theis_wells.py dynearthsol3d OUT '
                           '--spacings 5 2.5']},
 'open_flow': {'run': [('dynearthsol3d input.cfg',
                        'Restarts from the shipped seed.save.000000, seed.chkpt.000000 and seed.info, whose frame 0 '
                        'holds the initial sine pulse, and runs 200 steps of 0.05 s to 10 s. DES writes a '
                        'pulse.save.* frame at every step. To build the seed yourself, run the study command: it '
                        'runs one step, then writes the pulse into frame 0.')],
               'keys': [('sim.is_restarting', 'starts from a saved state instead of the built-in initial condition'),
                        ('sim.restarting_from_modelname', 'the shipped seed files, which hold the sine pulse'),
                        ('sim.restarting_from_frame', 'frame 0 of the seed, the initial pulse'),
                        ('sim.max_steps', '200 steps of 0.05 s, to 10 s'),
                        ('sim.output_step_interval',
                         'a frame at every step, so every step enters the pressure, storage and budget checks'),
                        ('control.has_hydraulic_diffusion', 'turns on Darcy flow of pore pressure'),
                        ('control.has_poroelastic_pressure_feedback',
                         'off: deformation does not change the pressure, so storage stays constant'),
                        ('control.gravity', 'no gravity, so the pulse drains by diffusion alone'),
                        ('control.fixed_dt',
                         'outer step of 0.05 s; DES splits it into stable flow substeps by itself'),
                        ('bc.hbc_x0',
                         'type 2: the end face x = 0 is held at the pressure in bc.hbc_pressure, here 0'),
                        ('bc.hbc_x1',
                         'type 2: the end face x = 10 m is also held at 0; the other faces are sealed (type 0)'),
                        ('mat.porosity', 'with the fluid bulk modulus, sets the constant storage'),
                        ('mat.hydraulic_perm', 'chosen so that the diffusivity is 1 m²/s'),
                        ('mat.bulk_modulus', 'very stiff skeleton (1e12 Pa), so the solid barely moves')],
               'fields': [('pore pressure',
                           'pressure along the column axis, compared with the decaying sine, and summed with element '
                           'volumes for the stored fluid'),
                          ('coordinate', 'x position of each node and element volumes for the stored fluid'),
                          ('fluid balance',
                           'five numbers per frame for the last step, in m³: storage, mechanics, wells, boundary '
                           'exchange and residual; the boundary entry, summed over steps, is the outflow in (b)')],
               'caption': '(a) Pressure p/P along the column axis at 0, 1, 5 and 10 s: DES nodes (circles) and the '
                          'decaying sine (lines). (b) Stored fluid and cumulative outflow, each divided by the '
                          'initial analytical storage M_0: DES (circles; stored fluid from pressure and element '
                          'volumes, outflow from the boundary term of the fluid balance) against the exponential '
                          'reference (lines). (c) Maximum pressure error relative to P and stored-volume error '
                          'relative to M_0 against mesh spacing.',
               'study': ['python3 tests/functional/hydromechanics/check_open_flow.py dynearthsol3d OUT --spacings 1 '
                         '0.5 0.25'],
               'replace': {'meaning': [('This checks open-boundary drainage and the discrete global fluid ledger at '
                                        'constant capacity. Unsaturated coupling, changing storage and remapping '
                                        'still require separate budgets.',
                                        'This checks drainage through open boundaries and the global water budget: '
                                        'stored water plus outflow stays equal to the initial volume. Storage is '
                                        'constant here; unsaturated coupling, changing storage and remeshing are not '
                                        'covered.')]}},
 'anisotropic_time': {'run': [('dynearthsol2d input.cfg',
                               'First phase: restarts from the shipped seed.save.000000, seed.chkpt.000000 and '
                               'seed.info, whose frame 0 holds the two-mode initial pressure, and runs 40,000 steps '
                               'of 1.25e-5 s to t_s = 0.5 s with the first permeability tensor. DES writes first.* '
                               'frames and checkpoints at 0.25 and 0.5 s.'),
                              ('dynearthsol2d continuation.cfg',
                               'Second phase: restarts from frame 2 of the first run (0.5 s) with doubled scalar '
                               'permeability and the new tensor, and runs to 1 s. DES writes second.* frames at 0.75 '
                               'and 1 s. It needs the first.* checkpoint written by the previous command.')],
                      'keys': [('sim.is_restarting',
                                'starts from a saved state instead of the built-in initial condition'),
                               ('sim.restarting_from_modelname',
                                'the shipped seed files, which hold the two sine modes; continuation.cfg restarts '
                                'from first instead'),
                               ('sim.max_steps',
                                '40,000 steps of 1.25e-5 s end the first phase at 0.5 s; continuation.cfg counts on '
                                'to 80,000 steps (1 s)'),
                               ('sim.output_step_interval',
                                'a frame every 20,000 steps (0.25 s), each with a checkpoint for the restart'),
                               ('control.fixed_dt', 'physical step of 1.25e-5 s, inside the elastic stability limit'),
                               ('control.is_quasi_static',
                                'off: mechanics runs with its physical inertia, without mass scaling'),
                               ('control.has_poroelastic_pressure_feedback',
                                'off: deformation does not change the pressure, so pressure obeys pure diffusion'),
                               ('control.gravity', 'no gravity'),
                               ('bc.hbc_x0',
                                'type 2: the side x = 0 is held at the pressure in bc.hbc_pressure, here 0; x = 1 m, '
                                'top and bottom are the same, so every side drains'),
                               ('bc.hbc_z0', 'type 2: the bottom is held at p = 0'),
                               ('mat.hydraulic_perm',
                                'scalar permeability k_0 = 1e-15 m² in the first phase; continuation.cfg doubles it'),
                               ('mat.permeability_tensor',
                                'dimensionless tensor multiplying hydraulic_perm, [xx, yy, zz, xy, xz, yz]; in 2D xx '
                                'is horizontal and zz vertical. [4, 1, 1] gives k_x = 4k_0, k_z = k_0; '
                                'continuation.cfg uses [1, 1, 9] with 2k_0, i.e. k_x = 2k_0, k_z = 18k_0'),
                               ('mat.porosity', 'with the fluid and grain moduli, sets the storage S_0'),
                               ('mat.bulk_modulus',
                                'soft skeleton (1e7 Pa, as the shear modulus), which sets the elastic step limit')],
                      'fields': [('pore pressure',
                                  'pressure along mid-depth, compared with the two-mode solution; the two mode '
                                  'amplitudes are projected from it'),
                                 ('coordinate', 'node positions for the mid-depth profile')],
                      'caption': '(a) Pressure p/P_0 along mid-depth (z = −0.5 m) at the four saved times of the two '
                                 'phases: DES nodes (circles) and the two-mode solution (lines). (b) Amplitudes of '
                                 'the two sine modes, projected from the DES pressure (circles), against their '
                                 'exponential decay (lines); the dotted line marks the permeability change at 0.5 s. '
                                 '(c) Maximum pressure error relative to P_0, in percent, against mesh spacing for '
                                 'the mesh study at 20,000 steps, the finest mesh at 80,000 steps (diamond), and a '
                                 'slope-2 guide.',
                      'study': ['python3 tests/functional/hydromechanics/check_anisotropic_flow.py dynearthsol2d '
                                'OUT-mesh',
                                'python3 tests/functional/hydromechanics/check_anisotropic_flow.py dynearthsol2d '
                                'OUT-steps40k --spacings .025 --steps 40000',
                                'python3 tests/functional/hydromechanics/check_anisotropic_flow.py dynearthsol2d '
                                'OUT-steps80k --spacings .025 --steps 80000'],
                      'replace': {'meaning': [('This checks directional diffusion, a prescribed coefficient jump and '
                                               'operator reconstruction at restart. Continuous state-dependent '
                                               'permeability, off-axis time changes and coupled fluid-content '
                                               'budgets remain separate requirements.',
                                               'This checks directional diffusion, a prescribed change of '
                                               'permeability during a run, and a restart after that change. '
                                               'Permeability that evolves with the state, rotated tensors that '
                                               'change in time and coupled fluid-content budgets are not '
                                               'covered.')]}},
 'point_sources': {'run': [('dynearthsol3d input.cfg',
                            'Runs the closed 3D box for 4 steps of 1 s with the two scheduled water and heat '
                            'sources. DES writes fresh.save.* frames at every step.')],
                   'keys': [('control.has_hydraulic_diffusion',
                             'turns on the pore-pressure equation that receives the water sources'),
                            ('control.has_poroelastic_pressure_feedback',
                             'off: deformation does not change the pressure, so stored water is storage times '
                             'pressure change'),
                            ('control.has_fluid_thermal_advection',
                             'turns on groundwater heat transport with the combined rock and water heat capacity, '
                             'which receives the heat sources'),
                            ('control.fixed_dt', 'steps of 1 s, so the source windows start and end inside steps'),
                            ('sim.max_steps', '4 steps of 1 s, to 4 s'),
                            ('mat.hydraulic_perm',
                             'zero permeability: water stays where it is supplied, so only the totals are tested'),
                            ('mat.therm_cond', 'zero conductivity: heat stays where it is supplied'),
                            ('mat.alpha',
                             'no thermal expansion of the solid (mat.fluid_alpha is also 0), so heat does not change '
                             'the pressure'),
                            ('mat.fluid_heat_capacity',
                             '4000 J/(kg K); with the rock values it gives C = 2.95e6 J/(m³ K)'),
                            ('injection.points_x',
                             'x of the two source points; points_y and points_z complete their positions '
                             '(heat_injection uses the same points)'),
                            ('injection.rate',
                             'water rates +8e-8 and −3e-8 m³/s: the first point supplies, the second extracts'),
                            ('injection.start_time_in_yr',
                             'start times 0.25 s and 0.8 s, in years (also used by heat_injection)'),
                            ('injection.end_time_in_yr', 'end times 2.75 s and 2.1 s, in years'),
                            ('heat_injection.rate', 'heat rates +4000 W and −1500 W at the same points and times')],
                   'fields': [('pore pressure',
                               'change since frame 0, times the nodal storage, gives the stored water'),
                              ('temperature',
                               'change since frame 0, times the nodal heat capacity, gives the stored heat'),
                              ('coordinate', 'with connectivity, gives the element volumes shared out to the nodes'),
                              ('connectivity', 'element-to-node map for the nodal volumes')],
                   'caption': '(a) Net water supplied to the closed box, in µL: the DES storage change from the '
                              'saved pore pressure (circles) against the time integral of the two signed source '
                              'schedules (line). (b) The same for heat, in kJ, from the saved temperature and the '
                              'combined heat capacity.',
                   'study': ['python3 tests/functional/hydromechanics/check_point_sources.py dynearthsol3d OUT3D '
                             '--dimension 3',
                             'python3 tests/functional/hydromechanics/check_point_sources.py dynearthsol2d OUT2D '
                             '--dimension 2']},
 'van_genuchten': {'run': [('dynearthsol2d input.cfg',
                            'Runs the 2D column with fixed 300 s steps for about 23 days (2e6 s), long enough to '
                            'reach the steady profile. DES writes column.save.* at the start, halfway and at the '
                            'end; the page compares the last frame.')],
                   'keys': [('control.has_hydraulic_diffusion',
                             'turns on Darcy flow, here partly saturated (Richards) flow'),
                            ('control.has_PT',
                             'PT solver: mechanical equilibrium is reached by iteration within each step'),
                            ('control.gravity', 'gravity drives the downward flow and sets the hydrostatic start'),
                            ('sim.max_time_in_yr',
                             'about 2e6 s in fixed 300 s steps, long enough for the column to reach steady state'),
                            ('sim.output_time_interval_in_yr', 'a frame halfway and one at the end'),
                            ('bc.hbc_z0',
                             'type 2: the base is held at the pressure in bc.hbc_pressure, here 0, the water table'),
                            ('bc.hbc_z1', 'type 3: the top has the prescribed flux in bc.hbc_flux'),
                            ('bc.hbc_flux',
                             'outward Darcy flux per face; the top value −1.389e-6 m/s is an inflow of half the '
                             'saturated conductivity'),
                            ('ic.water_table_z', 'starts hydrostatic with the water table at the base (z = −1 m)'),
                            ('mat.retention_model',
                             '1: van Genuchten retention with Mualem relative permeability (n is retention_n, '
                             'default 2)'),
                            ('mat.retention_alpha', 'pressure-based retention parameter α_p = 5e-4 1/Pa'),
                            ('mat.residual_saturation', 'residual saturation S_r = 0.1'),
                            ('mat.fluid_bulk_modulus',
                             'deliberately small (4e4 Pa) to allow a larger stable flow step; storage does not '
                             'affect the steady profile'),
                            ('mat.hydraulic_perm', 'gives saturated conductivity K_s = 2.78e-6 m/s (1 cm/h)')],
                   'fields': [('pore pressure',
                               'the steady pressure profile at the last frame, compared with the reference'),
                              ('coordinate', 'the height of each node above the base')],
                   'caption': 'Steady pressure (kPa) against height above the water table at the end of the run: DES '
                              'nodes (circles) and the steady-infiltration reference for the van Genuchten–Mualem '
                              'law (line).',
                   'study': ['python3 tests/functional/hydromechanics/check_unsaturated_column.py dynearthsol2d '
                             'OUT']},
 'brooks_corey': {'run': [('dynearthsol2d input.cfg',
                           'Runs the 2D column with fixed 300 s steps for about 23 days (2e6 s), long enough to '
                           'reach the steady profile. DES writes column.save.* at the start, halfway and at the end; '
                           'the page compares the last frame.')],
                  'keys': [('control.has_hydraulic_diffusion',
                            'turns on Darcy flow, here partly saturated (Richards) flow'),
                           ('control.has_PT',
                            'PT solver: mechanical equilibrium is reached by iteration within each step'),
                           ('control.gravity', 'gravity drives the downward flow and sets the hydrostatic start'),
                           ('sim.max_time_in_yr',
                            'about 2e6 s in fixed 300 s steps, long enough for the column to reach steady state'),
                           ('bc.hbc_z0',
                            'type 2: the base is held at the pressure in bc.hbc_pressure, here 0, the water table'),
                           ('bc.hbc_z1', 'type 3: the top has the prescribed flux in bc.hbc_flux'),
                           ('bc.hbc_flux',
                            'outward Darcy flux per face; the top value −1.389e-6 m/s is an inflow of half the '
                            'saturated conductivity'),
                           ('ic.water_table_z', 'starts hydrostatic with the water table at the base (z = −1 m)'),
                           ('mat.retention_model', '2: Brooks–Corey retention with Burdine relative permeability'),
                           ('mat.retention_alpha',
                            'α_p = 1e-3 1/Pa: the soil stays saturated down to the air-entry pressure −1/α_p = −1000 '
                            'Pa'),
                           ('mat.retention_n', 'Brooks–Corey pore-size index 1.5'),
                           ('mat.residual_saturation', 'residual saturation S_r = 0.05'),
                           ('mat.fluid_bulk_modulus',
                            'deliberately small (4e4 Pa) to allow a larger stable flow step; storage does not affect '
                            'the steady profile'),
                           ('mat.hydraulic_perm', 'gives saturated conductivity K_s = 2.78e-6 m/s (1 cm/h)')],
                  'fields': [('pore pressure',
                              'the steady pressure profile at the last frame, compared with the reference'),
                             ('coordinate', 'the height of each node above the base')],
                  'caption': 'Steady pressure (kPa) against height above the water table at the end of the run: DES '
                             'nodes (circles) and the steady-infiltration reference for the Brooks–Corey law (line). '
                             'The profile bends where the soil starts to desaturate at the air-entry pressure.',
                  'study': ['python3 tests/functional/hydromechanics/check_unsaturated_column.py dynearthsol2d OUT']},
 'gardner': {'run': [('dynearthsol3d input.cfg',
                      'Runs the 3D column for 2 h (7200 s) in fixed 30 s steps after the top inflow starts. DES '
                      'writes column.save.* every 30 min; the page shows the last frame and checks all of them.')],
             'keys': [('control.has_hydraulic_diffusion',
                       'turns on Darcy flow, here partly saturated (Richards) flow'),
                      ('control.has_PT',
                       'PT solver: mechanical equilibrium is reached by iteration within each step'),
                      ('control.gravity', 'gravity drives the downward flow and sets the hydrostatic start'),
                      ('control.fixed_dt', 'fixed 30 s steps'),
                      ('sim.max_time_in_yr', '7200 s (2 h) of infiltration'),
                      ('sim.output_time_interval_in_yr', 'a frame every 1800 s (30 min)'),
                      ('bc.hbc_z0',
                       'type 2: the base is held at the pressure in bc.hbc_pressure, here 0, the water table'),
                      ('bc.hbc_z1', 'type 3: the top has the prescribed flux in bc.hbc_flux'),
                      ('bc.hbc_flux',
                       'outward Darcy flux per face; the top value −2.5e-6 m/s is an inflow of 0.9 K_s, starting at '
                       't = 0'),
                      ('ic.water_table_z', 'starts hydrostatic with the water table at the base (z = −1 m)'),
                      ('mat.retention_model', '3: Gardner exponential law, S_e = k_r = exp(α_p p)'),
                      ('mat.retention_alpha', 'α_p = 1.019e-3 1/Pa, i.e. 0.1 per cm of head'),
                      ('mat.porosity', '0.4; with zero residual saturation the water content can change by Δθ = 0.4'),
                      ('mat.hydraulic_perm', 'gives saturated conductivity K_s = 2.78e-6 m/s (1 cm/h)')],
             'fields': [('pore pressure',
                         'converted to K/K_s = exp(α_p p) at every node and compared with the series solution'),
                        ('coordinate', 'the height of each node above the base')],
             'caption': 'Relative conductivity K/K_s = exp(α_p p), computed from the DES pore pressure (circles), '
                        'against the Srivastava–Yeh series (line) along the column at the end of the run (7200 s).',
             'study': ['python3 tests/functional/hydromechanics/check_unsaturated_column.py dynearthsol3d OUT']},
 'bishop': {'run': [('dynearthsol2d input.cfg',
                     'Runs the 10 m column with the van Genuchten law: 40 steps of 1e4 s. The base pressure is '
                     'lowered between steps 4 and 5 and the column then drains. DES writes bishop.save.* every 4 '
                     'steps; the page compares frame 1 (drained start) with the last frame. For the other laws set '
                     'mat.retention_model to 2 (retention_alpha 1e-4, retention_n 0.8, residual_saturation 0.05) or '
                     '3 (retention_alpha 5.097e-5, residual_saturation 0).')],
            'keys': [('control.has_PT',
                      'PT solver: mechanical equilibrium is reached by iteration within each step, with pressure '
                      'feedback on'),
                     ('control.gravity',
                      'gravity acts on the partly saturated bulk weight and sets the hydrostatic water table'),
                     ('control.fixed_dt',
                      'steps of 1e4 s, long enough for the column to drain to the new state well within the run'),
                     ('sim.output_step_interval',
                      'a frame every 4 steps; frame 1 is the drained start just before the base pressure drops'),
                     ('bc.hbc_z0',
                      'type 2: the base pressure is prescribed by bc.hbc_pressure and the schedule; the other faces '
                      'are sealed'),
                     ('bc.hbc_pressure', 'base pressure 88,290 Pa = ρ_f g × 9 m, a water table 9 m above the base'),
                     ('bc.hbc_schedule_time_in_yr',
                      'the schedule holds the base pressure until step 4 and changes it by step 5'),
                     ('bc.hbc_schedule_ratio',
                      'multipliers 1, 1, 5/9: the base pressure falls to that of a water table at 5 m'),
                     ('ic.water_table_z',
                      'initial hydrostatic water table 1 m below the top (z = −1 m), i.e. 9 m above the base'),
                     ('mat.bulk_modulus',
                      'drained bulk modulus K = 2e7 Pa (shear modulus 1e7 Pa) of the elastic skeleton'),
                     ('mat.biot_coeff', 'Biot coefficient α = 0.9'),
                     ('mat.retention_model',
                      '1: van Genuchten with retention_alpha 5e-5 1/Pa, retention_n 1.6 and residual_saturation 0.1; '
                      'Bishop χ = S_e'),
                     ('mat.fluid_bulk_modulus',
                      'deliberately small (1e5 Pa) to allow a larger stable flow step; the drained end states do not '
                      'depend on it'),
                     ('mat.hydraulic_perm', 'high permeability (1e-10 m²) so that the column drains quickly')],
            'fields': [('strain',
                        'vertical strain of each element; its change between the two drained frames, integrated up '
                        'from the fixed base, gives the settlement u_z'),
                       ('pore pressure', 'checks that the final state is drained to the new hydrostatic profile'),
                       ('coordinate', 'node heights for the profiles'),
                       ('connectivity', 'assigns element strains to the layers of the column')],
            'caption': '(a) Settlement u_z (mm) between the two drained states against height z, for the van '
                       'Genuchten, Brooks–Corey and Gardner laws and the column with two laws: DES (circles, '
                       'integrated from the element strains) and the exact solution (lines). The dotted lines mark '
                       'the initial and final water tables. (b) Change of vertical strain against height, same key. '
                       '(c) Maximum settlement error relative to the largest settlement, against mesh spacing: 2D as '
                       'lines with circles, 3D as open squares, and a slope-2 guide.',
            'study': ['python3 tests/functional/hydromechanics/check_bishop_column.py dynearthsol2d OUT --spacings 1 '
                      '0.5 0.25',
                      'python3 tests/functional/hydromechanics/check_bishop_column.py dynearthsol3d OUT3D --ndims 3 '
                      '--spacings 1 0.5 0.25']},
 'bishop_cell': {'run': [('dynearthsol2d input.cfg',
                          'Runs the moderate-suction cell: 320 PT steps of 3.125 s compress the top by 1 mm over '
                          '1000 s. DES writes cell.save.* every 32 steps. For the strong-suction state set '
                          'ic.pore_pressure_value = -2000 and mat.retention_alpha = [1e-3].')],
                 'keys': [('control.has_PT',
                           'PT solver with the fixed-stress coupling whose Bishop correction is tested; pressure '
                           'feedback is on'),
                          ('control.gravity', 'no gravity, so the cell is uniform'),
                          ('control.fixed_dt', '3.125 s steps: 320 steps over the 1000 s loading'),
                          ('sim.max_steps', '320 steps, the finest step count of the study'),
                          ('sim.output_step_interval', 'a frame every 32 steps (100 s)'),
                          ('bc.vbc_z1', 'type 1: the top moves vertically at the velocity bc.vbc_val_z1'),
                          ('bc.vbc_val_z1',
                           '−1e-6 m/s: the top moves down 1 mm in 1000 s, a vertical strain of −1e-3; the other '
                           'faces are held in their normal direction'),
                          ('ic.pore_pressure_option', '2: uniform initial pressure'),
                          ('ic.pore_pressure_value', 'initial pressure −5000 Pa (5 kPa suction)'),
                          ('mat.retention_model',
                           '3: Gardner exponential law (retention_alpha 1e-4 1/Pa, the default); Bishop χ = S_e'),
                          ('mat.residual_saturation', 'residual saturation S_r = 0.1'),
                          ('mat.hydraulic_perm',
                           'negligible permeability (1e-30 m²): no water moves, so the cell is undrained'),
                          ('mat.bulk_modulus',
                           'soft skeleton, K = 1e4 Pa, so the fixed-stress stabilising storage is large'),
                          ('mat.shear_modulus', 'G = 1e4 Pa')],
                 'fields': [('pore pressure',
                             'averaged over the cell at each frame and compared with the exact rate equation'),
                            ('strain',
                             'vertical strain averaged over the cell, the path along which the exact equation is '
                             'integrated'),
                            ('coordinate', 'identifies the vertical axis')],
                 'caption': '(a) Moderate and (b) strong suction: change of mean pore pressure p − p_0 against '
                            "compression −ε_zz. Line: the exact rate equation along the run's own strain path. "
                            'Circles: PT at the finest step, with the corrected fixed-stress source (coloured) and '
                            'without the Bishop correction (grey). Squares: the explicit strain form. (c) Maximum '
                            'error relative to the pressure change against the number of PT steps over the loading, '
                            'for both states in 2D, in 3D and in 2D without the correction, with a slope −1 guide.',
                 'study': ['python3 tests/functional/hydromechanics/check_bishop_cell.py dynearthsol2d OUT --steps '
                           '10 20 40 80 160 320 --explicit',
                           'python3 tests/functional/hydromechanics/check_bishop_cell.py dynearthsol3d OUT3D --ndims '
                           '3 --steps 10 20 40 80 160',
                           'python3 tests/functional/hydromechanics/check_bishop_cell.py dynearthsol2d OUT-restart '
                           '--restart'],
                 'replace': {'problem': [('A sealed, partly saturated soil cell is compressed slowly. Without '
                                          'drainage, the pore pressure follows the strain through the Bishop rate '
                                          'equation; the stabilised fixed-stress coupling used with the PT solver '
                                          'must recover it as the step shrinks, also where the Bishop pressure Se p '
                                          'decreases with suction.',
                                          'A sealed, partly saturated soil cell is compressed slowly. Without '
                                          'drainage, the pore pressure follows the strain through the Bishop rate '
                                          'equation. The fixed-stress coupling of the PT solver must recover it as '
                                          'the time step shrinks, including where the Bishop pressure $S_ep$ '
                                          'decreases with suction.')],
                             'result': [('Largest pressure error at 320 steps: 1.43% (moderate suction) and 0.64% '
                                         '(strong suction) in 2D; the error halves with each halving of the step. '
                                         'Without the Bishop correction of the fixed-stress source the error stayed '
                                         'at 78.89% and 67.62%.',
                                         'Largest pressure error at 320 steps: 1.43% (moderate suction) and 0.64% '
                                         '(strong suction) in 2D; the error halves with each halving of the step. '
                                         'Without the Bishop correction of the fixed-stress source the scheme does '
                                         'not converge: the error remains 78.89% and 67.62%.')],
                             'meaning': [('This checks the transient Bishop coupling of the PT scheme: the '
                                          'fixed-stress source corrected by $L(\\chi E-D)$, with the increments the '
                                          'mechanics consumed. It does not test flow, which the column checks cover, '
                                          'or large strain.',
                                          'This checks the transient Bishop coupling under PT: the fixed-stress '
                                          'source with its correction $L(\\chi E-D)$, built from the increments of '
                                          'the last mechanical step. It does not test flow, which the column pages '
                                          'cover, or large strain.')]}},
 'drainage': {'run': [('dynearthsol2d input.cfg',
                       'Runs the 1 m column on the 0.0125 m mesh: 480 PT steps of 15 s, 2 h in all. The base '
                       'pressure drops to 0 within the first millisecond and the column drains. DES writes '
                       'drain.save.* every 20 steps (5 min).')],
              'keys': [('control.has_PT',
                        'PT solver: mechanical equilibrium within each step, with two-way pressure feedback'),
                       ('control.gravity', 'gravity drives the drainage and sets the hydrostatic start'),
                       ('control.fixed_dt', '15 s steps: 480 steps over 2 h'),
                       ('sim.output_step_interval',
                        'a frame every 20 steps (5 min), which includes the compared times 5, 20, 60 and 120 min'),
                       ('bc.hbc_z0',
                        'type 2: the base pressure is prescribed by bc.hbc_pressure and the schedule; the top and '
                        'sides are sealed'),
                       ('bc.hbc_pressure',
                        'base pressure 9810 Pa, the hydrostatic value with the water table at the top '
                        '(ic.water_table_z = 0, the default)'),
                       ('bc.hbc_schedule_time_in_yr', 'schedule times 0 and 1 ms, in years'),
                       ('bc.hbc_schedule_ratio',
                        'multipliers 1 and 0: the base pressure ramps to 0 over the first millisecond, opening the '
                        'base'),
                       ('ic.has_body_force_adjustment',
                        'brings the column into mechanical equilibrium under gravity before the run'),
                       ('mat.bulk_modulus', 'K = 2.17e6 Pa, from E = 1.3 MPa and ν = 0.4'),
                       ('mat.shear_modulus', 'G = 4.64e5 Pa, from E = 1.3 MPa and ν = 0.4'),
                       ('mat.porosity', 'φ = 0.2975, the Liakopoulos sand'),
                       ('mat.hydraulic_perm', 'intrinsic permeability 4.5e-13 m²'),
                       ('mat.retention_model',
                        '1: van Genuchten–Mualem, with retention_alpha 1/3000 1/Pa, retention_n 1.8 and '
                        'residual_saturation 0.1; Bishop χ = S_e')],
              'fields': [('pore pressure', 'pressure profiles at 5, 20, 60 and 120 min, compared with the reference'),
                         ('strain',
                          'vertical element strain; its change since frame 0, summed up the column, gives the '
                          'settlement of the top'),
                         ('coordinate', 'node heights'),
                         ('connectivity', 'assigns element strains to the layers of the column')],
              'caption': '(a) Pore pressure (kPa) against height at 5, 20, 60 and 120 min: DES nodes on one side of '
                         'the column (circles) and the numerical reference (lines); the dotted line is the initial '
                         'hydrostatic profile. (b) Settlement of the top against time on a log scale: DES at each '
                         'saved frame (circles) and the reference (line). (c) Maximum pressure error relative to ρ_f '
                         'g H against mesh spacing at the four times: circles for 480 steps, diamonds for 960 steps '
                         'on the finest mesh, open squares for the explicit form.',
              'study': ['python3 benchmarks/poromechanics/drainage_reference.py',
                        'python3 tests/functional/hydromechanics/check_unsaturated_drainage.py dynearthsol2d OUT '
                        '--spacings 0.05 0.025 0.0125',
                        'python3 tests/functional/hydromechanics/check_unsaturated_drainage.py dynearthsol2d '
                        'OUT-explicit --spacings 0.05 --explicit',
                        'python3 tests/functional/hydromechanics/check_unsaturated_drainage.py dynearthsol2d '
                        'OUT-steps960 --spacings 0.0125 --steps 960']},
 'content': {'run': [('dynearthsol2d input.cfg',
                      'Runs the 2D column on a 0.05 m mesh for 800 steps of 0.005 s (4 s) and writes column.save.* '
                      'and column.chkpt.* at every step, with column.info listing the frames.')],
             'keys': [('control.has_hydraulic_diffusion',
                       'turns on Darcy flow, so pore pressure diffuses through the column'),
                      ('control.has_PT',
                       'solves the mechanics to equilibrium at every step with the pseudo-transient (PT) iteration'),
                      ('control.fixed_dt', 'physical time step; with max_steps it sets the 4 s simulated'),
                      ('control.gravity', 'no gravity: all pressure comes from the water pumped in at the base'),
                      ('bc.hbc_z0', 'type 3: the base is a prescribed-flux boundary, taking its value from hbc_flux'),
                      ('bc.hbc_flux',
                       'the bottom (z0) entry is an outward flux; its negative value pumps water in through the '
                       'base'),
                      ('bc.hbc_z1',
                       'type 2: the top is drained, pore pressure held at its hbc_pressure entry (0 Pa)'),
                      ('bc.stress_bc_z1',
                       'type 3: vertical traction on the top; with stress_val_z1 = 0 the top is traction-free and '
                       'can rise as the column swells'),
                      ('sim.output_step_interval',
                       'saves every step, so stored and drained water can be summed step by step'),
                      ('mat.porosity',
                       'pore fraction; with the fluid and grain moduli it sets the pressure storage S'),
                      ('mat.biot_coeff',
                       'Biot coefficient alpha: how strongly pore pressure and volume strain are coupled'),
                      ('mat.bulk_modulus_s', 'grain (solid) bulk modulus K_s, part of the storage S'),
                      ('mat.fluid_bulk_modulus', 'water bulk modulus K_f, part of the storage S'),
                      ('mat.hydraulic_perm',
                       'intrinsic permeability k; with the viscosity it sets the mobility k/mu and the 100 Pa steady '
                       'base pressure')],
             'fields': [('pore pressure',
                         'nodal pressure; the vertical profile on the column centre line is compared with the series '
                         'in panel (a), and its volume integral gives the pressure part of the stored water'),
                        ('strain',
                         'element strain; its volumetric part times the Biot coefficient gives the deformation part '
                         'of the stored water in panel (b)'),
                        ('fluid balance',
                         'five numbers per frame for the last step: storage change, mechanical source, wells, '
                         'boundary exchange and residual; the summed boundary exchange gives the drained water in '
                         'panel (b)')],
             'caption': '(a) Pressure down the centre of the column at several times while it fills: lines are the '
                        'series solution, symbols are DES. (b) Water stored in the column and water drained through '
                        'the top, per unit cross-section, against time: lines are the reference, symbols are DES '
                        '(stored water from pressure and volume strain, drained water from the boundary exchange).',
             'study': ['python3 benchmarks/poromechanics/fluid_content.py dynearthsol2d OUT/spatial-pt --spacings '
                       '.2 .1 .05 --dt .005 --duration 4 --figures',
                       'python3 benchmarks/poromechanics/fluid_content.py dynearthsol2d OUT/time-coarse --spacings '
                       '.05 --dt 0.02 --duration 1',
                       'python3 benchmarks/poromechanics/fluid_content.py dynearthsol2d OUT/time-medium --spacings '
                       '.05 --dt 0.01 --duration 1',
                       'python3 benchmarks/poromechanics/fluid_content.py dynearthsol2d OUT/time-fine --spacings '
                       '.05 --dt 0.005 --duration 1'],
             'replace': {'result': [('Late pressure error is 0.0286% of the steady scale. The coupled split budget '
                                     'closes within 2.4e-15 of injected volume across 19 cases.',
                                     'Late pressure error is 0.0286% of the steady base pressure. In all 19 runs the '
                                     'stored water, including the small increment the PT solver carries into the '
                                     'next step, equals the net inflow to within 2.4e-15 of the injected volume.')],
                         'meaning': [('Independent mesh-integrated pressure and strain, boundary reaction and Darcy '
                                      'outflow distinguish a physical budget from the pressure-equation residual. '
                                      'Temporal controls measure first-order PT lag; ordinary explicit, adaptive '
                                      'GVS, tiny 3D and restart controls are supporting evidence.',
                                      'Stored water is computed from the pressure and volume strain integrated over '
                                      'the mesh, and drained water both from the boundary exchange and from the '
                                      'Darcy flux at the top, so the check does not reuse the pressure equation. '
                                      'With PT the stored water lags by about one step; the lag halves when the time '
                                      'step is halved. Explicit, GVS, 3D and restarted runs close the same '
                                      'budget.')]}},
 'terzaghi': {'run': [('dynearthsol2d input.cfg',
                       'Runs the 10 m column on a 0.25 m mesh for 33600 steps of 0.25 s (8400 s, about 99% '
                       'consolidated) and writes terzaghi.save.000000 to terzaghi.save.000007, one frame every 1200 '
                       's; frame 0 is the state before the load.')],
              'keys': [('control.has_hydraulic_diffusion',
                        'turns on Darcy flow, so the excess pressure drains out through the top'),
                       ('control.has_PT',
                        'solves the mechanics to equilibrium at every step with PT, so the sudden load is carried '
                        'without elastic waves'),
                       ('control.has_moving_mesh',
                        'keeps the mesh fixed (small strain); settlement is computed from the vertical strain, not '
                        'from moved nodes'),
                       ('control.gravity', 'no gravity: only the applied load creates pressure'),
                       ('control.fixed_dt', 'physical time step'),
                       ('sim.max_steps', 'number of steps; with fixed_dt it sets the 8400 s simulated'),
                       ('sim.output_step_interval', 'writes a frame every 4800 steps (1200 s)'),
                       ('bc.hbc_z1', 'type 1: the top is drained, pore pressure held at its initial value (0)'),
                       ('bc.stress_bc_z1', 'type 3: a vertical traction is applied on the top'),
                       ('bc.stress_val_z1', 'the top load; negative is compression (100 kPa)'),
                       ('bc.vbc_z0', 'fixes vertical motion of the base'),
                       ('mat.bulk_modulus', 'drained bulk modulus K of the skeleton'),
                       ('mat.shear_modulus', 'shear modulus G; with K it sets the constrained modulus K + 4G/3'),
                       ('mat.hydraulic_perm',
                        'intrinsic permeability; with viscosity and storage it sets the consolidation coefficient '
                        'c_v')],
              'fields': [('pore pressure',
                          'nodal pressure; profiles against depth below the drained top are compared with the series '
                          'in panel (a)'),
                         ('strain',
                          'element strain; the vertical component integrated over the column gives the top '
                          'settlement in panel (b)')],
              'caption': '(a) Pressure profiles against depth below the drained top at several dimensionless times '
                         'T_v: lines are the series solution, symbols are DES. (b) Settlement of the top against '
                         'time over the full drainage: the line is the series solution, symbols are DES, and the '
                         'dashed line is the final drained settlement.',
              'study': ['python3 benchmarks/poromechanics/consolidation_refinement.py --output OUT/history-2d '
                        '--dimensions 2 --spacings 1 0.5 0.25 --timesteps 1 0.5 0.25 --spatial-dt 0.25 '
                        '--temporal-spacing 0.5 --duration 8400 --sample-interval 1200 --timeout 30000',
                        'python3 benchmarks/poromechanics/consolidation_refinement.py --output OUT/sanity-3d '
                        '--dimensions 3 --spacings 1 0.5 0.25 --timesteps 1 0.5 0.25 --spatial-dt 0.25 '
                        '--temporal-spacing 1 --duration 8400 --sample-interval 1200 --timeout 30000']},
 'mandel': {'run': [('dynearthsol2d input.cfg',
                     'Runs the quarter strip (a = b = 1 m) on a 0.05 m mesh for 800 steps over one consolidation '
                     'time a^2/c (41.8 s) and writes mandel.save.* every 20 steps, with a checkpoint every 10th '
                     'frame.')],
            'keys': [('control.has_hydraulic_diffusion',
                      'turns on Darcy flow, so water drains out through the side x = a'),
                     ('control.has_PT',
                      'solves the mechanics to equilibrium at every step with PT, so the sudden plate load is '
                      'carried without elastic waves'),
                     ('control.gravity', 'no gravity: only the plate load creates pressure'),
                     ('control.fixed_dt', 'physical time step, 1/800 of the consolidation time a^2/c'),
                     ('sim.max_steps', '800 steps; this, not max_time_in_yr, ends the run'),
                     ('sim.output_step_interval', 'writes a frame every 20 steps'),
                     ('bc.vbc_x1', 'type 0: the drained side x = a is free to move (traction-free)'),
                     ('bc.vbc_z0', 'fixes vertical motion on the bottom edge, the mid-plane of the strip (symmetry)'),
                     ('bc.vbc_z1',
                      'type 8: rigid frictionless platen; all top nodes share one vertical velocity and slide freely '
                      'sideways'),
                     ('bc.stress_bc_z1', 'type 3: the platen is pushed by a vertical traction'),
                     ('bc.stress_val_z1',
                      'the plate pressure; negative is compression (100 kPa, total force F = Pa)'),
                     ('bc.hbc_x1',
                      'type 2: the side x = a is drained, pore pressure held at its hbc_pressure entry (0 Pa); the '
                      'other sides stay sealed'),
                     ('mat.bulk_modulus', 'drained bulk modulus K (with G it sets the Poisson ratio and p0)'),
                     ('mat.hydraulic_perm',
                      'intrinsic permeability; with storage and stiffness it sets the consolidation coefficient c')],
            'fields': [('pore pressure',
                        'nodal pressure; centre and x = a/2 histories in panel (a), profiles across the strip in '
                        'panel (b), and the error at every node in panel (c)'),
                       ('strain',
                        'element strain; the vertical component averaged over the quarter gives the plate '
                        'settlement, whose error is in panel (c)'),
                       ('velocity',
                        'nodal velocity; the vertical velocity of the top row should be the same at every node under '
                        'the rigid platen')],
            'caption': '(a) Pressure history at the centre x = 0 and at x = a/2, scaled by the undrained pressure '
                       'p0, against ct/a^2: lines are the series solution, symbols are DES; the dotted line marks '
                       'p/p0 = 1, which the centre pressure first exceeds. (b) Pressure profiles across the strip at '
                       'four times. (c) Largest pressure error (over all nodes and times, divided by p0) and largest '
                       'relative settlement error against the time step, on log axes, with a slope-1 guide.',
            'study': ['python3 tests/functional/hydromechanics/check_mandel.py dynearthsol2d OUT/mandel-h --spacings '
                      '0.1 0.05 0.025 --steps 200',
                      'python3 tests/functional/hydromechanics/check_mandel.py dynearthsol2d OUT/mandel-t100 '
                      '--spacings 0.05 --steps 100',
                      'python3 tests/functional/hydromechanics/check_mandel.py dynearthsol2d OUT/mandel-t400 '
                      '--spacings 0.05 --steps 400',
                      'python3 tests/functional/hydromechanics/check_mandel.py dynearthsol2d OUT/mandel-t800 '
                      '--spacings 0.05 --steps 800']},
 'cryer': {'run': [('dynearthsol3d input.cfg',
                    'Meshes the sphere octant from octant.poly (shipped in the same folder) at 0.07 m, runs 800 '
                    'steps over 0.5 R^2/c (20.9 s) and writes cryer.save.* every 20 steps, with a checkpoint every '
                    '10th frame. This is a 3D run and takes much longer than the 2D examples.')],
           'keys': [('mesh.meshing_option', '90: the domain boundary is read from a .poly file instead of a box'),
                    ('mesh.poly_filename',
                     'octant.poly: one octant of a unit sphere; the curved surface carries the top (z1) boundary '
                     'flag, the flat faces x = 0, y = 0 and z = 0 the x0, y0 and z0 flags'),
                    ('mesh.resolution', 'target element size used when meshing the octant'),
                    ('control.has_hydraulic_diffusion',
                     'turns on Darcy flow, so water drains out through the curved surface'),
                    ('control.has_PT',
                     'solves the mechanics to equilibrium at every step with PT, so the sudden surface load is '
                     'carried without elastic waves'),
                    ('control.gravity', 'no gravity: only the surface load creates pressure'),
                    ('control.fixed_dt', 'physical time step, 1/800 of 0.5 R^2/c'),
                    ('sim.max_steps', '800 steps; this, not max_time_in_yr, ends the run'),
                    ('bc.vbc_y0',
                     'fixes normal motion on the symmetry plane y = 0 (x = 0 is fixed by the default vbc_x0)'),
                    ('bc.vbc_z0',
                     'fixes normal motion on the symmetry plane z = 0, which carries the z0 flag in octant.poly'),
                    ('bc.hbc_z1',
                     'type 2: the curved surface is drained, pore pressure held at its hbc_pressure entry (0 Pa); '
                     'the symmetry planes stay sealed'),
                    ('bc.stress_bc_z1',
                     'type 4: a normal stress on facets of any orientation, here the curved sphere surface'),
                    ('bc.stress_val_z1', 'the surface pressure; negative is compression (100 kPa)'),
                    ('mat.hydraulic_perm',
                     'intrinsic permeability; with storage and stiffness it sets the consolidation coefficient c')],
           'fields': [('pore pressure',
                       'nodal pressure; the node nearest the centre gives panel (a), all nodes binned by radius give '
                       'the profiles in panel (b), and the largest nodal error gives panel (c)')],
           'caption': '(a) Pressure at the centre of the sphere, scaled by the undrained pressure p0, against '
                      'ct/R^2: the line is the series solution, symbols are DES; the dotted line marks p/p0 = 1, '
                      'which the centre pressure first exceeds. (b) Radial pressure profiles at four times: lines '
                      'are the series solution, symbols are the mean of the DES nodes in each radius bin and the '
                      'shaded band their spread. (c) Largest pressure error after ct/R^2 = 0.05 against mesh '
                      'spacing, one line per number of time steps.',
           'study': ['python3 tests/functional/hydromechanics/check_cryer.py dynearthsol3d OUT/cryer-h --spacings '
                     '0.2 0.1 0.07 --steps 200',
                     'python3 tests/functional/hydromechanics/check_cryer.py dynearthsol3d OUT/cryer-t400 --spacings '
                     '0.1 --steps 400',
                     'python3 tests/functional/hydromechanics/check_cryer.py dynearthsol3d OUT/cryer-h0.1-t800 '
                     '--spacings 0.1 --steps 800',
                     'python3 tests/functional/hydromechanics/check_cryer.py dynearthsol3d OUT/cryer-h0.07-t400 '
                     '--spacings 0.07 --steps 400',
                     'python3 tests/functional/hydromechanics/check_cryer.py dynearthsol3d OUT/cryer-h0.07-t800 '
                     '--spacings 0.07 --steps 800']},
 'averaged_aquifer': {'run': [('dynearthsol3d input.cfg',
                               'Runs the 3D quarter slab on a 5 m hexahedral mesh: the corner well pumps from t = 0 '
                               'and DES writes frames slab.* at 0, 2500, 5000, 7500 and 10000 s. No mesh or restart '
                               'file is needed.')],
                      'keys': [('control.has_hydraulic_diffusion',
                                'turns on Darcy flow, so the drawdown spreads from the well by diffusion'),
                               ('control.has_PT',
                                'iterates each time step to mechanical equilibrium, so the slab deforms '
                                'quasi-statically with the pressure (two-way Biot coupling)'),
                               ('control.gravity',
                                'no gravity: every field is a change from the unstressed, zero-pressure initial '
                                'state'),
                               ('control.has_moving_mesh',
                                'keeps the mesh fixed (small strain), as the reference assumes'),
                               ('mesh.meshing_elem_shape',
                                'structured hexahedra, so the mesh spacing follows mesh.resolution exactly'),
                               ('bc.vbc_y0',
                                'symmetry plane y = 0: normal displacement fixed (the x = 0 plane is a roller by '
                                'default)'),
                               ('bc.vbc_z0',
                                'roller base: vertical displacement fixed, horizontal free; the top is '
                                'traction-free'),
                               ('bc.hbc_x1', 'far side x = 400 m holds its initial pore pressure (zero change)'),
                               ('bc.hbc_y1', 'far side y = 400 m holds its initial pore pressure (zero change)'),
                               ('bc.num_wells',
                                'one rate-controlled vertical well, placed on the corner edge x = y = 0'),
                               ('bc.well_rate',
                                'extraction rate of the quarter well: one quarter of the full rate Q, because the '
                                'model is one quarter of the slab'),
                               ('bc.well_z_bottom',
                                'screens the well over the whole 10 m thickness (fully penetrating)'),
                               ('mat.hydraulic_perm',
                                'permeability; with fluid viscosity and storage it sets the hydraulic diffusivity'),
                               ('mat.biot_coeff',
                                'Biot coefficient, the strength of the pressure–deformation coupling')],
                      'fields': [('pore pressure',
                                  'nodal pressure change against radius from the well at each saved time; panel (a) '
                                  'and the pressure error'),
                                 ('strain',
                                  'element strain tensor: the areal strain is the xx + yy components and the '
                                  'vertical strain the zz component; panel (b) and the strain errors'),
                                 ('coordinate',
                                  'node positions, used for the radius of nodes and element centres from the well')],
                      'caption': '(a) Pressure drop against distance from the well at the first and last saved '
                                 'times: lines are the averaged-aquifer line-well solution, symbols the mean DES '
                                 'value in radius bins with the spread shaded. (b) Areal and vertical contraction at '
                                 'the last time against the analytical solution, which predicts the two to be equal. '
                                 '(c) RMS errors in pressure, areal strain and vertical strain against mesh spacing, '
                                 'with a first-order slope for reference.',
                      'study': ['python3 benchmarks/poromechanics/slab_refinement.py dynearthsol3d OUT '
                                '--resolutions 20 10 5']},
 'regional': {'run': [('dynearthsol3d input.cfg',
                       'Runs the 3D quarter basin on a 25 m hexahedral mesh for two years (416 steps): pumping, '
                       'rest, pumping, rest. DES writes frames basin.* every 26 steps (17 times) and a restart '
                       'checkpoint every fourth frame. No mesh or restart file is needed; the run is long on a '
                       'laptop.')],
              'keys': [('sim.max_steps',
                        'ends the run after 416 steps of control.fixed_dt, i.e. two years; sim.max_time_in_yr is set '
                        'longer and is not reached'),
                       ('control.fixed_dt',
                        'time step of about 1.76 days, so each six-month season takes 104 steps and a rate change '
                        'falls on a step boundary'),
                       ('sim.output_step_interval', 'writes a frame every 26 steps, about every month and a half'),
                       ('control.has_hydraulic_diffusion', 'turns on Darcy flow of the pore pressure'),
                       ('control.has_PT',
                        'iterates each step to mechanical equilibrium, so the basin deforms quasi-statically with '
                        'the pressure'),
                       ('control.gravity', 'no gravity: pressure and subsidence are changes from the initial state'),
                       ('bc.hbc_x1', 'outer side x = 3 km held at the pressure in bc.hbc_pressure (zero change)'),
                       ('bc.hbc_y1', 'outer side y = 3 km held at the pressure in bc.hbc_pressure (zero change)'),
                       ('bc.well_rate',
                        'rate of the quarter well field at schedule ratio 1: one quarter of the full Q'),
                       ('bc.well_z_bottom', 'screens the well over the whole 50 m aquifer'),
                       ('bc.well_schedule_time_in_yr',
                        'times of the pumping schedule; a repeated time is an instant on/off switch'),
                       ('bc.well_schedule_ratio',
                        'rate multiplier at each schedule time: 1 while pumping, 0 while resting'),
                       ('mat.bulk_modulus', 'drained bulk modulus of the soft basin sediment'),
                       ('mat.hydraulic_perm',
                        'aquifer permeability; with the storage it sets how far the drawdown spreads in a season')],
              'fields': [('pore pressure',
                          'nodal pressure change, compared with the closed-basin reference for 250–1500 m; the '
                          'pressure errors in panel (c)'),
                         ('strain',
                          'the zz component times the 50 m thickness is the subsidence (positive down); panels (a), '
                          '(b) and the subsidence errors'),
                         ('coordinate',
                          'node positions, used for the distance of nodes and element centres from the well field')],
              'caption': '(a) Subsidence against distance from the well field at the ends of the four seasons: lines '
                         'are the closed-basin averaged-aquifer reference, symbols DES elements along the x axis; '
                         'the shaded band is the compared range. (b) Subsidence in time at three observation wells, '
                         'reference lines and DES symbols; dotted lines mark the season changes. (c) Largest '
                         'pressure and subsidence errors after the first frame against mesh spacing: open markers '
                         'with the longer time step, solid markers with four times as many steps.',
              'study': ['python3 benchmarks/poromechanics/regional_subsidence.py dynearthsol3d OUT-mesh --spacings 100 '
                        '50 25',
                        'python3 benchmarks/poromechanics/regional_subsidence.py dynearthsol3d OUT-steps416 --spacings 50 '
                        '--steps 416',
                        'python3 benchmarks/poromechanics/regional_subsidence.py dynearthsol3d OUT-h25-s416 --spacings 25 '
                        '--steps 416']},
 'embankment': {'run': [('dynearthsol2d seed.cfg',
                         'Optional: meshes slope.poly at 0.25 m and writes the unstressed, gravity-free state as '
                         'frame seed 0 (zero steps).'),
                        ('dynearthsol2d gravity.cfg',
                         'Optional: restarts from seed frame 0 and applies gravity to the elastic slope in two '
                         'equilibrium (PT) steps; writes gravity frames with checkpoints. The shipped '
                         'gravity.chkpt.000002, gravity.save.000002 and gravity.info are this result, so both '
                         'optional stages can be skipped.'),
                        ('dynearthsol2d input.cfg',
                         'One strength-reduction trial at F = 1.36: restarts from gravity frame 2 as Mohr–Coulomb '
                         'with c/F and atan(tan φ/F), runs explicitly for 2000 s of mass-scaled time and saves a '
                         'frame every 10 s as srf1.360.*. At this F the slope keeps sliding.'),
                        ('cp input.cfg srf1.300.cfg',
                         'To try another factor F (here 1.30), edit the copy: set sim.modelname to a new name '
                         '(srf1.300) so the outputs do not overwrite each other, cohesion0 and cohesion1 to 10000/F '
                         'Pa, and friction_angle0 and friction_angle1 to atan(tan 20°/F) in degrees. Everything '
                         'else, including the restart from gravity frame 2, stays.'),
                        ('dynearthsol2d srf1.300.cfg',
                         'Runs that trial. Judge it from the velocity history: if the peak nodal speed over the last '
                         'frames is below that about three quarters through the run, the slope is coming to rest and '
                         'F is below the factor of safety; otherwise it is sliding. The factor of safety lies '
                         'between the largest F at rest and the smallest F sliding.')],
                'keys': [('sim.is_restarting', 'starts from a saved state instead of a new mesh'),
                         ('sim.restarting_from_modelname',
                          'restarts from the elastic gravity state written by gravity.cfg'),
                         ('sim.restarting_from_frame', 'uses gravity frame 2, the slope at rest under gravity'),
                         ('mat.rheology_type',
                          'Mohr–Coulomb elasto-plasticity for the trial (the seed and gravity stages are elastic)'),
                         ('mat.cohesion0', 'reduced cohesion c/F; cohesion1 is the same, so there is no softening'),
                         ('mat.friction_angle0',
                          'reduced friction angle atan(tan φ/F); friction_angle1 is the same, so there is no '
                          'softening'),
                         ('mat.is_plane_strain', 'plane-strain plasticity, as in the 2D slope-stability problem'),
                         ('mat.rho0', 'density giving a unit weight of 20 kN/m³ with control.gravity'),
                         ('control.has_moving_mesh',
                          'keeps the original geometry (small strain): the factor of safety refers to the undeformed '
                          'slope'),
                         ('control.has_poroelastic_pressure_feedback', 'no pore-pressure coupling: the slope is dry'),
                         ('bc.vbc_z0', 'fixes the base in both directions; the sides are rollers'),
                         ('control.characteristic_speed',
                          'speed used to size the stable explicit time step, since no boundary moves'),
                         ('sim.max_time_in_yr',
                          'runs each trial for 2000 s of mass-scaled time, long enough to tell rest from sliding'),
                         ('sim.output_time_interval_in_yr',
                          'saves a frame every 10 s, the velocity history used to judge the trial')],
                'fields': [('velocity',
                            'nodal velocity at each saved frame; integrated in time (trapezoid rule) it gives the '
                            'displacement in panels (a)–(c), and its peak over the mesh decides whether a trial is '
                            'coming to rest or sliding'),
                           ('plastic strain',
                            'shows the shear band of the mechanism and how much of the slope yields'),
                           ('coordinate', 'node positions of the fixed mesh, for plotting the mechanism')],
                'caption': "(a) Displacement on the finest mesh in the first trial that keeps sliding, with Bishop's "
                           'critical circle dashed: the mechanism passes through the toe. (b) Largest displacement '
                           'against mass-scaled time for each strength-reduction factor on the finest mesh; dashed '
                           'lines are trials still sliding at the end. (c) Final largest displacement against the '
                           "factor F on each mesh, open circles at rest and filled squares sliding, with Bishop's "
                           'factor of safety as a vertical dashed line.',
                'study': ['python3 benchmarks/poromechanics/bishop_reference.py',
                          'python3 benchmarks/poromechanics/strength_reduction_slope.py dynearthsol2d OUT '
                          '--spacings 1 0.5 0.25']},
 'drawdown': {'run': [('dynearthsol2d input.cfg',
                       'Restarts from the shipped filled-reservoir checkpoint, frame 100. It removes the water '
                       'load, drains the face and toe platform at p = 0 and seals the crest and base, then runs '
                       '20000 s of drainage. Writes drawdown frames 101–180, one every 250 s, each with a checkpoint '
                       'for the trials.'),
                      ('cp input.cfg trial.cfg',
                       'A strength-reduction trial at one time after drawdown, here F = 0.85 at 250 s. In the copy '
                       'set: sim.modelname to a new name (srf250-0.850); restarting_from_modelname = drawdown and '
                       'restarting_from_frame = 101 (250 s; frames 104, 112, 124, 148 and 180 are 1000, 3000, 6000, '
                       "12000 and 20000 s); max_time_in_yr to that frame's time in drawdown.info (third column, s) "
                       'plus 2000 s, divided by 31556926 (DES’s year, 365.2422 days); output_time_interval_in_yr = 3.168876453507857e-07 (10 s); '
                       'checkpoint_frame_interval = 100000; rheology_type = elasto-plastic; hydraulic_perm = [1e-25] '
                       'and fluid_bulk_modulus = [1000.0], which hold the pore pressure during the trial; cohesion0 '
                       'and cohesion1 = 10000/F Pa; friction_angle0 and friction_angle1 = atan(tan 20°/F) in '
                       'degrees.'),
                      ('dynearthsol2d trial.cfg',
                       'Runs the trial for 2000 s of mass-scaled time. Judge it from the velocity history: if the '
                       'peak nodal speed over the last frames is below that about three quarters through the run, '
                       'the slope is coming to rest and F is below the factor of safety; otherwise it is sliding. '
                       'Repeat with other F to bracket the factor of safety at that time.')],
              'keys': [('sim.restarting_from_modelname',
                        'starts from the shipped filled-reservoir state'),
                       ('sim.restarting_from_frame',
                        'uses reservoir frame 100, the slope with hydrostatic pore pressure under the full '
                        'reservoir'),
                       ('sim.max_time_in_yr',
                        'end time counted from the start of the reservoir stage: 25000 s of filling plus 20000 s of '
                        'drainage'),
                       ('sim.checkpoint_frame_interval',
                        'writes a restart checkpoint with every frame, so a strength-reduction trial can start at '
                        'any saved time'),
                       ('control.has_hydraulic_diffusion',
                        'turns on Darcy flow, so the slope drains through the face after drawdown'),
                       ('bc.has_water_loading',
                        "removes the reservoir's pressure on the face and toe platform"),
                       ('bc.hbc_z1',
                        'drains the face and toe platform: pore pressure held at 0 (the crest is marked as a '
                        'separate boundary in slope.poly and stays sealed)'),
                       ('mat.retention_model',
                        'Gardner exponential retention, so the slope desaturates where the pressure turns negative'),
                       ('mat.retention_alpha', 'Gardner parameter: how quickly saturation falls with suction'),
                       ('mat.residual_saturation', 'saturation left at high suction'),
                       ('mat.hydraulic_perm',
                        'permeability; drainage takes thousands of seconds, slow against the mass-scaled mechanics'),
                       ('mat.bulk_modulus_s', 'incompressible grains (with Biot coefficient 1)'),
                       ('mat.rho0', 'grain density chosen so the saturated unit weight is 20 kN/m³'),
                       ('mat.rheology_type',
                        'elastic during drainage; each trial switches to Mohr–Coulomb elasto-plasticity')],
              'fields': [('pore pressure',
                          'pressure (and suction above the zero contour) at each saved time after drawdown; panel '
                          "(a) and the input to Bishop's method"),
                         ('velocity',
                          'in a trial, integrated in time it gives the displacement of the mechanism in panel (c); '
                          'its peak over the mesh decides rest or sliding'),
                         ('coordinate',
                          'node positions of the fixed mesh, for the maps and for interpolating the pore pressure on '
                          "Bishop's slices")],
              'caption': '(a) Pore pressure at the first comparison time after drawdown, with the zero-pressure '
                         "contour and Bishop's critical circle on the same pore pressure. (b) Factor of safety "
                         "against time after drawdown: Bishop's method on the DES pore pressure (line and circles) "
                         'and the DES strength-reduction bracket, from the last F at rest to the first F sliding '
                         '(bars). (c) Displacement in the first sliding trial at the first comparison time, with '
                         "Bishop's circle.",
              'study': ['python3 benchmarks/poromechanics/drawdown_slope.py dynearthsol2d OUT --spacings .5 --workers 12 --flow-threads 12']},
 'slope': {'run': [('dynearthsol2d seed.cfg',
                    'Optional. Meshes the slope from slope.poly and writes the checkpoint seed at frame 0 without '
                    'gravity and without stress (zero steps). The shipped seed.chkpt.000000, seed.save.000000 and '
                    'seed.info are this output, so you can skip this step.'),
                   ('dynearthsol2d input.cfg',
                    'Restarts from seed frame 0 with gravity switched on, so the first step finds the gravitational '
                    'equilibrium of the slope; after 20 s the bedrock pressure starts to rise. Saves every 5 s step '
                    'as slope.*. The run ends with exit code 54 (mechanical solve did not converge) once the '
                    'accelerating basal slip outruns PT, before max_steps is reached; that is the expected end of '
                    'this example.')],
           'keys': [('sim.is_restarting',
                     'starts from the unstressed seed checkpoint instead of building a new initial state, so gravity '
                     'is applied as a load'),
                    ('control.gravity', 'gravity is on in this run; seed.cfg has it at 0'),
                    ('mesh.meshing_option',
                     'reads the slope outline, and the soil and strong-end regions, from slope.poly'),
                    ('mat.num_materials',
                     'two materials: the soil (0) and the strong 4 m end blocks (1), assigned by the regions of '
                     'slope.poly'),
                    ('mat.rheology_type', 'Mohr–Coulomb elasto-plasticity on the effective stress'),
                    ('mat.cohesion0',
                     'soil cohesion and the high cohesion of the strong ends, which keep the corners elastic; '
                     'cohesion1 is the same, so there is no softening'),
                    ('mat.friction_angle0',
                     'friction angle of the soil; dilation stays at its default of zero (non-associated flow)'),
                    ('mat.shear_modulus',
                     "with bulk_modulus sets Poisson's ratio 0.45, which fixes the slope-parallel stress K0 in the "
                     'pointwise estimate'),
                    ('mat.hydraulic_perm',
                     'high permeability, so the seepage stays quasi-steady while the pressure rises'),
                    ('bc.vbc_z0', 'fixes the bedrock (base): normal and tangential velocity zero, no slip'),
                    ('bc.hbc_z0', 'bedrock pressure prescribed from hbc_pressure: water is supplied from below'),
                    ('bc.hbc_z1', 'drained surface: pore pressure held at 0'),
                    ('bc.hbc_pressure', 'the bedrock pressure at schedule ratio 1 (fifth entry, z0)'),
                    ('bc.hbc_schedule_ratio',
                     'with hbc_schedule_time_in_yr: no bedrock pressure for the first 20 s, then a linear rise to '
                     'the full value at the end time')],
           'fields': [('plastic strain',
                       'panel (a) maps it; its mean over the basal 0.6 m in the middle third of the slope is the '
                       'slip measure in (b), and its first nonzero value marks the onset of basal yield in (c)'),
                      ('pore pressure',
                       'averaged over the bedrock nodes in the middle third: the base pressure on the horizontal '
                       'axis of (b) and the vertical axis of (c)'),
                      ('velocity',
                       'downslope velocity of the surface in the middle third, recorded in summary.json to show the '
                       'slip accelerating')],
           'caption': '(a) Plastic strain, on a logarithmic scale, in the middle of the slope at the last time step '
                      'for which PT finds equilibrium on the finest mesh; plastic flow is concentrated along the '
                      'bedrock. (b) Mean plastic strain of the basal 0.6 m against the bedrock pore pressure for '
                      'each mesh (solid) and for the run with twice as many steps (dashed); crosses mark the last PT '
                      'equilibrium, vertical lines the pointwise-yield, Davis non-associated and classical '
                      'factor-of-safety estimates. (c) Bedrock pressure at the onset of basal yield, at two fixed '
                      'basal plastic strains and at the last PT equilibrium, against mesh spacing; filled symbols '
                      'are the run with twice as many steps, horizontal lines the three estimates.',
           'study': ['python3 benchmarks/poromechanics/slope_application.py dynearthsol2d OUT --spacings 0.5 0.25 '
                     '0.125',
                     'python3 benchmarks/poromechanics/slope_application.py dynearthsol2d OUT-808 --spacings 0.25 '
                     '--steps 808']},
 'rainfall': {'run': [('dynearthsol2d seed.cfg',
                       'Optional. Meshes the slope from slope.poly and writes the checkpoint seed at frame 0: no '
                       'gravity, no stress, uniform initial suction (zero steps). The shipped seed.chkpt.000000, '
                       'seed.save.000000 and seed.info are this output, so you can skip this step.'),
                      ('dynearthsol2d input.cfg',
                       'Restarts from seed frame 0 with gravity on and rain entering through the surface; 125 s '
                       'steps, saved every fourth step as slope.*. The run ends with exit code 54 (mechanical solve '
                       'did not converge) once the accelerating basal slip outruns PT, well before max_steps; that '
                       'is the expected end of this example.')],
              'keys': [('sim.is_restarting',
                        'starts from the seed checkpoint, so gravity is applied as a load to the unstressed slope'),
                       ('control.gravity', 'gravity is on in this run; seed.cfg has it at 0'),
                       ('mat.retention_model',
                        'Gardner exponential retention: the soil is partly saturated where the pore pressure is '
                        'negative, with Bishop effective stress chi = S_e'),
                       ('mat.retention_alpha',
                        'Gardner alpha: how quickly saturation and relative permeability fall with suction'),
                       ('ic.pore_pressure_option', 'uniform initial pore pressure, given by pore_pressure_value'),
                       ('ic.pore_pressure_value', 'initial suction of 5 kPa (negative pressure) everywhere'),
                       ('bc.hbc_z1', 'the surface takes a prescribed Darcy flux (hbc_flux): the rain'),
                       ('bc.hbc_flux',
                        'rain rate: the z1 entry is the outward flux, negative for inflow (1e-5 m/s, 36 mm/h)'),
                       ('bc.hbc_x1',
                        'the downslope end (toe) drains: pore pressure held at 0; the bedrock and crest end keep the '
                        'default, sealed'),
                       ('bc.vbc_z0', 'fixes the bedrock (base): normal and tangential velocity zero, no slip'),
                       ('mat.cohesion0',
                        'soil cohesion (weaker than the pressure-driven slope) and the high cohesion of the strong '
                        'end blocks'),
                       ('mat.friction_angle0', 'friction angle of the soil; dilation stays zero'),
                       ('mat.num_materials',
                        'soil (0) and strong end blocks (1), assigned by the regions of slope.poly'),
                       ('control.fixed_dt', '125 s time step')],
              'fields': [('pore pressure',
                          'panel (a) maps it with the p = 0 contour as the perched water table; along the bedrock it '
                          'gives the base pressure of (b) and the profiles of (c); its maximum on the surface is '
                          'checked to stay negative, since DES has no ponding'),
                         ('plastic strain',
                          'mean over the basal 0.6 m in the middle third of the slope: the slip measure of (b) and '
                          'the onset and failure times'),
                         ('velocity', 'downslope surface velocity in the middle third, recorded in summary.json')],
              'caption': '(a) Pore pressure in the middle of the slope at the last time step for which PT finds '
                         'equilibrium; the black line p = 0 is the perched water table, with suction above it. (b) '
                         'Mean plastic strain of the basal 0.6 m against the bedrock pore pressure under rain '
                         '(solid) and for the same soil loaded by pressure from the bedrock (dashed), on each mesh, '
                         'with the pointwise-yield and Davis non-associated estimates. (c) Pore pressure along the '
                         'bedrock at five times; the shaded ends are the strong blocks.',
              'study': ['python3 benchmarks/poromechanics/rainfall_slope.py dynearthsol2d OUT --spacings 0.5 0.25',
                        'python3 benchmarks/poromechanics/slope_application.py dynearthsol2d OUT-driven --spacings '
                        '0.5 0.25 0.125 --cohesion 1000 --friction 30 --elastic-gravity']},
 'rudnicki': {'run': [('dynearthsol3d input.cfg',
                       '3D. DES builds the box mesh itself, refined near the source (no mesh file needed), injects '
                       'at the corner node and runs 48 steps of 5 s to 240 s, saving every fourth step as '
                       'rudnicki.*.')],
              'keys': [('injection.enabled', 'turns on scheduled point sources'),
                       ('injection.points_x',
                        'with points_y and points_z: the source at the origin, the corner node shared by the three '
                        'symmetry planes'),
                       ('injection.rate',
                        'injected volume rate: Q/8 for one octant of the full-space source, in m^3/s'),
                       ('mesh.meshing_option', 'box mesh with a refined zone'),
                       ('mesh.refined_zonex',
                        'with refined_zoney and refined_zonez: the refined zone is the first 10 m from the source '
                        'along each axis (refined_zonez is measured downward from the top, z = 0)'),
                       ('mesh.largest_size',
                        'lets the elements grow large away from the refined zone, toward the far faces of the 60 m '
                        'box'),
                       ('bc.vbc_z0',
                        'with vbc_z1, vbc_y0 and vbc_y1 (and the default vbc_x0, vbc_x1): zero normal displacement '
                        'on every face; symmetry on x = 0, y = 0, z = 0'),
                       ('bc.hbc_x1',
                        'with hbc_y1 and hbc_z0: the far faces hold p = 0; the symmetry planes keep the default, '
                        'sealed'),
                       ('control.gravity',
                        'no gravity: pressure and stress are changes caused by the injection only'),
                       ('control.has_PT',
                        'PT keeps mechanical equilibrium at every step (quasi-static poroelasticity)'),
                       ('mat.hydraulic_perm',
                        'intrinsic permeability; with the storage it sets the coupled diffusivity c'),
                       ('mat.biot_coeff', 'Biot coefficient 0.8: couples swelling to pressure'),
                       ('mat.bulk_modulus_s', 'grain bulk modulus, entering the storage 1/M'),
                       ('control.fixed_dt', '5 s time step')],
              'fields': [('pore pressure',
                          'nodal pressure against distance from the source, averaged in radial bins, in (a); RMS '
                          'error over 2.5 ≤ r ≤ 8 m in (c)'),
                         ('strain',
                          'element strain tensor, rotated to radial and hoop components about the source, in (b); '
                          'RMS errors in (c)')],
              'caption': '(a) Pore pressure against distance from the source at four times: lines are the Rudnicki '
                         'solution, symbols DES nodal pressures averaged in radial bins; the shaded band is the '
                         'compared range 2.5 ≤ r ≤ 8 m. (b) Radial and hoop strain at the last saved time, DES '
                         'element values binned by radius, against the solution. (c) RMS error of pressure, radial '
                         'and hoop strain after the front reaches the band, against mesh spacing, with a first-order '
                         'slope for reference.',
              'study': ['python3 tests/functional/hydromechanics/check_rudnicki.py dynearthsol3d OUT --spacings 1 '
                        '0.7 0.5 --steps 48',
                        'python3 tests/functional/hydromechanics/check_rudnicki.py dynearthsol3d OUT-96 --spacings '
                        '0.7 --steps 96']},
 'undrained_wave': {'run': [('dynearthsol3d input.cfg',
                             '3D. Restarts from the shipped checkpoint seed, frame 0: the column mesh with zero '
                             'stress, strain and pressure and the initial sine velocity (listed in '
                             'initial-velocity.csv). Runs 8192 explicit steps, one wave period, saving every 128th '
                             'step as wave.* and writing a checkpoint at every saved frame.')],
                    'keys': [('sim.is_restarting',
                              'starts from the seed checkpoint, which holds the initial velocity; an input file '
                              'cannot set an initial velocity, so undrained_wave.py writes it into the seed frame'),
                             ('sim.restarting_from_modelname', 'the checkpoint to start from: seed (shipped)'),
                             ('control.is_quasi_static',
                              'dynamic run with physical (unscaled) mass, so inertial_scaling has no effect'),
                             ('control.damping_option', 'no damping: the wave keeps its amplitude'),
                             ('control.fixed_dt', 'one wave period divided by 8192'),
                             ('sim.max_steps', '8192 steps: exactly one period'),
                             ('sim.output_step_interval', 'a frame every 128 steps: 65 saved times over the period'),
                             ('mat.hydraulic_perm', 'zero permeability: no flow, so each point stays undrained'),
                             ('control.has_hydraulic_diffusion',
                              'solves the pore-pressure equation; with zero permeability the pressure follows the '
                              'volume change only'),
                             ('bc.vbc_z0',
                              'with vbc_z1: the ends cannot move vertically; vbc_y0, vbc_y1 and the default vbc_x0, '
                              'vbc_x1 fix the side-normal motion'),
                             ('control.gravity', 'no gravity: the column is unstressed apart from the wave'),
                             ('mesh.meshing_elem_shape', 'regular structured mesh of the 0.025 × 0.025 × 1 m column'),
                             ('mat.bulk_modulus',
                              'with shear_modulus, porosity and the fluid and grain moduli, sets the undrained '
                              'modulus and the wave speed'),
                             ('control.has_moving_mesh', 'fixed geometry: small-strain wave')],
                    'fields': [('velocity',
                                'vertical component at the node nearest mid-height, normalized by the initial '
                                'amplitude V, in (a); maximum error over all nodes in the result'),
                               ('pore pressure',
                                'at the base node, normalized by the pressure amplitude, in (b); maximum error over '
                                'all nodes in the result'),
                               ('strain',
                                'volumetric strain, combined with the pressure to check that the local fluid content '
                                'stays zero')],
                    'caption': '(a) Vertical velocity at mid-height over one period, normalized by the initial '
                               'amplitude V: line, eigenmode solution; symbols, DES saved frames. (b) Pore pressure '
                               'at the sealed base over the same period, normalized by the pressure amplitude.',
                    'study': ['python3 benchmarks/poromechanics/undrained_wave.py dynearthsol3d OUT --spacing .0125 '
                              '--steps 1024 2048 4096 8192 --samples 64 --restart-check',
                              'python3 benchmarks/poromechanics/undrained_wave.py dynearthsol3d OUT-h005 --spacing '
                              '.05 --steps 4096 --samples 32',
                              'python3 benchmarks/poromechanics/undrained_wave.py dynearthsol3d OUT-h0025 --spacing '
                              '.025 --steps 4096 --samples 32']},
 'strength': {'run': [('dynearthsol3d input.cfg',
                       'Runs the 3D unit cell for 64 steps of 1 s (the input is 3D: six initial-stress components) '
                       'and writes cell.save.* and cell.chkpt.* every 4 steps, with cell.info listing the frames.')],
              'keys': [('mat.rheology_type',
                        'elasto-plastic: Mohr–Coulomb plasticity, whose strength DES evaluates with the effective '
                        'stress'),
                       ('mat.cohesion0', 'cohesion c = 100 kPa; equal to cohesion1, so there is no strain weakening'),
                       ('mat.cohesion1', 'same as cohesion0: no weakening'),
                       ('mat.friction_angle1',
                        '30°, equal to friction_angle0, so the friction angle stays fixed (dilation stays at its '
                        'default of zero)'),
                       ('mat.bulk_modulus', 'K = 1 GPa'),
                       ('mat.shear_modulus', 'G = 1 GPa; the elastic branch of the load path is t = 2Ge'),
                       ('ic.initial_stress_option',
                        '1: start from the homogeneous total stress given in initial_stress'),
                       ('ic.initial_stress', 'isotropic total stress −C with C = 1 MPa (compression negative)'),
                       ('ic.pore_pressure_option', '2: uniform initial pore pressure'),
                       ('ic.pore_pressure_value',
                        'pore pressure p = 0.8 MPa; with Biot coefficient 1 it lowers the effective confinement to '
                        '0.2 MPa'),
                       ('control.has_poroelastic_pressure_feedback',
                        'off: deformation does not change the pore pressure, so p stays fixed'),
                       ('bc.vbc_val_x1',
                        'the right face moves outward, extending the cell in x at strain rate e-dot'),
                       ('bc.vbc_val_z1',
                        'the top moves down at the same rate, shortening the cell in z; the y faces are held '
                        '(vbc_y0, vbc_y1 = 1), so y strain is zero'),
                       ('control.gravity', 'no gravity: stress and pressure are uniform')],
              'fields': [('strain',
                          'element strain; its xx component, averaged over the cell, is the imposed normal strain e '
                          'on the horizontal axis'),
                         ('stress',
                          'element total stress; half of (σxx − σzz), averaged over the cell, is the principal shear '
                          'load t on the vertical axis')],
              'caption': 'Principal shear load t (half the difference of the extreme principal stresses) against the '
                         'imposed normal strain e, for the 3D elasto-plastic cell at 1 MPa total confinement and 0.8 '
                         'MPa pore pressure. Line: the elastic–perfectly plastic analytical path, which yields at '
                         'the effective-stress load t_y. Symbols: DES at each saved frame.',
              'study': ['python3 tests/functional/hydromechanics/check_pore_strength.py dynearthsol2d OUT/2d --ndims '
                        '2',
                        'python3 tests/functional/hydromechanics/check_pore_strength.py dynearthsol3d OUT/3d --ndims '
                        '3'],
              'replace': {'meaning': [('This isolates pressure-dependent skeleton strength and its plastic strain. '
                                       'Constant-friction EP-RSF is checked; evolving RSF is checked on the '
                                       'rate-and-state page; flow, dilation, tensile returns and slope/injection '
                                       'applications remain separate requirements.',
                                       'This isolates the pressure-dependent strength of the skeleton and its '
                                       'plastic strain. EP-RSF is checked at constant friction here and with '
                                       'evolving friction on the rate-and-state page. Flow, dilation, tensile '
                                       'failure and slope or injection applications are not tested here.')]}},
 'rsf': {'run': [('dynearthsol2d input.cfg',
                  'Loads the 2D cell at the slow rate for 512 steps of 0.125 s (64 s) and writes cell.save.* and '
                  'cell.chkpt.* at every step; the last checkpoint (frame 512) is the starting point of the second '
                  'stage.'),
                 ('dynearthsol2d step.cfg',
                  'Restarts from cell frame 512 with ten times the boundary velocity and runs to step 1024 (another '
                  '64 s), writing step.save.*; no checkpoint is shipped, so the first stage must be run first.')],
         'keys': [('mat.rheology_type',
                   'elasto-plastic-rate-state-friction: Mohr–Coulomb plasticity whose friction coefficient follows '
                   'rate-and-state friction'),
                  ('mat.state_var_model', '1: aging law for the state variable θ'),
                  ('mat.direct_a',
                   'direct-effect parameter a = 0.01: the instant friction jump at the rate step is a ln 10'),
                  ('mat.evolution_b',
                   'evolution parameter b = 0.015 > a, so the material is velocity weakening and the new steady '
                   'friction is lower by (b − a) ln 10'),
                  ('mat.characteristic_distance',
                   'slip distance D_c = 1 mm over which friction relaxes to the new steady state'),
                  ('mat.characteristic_velocity',
                   'reference slip rate V0, set equal to the initial slip rate V1 so the cell starts at steady state '
                   'with μ = tan 30°'),
                  ('control.rsf_slip_rate_projection_option',
                   "1: slip rate V = 2w ε̇_II from the total deviatoric strain rate (w is the element's minimum "
                   'altitude)'),
                  ('monitor.output_state_variable',
                   'writes the state variable and the dynamic friction coefficient to every frame'),
                  ('bc.vbc_val_x1', 'extension rate of the cell (right face velocity); step.cfg multiplies it by 10'),
                  ('bc.vbc_val_z1', 'equal shortening rate on the top; step.cfg multiplies it by 10'),
                  ('ic.pore_pressure_value', 'fixed pore pressure 0.8 MPa (uniform, pore_pressure_option 2)'),
                  ('ic.initial_stress', 'isotropic total stress of −1 MPa'),
                  ('control.fixed_dt', 'time step 0.125 s, the finest step of the study'),
                  ('sim.checkpoint_frame_interval',
                   '1: a checkpoint at every frame, so the rate step can restart from the last one')],
         'fields': [('dynamic friction coefficient',
                     'element friction coefficient μ, plotted across the rate step in panel (a)'),
                    ('friction state variable',
                     'element state θ; compared with the exact forward-Euler recursion at every step'),
                    ('stress',
                     'element total stress; half the principal-stress difference is the shear stress at yield in '
                     'panel (b)')],
         'caption': '(a) Friction coefficient across the tenfold step in loading rate, against time from the step: '
                    'solid line, continuum solution; thin stepped line, the forward-Euler recursion at the same time '
                    'step; symbols, DES at the finest time step. (b) Principal shear stress at yield, which follows '
                    'the friction coefficient through the Mohr–Coulomb load. (c) Largest friction error against the '
                    'continuum solution after the step, for each time step, with a first-order reference slope.',
         'study': ['python3 tests/functional/hydromechanics/check_rsf_step.py dynearthsol2d OUT --refinements 1 2 4 '
                   '8']},
 'plastic_content': {'run': [('dynearthsol2d input.cfg',
                              'Compresses the 2D column for 100 steps of 0.01 s (1 s) and writes column.save.* and '
                              'column.chkpt.* at every step, with column.info listing the frames.')],
                     'keys': [('mat.rheology_type', 'elasto-plastic: Mohr–Coulomb plasticity'),
                              ('mat.cohesion0',
                               'cohesion 100 Pa (equal to cohesion1: no weakening), small so the column yields '
                               'early'),
                              ('mat.friction_angle0', 'friction angle 20° (friction_angle1 is the same)'),
                              ('mat.dilation_angle0',
                               'dilation angle 10° (dilation_angle1 is the same): yielding makes the skeleton '
                               'dilate'),
                              ('bc.vbc_val_z1',
                               'the top moves down at 1 mm/s; base and sides are held normally, so the volume loss '
                               'is fixed by this velocity'),
                              ('mat.hydraulic_perm',
                               'zero permeability, and every hbc is 0 (no flow): the column is sealed'),
                              ('mat.fluid_bulk_modulus',
                               'K_f = 1 MPa, a soft fluid, so the pressure stays small (4 kPa at 1 s)'),
                              ('mat.porosity', 'φ = 0.25; with K_f it gives the storage φ/K_f'),
                              ('mat.bulk_modulus_s',
                               'grain modulus 1e20 Pa: incompressible grains, Biot coefficient 1'),
                              ('mat.bulk_modulus', 'K = 20 MPa'),
                              ('mat.shear_modulus', 'G = 10 MPa'),
                              ('control.has_PT',
                               'PT solves the mechanics to equilibrium every step; the pressure increment it carries '
                               'to the next step is counted separately'),
                              ('control.fixed_dt', '0.01 s step; 100 steps to 1 s'),
                              ('mat.is_plane_strain', '2D plane strain, including the out-of-plane plastic strain')],
                     'fields': [('pore pressure',
                                 'nodal pressure; its volume average is the mean pressure in panel (a) and gives the '
                                 'pressure part of the stored water'),
                                ('strain',
                                 'element strain; its volumetric part integrated over the mesh is the imposed volume '
                                 'loss'),
                                ('plastic strain', 'shows that the column has yielded before the balance is checked'),
                                ('fluid balance',
                                 'per-step water ledger; its boundary entry stays zero in this sealed column')],
                     'caption': '(a) Mean pore pressure of the sealed column against time, through yield: line, the '
                                'undrained reference set by the imposed volume loss; symbols, DES. (b) Volume '
                                'balance: line, the imposed pore-volume loss; symbols, DES water stored by pressure '
                                'plus the pressure increment still pending in the PT step.',
                     'study': ['python3 tests/functional/hydromechanics/check_plastic_content.py dynearthsol2d '
                               'OUT/2d',
                               'python3 tests/functional/hydromechanics/check_plastic_content.py dynearthsol3d '
                               'OUT/3d --ndims 3'],
                     'replace': {'result': [('The corrected source conserves fluid content plus the pending PT '
                                             'increment to round-off. The parent drift was 1.24% in the selected 2D '
                                             'case and 4.49% in the 3D dilation control.',
                                             'Through yield, the stored fluid content plus the pressure increment '
                                             'still pending in the PT step equals the imposed volume loss to '
                                             'round-off, and the mean pressure follows its reference recurrence to '
                                             'round-off.')],
                                 'meaning': [('Forty-eight primary cases cover both dimensions and containers, '
                                              'negative/zero/positive dilation, hydraulic substeps, open drainage, '
                                              'Gardner suction, evolving-pressure EP-RSF and thermal expansion. Old '
                                              'checkpoints preserve earlier offsets. Conservative remapping, '
                                              'finite-grain plastic porosity and general energy remain open.',
                                              'This checks that plastic yielding neither creates nor destroys pore '
                                              'water: the fluid source uses the accepted volume change, not an '
                                              'elastic estimate. Forty-eight cases cover 2D and 3D, negative, zero '
                                              'and positive dilation, hydraulic substeps, open drainage, Gardner '
                                              'suction, EP-RSF with evolving pressure and thermal expansion. '
                                              'Remeshing, porosity change from plastic strain and the general energy '
                                              'balance are not covered.')]}},
 'pressure_yield': {'run': [('dynearthsol2d input.cfg',
                             'Runs the 2D layer (one 0.025 m element wide, 2 m deep) for 6400 steps of 0.3125 s '
                             '(2000 s) and saves every step to layer.save.*, with layer.info listing the frames; it '
                             'first brings the layer to gravitational equilibrium.')],
                    'keys': [('mat.rheology_type',
                              'elasto-plastic: elastic–perfectly plastic Mohr–Coulomb (dilation stays at its default '
                              'of zero)'),
                             ('mat.cohesion0', 'cohesion c = 10 kPa (cohesion1 equal: no weakening)'),
                             ('mat.friction_angle0', 'friction angle φ = 20° (friction_angle1 equal)'),
                             ('bc.stress_bc_z1', 'type 3: vertical traction on the top surface'),
                             ('bc.stress_val_z1',
                              'the surcharge q = 20 kPa (compression negative), which keeps the yielded top in '
                              'compression'),
                             ('bc.hbc_z1', 'type 2: the top pore pressure is prescribed from hbc_pressure'),
                             ('bc.hbc_pressure',
                              'the z1 (top) entry is the final top pore pressure, 95% of the pressure that would '
                              'bring the surface to the apex of the yield cone'),
                             ('bc.num_hbc_schedule', 'two schedule points for the prescribed pressure'),
                             ('bc.hbc_schedule_time_in_yr',
                              '0 and 2000 s: the top pressure ramps linearly over the run'),
                             ('bc.hbc_schedule_ratio', '0 to 1: the ramp starts from the hydrostatic state'),
                             ('mat.hydraulic_perm',
                              'high permeability, so the whole layer follows the top pressure rise almost at once'),
                             ('control.gravity',
                              '9.81 m/s²: lithostatic stress and hydrostatic pressure set the depth dependence of '
                              'the yield pressure'),
                             ('ic.has_body_force_adjustment',
                              'brings the initial stress into equilibrium with gravity before loading'),
                             ('control.fixed_dt', '0.3125 s; 6400 steps reach 2000 s')],
                    'fields': [('pore pressure',
                                'element-averaged pressure; its rise since the start, read at the first saved step '
                                "with plastic strain, is each element's yield pressure"),
                               ('plastic strain', 'marks the first saved step at which an element has yielded'),
                               ('stress',
                                'element total stress; with the pore pressure it gives the effective-stress paths in '
                                'panel (b)')],
                    'caption': '(a) Depth of the yield front against the pressure rise: line, the exact onset depth; '
                               'stepped lines, DES on the coarsest mesh and on the finest mesh with the most steps. '
                               '(b) Effective-stress paths to yield at two depths: lines, the exact paths up to the '
                               'Mohr–Coulomb line; symbols, DES elements. (c) Mean error of the yield pressure for '
                               'each mesh and number of steps.',
                    'study': ['python3 benchmarks/poromechanics/pressure_yield_reference.py',
                              'python3 tests/functional/hydromechanics/check_pressure_yield.py dynearthsol2d '
                              'OUT/mesh --spacings 0.2 0.1 0.05 --steps 400',
                              'python3 tests/functional/hydromechanics/check_pressure_yield.py dynearthsol2d '
                              'OUT/steps1600 --spacings 0.05 --steps 1600',
                              'python3 tests/functional/hydromechanics/check_pressure_yield.py dynearthsol2d '
                              'OUT/pt1-steps1600 --spacings 0.05 0.025 --steps 1600 --pt-option 1',
                              'python3 tests/functional/hydromechanics/check_pressure_yield.py dynearthsol2d '
                              'OUT/pt1-steps6400 --spacings 0.025 --steps 6400 --pt-option 1']},
 'kirsch': {'run': [('dynearthsol2d input.cfg',
                     'Meshes the quarter annulus from annulus.poly (shipped alongside) at 0.025 m inner spacing, '
                     'solves one quasi-static step to equilibrium and writes kirsch.save.* and kirsch.chkpt.*, with '
                     'kirsch.info listing the frames.')],
            'keys': [('mesh.meshing_option', '90: mesh the polygon read from poly_filename'),
                     ('mesh.poly_filename',
                      'annulus.poly: quarter annulus a = 1 m to R = 20 m; its boundary flags make the vertical cut '
                      'x0, the opening x1, the horizontal cut z0 and the outer arc z1'),
                     ('mesh.resolution', '0.025 m element size, the finest mesh of the study'),
                     ('mesh.min_angle', '28°, keeps triangles well shaped along the curved boundaries'),
                     ('bc.stress_bc_z1',
                      'type 4: normal stress on facets of any orientation, here the curved outer arc'),
                     ('bc.stress_val_z1',
                      'P = 100 kPa isotropic compression on the outer arc (tension positive, so negative)'),
                     ('bc.vbc_x1', '0: the opening (flag x1) is free'),
                     ('bc.vbc_z0',
                      '1: zero normal displacement on the horizontal symmetry cut (the vertical cut x0 uses the '
                      'default 1)'),
                     ('mat.rheology_type', 'elastic'),
                     ('mat.bulk_modulus', 'K = 1 GPa'),
                     ('mat.shear_modulus', 'G = 1 GPa'),
                     ('mat.is_plane_strain', '2D plane strain; the out-of-plane stress σyy is not zero'),
                     ('control.gravity', 'no gravity: the outer traction is the only load'),
                     ('sim.max_steps', '1: one step solved to equilibrium by PT is the static answer')],
            'fields': [('stress',
                        'element stress, rotated to radial and hoop components and compared over 1.1a ≤ r ≤ 5a in '
                        'panel (a)'),
                       ('strain', 'element strain, rotated to radial and hoop components in panel (b)')],
            'caption': '(a) Radial and hoop compression, scaled by P, against distance from the opening centre: '
                       'lines, the finite-annulus solution; symbols, DES elements (thinned for display). (b) Radial '
                       'and hoop elastic strain, scaled by P/G. (c) RMS stress error, scaled by P, against the inner '
                       'mesh spacing.',
            'study': ['python3 benchmarks/poromechanics/kirsch_reference.py',
                      'python3 tests/functional/hydromechanics/check_kirsch.py dynearthsol2d OUT --spacings .2 .1 '
                      '.05 .025']},
 'borehole': {'run': [('dynearthsol2d input.cfg',
                       'Meshes the quarter annulus from annulus.poly (shipped alongside) at 0.05 m inner spacing and '
                       'runs 144 steps of 0.25 s (36 s), saving hole.save.* every 36 steps (9, 18, 27 and 36 s) with '
                       'hole.info listing the frames.')],
              'keys': [('mesh.meshing_option', '90: mesh the polygon read from poly_filename'),
                       ('mesh.poly_filename',
                        'annulus.poly: quarter annulus a = 1 m to R = 20 m with a refined ring to 5a; flags: '
                        'vertical cut x0, borehole wall x1, horizontal cut z0, outer arc z1'),
                       ('mesh.resolution', '0.05 m element size near the wall, the finest mesh of the study'),
                       ('bc.hbc_x1', 'type 2: pore pressure prescribed on the borehole wall from hbc_pressure'),
                       ('bc.hbc_pressure', 'the x1 (wall) entry is the fluid pressure p_w = 1 MPa'),
                       ('bc.stress_bc_x1', 'type 4: normal stress on the curved wall'),
                       ('bc.stress_val_x1', '−p_w: the wall fluid pushes on the rock (tension positive)'),
                       ('bc.num_hbc_schedule', 'two schedule points for the wall pressure'),
                       ('bc.hbc_schedule_time_in_yr', '0 and 1 ms: the wall pressure ramps on almost instantly'),
                       ('bc.hbc_schedule_ratio', '0 to 1, then held at p_w'),
                       ('mat.hydraulic_perm',
                        'k = 1e-12 m²; with water viscosity 1e-3 Pa s the mobility k/μ is 1e-9 m²/(Pa s)'),
                       ('mat.fluid_bulk_modulus',
                        'K_f = 200 MPa; with porosity 0.2 the storage 1/M = φ/K_f is 1e-9 /Pa'),
                       ('mat.bulk_modulus',
                        'K = 200 MPa (shear modulus 100 MPa); together they give diffusivity c = 0.25 m²/s'),
                       ('sim.output_step_interval',
                        '36 steps: frames at 9, 18, 27 and 36 s, the four compared times')],
              'fields': [('pore pressure',
                          'nodal pressure, binned by radius and compared with the reference in panel (a)'),
                         ('stress', 'element total stress, rotated to radial and hoop components in panel (b)')],
              'caption': '(a) Pore pressure, scaled by the wall pressure p_w, against distance from the borehole '
                         'centre at the four saved times: lines, the finite-annulus poroelastic solution; symbols, '
                         'DES. (b) Total radial and hoop stress, scaled by p_w, at the last saved time. (c) Largest '
                         'pressure and stress errors against mesh spacing.',
              'study': ['python3 benchmarks/poromechanics/borehole_reference.py',
                        'python3 tests/functional/hydromechanics/check_borehole.py dynearthsol2d OUT/mesh --spacings '
                        '0.2 0.1 0.05',
                        'python3 tests/functional/hydromechanics/check_borehole.py dynearthsol2d OUT/steps '
                        '--spacings 0.1 --steps 576'],
              'replace': {'meaning': [('This checks the coupled poroelastic response around an opening: diffusion '
                                       'from a curved boundary, the pressure-driven hoop stress and the wall load '
                                       'together. It precedes the tunnel-excavation problems.',
                                       'This checks the coupled poroelastic response around an opening: diffusion '
                                       'from a curved boundary, the pressure-driven hoop stress and the wall load '
                                       'together.')]}},
 'ogata_banks': {'run': [('dynearthsol3d input.cfg',
                          'Runs the 3D column (1 m long, 2×2 elements across) for 8192 steps to 0.47 s and writes '
                          'column.save.* and column.chkpt.* every 2048 steps, with column.info listing the five '
                          'frames.')],
                 'keys': [('control.has_fluid_thermal_advection',
                           'turns on heat transport by the Darcy flux, with the combined rock–water heat capacity'),
                          ('bc.hbc_x0', 'type 2: pressure prescribed at the inlet end'),
                          ('bc.hbc_x1', 'type 2: pressure prescribed at the outlet end'),
                          ('bc.hbc_pressure',
                           '2 Pa at x = 0 and 1 Pa at x = 1 m: a steady pressure drop that drives q_x = 1 m/s'),
                          ('ic.pore_pressure_option',
                           '1: start from the linear pressure profile between the two prescribed ends, so the flow '
                           'is steady from the start'),
                          ('bc.tbc_x0', '1: fixed temperature at the inlet'),
                          ('bc.tbc_val_x0', 'inlet temperature T_in = T_0 + 1 K'),
                          ('bc.tbc_z1', '0: the top is insulated, like the other sides and the outlet'),
                          ('ic.temp_background', 'initial temperature T_0 = 273.15 K'),
                          ('mat.hydraulic_perm', '1 m², with fluid_visc = 1 Pa s: mobility κ/μ = 1 m²/(Pa s)'),
                          ('mat.fluid_bulk_modulus',
                           'deliberately tiny, which makes the hydraulic storage large without changing the steady '
                           'flux'),
                          ('mat.therm_cond',
                           'λ chosen with the capacities to give thermal diffusivity D = 0.005 m²/s'),
                          ('mat.heat_capacity',
                           'solid specific heat; with rho0, porosity and fluid_heat_capacity it sets the bulk '
                           'capacity C = 4.7 MJ/(m³ K)'),
                          ('mat.fluid_heat_capacity',
                           'water specific heat; ρ_f c_f q_x / C sets the front speed v_T')],
                 'fields': [('temperature',
                             'nodal temperature, normalized as (T − T_0)/(T_in − T_0) and plotted along the column'),
                            ('pore pressure', 'checked to stay on the steady linear profile at every saved time')],
                 'caption': 'Normalized temperature along the column at the last saved time: line, the Ogata–Banks '
                            'solution moving at the bulk-capacity speed v_T; symbols, DES nodes (about 40 shown).',
                 'study': ['python3 tests/functional/hydromechanics/check_heat_transport.py dynearthsol3d OUT/3d '
                           '--spacings 0.0125 0.00625 0.003125 --steps 8192 --restart',
                           'python3 tests/functional/hydromechanics/check_heat_transport.py dynearthsol2d OUT/2d '
                           '--spacings 0.025 0.0125 0.00625 0.003125 --steps 8192 --restart',
                           'python3 tests/functional/hydromechanics/check_heat_transport.py dynearthsol2d '
                           'OUT/2d-steps --spacings 0.0125 --steps 512 1024 2048 4096 --restart']},
 'thermal_cell': {'run': [('dynearthsol2d input.cfg',
                           'Heats the 2D cell for 100 steps of 0.01 s (1 s) with ordinary explicit mechanics and '
                           'writes heated.save.* and heated.chkpt.* at every step, with heated.info listing the '
                           'frames.')],
                  'keys': [('control.has_thermal_pressurization',
                            'couples thermal expansion of water and solid to the fluid content and stress'),
                           ('control.has_fluid_thermal_advection',
                            'uses the combined rock–water heat capacity C_0 (there is no flow to advect heat here)'),
                           ('heat_injection.enabled', 'turns on scheduled point heat sources'),
                           ('heat_injection.num_points', 'nine sources, one at each node of the 3×3-node mesh'),
                           ('heat_injection.points_x', 'x coordinates of the nine nodes'),
                           ('heat_injection.points_z', 'z coordinates of the nine nodes'),
                           ('heat_injection.rate',
                            'heat supplied at each node (W per metre of thickness in 2D), proportional to the volume '
                            'the node represents, so the cell warms uniformly at 40 K/s'),
                           ('heat_injection.duration_in_yr',
                            '0.25 s of heating (in years), followed by a hold to 1 s'),
                           ('bc.vbc_x0',
                            'type 3: every velocity component fixed; with the other faces also 3 the cell is rigid'),
                           ('mat.therm_cond', '0: no conduction, so the temperature stays uniform'),
                           ('mat.fluid_alpha', 'water thermal expansion α_f = 2e-4 /K, the source of the pressure'),
                           ('mat.alpha', '0: no solid thermal expansion'),
                           ('mat.biot_coeff', 'α = 0.98 = 1 − K/K_s'),
                           ('mat.fluid_bulk_modulus', 'K_f = 100 MPa; with porosity and K_s it sets the storage S')],
                  'fields': [('temperature',
                              'nodal temperature; its volume-weighted mean is the temperature rise in panel (a) and '
                              'gives the stored heat in panel (b)'),
                             ('pore pressure',
                              'nodal pressure; its volume average is the sealed pressure in panel (a)')],
                  'caption': '(a) Pore pressure of the sealed, rigid cell against its temperature rise: line, the '
                             'sealed linear reference βθ/S; symbols, DES. (b) Heat stored in the cell (symbols, DES) '
                             'and heat supplied by the sources (line) against time.',
                  'study': ['python3 benchmarks/poromechanics/thermal_cell.py dynearthsol2d OUT/2d --dt .01 '
                            '--figures',
                            'python3 benchmarks/poromechanics/thermal_cell.py dynearthsol2d OUT/2d-fine --methods '
                            'pt --boundaries fixed --dt .000625',
                            'python3 benchmarks/poromechanics/thermal_cell.py dynearthsol3d OUT/3d --ndims 3 --dt '
                            '.01'],
                  'replace': {'result': [('The 10 K pulse gives 154.799 kPa. Sealed pressure and fixed-capacity heat '
                                          'match their references to round-off; general open-system energy remains '
                                          'unqualified.',
                                          'The 10 K pulse raises the pore pressure to 154.799 kPa. Pressure and '
                                          'stored heat match their references to round-off.')],
                              'meaning': [('Thirty-one cases cover 2D/3D, native/HDF5 and timestep controls; eight '
                                           'thermal restarts are exact. PT has a first-order pending-pressure lag '
                                           'during heating. The drained retained-fluid ledger has a persistent '
                                           '10,450 J deficit (0.03446% of supplied heat), so the sealed comparison '
                                           'does not qualify thermal-plume energy.',
                                           'This checks thermal pressurization of pore water and the heat equation '
                                           'with a supplied heat source. The same comparison holds in 2D and 3D, for '
                                           'several time steps and across restarts. With PT the pressure lags the '
                                           'heating by a first-order amount. A drained version of the cell misses '
                                           '10,450 J (0.03446% of the supplied heat) because the heat capacity stays '
                                           'fixed while water leaves, so the energy balance of open systems such as '
                                           'thermal plumes is not covered.')]}},
 'mctigue': {'run': [('dynearthsol2d start.cfg',
                      'Zero-step run of the column at the uniform initial temperature T_0; it writes '
                      'start.save.000000, start.chkpt.000000 and start.info. These three files are shipped, so this '
                      'stage can be skipped.'),
                     ('dynearthsol2d input.cfg',
                      'Restarts from start frame 0 with the surface held at T_0 + 10 K and runs 400 steps of 100 s '
                      '(40 ks), saving mctigue.save.* every 50 steps (eight times).')],
             'keys': [('sim.is_restarting',
                       'starts from the checkpoint written by start.cfg rather than from a fresh initial state'),
                      ('sim.restarting_from_modelname', 'start: the zero-step run at T_0'),
                      ('bc.tbc_val_z1',
                       'surface temperature T_0 + ΔT = 283 K (273 K in start.cfg). A boundary temperature present at '
                       'initialization becomes part of the initial state, so the step is applied through the '
                       'restart'),
                      ('control.has_thermal_pressurization',
                       'couples thermal expansion of water and solid to the fluid content and stress'),
                      ('bc.hbc_z1',
                       'type 2: the top is drained, pore pressure held at its hbc_pressure entry (0 Pa)'),
                      ('mat.hydraulic_perm',
                       'k = 1e-19 m²: a tight rock, so pressure diffuses about five times slower than heat (c/κ = '
                       '0.21)'),
                      ('mat.therm_cond',
                       'λ = 2.65 W/(m K); with the default rock heat capacity and density it gives κ = 1e-6 m²/s'),
                      ('mat.biot_coeff', 'α = 0.8'),
                      ('mat.bulk_modulus_s', 'grain modulus K_s = 5 GPa, part of the storage S_0'),
                      ('mat.porosity', 'φ = 0.2'),
                      ('mat.bulk_modulus', 'K = 1 GPa (shear modulus also 1 GPa)'),
                      ('control.fixed_dt', '100 s; 400 steps reach 40 ks'),
                      ('sim.output_step_interval', '50 steps: eight saved times compared with the solution'),
                      ('control.gravity', 'no gravity: only heating drives pressure and strain')],
             'fields': [('pore pressure', 'nodal pressure along the column, panel (a)'),
                        ('temperature',
                         'nodal temperature along the column, panel (b); also integrated for the heat ledger'),
                        ('strain',
                         'element vertical strain by layer; integrated over depth it gives the surface heave')],
             'caption': "(a) Pore pressure against depth at the saved times: lines, McTigue's solution; symbols, "
                        'DES. (b) Temperature rise against depth. (c) Largest relative errors of pressure, '
                        'temperature, vertical strain and surface heave against mesh spacing.',
             'study': ['python3 benchmarks/poromechanics/mctigue_reference.py',
                       'python3 tests/functional/hydromechanics/check_mctigue.py dynearthsol2d OUT/mesh --spacings '
                       '0.05 0.025 0.0125 0.00625',
                       'python3 tests/functional/hydromechanics/check_mctigue.py dynearthsol2d OUT/steps400 '
                       '--spacings 0.00625 --steps 400',
                       'python3 tests/functional/hydromechanics/check_mctigue.py dynearthsol3d OUT/3d --ndims 3 '
                       '--spacings 0.025 0.0125']},
 'convection': {'run': [('dynearthsol2d input.cfg',
                         'Runs the seeded 2D box at Ra = 2 Ra_c for 400 steps of 1000 s and writes box.save.* and '
                         'box.chkpt.* every 10 steps, with box.info listing the frames.')],
                'keys': [('control.has_fluid_thermal_advection', 'turns on heat transport by the Darcy flux'),
                         ('mat.fluid_alpha',
                          'water thermal expansion α_f = 1e-3 /K; with the default gravity it makes warm water '
                          'buoyant'),
                         ('mat.hydraulic_perm',
                          'permeability κ, chosen to give Ra = 2 Ra_c in this input; the study varies it from 0.5 to '
                          '2 Ra_c'),
                         ('bc.tbc_z0', '1: fixed temperature at the base'),
                         ('bc.tbc_val_z0', 'base temperature T_0 + ΔT = 283 K'),
                         ('bc.tbc_val_z1',
                          'top temperature T_0 = 273 K; the sides keep the default insulated condition'),
                         ('ic.temperature_option',
                          '5: linear conduction profile between the top and base, plus a Gaussian seed'),
                         ('ic.temp_anomaly_amplitude', 'seed amplitude 1 mK, small enough for linear growth'),
                         ('ic.temp_anomaly_xcenter',
                          'seed centred at x = 0.25 m, off the box centre, so it projects on the first mode'),
                         ('ic.temp_anomaly_zcenter', 'seed at mid-depth'),
                         ('mat.fluid_bulk_modulus',
                          'a compressible pore fluid; pressure still diffuses about 1000 times faster than heat'),
                         ('control.has_poroelastic_pressure_feedback',
                          'off: a rigid matrix, as the reference assumes'),
                         ('mat.therm_cond', 'λ = 1 W/(m K); with C = 1 MJ/(m³ K) the conductive time CH²/λ is 1e6 s'),
                         ('control.fixed_dt', '1000 s; 400 steps cover 0.4 conductive times')],
                'fields': [('temperature',
                            'nodal temperature; its departure from the conductive profile, projected on cos(πx/W) '
                            'sin(πz/H), is the mode amplitude in panel (b), and is mapped in panel (c)')],
                'caption': '(a) Fitted growth rate of the first convective mode against Ra/Ra_c: line, linear '
                           'stability theory; symbols, DES. (b) Mode amplitude against time for each Rayleigh '
                           'number: symbols, DES; dashed lines, the theoretical rate over the second half of each '
                           'run. (c) Temperature departure from conduction at the end of the most supercritical run, '
                           'showing the single convection cell.',
                'study': ['python3 tests/functional/hydromechanics/check_convection.py dynearthsol2d OUT/h0.05 '
                          '--spacing 0.05 --ratios 0.5 0.9 1.1 1.5 2',
                          'python3 tests/functional/hydromechanics/check_convection.py dynearthsol2d OUT/h0.025 '
                          '--spacing 0.025 --ratios 0.5 0.9 1.1 1.5 2']},
 'nusselt': {'run': [('dynearthsol2d input.cfg',
                      'Runs the 2D cell at Ra = 200 on the 0.0125 m mesh for 1600 steps of 1000 s, until the heat '
                      'transport is steady, saving box.save.* every 50 steps with box.info listing the frames.')],
             'keys': [('control.has_fluid_thermal_advection', 'turns on heat transport by the Darcy flux'),
                      ('mat.fluid_alpha',
                       'water thermal expansion α_f = 1e-3 /K; with the default gravity it drives buoyancy'),
                      ('mat.hydraulic_perm', 'permeability chosen to give Ra = 200 (the study also runs Ra = 100)'),
                      ('bc.tbc_z0', '1: fixed temperature at the base'),
                      ('bc.tbc_val_z0', 'base temperature T_0 + ΔT = 283 K'),
                      ('bc.tbc_val_z1', 'top temperature T_0 = 273 K; the sides are insulated'),
                      ('ic.temperature_option', '5: linear conduction profile plus a Gaussian seed'),
                      ('ic.temp_anomaly_amplitude', '0.5 K seed, large enough to start a single roll quickly'),
                      ('mat.fluid_bulk_modulus',
                       'a very compressible pore fluid; storage only sets the transient, so this shortens the run '
                       'without changing the steady state'),
                      ('control.has_poroelastic_pressure_feedback', 'off: a rigid matrix'),
                      ('mesh.resolution', '0.0125 m, the finest mesh of the study'),
                      ('sim.max_steps', '1600 steps of 1000 s, long enough for Nu to stop changing'),
                      ('sim.output_step_interval', '50 steps: the Nu history in panel (b)')],
             'fields': [('temperature',
                         'nodal temperature; isotherms in panel (a) and the heat transport in the Nusselt number'),
                        ('pore pressure',
                         'with the temperature, rebuilds the Darcy flux used in the volume form of Nu')],
             'caption': '(a) Isotherms of the steady roll at Ra = 200: solid, DES; dashed, the independent fine-grid '
                        'solution. (b) Nusselt number (volume form) against time as the cell approaches steady '
                        'transport. (c) Nusselt-number error against mesh spacing for both Rayleigh numbers.',
             'study': ['python3 benchmarks/poromechanics/convection_steady.py',
                       'python3 tests/functional/hydromechanics/check_convection_nusselt.py dynearthsol2d OUT '
                       '--spacings 0.05 0.025 0.0125']},
 'injection_temperature': {'run': [('dynearthsol2d input.cfg',
                                    'Avdonin case: injects hot water at the corner of the 8 m quarter plane (0.1 m '
                                    'mesh) for 500 steps of 1000 s (5.8 d), saving hot.save.* every 125 steps (four '
                                    'times) with hot.info listing the frames.'),
                                   ('mkdir -p doublet && cp input-doublet.cfg doublet/ && (cd doublet && dynearthsol2d '
                                    'input-doublet.cfg)',
                                    'Doublet case, run in its own directory because both inputs write files named '
                                    'hot.*: 400 steps of 1500 s in the sealed 8 m box, saving every step.')],
                           'keys': [('injection.enabled', 'turns on scheduled point sources of fluid'),
                                    ('injection.points_x',
                                     'source at x = 0, the corner of the quarter plane (points_z = 0, the top)'),
                                    ('injection.rate',
                                     'a quarter of the line-source rate q = 2e-5 m²/s, since the model is a quarter '
                                     'plane (input-doublet.cfg: +1e-4 injector, −1e-4 producer)'),
                                    ('injection.temperature',
                                     'T_in = 333 K, the temperature of the injected water; the injecting node gains '
                                     'ρ_f c_f Q (T_in − T)'),
                                    ('control.has_fluid_thermal_advection',
                                     'turns on heat transport by the Darcy flux, required by injection.temperature'),
                                    ('bc.hbc_x1',
                                     'type 2: the far side x = 8 m is drained (pressure from hbc_pressure, 0 Pa), so '
                                     'water and heat can leave'),
                                    ('bc.hbc_z0',
                                     'type 2: the bottom is drained in the same way; the two symmetry sides stay '
                                     'sealed (input-doublet.cfg seals all four)'),
                                    ('bc.tbc_z1',
                                     '0: insulated top; all boundaries are insulated, and heat leaves only with the '
                                     'outflowing water'),
                                    ('ic.temp_background', 'initial temperature T_0 = 283 K'),
                                    ('mat.hydraulic_perm', 'k = 1e-13 m²'),
                                    ('mat.bulk_modulus', '1e12 Pa (shear modulus too): a nearly rigid skeleton'),
                                    ('control.has_poroelastic_pressure_feedback',
                                     'off: deformation does not change the pressure'),
                                    ('mat.heat_capacity',
                                     "rock specific heat; with the water's it sets the bulk capacity C"),
                                    ('mat.therm_cond', 'λ = 2 W/(m K)')],
                           'fields': [('temperature',
                                       'nodal temperature: along the symmetry edge against the Avdonin solution in '
                                       'panel (a), and summed with the heat capacity into the stored heat in panel '
                                       '(b)')],
                           'caption': '(a) Temperature rise against distance from the injection point along a '
                                      "symmetry side at four times: lines, Avdonin's solution; symbols, DES. (b) Hot "
                                      'doublet: heat gained by the box (DES) against the ledger of injected minus '
                                      'extracted heat. (c) Errors against mesh spacing: largest Avdonin temperature '
                                      'error, as a percentage of T_in − T_0, and doublet heat-ledger error, as a '
                                      'percentage of injected heat.',
                           'study': ['python3 benchmarks/poromechanics/avdonin_reference.py',
                                     'python3 tests/functional/hydromechanics/check_injection_temperature.py '
                                     'dynearthsol2d OUT/mesh --spacings 0.2 0.1',
                                     'python3 tests/functional/hydromechanics/check_injection_temperature.py '
                                     'dynearthsol2d OUT/mesh-coarse --spacings 0.4',
                                     'python3 tests/functional/hydromechanics/check_injection_temperature.py '
                                     'dynearthsol2d OUT/steps200 --spacings 0.2 --steps 200',
                                     'python3 tests/functional/hydromechanics/check_injection_temperature.py '
                                     'dynearthsol2d OUT/steps1000 --spacings 0.2 --steps 1000',
                                     'python3 tests/functional/hydromechanics/check_injection_energy.py '
                                     'dynearthsol2d OUT/doublet --spacings 0.2 0.1 --steps 400']}}

MANUAL['plastic_heat'] = dict(
    run=[('dynearthsol2d input.cfg', 'Runs 1000 steps to 0.2 s; every physical step is saved under cell.')],
    keys=[('control.has_shear_heating', 'converts accepted nonelastic work into heat; unit heat fraction'),
          ('control.has_thermal_pressurization', 'couples fluid thermal expansion to sealed pore pressure'),
          ('control.has_fluid_thermal_advection', 'uses the combined solid and water heat capacity; zero flux here'),
          ('control.has_hydraulic_diffusion', 'enables pressure evolution; zero permeability prevents drainage'),
          ('control.fixed_dt', '0.0002 s, the finest of three separately checked physical steps'),
          ('control.is_using_mixed_stress', 'off, so this reference checks the constitutive work before mixed-stress projection'),
          ('mat.rheology_type', 'plasticity with evolving aging rate-and-state friction'),
          ('mat.direct_a', 'direct velocity dependence; b is larger, giving velocity weakening'),
          ('mat.evolution_b', 'state-dependent friction evolution'),
          ('mat.characteristic_distance', 'sets the state relaxation distance'),
          ('bc.vbc_val_x1', 'right-face extension velocity; bottom velocity provides matching compression')],
    fields=[('temperature', 'physical temperature rise from 300 K, compared with integrated plastic work divided by mixture capacity'),
            ('pore pressure', 'sealed thermal pressurization; compared with βΔT/S0'),
            ('stress', 'principal half-span τ=(σxx−σzz)/2'),
            ('shear heating', 'last accepted nonelastic power density in W/m³'),
            ('pending shear heat', 'accepted energy awaiting the next thermal update, J/m per element in 2D'),
            ('friction state variable', 'aging state θ, advanced once per physical step')],
    caption='(a) Principal shear stress during elastic loading and aging plastic flow. (b) Physical temperature rise. (c) Sealed pore pressure. Lines evaluate the continuum aging solution and work integral; circles are the finest DES run. The physical fields retain the measured first-order lag.',
    study=['python3 tests/functional/hydromechanics/check_heat_restart.py dynearthsol2d OUT --ndims 2 --study --coupled --dt 0.0002 --steps 1000',
           'python3 benchmarks/poromechanics/plastic_heat.py --exe dynearthsol2d --dimensions 2 --output OUT-ep'])

MANUAL['fault_heat'] = dict(
    run=[('dynearthsol2d input.cfg', 'Slides at 1 m/s to 0.1 s (100 mm). Saves 600 time intervals; physical timestep is selected from the mixture mass.')],
    keys=[('control.has_shear_heating', 'converts accepted plastic work to heat, with unit heat fraction'),
          ('control.has_thermal_pressurization', 'expanding sealed water raises pressure'),
          ('control.has_fluid_thermal_advection', 'selects mixture heat capacity; this sealed case has zero Darcy flux'),
          ('control.is_quasi_static', 'false: use physical saturated mixture mass'),
          ('control.has_PT', 'false: ordinary explicit mechanics, with no PT solve'),
          ('control.damping_option', 'zero: no damping in the physical case'),
          ('control.fixed_dt', 'zero: automatic physical wave timestep'),
          ('control.dt_fraction', '0.2: fraction of the mechanical stability bound'),
          ('bc.vbc_z1', '8: the top has a uniform free normal velocity under its applied load'),
          ('bc.stress_val_z1', '−50 MPa total normal traction'),
          ('bc.vbc_val_x1', 'top horizontal velocity; the configured ratios give a linear profile on both sides'),
          ('mat.rheology_type', 'constant-friction Mohr–Coulomb plasticity, zero dilation')],
    fields=[('stress', 'σxz is shear resistance; σzz checks departure from imposed normal equilibrium'),
            ('temperature', 'mixture temperature minus the initial 300 K'),
            ('pore pressure', 'sealed thermal pressurization with normal compliance'),
            ('shear heating', 'last accepted plastic power density, W/m³'),
            ('pending shear heat', 'accepted heat awaiting its thermal update, J/m per element in 2D'),
            ('strain', 'trace gives the normal volume change in the fluid-content ledger')],
    caption='(a) Shear resistance under fixed total normal traction; the dashed line is the zero-water-expansion reference. (b) Temperature rise from plastic work. (c) Sealed pore pressure. Black curves solve the independent continuous MC/Jaumann equations; markers are ordinary explicit DES with physical mass and zero damping.',
    study=['python3 benchmarks/poromechanics/fault_heat.py --exe ./dynearthsol2d --output OUT --mode physical --strain 10 --samples 600',
           'python3 benchmarks/poromechanics/fault_heat.py --exe ./dynearthsol2d --output OUT-switch --mode gvs --scaling 1e5 --velocity 1e-5 --strain .02 --physical-transition'])

MANUAL['salencon'] = dict(
    run=[('dynearthsol2d seed.cfg',
          'Meshes annulus.poly (0.05 m near the opening, coarser outward) and writes the initial state: uniform 500 kPa '
          'compression, the wall still supported at 500 kPa, frame 0 with a checkpoint.'),
         ('ls load-0*.cfg | xargs -n 1 dynearthsol2d',
          'Fifteen unloading steps, load-001.cfg to load-015.cfg in order. Each restarts from the previous one and lowers the wall support pressure by P/16 '
          '(load-path.json lists them); one PT step per file brings the rock to equilibrium.'),
         ('dynearthsol2d salencon.cfg',
          'The last step removes the support completely; this is the state compared with Salençon. input.cfg is the same '
          'file, and the shipped load-015 frame lets you run it without the earlier steps (but the displacement in the '
          'figure sums the velocities of every step).')],
    keys=[('ic.initial_stress', 'uniform in-situ compression of 500 kPa in x and z (the out-of-plane stress starts equal)'),
          ('bc.stress_bc_x1', 'normal traction on the opening wall: the support pressure'),
          ('bc.stress_val_x1', 'support pressure in this step, negative in compression; it is −500 kPa in seed.cfg, '
                               'steps down by 31.25 kPa per load file, and is zero here'),
          ('bc.stress_bc_z1', 'normal traction on the outer boundary'),
          ('bc.stress_val_z1', 'the outer boundary keeps the in-situ 500 kPa compression throughout'),
          ('bc.vbc_z0', 'the straight cut along x has zero normal velocity (symmetry); vbc_x0 does the same for the other cut'),
          ('mat.rheology_type', 'Mohr–Coulomb elasto-plasticity'),
          ('mat.cohesion0', 'cohesion 100 kPa; cohesion1 is the same, so there is no softening'),
          ('mat.dilation_angle0', 'dilation 30°, equal to the friction angle: associated flow'),
          ('control.has_PT', 'each step is solved to equilibrium by PT; the material is rate independent, so the 1 s step '
                             'has no physical meaning'),
          ('control.PT_relative_tolerance', 'tight equilibrium tolerance for every unloading step'),
          ('control.is_using_mixed_stress', 'smooths element pressure after the plastic return; needed here, because '
                                            'without it the plastic ring next to the wall develops an alternating pressure '
                                            'between neighbouring elements and diverges'),
          ('control.has_moving_mesh', 'small strain: the mesh stays in its original position'),
          ('mesh.poly_filename', 'the quarter-annulus outline, with a finer region within 3a of the centre')],
    fields=[('stress', 'radial and hoop stress at element centres, rotated to polar components (panel a)'),
            ('stressyy', 'out-of-plane stress; written in the checkpoint (salencon.chkpt.*), not the .save frame (panel a)'),
            ('velocity', 'times the 1 s step and summed over every step, it gives the radial displacement (panel b)'),
            ('strain', 'total strain from the initially compressed state, compared as radial and hoop strain'),
            ('plastic strain', 'marks the yielded ring; its outer edge is the plastic radius')],
    caption=('(a) Radial, hoop and out-of-plane stress against $r/a$ after the support is removed, finest mesh: DES at element '
             'centres (dots) and Salençon (lines); the dashed line marks the analytical plastic radius. (b) Inward radial '
             'displacement of the nodes, scaled by the analytical wall convergence, with the analytical curve. (c) RMS errors '
             'of hoop stress, displacement and radial strain over $1.05a\\le r\\le5a$ against the inner mesh spacing.'),
    study=['python3 tests/functional/hydromechanics/check_salencon.py dynearthsol2d OUT --spacings .2 .1 .05 '
           '--dilation 30 --mixed --refined --load-steps 16',
           'python3 benchmarks/poromechanics/qualify_salencon.py OUT FIGURES'])

MANUAL['elder'] = dict(
    run=[('dynearthsol2d input.cfg', 'Runs the selected 2D half box for 100,000 s in 2,000 steps of 50 s. elder.poly separates the hot and insulated base segments.')],
    keys=[('control.has_hydraulic_diffusion', 'enables Darcy pressure transport'),
          ('control.has_fluid_thermal_advection', 'carries heat with the Darcy flux'),
          ('control.has_poroelastic_pressure_feedback', 'off: the matrix is rigid for the flow comparison'),
          ('control.fixed_dt', '50 s physical interval; stable transport subcycles are internal'),
          ('mat.fluid_bulk_modulus', '100 Pa: deliberately scaled finite pressure storage'),
          ('mat.hydraulic_perm', 'sets Ra = 400'),
          ('mat.fluid_alpha', 'Boussinesq fluid expansion gives buoyancy'),
          ('mesh.poly_filename', 'splits the base into a hot z0 segment and an insulated slant segment'),
          ('bc.tbc_z0', 'heater temperature T0 + 10 K; tbc_z1 fixes the cold top'),
          ('bc.vbc_n0', 'zero normal velocity on the insulated base segment')],
    fields=[('temperature', 'scaled temperature and the retained six-time heat-content comparison'),
            ('coordinate', 'positions for bilinear reference interpolation and the selected final field'),
            ('connectivity', 'triangle areas for lumped integration of heat content')],
    caption='(a) Independent 256 × 512 reference and (b) selected DES mesh at dimensionless time 0.1, with the same temperature scale and white 0.2 and 0.6 contours. The red bar marks the heater. (c) Heat content at every comparison time, including the early discrepancy. Complete field errors and reference grid checks are in the tables.',
    study=['python3 benchmarks/poromechanics/elder_reference.py reference.npz --grids 32 64 128 256',
           'python3 tests/functional/hydromechanics/check_elder.py dynearthsol2d reference.npz OUT --spacings .0625 .03125 .015625 --dt 50 --fluid-bulk 100 --threads 8'])

MANUAL['convection_3d'] = dict(
    run=[('dynearthsol3d input.cfg', 'Runs the selected cube for 200,000 s in 400 explicit physical steps. DES subcycles transport within those intervals.')],
    keys=[('control.has_hydraulic_diffusion', 'enables Darcy pressure transport'),
          ('control.has_fluid_thermal_advection', 'carries heat with the Darcy flux'),
          ('control.has_poroelastic_pressure_feedback', 'off: the stiff matrix supplies no deformation feedback'),
          ('control.fixed_dt', '500 s physical interval; the temporal control uses 250 s'),
          ('mat.fluid_bulk_modulus', '10 Pa: deliberately scaled storage, retained in the analytical reference'),
          ('mat.hydraulic_perm', 'sets the Rayleigh number; held fixed through each run'),
          ('mat.fluid_alpha', 'Boussinesq buoyancy coefficient'),
          ('ic.temp_anomaly_amplitude', '0.1 mK seed, kept small for the linear reference'),
          ('bc.tbc_z0', 'fixed hot base; tbc_z1 fixes the cold top'),
          ('mesh.meshing_elem_shape', 'structured tetrahedral cube')],
    fields=[('temperature', 'departure from the conductive profile, projected onto the (1,1,1) mode'),
            ('pore pressure', 'pressure disturbance driving Darcy flux'),
            ('coordinate', 'positions for the three-dimensional mode and y/H=0.25 section'),
            ('connectivity', 'tetrahedral volumes for the mode projection')],
    caption='(a) Late growth rate against the finite-storage analytical eigenvalue. (b) DES amplitudes normalized at dimensionless time 0.1; lines use the analytical slope with that one amplitude anchor. (c) Temperature departure at y/H=0.25 and dimensionless time 0.2 on the selected mesh. Rates and all refinement controls are in the accompanying table.',
    study=['python3 tests/functional/hydromechanics/check_convection_3d.py dynearthsol3d OUT --spacing .03125 --ratio 1.4 --fluid-bulk 10 --steps 400 --dt 500',
           'Repeat with spacings .125 and .0625, ratios .7 and 1.2, and --seed 0 for the rest control. Use 800 steps of 250 s for the separate time control.'])

for _key, _description, _cases in (
        ('mass', 'Compresses the sealed unit cell for 0.1 s in 4,000 steps and saves its pressure, physical water inventory and numerical buffer.', 'mass lag'),
        ('drained_heat', 'Supplies 30 MW of heat for 5 s in 200 steps while the pressure-controlled boundary exports expanding water.', 'heat'),
        ('mixing', 'Injects 323 K water for 5 s in 200 steps into a drained, initially 273 K unit cell.', 'heat')):
    MANUAL[_key] = dict(
        run=[('dynearthsol2d input.cfg', _description)],
        keys=[('control.hydraulic_content_model', 'selects the finite-density saturated inventory law'),
              ('mat.porosity_model', 'evolves Eulerian Biot porosity from the accepted strain and pressure'),
              ('control.fixed_dt', 'physical coupling interval; transport can subcycle within it'),
              ('mat.fluid_bulk_modulus', 'constant fluid compressibility in the exponential density law'),
              ('mat.bulk_modulus_s', 'grain compressibility in the porosity law'),
              ('mat.fluid_alpha', 'volumetric fluid thermal expansivity'),
              ('sim.output_step_interval', 'saved-frame interval; ledger entries belong to individual physical steps')],
        fields=[('pore pressure', 'physical nodal pressure in Pa'),
                ('fluid mass', 'retained physical nodal water mass in kg'),
                ('hydraulic split buffer', 'numerical fixed-stress mass increment, excluded from density and heat'),
                ('solid heat capacity', 'material heat-capacity inventory, unchanged by fluid drainage'),
                ('temperature', 'temperature in K'),
                ('fluid mass balance', 'physical and numerical mass changes, measured boundary/source exchange and residual')],
        study=['python3 tests/functional/hydromechanics/check_mass.py dynearthsol2d OUT --cases '+_cases])
