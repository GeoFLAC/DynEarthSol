#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def terrig(xi,zi,S0,C0,C1,is_strip=False):
    ind = zi < 0
    x, z = xi[ind], zi[ind]
    dxi, dzi = np.diff(x), np.diff(z)

    si = dzi/dxi
    si = np.concatenate(([0], si, [0]))
    
    si[0] = si[1] - S0 / C0
    si[-1] = si[-2]

    dh = C0 * np.exp(-C1*z) * np.diff(si) / np.gradient(x)
    if is_strip:    
        is_terrig = True
        for i, h in enumerate(dh):
            if is_terrig == True:
                if h <= 0:
                    is_terrig = False
                if i > 0:
                    if h > dh[i-1]:
                        dh[i] = dh[i-1]
            else:
                dh[i] = 0

    dh_all = np.zeros(len(xi))
    dh_all[ind] = dh
        
    return dh_all

def main(vars):
    xmin, xmax = 0, 100e3

    x_surface = [0, 20e3, 50e3, 100e3]
    z_surface = [0, -0.5e3, -0.5e3, 0]
    dx_max = 10e3
    dx_min = 1e3
    nx = int(xmax / vars.resolution)
    dx_arr = np.random.uniform(dx_min, dx_max,nx)
    xi = np.cumsum(dx_arr)
    xi = np.insert(xi, 0, 0)/xi[-1]*xmax

    zi = np.interp(xi, x_surface, z_surface)

    year2sec =  60*60*24*365.25
    
    print(f'S0 = {vars.S0:.2e} C0 = {vars.C0:.2e} C1 = {vars.C1:.2e} dt = {vars.dt/year2sec:.2f} yr\n')
    print(f'params: flux = {vars.S0:.2e} m^2/s ({vars.S0*year2sec:.2e} m^2/yr)')
    print('time     flux (m^2/s)   flux (m^2/yr)')
    t = 0
    nstep = 0
    ninterval = vars.dt_print // vars.dt

    fig, axes = plt.subplots(2, 1, figsize=vars.figsize)

    while True:
        
        dh_rate = terrig(xi,zi,vars.S0,vars.C0,vars.C1,is_strip=vars.is_strip)

        flux = np.sum(dh_rate * np.gradient(xi))

        zi += dh_rate * vars.dt
        zi[zi > 0] = 0
        # zi = np.min(zi,0)
        if nstep % ninterval == 0:
            flux_m_yr = flux * year2sec
            print(f'{t/year2sec:4.2e} {flux:.2e} m^2/s ({flux_m_yr:.2e} m^2/yr) {flux/vars.S0*100:.1f}%')

            dh = dh_rate * vars.dt
            ax = axes[0]
            ax.plot(xi, dh, lw=1)

            ax = axes[1]
            ax.plot(xi, zi, '-o', lw=1,ms=4)


        t += vars.dt
        nstep += 1
        if t >= vars.tmax:
            break
    
    for ax in axes:
        ax.grid('--')
        ax.set_xlim(xmin, xmax)
        
    filename = f'terrig-{vars.dt/year2sec:.0f}-{vars.S0:.2e}-{vars.C0:.2e}-{vars.C1:.2e}-{vars.is_strip}.png'
    print(f'save {filename}')
        
    plt.tight_layout()
    plt.savefig(filename,dpi=300)
    plt.close()


from types import SimpleNamespace

year2sec = 60*60*24*365.25
if __name__ == '__main__':
    vars = SimpleNamespace()
    vars.S0 = 1e-7 # terrig_sediment_area
    vars.C0 = 1.e-6 # terrig_sediment_diffusivity
    vars.C1 = 5e-4 # terrig_depth_coefficient
    vars.dt = 5 * year2sec # yr
    vars.dt_print = 1e4 * year2sec # yr
    vars.tmax = 1e5 * year2sec # yr
    vars.is_strip = False
    vars.resolution = 1e3 # m
    vars.figsize = (9, 3)

    main(vars)