#!/usr/bin/env python
import numpy as np
from math import sqrt
from scipy.special import erf

import matplotlib.pyplot as plt

def half_space_cooling_T(z, T0, Tm,  age_in_myrs, alpha):
    myrs2sec = 86400 * 365.2425e6

    T = T0 + (Tm - T0) * erf(z /
            sqrt(4 * alpha * age_in_myrs * myrs2sec) )
    return T

def continental_radiogenic_T(z,T0,hr,k,qm,rhoH0):
    T = T0 + qm/k*z + rhoH0*hr**2/k*(1-np.exp(-z/hr))
    return T

def continental_radiogenic_T2(z,T0,Tm,bdy_arr,hr,k_arr,rhoH0_arr,is_hr_from_layer_top=False):
    layer_top = bdy_arr[:-1]
    layer_bot = bdy_arr[1:]
    thick_arr = layer_bot - layer_top

    dT = sum_dt_of_layers(bdy_arr,k_arr,rhoH0_arr,hr,is_hr_from_layer_top)
    qm = (Tm - T0 - dT) / sum(thick_arr) * avg_k(k_arr,thick_arr)

    dT_h_layer = dT_h(layer_top,layer_bot,k_arr,rhoH0_arr,hr,is_hr_from_layer_top)
    dtlayer = dT_h_layer + thick_arr * qm / k_arr
    
    T = np.zeros_like(z)
    for i, zt in enumerate(z):
        is_value = False
        for j in range(len(k_arr)):
            if zt <= layer_bot[j]:
                dTh = dT_h(layer_top[j],zt,k_arr[j],rhoH0_arr[j],hr,is_hr_from_layer_top)
                T[i] = T0 + np.sum(dtlayer[:j]) + qm/k_arr[j]*(zt - layer_top[j]) + dTh
                is_value = True
                break
        if not is_value:
            raise ValueError(f'z={zt} is out of range')
    return T, qm

def continental_radiogenic_T3(z,dTlayer_acc,bdy_top,rhoH0,k,T0,hr,qm,is_hr_from_layer_top=False):
    dTh = dT_h(bdy_top,z,k,rhoH0,hr,is_hr_from_layer_top)
    T = T0 + dTlayer_acc + qm * (z - bdy_top) / k + dTh
    return T

def get_h_value(z,hr,rhoH0):
    return rhoH0*np.exp(-z/hr)
def get_h_integral(z,hr,rhoH0):
    return rhoH0*hr**2*(1-np.exp(-z/hr))

def dT_h(zstart,zend,k,rhoH0,hr,is_hr_from_layer_top):
    if is_hr_from_layer_top:
        zend = zend - zstart
        zstart = zstart - zstart
    hstart = get_h_integral(zstart,hr,rhoH0)
    hend = get_h_integral(zend,hr,rhoH0)
    return (hend - hstart) / k

def sum_dt_of_layers(bdy_arr,k_arr,rhoH0_arr,hr,is_hr_from_layer_top):
    dt = dT_h(bdy_arr[:-1],bdy_arr[1:],k_arr,rhoH0_arr,hr,is_hr_from_layer_top)
    return np.sum(dt)

def avg_k(k_arr,thick_arr):
    return np.sum(thick_arr)/np.sum(thick_arr/k_arr)

def get_mat_profile(z,bdy_arr,phase_arr):
    mat = np.zeros_like(z,dtype=int)
    for i, zt in enumerate(z):
        is_value = False
        for j in range(len(phase_arr)):
            if zt <= bdy_arr[j+1]:
                mat[i] = phase_arr[j]
                is_value = True
                break
        if not is_value:
            raise ValueError(f'z={zt} is out of range')
    return mat

def get_visc(edot, T,  n, A, E, P=0.0, V=0.0):
    '''
    from: https://github.com/tan2/geoflac/blob/master/util/visc_profile.py
    edot: second invariant of strain rate
    T: temperature in Celsius
    n, A, E: viscosity parameters
    P: pressure in MPa (compression -> negative)
    V: activation volume in cm^3/mol

    return viscosity in Pascal.s
    '''
    R = 8.31448  # gas constant
    pow = 1.0/n - 1
    pow1 = -1.0/n
    visc = 0.25 * (edot**pow) * (0.75*A)**pow1 * np.exp((E-P*V)/ (n * R * (T + 273))) * 1e6
    return visc

def get_dTlayer_acc(bdy_arr,k_arr,rhoH0_arr,qm,hr,is_hr_from_layer_top):
    layer_top = bdy_arr[:-1]
    layer_bot = bdy_arr[1:]
    thick_arr = layer_bot - layer_top

    dT_h_layer = dT_h(layer_top,layer_bot,k_arr,rhoH0_arr,hr,is_hr_from_layer_top)
    dtlayer_arr = dT_h_layer + thick_arr * qm / k_arr
    dtlayer_acc = np.cumsum(dtlayer_arr)
    dtlayer_acc[1:] = dtlayer_acc[:-1]
    dtlayer_acc[0] = 0
    return dtlayer_acc

def main():
    # material properties
    # Brune 2014
    mat_rho = [3280, 2850, 2700] # kg/m^3 density
    mat_k = [3.3, 2.5, 2.5] # W/mK thermal conductivity
    mat_cp = [1200, 1200, 1200] # J/kgK specific heat capacity
    mat_H0 = [0, 1e-10,6.0e-10] # W/kg heat production
    mat_n = np.array([3.5, 3.0, 4.0])
    mat_B = [-15.56,-15.40,-28.00]

    # mat_A = np.array([2.87E+03, 8.29E+00, 5.21E-07])
    mat_A = np.power(4,1-mat_n) / 3 * np.power(10,mat_B) * np.power(1e6,mat_n)
    # mat_A = np.power(4,1-mat_n) / 3 * np.exp(mat_B) * np.power(1e6,mat_n)
    # mat_A = np.power(4,1-mat_n) / 3 * np.exp(mat_B) * np.power(1e6,1)

    mat_E = [5.30E+05, 3.56E+05, 2.23E+05]
    mat_V = [13.0, 0.0, 0.0] # cm^3/mol
    # Andres 2019
    # mat_rho = [3300,3300,2850,2700] # kg/m^3 density
    # mat_k = [3.3, 3.3, 2.5, 2.1]
    # mat_cp = [1200, 1200, 1200, 1200]
    # mat_H0 = [0.0, 0.0, 7.02E-11, 4.81e-10]
    # mat_A = [1.17E+03, 5.58E+01, 8.29E+00, 5.21E-07]
    # mat_n = [3.5, 4.2, 3.0, 4.0]
    # mat_E = [5.30E+05, 4.45E+05, 3.56E+05, 2.23E+05]
    # mat_V = [13.0, 0.0, 0.0, 0.0] # cm^3/mol


    # 9.6e-10 W/kg for granite
    hr = 20e3 # km length_scale_for_the_decrease_of_heat_production
    is_hr_from_layer_top = True

    mat_k = np.array(mat_k)
    mat_rho = np.array(mat_rho)
    mat_cp = np.array(mat_cp)
    mat_H0 = np.array(mat_H0)
    mat_E = np.array(mat_E)
    mat_V = np.array(mat_V)
    mat_rhoH0 = mat_rho * mat_H0
    # Brune 2014
    mat_rhoH0[2] = 1.5e-6 # W/m^3
    mat_rhoH0[1] = 0.2e-6 # W/m^3

    T0 = 273
    Tm = 1300 + 273

    Zulb = 22e3 # km
    Zmoho = 33e3 # km
    Zbot = 90e3 # km
    edot = 1e-15 # 1/s
    
    phase_arr = [2, 1, 0]
    bdy_arr = [0, Zulb, Zmoho, Zbot]
    bdy_arr = np.array(bdy_arr)
    nlayer = len(phase_arr)    
    k_arr = mat_k[phase_arr]
    rhoH0_arr = mat_rhoH0[phase_arr]

    layer_top = bdy_arr[:-1]
    layer_bot = bdy_arr[1:]
    thick_arr = layer_bot - layer_top
    
    z = np.linspace(0, Zbot, int(Zbot*0.1)+1)
    dz = z[1]-z[0]
    
    mat = get_mat_profile(z,bdy_arr,phase_arr)
    visc_A = mat_A[mat]
    visc_n = mat_n[mat]
    visc_E = mat_E[mat]
    visc_V = mat_V[mat]
    rho = mat_rho[mat]
    rhoH0 = mat_rhoH0[mat]
    k = mat_k[mat]

    dT = sum_dt_of_layers(bdy_arr,k_arr,rhoH0_arr,hr,is_hr_from_layer_top)
    qm = (Tm - T0 - dT) / sum(thick_arr) * avg_k(k_arr,thick_arr)
    dtlayer_acc_arr = get_dTlayer_acc(bdy_arr,k_arr,rhoH0_arr,qm,hr,is_hr_from_layer_top)
    
    bdy_top = np.zeros_like(z)
    dtlayer = np.zeros_like(z)
    for i, ph in enumerate(phase_arr):
        ind = mat==ph
        bdy_top[ind] = bdy_arr[i]
        dtlayer[ind] = dtlayer_acc_arr[i]

    # calculate temperature increase by radiogenic heat production in the crust    
    # temp, qm = continental_radiogenic_T2(z,T0,Tm,bdy_arr,hr,k_arr,rhoH0_arr,is_hr_from_layer_top=is_hr_from_layer_top)
    temp = continental_radiogenic_T3(z,dtlayer,bdy_top,rhoH0,k,T0,hr,qm,is_hr_from_layer_top=is_hr_from_layer_top)

    pres = np.cumsum(rho)*9.8*dz
    pres -= pres[0]
    stren_fric = pres * np.tan(np.pi/6.)

    visc = get_visc(edot, temp-273, visc_n, visc_A, visc_E, pres/1e6, visc_V)
    stren_visc = 2. * edot * visc
    
    stren = np.minimum(stren_fric,stren_visc)

    if is_hr_from_layer_top:
        # use z in layer to calculate heat production
        z_layer = z - bdy_top
    else:
        z_layer = z
    hp = get_h_value(z_layer,hr,rhoH0)

    print(f'temperature at surface: {temp[0]-273:.0f} degC')
    print(f'temperature at moho: {np.interp(Zmoho,z,temp)-273:.0f} degC')
    print(f'temperature at bottom: {temp[-1]-273:.0f} degC')

    # plot
    fig, axes = plt.subplots(1,2,figsize=(8,8))
    
    ax = axes[1]
    ax.plot(stren_visc/1e6,z/1e3,'--',color='r',label="Viscosity")
    ax.plot(stren_fric/1e6,z/1e3,'--',color='b',label="Friction")
    ax.plot(stren/1e6,z/1e3,'-',color='k',label="Viscosity")
    ax.set_xlabel(r"Strength (MPa)")
    ax.set_xlim(0,800)

    ax = axes[0]
    ax.vlines(Tm-273,0,Zbot/1e3,ls='--',color='grey',lw=1)
    ax.text(Tm-273+0.1,0.1,"Tm",fontsize=12,color='grey')
    
    for bdy in bdy_arr[1:-1]:
        ax.hlines(bdy*1e-3,0,1500,ls='--',color='grey',lw=1)

    color = 'k'
    ax.plot(temp-273, z/1e3,'-',color=color,label="Radiogenic")
    for bdy in bdy_arr[1:-1]:
        tmp = np.interp(bdy,z,temp)
        ax.plot([tmp-273],[bdy*1e-3],'o',color=color,ms=4,mfc='w')
        ax.text(tmp-273+10,bdy*1e-3-1,f"{tmp-273:.0f}"+'$^\circ{C}$',fontsize=12,ha='left',va='bottom',color=color)

    ax.plot([temp[0]-273],[0.],'o',color=color,ms=6,mfc='w')
    ax.text(100,0.5,f"{temp[0]-273:.0f}"+'$^\circ{C}$',fontsize=12,ha='left',va='top',color=color)
    ax.plot([temp[-1]-273],[Zbot*1e-3],'o',color=color,ms=4,mfc='w')
    ax.text(900,Zbot*1e-3-0.3,f"{temp[-1]-273:.0f}"+'$^\circ{C}$',fontsize=12,ha='left',va='bottom',color=color)

    ax.text(300,93,'q$_m$='+f'{qm*1e3:.0f} mW/m$^2$',fontsize=12,ha='left',va='top',color=color)


    ax.text(100,hr/1e3,'hr='+f'{hr*1e-3:.0f} km',fontsize=12,ha='left',va='bottom',color=color)

    ax2 = ax.twiny()
    
    ax2.plot(hp*1e6,z/1e3,'--',color='r',label="Heat production",lw=1)
    ax2.set_xlabel(r"Heat production ($\mu$W/m3)")
    ax2.legend(loc='upper right')
    ax2.set_xlim(0,4)

    ax.legend(loc='lower left')
    ax.set_xlabel(r"Temperature ($^\circ{C}$)")
    ax.set_ylabel("Depth (km)")

    ax.set_xlim(0,1500)
    
    for ax in axes:
        ax.set_ylim(Zbot/1e3,0)
        ax.grid(ls='--')
    fig.tight_layout()
    filename = f'geo-{nlayer:d}-{Zmoho/1e3:.0f}-{rhoH0_arr[0]:.1e}-{hr/1e3:.0f}km.png'
    filename = 'test.png'
    
    print(f'save figure to {filename}')
    fig.savefig(filename)
    
    return

if __name__ == "__main__":
    main()