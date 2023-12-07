#!/usr/bin/env python
import numpy as np
from math import sqrt
from scipy.special import erf
import sympy as sp

import matplotlib.pyplot as plt

def half_space_cooling_T(z, T0, Tm,  age_in_myrs, alpha):
    myrs2sec = 86400 * 365.2425e6

    T = T0 + (Tm - T0) * erf(z /
            sqrt(4 * alpha * age_in_myrs * myrs2sec) )
    return T

def continental_radiogenic_T(z,T0,hr,k,qm,rhoH0):
    T = T0 + qm/k*z + rhoH0*hr**2/k*(1-np.exp(-z/hr))
    return T

def continental_radiogenic_T2(z,T0,thick_arr,hr,k_arr,qm,rhoH0_arr):
    '''
    No heat production below the moho
    '''
    dtlayer = np.array([dt_of_heat_production(thick_arr[i],k_arr[i],rhoH0_arr[i],hr) for i in range(len(thick_arr))])

    dtlayer = dtlayer + qm*thick_arr/k_arr
    
    T = np.zeros_like(z)
    t_arr_inv = thick_arr.copy()[::-1]
    k_arr_inv = k_arr.copy()[::-1]
    rhoH0_arr_inv = rhoH0_arr.copy()[::-1]
    dtlayer_inv = dtlayer.copy()[::-1]

    for i, zt in enumerate(z):
        for j in range(len(t_arr_inv)):
            if zt <= np.sum(t_arr_inv[:j+1]):
                dz = zt - np.sum(t_arr_inv[:j])
                T[i] = T0 + np.sum(dtlayer_inv[:j]) + qm/k_arr_inv[j]*dz + rhoH0_arr_inv[j]*hr**2/k_arr_inv[j]*(1-np.exp(-dz/hr))
                break
    return T

def dt_of_heat_production(thickness,k,rhoH0,hr):
    return rhoH0*hr**2*(1-np.exp(-thickness/hr))/k

def sum_dt_of_layers(thick_arr,k_arr,rhoH0_arr,hr):
    dt = 0
    for i in range(len(thick_arr)):
        dt += dt_of_heat_production(thick_arr[i],k_arr[i],rhoH0_arr[i],hr)
    return dt

def avg_k(k_arr,thick_arr):
    k = 0
    for i in range(len(thick_arr)):
        k += k_arr[i]*thick_arr[i]
    return k/sum(thick_arr)

def main():
    k_arr = [3.3, 2.5] # W/mK thermal conductivity
    rho_arr = [3300, 2800] # kg/m^3 density
    cp_arr = [1000, 1000] # J/kgK specific heat capacity
    H0_arr = [0, 6.0e-10] # W/kg heat production
    H0 = 6.0e-10 # W/kg heat production
    # 9.6e-10 W/kg for granite

    k_arr = np.array(k_arr)
    rho_arr = np.array(rho_arr)
    cp_arr = np.array(cp_arr)
    H0_arr = np.array(H0_arr)

    alpha_arr = k_arr/rho_arr/cp_arr
    rhoH0_arr = rho_arr * H0_arr
    
    Zbot = 120e3 # km
    Zmoho = 33e3 # km
    thick_arr = [Zbot-Zmoho, Zmoho]
    thick_arr = np.array(thick_arr)
    hr = 50e3 # km length_scale_for_the_decrease_of_heat_production

    # Brune 2014
    rhoH0_arr[-1] = 1.5e-6 # W/m^3

    T0 = 273
    Tm = 1300 + 273
    
    age = 26 # Myrs

    z = np.linspace(0, Zbot, 1000)

    th = half_space_cooling_T(z, T0, Tm, age,alpha_arr[0])

    hz = rhoH0_arr[-1]*np.exp(-z/hr)

    thick_arr_inv = thick_arr.copy()[::-1]
    rhoH0_arr_inv = rhoH0_arr.copy()[::-1]
    hz3 = np.zeros_like(z)
    for i, zt in enumerate(z):
        for j in range(len(thick_arr_inv)):
            if zt < np.sum(thick_arr_inv[:j+1]):
                dz = zt - np.sum(thick_arr_inv[:j])
                hz3[i] = rhoH0_arr_inv[j]*np.exp(-dz/hr)
                break
    
    
    # heat flux generated by radiogenic heat production
    qh = rhoH0_arr[-1]*hr**2*(1-np.exp(-Zbot/hr))/Zbot # W/m^3
    q = (Tm-T0)/Zbot*k_arr[0]
    
    qm = 30e-3 # W/m^2 heat flux from mantle
    qm = q - qh

    tr = continental_radiogenic_T(z,T0,hr,k_arr[0],qm,rhoH0_arr[-1])
    
    # calculate temperature increase by radiogenic heat production in the crust    
    dT = sum_dt_of_layers(thick_arr,k_arr,rhoH0_arr,hr)
    k = avg_k(k_arr,thick_arr)
    qm3 = (Tm-T0-dT)/sum(thick_arr)*k
    
    tr3 = continental_radiogenic_T2(z,T0,thick_arr,hr,k_arr,qm3,rhoH0_arr)

    print(f'temperature at surface: {tr[0]-273:.0f} degC')
    print(f'temperature at moho: {np.interp(Zmoho,z,tr)-273:.0f} degC')
    print(f'temperature at bottom: {tr[-1]-273:.0f} degC')
    
    fig, ax = plt.subplots(figsize=(4,8))

    ax.hlines(hr*1e-3,0,1500,ls='--',color='grey')
    ax.text(0.1,hr*1e-3+0.1,"hr",fontsize=12,color='grey')
    ax.vlines(Tm-273,0,Zbot/1e3,ls='--',color='grey')
    ax.text(Tm-273+0.1,0.1,"Tm",fontsize=12,color='grey')
    ax.hlines(Zmoho*1e-3,0,1500,ls='--',color='grey')
    ax.text(0.1,Zmoho*1e-3+0.1,"Zmoho",fontsize=12,color='grey')

    ax.plot(th-273, z/1e3,'r',label="Half-space cooling")


    ax.plot(tr-273, z/1e3,'k',label="Radiogenic")
    tmoho = np.interp(Zmoho,z,tr)
    ax.plot([tmoho-273],[Zmoho*1e-3],'o',color='k')
    ax.text(tmoho-273-10,Zmoho*1e-3+1,f"{tmoho-273:.0f}",fontsize=12,ha='right',va='top')

    ax.plot(tr3-273, z/1e3,'--b',label="Radiogenic in Zmoho")
    tmoho3 = np.interp(Zmoho,z,tr3)
    ax.plot([tmoho3-273],[Zmoho*1e-3],'o',color='b')
    ax.text(tmoho3-273+10,Zmoho*1e-3-1,f"{tmoho3-273:.0f}",fontsize=12,ha='left',va='bottom',color='b')
    ax.plot([tr3[-1]-273],[Zbot*1e-3],'o',color='b')
    ax.text(tr3[-1]-273-40,Zbot*1e-3-0.3,f"{tr3[-1]-273:.0f}",fontsize=12,ha='right',va='bottom',color='b')

    ax.text(300,90,f'q={qm*1e3:.0f} mW/m$^2$',fontsize=12,ha='left',va='top')
    ax.text(300,93,f'q={qm3*1e3:.0f} mW/m$^2$',fontsize=12,ha='left',va='top',color='b')


    ax2 = ax.twiny()
    
    ax2.plot(hz*1e6,z/1e3,'--',color='k',label="Heat production",lw=1)
    ax2.plot(hz3*1e6,z/1e3,'--',color='b',label="H only in crust",lw=1)
    ax2.set_xlabel(r"Heat production ($\mu$W/m3)")
    ax2.legend(loc='upper right')
    ax2.set_xlim(0,4)

    ax.legend(loc='lower left')
    ax.set_xlabel(r"Temperature ($^\circ{C}$)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim(Zbot/1e3,0)
    ax.set_xlim(0,1500)
    ax.grid(ls='--')
    fig.tight_layout()
    filename = f'geo-{H0:e}-{hr/1e3:.0f}km.png'
    fig.savefig(filename)
    
    return

if __name__ == "__main__":
    main()