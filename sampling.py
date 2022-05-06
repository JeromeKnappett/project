#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:58:14 2022

@author: jerome
"""

import matplotlib.pyplot as plt
import numpy as np


def minDistOverSFF(w,dx,lam):
    """
    Only valid for far-field!
    
    w: width of object
    dx: detector pixel width
    lam: wavelength
    
    returns:
        z: minimum propagation distance for oversampling
    """
    z = (2*w*dx)/lam
    return z

def maxWidthOverSFF(z,lam,dx):
    """
    Only valid for far-field!
    
    z: distance from object to detector
    lam: wavelength
    dx: detector pixel width
    
    returns:
        w: maximum width of object
    """  
    w = (z*lam)/(2*dx)
    return w

def objectResFF(lam,z,N,dx):
    """
    Only valid for far-field!
    
    lam: wavelength
    z: distance from sample to detector
    N: number of pixels in detector
    dx: detector pixel width
    
    returns:
        dxs: resolution at sample
    """
    dxs = (lam*z)/(N*dx)
    return dxs

def minDistOverSNF(w,dXd,dXo,lam,N):
    """
    Valid for Fresnel propagation.
    
    w: width of object
    dXd: detector plane pixel width
    dXo: object plane pixel width
    lam: wavelength
    N: number of pixels in detector
    
    returns:
        z: minimum propagation distance for oversampling
    """
    D = N*dXd
#    z = ((w*dXd + D*dXo)/lam)
    z = ((dXd*w)/lam) - ((dXo*D)/lam)
#    u = (w*N*dXd)
#    d = ((lam*(N - (D/dXd))))
#    z = u/d #(w*N*dXd)/((lam*(N - (D/dXd))))
    return z
    
def maxWidthOverSNF(z,lam,dXd,dXo,N):
    """
    Valid for Fresnel propagation
    
    z: distance from object to detector
    lam: wavelength
    dXd: detector plane pixel width
    dXo: object plane pizel width
    N: number of 
    
    returns:
        w: maximum width of object
    """  
    D = N*dXd
    w = ((D*dXo)/dXd) + ((lam*z)/dXd)
    return w

def objectResNF(lam,z,N,dx,w):
    """
    Valid for Fresnel propagation.
    
    lam: wavelength
    z: distance from sample to detector
    N: number of pixels in detector
    dx: detector pixel width
    w = width of object
    
    returns:
        dxs: resolution at sample
    """
    D = N*dx
    dxs = ((dx*w)/D)-((lam*z)/D)
    return dxs

def propReq(lam,z,N,dx):
    print("Checking Sampling Requirements... ")
    
    rCon = (2*np.sqrt(2)*lam*z)/(N*(dx**2))
    rFFT = (2*np.sqrt(2)*N*(dx**2))/(lam*z)
    rAngSpec = (2*np.sqrt(2)*lam*z)/((N*(dx**2))*np.sqrt(1-((lam**2)/(2*(dx**2)))))
        
    print("|=============================================|")
    print("|Propagator       ||     Requirment Met (y/n)|")
    print("|-----------------||-------------------------|")
    if rCon < 1:
        print("|Convolution      ||                     yes |")
    else:
        print("|Convolution      ||                      no |")
    if rFFT < 1:
        print("|FFT              ||                     yes |")
    else:
        print("|FFT              ||                      no |")
    if rAngSpec < 1:
        print("|Angular Spectrum ||                     yes |")
    else:
        print("|Angular Spectrum ||                      no |")
    print("|=============================================|")
    

def testSamplingFF():
    
    lam1 = 6.7e-9
    lam2 = 13.5e-9
    px = 11e-6
    N = 2048
    zMax = 1      # maximum allowable propagation distance
    zMin = 21e-3  # minimum allowable propagation distance
    
    Gx = 40e-6    # size of grating
    Mx = 120e-6   # size of four-grating mask
    
    deltaZ = np.linspace(zMin,0.1,1000)
    deltaW = np.linspace(Gx,Mx*4,1000)
    
    
    fig, ax = plt.subplots(1,3)
    print("FAR-FIELD SAMPLING REQUIREMENTS")
    print(" ")
    print(f"Detector resolution:                    {px*1e6} um")
    print(f"Number of pixels in detector:           {N}")
    print(f"Minimum allowable propagation distance: {zMin*1e3} mm")
    print(f"Maximum allowable propagation distance: {zMax} m")
    print(f"Size of grating:                        {Gx*1e6} um")
    print(f"Size of four grating mask:              {Mx*1e6} um")
    
    
    for l in (lam1,lam2):
        W_zmax = maxWidthOverSFF(zMax,l,px) # maximum object width for maximum propagation distance 
        W_zmin = maxWidthOverSFF(zMin,l,px) # maximum object width for minimum propagation distance
        
        z_g = minDistOverSFF(Gx,px,l)   # minimum propagation distance to oversample grating
        z_m = minDistOverSFF(Mx,px,l)   # minimum propagation distance to oversample mask
        
        dx_zmin = objectResFF(l,zMin,N,px)  # object resolution for minimum propagation distance
        dx_zmax = objectResFF(l,zMax,N,px)  # object resolution for maximum propagation distance
        
        
        print("===========================================================================")
        print(f'|| Beam wavelength:       {l*1e9:.3g} nm')
        print("||")
        print(f'|| Maximum object width for minimum propagation distance:    || {W_zmin*1e6:.3g} um')
        print(f'|| Maximum object width for maximum propagation distance:    || {W_zmax*1e6:.3g} um')
        print('|| -')
        print(f'|| Minumum propagation distance to oversample grating:       || {z_g*100:.3g} cm')
        print(f'|| Minumum propagation distance to oversample mask:          || {z_m*100:.3g} cm')
        print('|| -')
        print(f'|| Object resolution for minimum propagation distance:       || {dx_zmin*1e9:.3g} nm')
        print(f'|| Object resolution for maximum propagation distance:       || {dx_zmax*1e9:.3g}  nm')
        
        print("===========================================================================")
        print(" ")
        
        oWidth = [maxWidthOverSFF(z,l,px)*1e3 for z in deltaZ]
        propD = [minDistOverSFF(w,px,l)*1e3 for w in deltaW]
        oRes = [objectResFF(l,z,N,px)*1e9 for z in deltaZ]
        
        ax[0].plot(deltaZ,oWidth, label=f'$\lambda =${l*1e9:.3g} nm')
        ax[0].set_xlabel('Propagation Distance [m]')
        ax[0].set_ylabel('Maximum Object Width [mm]')
        ax[1].plot(deltaZ,oRes, label=f'$\lambda =${l*1e9:.3g} nm')
        ax[1].set_xlabel('Propagation Distance [m]')
        ax[1].set_ylabel('Object Resolution [nm]')
        ax[1].set_title('Far field')
        ax[2].plot(deltaW*1e6,propD, label=f'$\lambda =${l*1e9:.3g} nm')
        ax[2].set_xlabel('Object Width [microns]')
        ax[2].set_ylabel('Minimum Propagation Distance [mm]')
    plt.legend()
    plt.tight_layout()
    plt.plot()
    
def testSamplingNF():
    
    lam1 = 6.7e-9
    lam2 = 13.5e-9
    px = 11e-6
    N = 2048
    zMax = 1      # maximum allowable propagation distance
    zMin = 21e-3  # minimum allowable propagation distance
    
    Gx = 40e-6    # size of grating
    Mx = 120e-6   # size of four-grating mask
    Bx = 2e-3     # maximum extend of beam at mask plane
    
    deltaZ = np.linspace(zMin,0.1,1000)
    deltaW = np.linspace(Gx,Mx*4,1000)
    
    
    fig, ax = plt.subplots(1,3)
    print("NEAR-FIELD SAMPLING REQUIREMENTS")
    print(" ")
    print(f"Detector resolution:                    {px*1e6} um")
    print(f"Number of pixels in detector:           {N}")
    print(f"Minimum allowable propagation distance: {zMin*1e3} mm")
    print(f"Maximum allowable propagation distance: {zMax} m")
    print(f"Size of grating:                        {Gx*1e6} um")
    print(f"Size of four grating mask:              {Mx*1e6} um")
    for l in (lam1,lam2):
        
        dx_zmin = objectResFF(l,zMin,N,px)  # object resolution for minimum propagation distance
        dx_zmax = objectResFF(l,zMax,N,px)  # object resolution for maximum propagation distance
        
        dX_zmin = testN(l,zMin,px,Bx,N)  # object resolution for minimum propagation distance
        dX_zmax = testN(l,zMin,px,Bx,N)  # object resolution for maximum propagation distance
        
        W_zmax = maxWidthOverSNF(zMax,l,px,dx_zmax,N) # maximum object width for maximum propagation distance 
        W_zmin = maxWidthOverSNF(zMin,l,px,dx_zmin,N) # maximum object width for minimum propagation distance

        z_g = minDistOverSNF(Gx,px,dx_zmin,l,N)   # minimum propagation distance to oversample grating
        z_m = minDistOverSNF(Mx,px,dx_zmax,l,N)   # minimum propagation distance to oversample mask
#        w,dXd,dXo,lam,N
        
        print("===========================================================================")
        print(f'|| Beam wavelength:       {l*1e9:.3g} nm')
        print("||")
        print(f'|| Maximum object width for minimum propagation distance:    || {W_zmin*1e6:.3g} um')
        print(f'|| Maximum object width for maximum propagation distance:    || {W_zmax*1e6:.3g} um')
        print('|| -')
        print(f'|| Minumum propagation distance to oversample grating:       || {z_g*100:.3g} cm')
        print(f'|| Minumum propagation distance to oversample mask:          || {z_m*100:.3g} cm')
        print('|| -')
        print(f'|| Object resolution for minimum propagation distance:       || {dx_zmin*1e9:.3g} nm')
        print(f'|| Object resolution for maximum propagation distance:       || {dx_zmax*1e9:.3g}  nm')
        
        print("===========================================================================")
        
        oRes = [objectResNF(l,z,N,px,Mx)*1e9 for z in deltaZ]
        oRes1 = [objectResFF(l,z,N,px)*1e9 for z in deltaZ]
        oRes2 = [testN(l,z,px,Mx,N)*1e9 for z in deltaZ]
        
        oWidth = [maxWidthOverSNF(z,l,px,d,N)*1e3 for z,d in zip(deltaZ,oRes2)]
        propD = [minDistOverSNF(w,px,d,l,N)*1e3 for w,d in zip(deltaW,oRes2)]
        
        
        ax[0].plot(deltaZ,oWidth, label=f'$\lambda =${l*1e9:.3g} nm')
        ax[0].set_xlabel('Propagation Distance [m]')
        ax[0].set_ylabel('Maximum Object Width [mm]')
#        ax[1].plot(deltaZ,oRes, label=f'$\lambda =${l*1e9:.3g} nm')
        ax[1].plot(deltaZ,oRes1, ':', label=f'$\lambda =${l*1e9:.3g} nm - FF')
        ax[1].plot(deltaZ,oRes2, label=f'$\lambda =${l*1e9:.3g} nm - NF')
        ax[1].set_xlabel('Propagation Distance [m]')
        ax[1].set_ylabel('Object Resolution [nm]')
        ax[1].set_title('Near field')
        ax[1].legend()
        ax[2].plot(deltaW*1e6,propD, label=f'$\lambda =${l*1e9:.3g} nm')
        ax[2].set_xlabel('Object Width [microns]')
        ax[2].set_ylabel('Minimum Propagation Distance [mm]')
    plt.legend()
    plt.tight_layout()
    plt.plot()
    
    testN(lam2,zMin,px,Mx,N)

def testSampling():
    
    lam1 = 6.7e-9
    lam2 = 13.5e-9
    px = 11e-6
    N = 2048
    zMax = 1      # maximum allowable propagation distance
    zMin = 21e-3  # minimum allowable propagation distance
    
    Gx = 40e-6    # size of grating
    Mx = 120e-6   # size of four-grating mask
    
    deltaZ = np.linspace(zMin,zMax,1000)
    deltaW = np.linspace(Gx,Mx,1000)
    
    propReq(lam1,zMin,N,px)
    propReq(lam1,zMax,N,px)
    propReq(lam2,zMin,N,px)
    propReq(lam2,zMax,N,px)
    
    NearField = True
    FarField = True
    
    if NearField:
        testSamplingNF()
    else:
        pass
    if FarField:
        testSamplingFF()
    else:
        pass
    
def testN(lam,z,dx,w,N):
    D = dx*N
#    dxs_n = ((lam*z)*(1 - np.sqrt((1 + ((4*D*w)/(lam*z*N))))))/(2*D)
    dxs_p = ((lam*z)*(1 + np.sqrt((1 + ((4*D*w)/(lam*z*N))))))/(2*D)
    
#    print(dxs_n)
#    print(dxs_p)
    return dxs_p

if __name__ == "__main__":

#    testSamplingFF()
    testSampling()
#    testSamplingNF()