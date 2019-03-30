## @ingroup Methods-Aeroelasticity
# dynamic_aeroelasticity.py
#
# Created:  Feb 2019, S. Karpuk
# Modified:

import Aerodynamics as weiss
from static_aeroelasticity import compute_influence_coeffcieints
from static_aeroelasticity import compute_mass_per_length
from static_aeroelasticity import compute_MOI_per_length
from scipy import special as spec

# package imports
import math
import numpy as np
from numpy import linalg as LA
import matplotlib
import matplotlib.pyplot as plt

def compute_flutter_speed(General,Dynamic,Wing):

    """ Compute flutter speed using a U-g method

    References:
    Bisplinghoff, Ashley "Principles of Aeroelasticity"

    """

    # ---------------
    #  Unpack input
    # ---------------
    rho      = General.aerodynamics.density
    sound_sp = 12*(1.4*1716*General.aerodynamics.temperature)**0.5
    Y_data   = Wing.sections.Y
    Iy_data  = Dynamic.structures.Iyy
    mg_data  = Dynamic.structures.mg
    Y        = Wing.sections.Y_centroid
    chords   = Wing.sections.chords
    span     = Wing.span
    C2_sweep = General.structures.C2sweep
    e        = Dynamic.structures.distances.to_AC
    d        = Dynamic.structures.distances.to_CG

    N = len(Y)

    # Convert input to inches
    rho = rho / 12**3

    e = e[::-1]
    d = d[::-1]
    
    # Compute the system's natural model and frequencies
    phi_h,omega_h,phi_a,omega_a = compute_natural_modes(General,Wing,Dynamic,N)

    # Generate necessary inputs for the flutter stability determinant evaluation
    #--------------------------------------------------------------------------------

    Sa      = np.zeros(N)
    COEFS   = np.zeros(3)
    a       = np.zeros(N)
    b       = np.zeros(N)
    y       = np.zeros(N)
    A1_coef = np.zeros(N, dtype=np.complex)
    A2_coef = np.zeros(N, dtype=np.complex)
    B1_coef = np.zeros(N, dtype=np.complex)
    B2_coef = np.zeros(N, dtype=np.complex)
    C1_coef = np.zeros(N, dtype=np.complex)
    C2_coef = np.zeros(N, dtype=np.complex)
    D1_coef = np.zeros(N, dtype=np.complex)
    D2_coef = np.zeros(N, dtype=np.complex)

    y = Y / (0.5 * span)

    # Find mass and MOI per unit length
    Nsect = len(Y)

    M     = compute_mass_per_length(Nsect,Y_data,mg_data)
    Iy    = compute_MOI_per_length(Nsect,Y_data,Iy_data)

    # Normalize the mass and inertia data and mass moment about the mid-chord
    for i in range(N):
        M[i][0]  = M[i][0]/32.2
        Iy[i][0] = Iy[i][0]/32.2
        Sa[i]    = -M[i] * d[i]

    # Set up the reference station (3/4 of the wing span)
    b       = 0.5 * chords
    a       = np.divide(e+0.25*chords-b, b)
    br      = 0.5 * np.interp(3/8*span,Y[::-1],chords[::-1])
    kr      = np.linspace(0.02, 1, 500)
    omega   = np.zeros((2,len(kr)))
    Mach    = np.zeros((2,len(kr)))
    g       = np.zeros((2,len(kr)))
    Uf      = np.zeros(1)
    Mf      = np.zeros(1)
    omega_f = np.zeros(1)

    # evaluate integral terms of the flutter determinant
    for j in range(len(kr)):
        K1_Lh, K2_Lh, K1_La, K2_La, K3_La, K1_Mh, K1_Ma, K2_Ma = Unsteady_Aero_coefs(kr[j])
        for i in range(N-1):
            ratio = br/b[i]
            Lh = K1_Lh+ratio*K2_Lh
            La = K1_La+ratio*K2_La+ratio**2*K3_La
            Mh = K1_Mh
            Ma = K1_Ma+ratio*K2_Ma
            A1_coef[i] = M[i]/(math.pi*rho*br**2)*(phi_h[i][0])**2
            A2_coef[i] = math.cos(math.radians(C2_sweep[i]))*(b[i]/br)**2*Lh*(phi_h[i][0])**2
            B1_coef[i] = - Sa[i]/(math.pi*rho*br**3)*(phi_h[i][0])*(phi_a[i][0])
            B2_coef[i] = - math.cos(math.radians(C2_sweep[i]))*(b[i]/br)**3 \
                         *(La-Lh*(0.5+a[i]))*phi_h[i][0]*phi_a[i][0]
            C1_coef[i] = B1_coef[i]
            C2_coef[i] = - math.cos(math.radians(C2_sweep[i]))*(b[i]/br)**3 \
                         *(Mh-Lh*(0.5+a[i]))*phi_h[i][0]*phi_a[i][0]
            D1_coef[i] = Iy[i]/(math.pi*rho*br**4)*phi_a[i][0]**2
            D2_coef[i] = math.cos(math.radians(C2_sweep[i]))*(b[i]/br)**4* \
                         (Ma-(Mh+La)*(0.5+a[i])+Lh*(0.5+a[i])**2)*(phi_a[i][0])**2

        A1 = -np.trapz(A1_coef,y)
        A2 = -np.trapz(A2_coef,y)
        B1 = -np.trapz(B1_coef,y)
        B2 = -np.trapz(B2_coef,y)
        C1 = -np.trapz(C1_coef,y)
        C2 = -np.trapz(C2_coef,y)
        D1 = -np.trapz(D1_coef,y)
        D2 = -np.trapz(D2_coef,y)

        COEFS       = [ A1*D1*(omega_h[0]/omega_a[0])**2 , \
                        -((1+(omega_h[0]/omega_a[0])**2)*A1*D1+(omega_h[0]/omega_a[0])**2*A1*D2+A2*D1), \
                        (A1+A2)*(D1+D2)-(C1+C2)*(B1+B2) ]
        
        X           = np.roots(COEFS)
        omega[0][j] = omega_a[0] / (X[0].real)**0.5
        omega[1][j] = omega_a[0] / (X[1].real)**0.5
        g[0][j]     = X[0].imag / (X[0].real)**0.5
        g[1][j]     = X[1].imag / (X[1].real)**0.5
        Mach[0][j]  = omega[0][j] * br / kr[j] / sound_sp
        Mach[1][j]  = omega[1][j] * br / kr[j] / sound_sp

    for i in range(len(kr)-1):
        if np.sign(g[0][i]) != np.sign(g[0][i+1]):
            if len(Uf) == 1:
                Mf[0] = Mach[0][i]
                Uf[0] = Mf[0] * sound_sp
                omega_f = omega[0][i]
            else:
                Mf      = np.append(Mf,Mach[0][i])
                omega_f = np.append(omega_f,omega[0][i])
                
        elif np.sign(g[1][i]) != np.sign(g[1][i+1]):
            if len(Uf) == 1:
                Mf[0] = Mach[1][i]
                Uf[0] = Mf[0] * sound_sp
                omega_f = omega[1][i]
            else:
                Mf      = np.append(Mf,Mach[1][i])
                omega_f = np.append(omega_f,omega[1][i])

    # ---------------
    #  Pack output
    # ---------------
    Dynamic.bending_modes      = phi_h
    Dynamic.torsion_modes      = phi_a
    Dynamic.bending_frequency  = omega_h
    Dynamic.torsion_frequency  = omega_a
    Dynamic.structural_damping = g
    Dynamic.Mach               = Mach
    
    if Uf[0] == 0:
        print('Flutter does not occur')
        Dynamic.Flutter_Mach = 'does not occur'
        Dynamic.Flutter_speed = 'does not occur'
        Dynamic.Flutter_natural_frequiency = 'does not occur'     
    else:
        Dynamic.Flutter_Mach = Mf
        Dynamic.Flutter_speed = Uf
        Dynamic.Flutter_natural_frequiency = omega_f
    
    return Dynamic


def compute_natural_modes(General,Wing,Dynamic,Nsect):

    """ Compute natural modes of the structure using
        an integral method for concentrated mass systems

    References:
    Bisplinghoff, Ashley "Principles of Aeroelasticity"

    """

    # ---------------
    #  Unpack input
    # ---------------
    Y_data    = General.structures.EIGJ.Y
    EI_data   = General.structures.EIGJ.EI
    GJ_data   = General.structures.EIGJ.GJ
    GK_data   = General.structures.EIGJ.GK
    EAsweep   = General.structures.EAsweep
    Y_struct  = Wing.sections.Y
    Y         = Wing.sections.Y_centroid
    Mass_data = Dynamic.structures.mg
    Iyy_data  = Dynamic.structures.Iyy

    EI   = np.zeros((1,Nsect))
    GJ   = np.zeros((1,Nsect))
    Mass = np.zeros((Nsect,Nsect))
    Iyy  = np.zeros((Nsect,Nsect))
    
    # Interpolate EIGJ data
    EI   = np.interp(Y, Y_data, EI_data)
    GJ   = np.interp(Y, Y_data, GJ_data)
    GK   = np.interp(Y, Y_data, GK_data)

    # Generate a diagonal mass and intertia matrices
    Mass_data = Mass_data[::-1]
    Iyy_data  = Iyy_data[::-1]
    for i in range(Nsect):
        Mass[i][i] = Mass_data[i] * 1/386.4
        Iyy[i][i]  = Iyy_data[i] * 1/386.4

    # Compute influence coefficeints
    Ctt,Ctz,Czz = compute_influence_coeffcieints(EAsweep,Y,EI,GJ,GK,Nsect)
    
    # Compute natural modes and frequencies
    phi_h, omega_h = nat_freq_loop(Czz,Mass,Nsect)
    phi_a, omega_a = nat_freq_loop(Ctt,Iyy,Nsect)    
    
    return phi_h, omega_h, phi_a, omega_a


def nat_freq_loop(C,M,Nsect):
    
    D       = np.ones((Nsect,Nsect))
    W       = np.ones((Nsect-1,1))
    W_old   = np.ones((Nsect-1,1))
    phi     = np.zeros((Nsect,1))

    EPS        = 1E-7
    error      = 1
    omega_old  = 1
    iter_count = 0
    
    D = np.matmul(C,M)

    while error >= EPS:
        W          = np.matmul(D[0:Nsect-1,0:Nsect-1],W_old)
        omega      = max(W)
        W          = W / omega
        error      = abs(omega - omega_old)
        W_old      = W
        omega_old  = omega
        iter_count += 1

    omega = (1/omega)**0.5
    
    # write a complete matrix of natural modes
    for i in range(Nsect-1):
        phi[i][0] = W[i][0]
    phi[Nsect-1][0] = 0

    return phi, omega


def Unsteady_Aero_coefs(k):
    """ Computes Required unsteady aerodynamic
        coefficients
    """

    F, G = Theodorsen_functions(k)

    K1_Lh = 1
    K2_Lh = -2j/k*(F+1j*G)
    K1_La = 0.5
    K2_La = -1j/k-2j/k*(F+1j*G)
    K3_La = -2/k**2*(F+1j*G)
    K1_Mh = 0.5
    K1_Ma = 3/8
    K2_Ma = -1j/k

    return K1_Lh, K2_Lh, K1_La, K2_La, K3_La, K1_Mh, K1_Ma, K2_Ma


def Theodorsen_functions(i):
    """ Computes Theodorsen's F and G functions
        for a given reduced frequency
    """
    
    F = (spec.j1(i)*(spec.j1(i)+spec.y0(i))+spec.y1(i)*(spec.y1(i)-spec.j0(i)))/ \
           ((spec.j1(i)+spec.y0(i))**2 + (spec.y1(i)-spec.j0(i))**2)
    
    G = -(spec.y1(i)*spec.y0(i)+spec.j1(i)*spec.j0(i))/  \
           ((spec.j1(i)+spec.y0(i))**2 + (spec.y1(i)-spec.j0(i))**2)
 
    return F, G
        
