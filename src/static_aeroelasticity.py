## @ingroup Methods-Aeroelasticity
# static_aeroelasticity.py
#
# Created:  Dec 2018, S. Karpuk
# Modified:

import Aerodynamics as weiss
from Wing_preprocessing import Wing_preprocessing as wing_pre

# package imports
import math
import numpy as np
from numpy import linalg as LA
import matplotlib
import matplotlib.pyplot as plt


def compute_divergence_speed(General,Aerodynamics,Static,Wing,Analysis):
    """ Compute divergence speed for an arbitrary
        slender wing planform at the steady level flight

    References:
    Bisplinghoff, Ashley "Principles of Aeroelasticity"

    """

    # ------------------
    #  Unpack the input
    # ------------------
    rho      = General.aerodynamics.density
    T        = General.aerodynamics.temperature
    M        = General.aerodynamics.Mach
    Y_data   = General.structures.EIGJ.Y
    EI_data  = General.structures.EIGJ.EI
    GJ_data  = General.structures.EIGJ.GJ
    GK_data  = General.structures.EIGJ.GK
    EAsweep  = General.structures.EAsweep
    Y_for_AC = Wing.sections.Y 
    AC_dist  = Static.structures.distances.to_AC
    span     = Wing.span
    AR       = Wing.aspect_ratio
    Y        = Wing.sections.Y_centroid 
    chords   = Wing.sections.chords 
    phi      = Static.aerodynamics.wing.phi

    Nsect = len(Y)

    # Initialize required arrays
    EI      = np.zeros(Nsect)
    GJ      = np.zeros(Nsect)
    GK      = np.zeros(Nsect)
    e       = np.zeros((Nsect,Nsect))
    Ebar    = np.zeros((Nsect,Nsect))
    
    # Interpolate EIGJ and AC data
    EI      = np.interp(Y, Y_data, EI_data)
    GJ      = np.interp(Y, Y_data, GJ_data)
    GK      = np.interp(Y, Y_data, GK_data)
    
    # Compute influence coefficeints    
    Ctt,Ctz,Czz = compute_influence_coeffcieints(EAsweep,Y,EI,GJ,GK,Nsect)
    
    # Compute weighting matrix coefficients using the Malthopp's method
    W = Malthopp_quad(Nsect,phi,span)
    
    # Compute required arrays of the equation
    e = np.diagflat(AC_dist)

    # Calculate the dive speed 
    EPS   = 10**(-4)
    M_old = 0
    count = 0
    while abs(M-M_old) > EPS and count < 1000:
        if Analysis == 2.1:
            count += 1
            Aerodynamics, Static = weiss.strip_theory(Wing, Aerodynamics, General, Static, Analysis)
            CLa  = Aerodynamics.CLa
            Ebar = np.matmul(np.diag(chords),np.matmul(np.add(Ctz,np.matmul(Ctt,e)),W))
        
            # Find minimum a eigenvalue 
            v, w = LA.eig(Ebar)
            v_max = max(v)

            if v_max <= 0 or isinstance(v_max, complex):
                print('The wing does not diverge\n')
                Ud = 'None'
                break
            else:
                Ud = (1 / (0.5 * rho * v_max * CLa))**0.5 * 12

            M_old = M
            M = Ud / (1.4 * 1716 * T)**0.5
           
            '''print('Divergence speed = ' + str(Ud) + ' ft/s')
            print('                   ' + str(Ud * 0.592484) + ' KTAS')
            print(' Divergence Mach = ' + str(M) + '\n')'''
    
            # ------------------
            #  Pack the output
            # ------------------
            Static.Divergence_speed = Ud
            Static.Divergence_Mach  = M
            Wing.sections.EAtoAC = AC_dist
           
    return Static


def modified_aerodynamics(General,Static,Aerodynamics,Wing,analysis):

    """ Compute modified lift distribuiton and aerodynamic coef-s
        under the influence of aeroelasticity

    References:
    Bisplinghoff, Ashley "Principles of Aeroelasticity"
    
    """

    # -----------
    #  Unpack
    # -----------
    Y_data    = General.structures.EIGJ.Y
    EI_data   = General.structures.EIGJ.EI
    GJ_data   = General.structures.EIGJ.GJ
    GK_data   = General.structures.EIGJ.GK
    mg_data   = Static.structures.mg
    EAsweep   = General.structures.EAsweep 
    Y_for_AC  = Wing.sections.Y 
    AC_dist   = Static.structures.distances.to_AC
    CG_dist  = Static.structures.distances.to_CG
    phi       = Static.aerodynamics.wing.phi  
    span      = Wing.span
    area      = Wing.area
    CmAC      = Static.aerodynamics.wing.CmAC
    Cmac      = Wing.chords.cmac
    AR        = Wing.aspect_ratio
    Y         = Wing.sections.Y_centroid   
    chords    = Wing.sections.chords 
    q         = General.aerodynamics.dynamic_pressure
    N         = Static.aerodynamics.load_factor
    Weight    = Wing.weight
        
    Nsect = len(Y)
    
    EI      = np.zeros((Nsect,1))
    GJ      = np.zeros((Nsect,1))
    GK      = np.zeros((Nsect,1))
    Ainv    = np.zeros((Nsect,Nsect))
    c       = np.zeros((Nsect,1))
    ccl     = np.zeros(Nsect)
    ccl_el  = np.zeros((Nsect+1,1))
    cl      = np.zeros(Nsect)
    cl_r    = np.zeros(Nsect)
    cl_el   = np.zeros(Nsect)
    alpha_e = np.zeros((Nsect,1))
    El_coef = np.zeros((Nsect,Nsect))
    cmac    = np.zeros((Nsect,1))
    csq     = np.zeros((Nsect,Nsect))
    e       = np.zeros((Nsect,Nsect))
    f       = np.zeros((Nsect,Nsect))
    f1      = np.zeros((Nsect,Nsect))
    RHS     = np.zeros((Nsect,1))
    W1      = np.zeros((Nsect,1))
    W3      = np.zeros((1,Nsect+1))
    qW      = np.zeros((Nsect,1))

    # Compute rigid body lift distribution
    Aerodynamics, Statics = weiss.strip_theory(Wing, Aerodynamics, General, Static, analysis)
    A    = Static.aerodynamics.wing.A
    ccl_r = Static.aerodynamics.wing.ccl      
    
    # Convert geometric parameters from ft to in
    q = q / 144

    # Find mass per unit length
    mg = compute_mass_per_length(Nsect,Y_for_AC,mg_data)

    # Interpolate EIGJ data
    EI      = np.interp(Y, Y_data, EI_data)
    GJ      = np.interp(Y, Y_data, GJ_data)
    GK      = np.interp(Y, Y_data, GK_data)

    # Compute influence coefficeints
    Ctt,Ctz,Czz = compute_influence_coeffcieints(EAsweep,Y,EI,GJ,GK,Nsect)

    # Compute weighting matrix coefficients using the Malthopp's method
    W = Malthopp_quad(Nsect,phi,span)
    
    # Compute required arrays of the equation
    e    = np.diagflat(AC_dist)
    d    = np.diagflat(CG_dist)
    mg   = row_reshape(mg,Nsect)
    c    = np.diag(chords)
    cmac = row_reshape(CmAC,Nsect)

    # Compute the modified lift distribution
    #------------------------------------------------------
    # Prepare the RHS and RHS matrices
    Ebar = np.matmul(np.add(Ctz,np.matmul(Ctt,e)),W)
    Fbar = np.matmul(np.matmul(Ctt,np.square(c)),W)
    Gbar = np.matmul(np.add(Ctz,np.matmul(Ctt,d)),W)
    
    W1    = np.diag(W)
    f     = np.matmul(q*Ebar,np.transpose(ccl_r))
    f     = np.add(f,np.matmul(q*Fbar,cmac))
    f     = np.subtract(f,np.matmul(N*Gbar,mg))
    f1    = f
    f1    = np.insert(f,Nsect,0,0)
    LHS   = np.insert(np.subtract(A,q*Ebar), Nsect, -1, axis=1)    

    for i in range(Nsect):
        W3[0,i] = W[i][i]

    LHS    = np.insert(LHS, Nsect, W3, 0)
    ccl_el = np.linalg.solve(LHS, f1)
    ccl    = np.add(ccl_r,np.transpose(ccl_el[0:Nsect][:]))

    for i in range(Nsect):
        cl[i]    = ccl[0][i]/chords[i] 
        cl_el[i] = ccl_el[i][0]/chords[i]
        cl_r[i]  = ccl_r[0][i]/chords[i]

    # Compute rigit and elastic stability coefficients
    El_coef = np.subtract(A,q*Ebar)
    El_coef = np.linalg.inv(np.add(El_coef,2*q/Weight * \
                                   np.matmul(np.matmul(np.matmul(Gbar,mg),np.ones((1,Nsect))),W)))

    CLa_r  = np.matmul(np.matmul(2/area*np.matmul(np.ones(Nsect),W),np.linalg.inv(A)),np.ones((Nsect,1)))
    CLa_el = np.matmul(np.matmul(2/area*np.matmul(np.ones(Nsect),W),El_coef),np.ones((Nsect,1)))

    # Calculate elastic angle-of-attack
    alpha_e = np.add(np.matmul(q*Ebar,ccl_el[:][0:Nsect]),f)
    for i in range(Nsect):
        alpha_e[:][i] = math.degrees(alpha_e[:][i])

    
    #Pack
    Static.aerodynamics.wing.ccl         = ccl
    Static.aerodynamics.wing.ccl_r       = ccl_r
    Static.aerodynamics.wing.elast_alpha = alpha_e
    Static.aerodynamics.wing.cl          = cl
    Static.aerodynamics.wing.cl_r        = cl_r
    Static.aerodynamics.wing.cl_el       = cl_el
    Static.aerodynamics.wing.ccl_el      = np.transpose(ccl_el)
    Static.aerodynamics.wing.CLa_r       = CLa_r
    Static.aerodynamics.wing.CLa_el      = CLa_el
    
    return Static


def elastic_control(General,Aerodynamics,Static,Wing,Analysis):
    """ Compute control reversal speed for an arbitrary
        slender wing planform at the steady level flight

    References:
    Bisplinghoff, Ashley "Principles of Aeroelasticity"

    """
    
    # ------------------
    #  Unpack the input
    # ------------------
    q        = General.aerodynamics.dynamic_pressure
    rho      = General.aerodynamics.density
    T        = General.aerodynamics.temperature
    M        = General.aerodynamics.Mach
    Y_data   = General.structures.EIGJ.Y
    EI_data  = General.structures.EIGJ.EI
    GJ_data  = General.structures.EIGJ.GJ
    GK_data  = General.structures.EIGJ.GK
    EAsweep  = General.structures.EAsweep
    Y_for_AC = Wing.sections.Y
    AC_dist  = Static.structures.distances.to_AC
    span     = Wing.span
    AR       = Wing.aspect_ratio
    Y        = Wing.sections.Y_centroid   
    chords   = Static.aerodynamics.wing.chords
    phi      = Static.aerodynamics.wing.phi
    Cla      = Aerodynamics.Cla

    Nsect = len(Y)

    # Initialize required arrays
    adel      = np.zeros((1,Nsect))
    Cm_beta   = np.zeros((1,Nsect))
    H         = np.zeros((Nsect,Nsect))
    y         = np.zeros((Nsect,Nsect))
    c         = np.zeros((Nsect,Nsect))
    numMat    = np.zeros((Nsect,Nsect))
    denumMat  = np.zeros((Nsect,Nsect))
    Gbar      = np.zeros((Nsect,Nsect))
    Ebar      = np.zeros((Nsect,Nsect))
    Fbar      = np.zeros((Nsect,Nsect))
    

    # Convert geometric parameters from ft to in
    q = q / 144
    
    # Calculate required stability and control coefficients
    Wing, cf_c = wing_pre.control_surfaces_calc(Wing,Static,Aerodynamics)

    for i in range(Nsect):
        theta = math.acos(2*cf_c[i] - 1)
        if theta == math.pi:
            miu = 0
        else:
            miu = 1.15 * 0.5 * (1 - cf_c[i]) * math.sin(theta)/(math.pi - (theta - math.sin(theta)))
        adel[0][i]    = 0.85 * (1 - (theta - math.sin(theta))/math.pi)
        Cm_beta[0][i] = -miu * adel[0][i] * Cla

    # Interpolate EIGJ and AC data
    EI  = np.interp(Y, Y_data, EI_data)
    GJ  = np.interp(Y, Y_data, GJ_data)
    GK  = np.interp(Y, Y_data, GK_data)
    
    # Compute influence coefficeints   
    Ctt,Ctz,Czz = compute_influence_coeffcieints(EAsweep,Y,EI,GJ,GK,Nsect)

    # Compute weighting matrix coefficients using the Malthopp's method
    W = Malthopp_quad(Nsect,phi,span)
    
    # Compute required arrays of the equation
    Aerodynamics, Static = weiss.strip_theory(Wing, Aerodynamics, General, Static, Analysis)
    A = Static.aerodynamics.wing.A
    
    e = np.diagflat(AC_dist)
    c = np.diagflat(chords)
    y = np.diagflat(Y)
    H = np.matmul(np.matmul(np.ones(Nsect),y),W)

    Fbar = np.matmul(np.matmul(Ctt,np.square(c)),W)
    Ebar = np.matmul(np.add(Ctz,np.matmul(Ctt,e)),W)
    Gbar = np.matmul(H,np.linalg.inv(A-q*Ebar))

    # Compute the elastic aileron effectiveness
    numMat    = np.matmul(q*Fbar,np.transpose(Cm_beta))
    numMat    = np.add(np.transpose(adel),numMat)
    num       = np.matmul(Gbar,numMat)
    denum     = np.matmul(Gbar,np.transpose(-Y/(0.5*span)))
    ailer_eff = -num/denum

    # Compute the control reversal speed        
    qr_old = 5000
    qr     = 0
    EPS    = 0.001
    count  = 0
    while abs(qr - qr_old) > EPS and count < 1000:
        count += 1
        qr_old = qr
        Aerodynamics, Static = weiss.strip_theory(Wing, Aerodynamics, General, Static, Analysis)
        A     = Static.aerodynamics.wing.A 
        Gbar  = np.matmul(H,np.linalg.inv(np.subtract(A,qr*Ebar)))
        num   = np.matmul(Gbar,np.transpose(adel))
        denum = np.matmul(Gbar,np.matmul(Fbar,np.transpose(Cm_beta)))
        qr    = -num/denum
        Ur    = (2*qr*144/rho)**0.5
        Mr     = Ur / (1.4 * 1716 * T)**0.5

        General.aerodynamics.Mach = M

    # ------------------
    #  Pack the output
    # ------------------
    Static.control_reversal_dynamic_pressure = 144 * qr
    Static.Control_reversal_speed            = Ur 
    Static.Control_reversal_Mach             = Mr
    Static.Aileron_effectiveness             = ailer_eff
    
    return Static

def compute_influence_coeffcieints(EAsweep,Y,EI,GJ,GK,Nsect):
    """ Calculate structural influence coefficients

    References:
    Bisplinghoff, Ashley "Principles of Aeroelasticity"

    """

    Ctt  = np.zeros((Nsect,Nsect))
    Ctz  = np.zeros((Nsect,Nsect))
    Czz  = np.zeros((Nsect,Nsect))
    
    for i in range(Nsect-1):
        for j in range(Nsect-1):
            if (Y[i] >= Y[j]):
                for k in range(j+1,Nsect):
                    Ctt[i][j] = Ctt[i][j] + 0.5 * ((math.cos(math.radians(EAsweep[k])))**2 / GJ[k] + \
                                                   (math.sin(math.radians(EAsweep[k])))**2 / EI[k] + \
                                                   (math.cos(math.radians(EAsweep[k-1])))**2 / GJ[k-1] + \
                                                   (math.sin(math.radians(EAsweep[k-1])))**2 / EI[k-1]) * abs(Y[k]/math.cos(math.radians(EAsweep[k]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))
                   
                    Ctz[i][j] = Ctz[i][j] + 0.5 * ((Y[j]/math.cos(math.radians(EAsweep[k]))-Y[k]/math.cos(math.radians(EAsweep[k])))/EI[k] + \
                                                   (Y[j]/math.cos(math.radians(EAsweep[k-1]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))/EI[k-1]) * \
                                                                                                           abs(Y[k]/math.cos(math.radians(EAsweep[k]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))

                    A1 = (Y[j]/math.cos(math.radians(EAsweep[k]))-Y[k]/math.cos(math.radians(EAsweep[k])))*(Y[i]/math.cos(math.radians(EAsweep[k]))-Y[k]/math.cos(math.radians(EAsweep[k])))/EI[k] + 1/GK[k]
                    A2 = (Y[j]/math.cos(math.radians(EAsweep[k-1]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))*(Y[i]/math.cos(math.radians(EAsweep[k-1]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))/EI[k-1] + 1/GK[k-1]
                    
                    Czz[i][j] = Czz[i][j] + 0.5 * (A1 + A2)*abs(Y[k]/math.cos(math.radians(EAsweep[k]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))

                Ctz[i][j] = - math.sin(math.radians(EAsweep[i])) * Ctz[i][j]
                
            else:
                for k in range(i+1,Nsect):
                    Ctt[i][j] = Ctt[i][j] + 0.5 * ((math.cos(math.radians(EAsweep[k])))**2 / GJ[k] + \
                                                   (math.sin(math.radians(EAsweep[k])))**2 / EI[k] + \
                                                   (math.cos(math.radians(EAsweep[k-1])))**2 / GJ[k-1] + \
                                                   (math.sin(math.radians(EAsweep[k-1])))**2 / EI[k-1]) * abs(Y[k]/math.cos(math.radians(EAsweep[k]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))

                    Ctz[i][j] = Ctz[i][j] + 0.5 * ((Y[j]/math.cos(math.radians(EAsweep[k]))-Y[k]/math.cos(math.radians(EAsweep[k])))/EI[k] + \
                                                   (Y[j]/math.cos(math.radians(EAsweep[k-1]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))/EI[k-1]) * \
                                                                                                           abs(Y[k]/math.cos(math.radians(EAsweep[k]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))
                    
                    A1 = (Y[j]/math.cos(math.radians(EAsweep[k]))-Y[k]/math.cos(math.radians(EAsweep[k])))*(Y[i]/math.cos(math.radians(EAsweep[k]))-Y[k]/math.cos(math.radians(EAsweep[k])))/EI[k] + 1/GK[k]
                    A2 = (Y[j]/math.cos(math.radians(EAsweep[k-1]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))*(Y[i]/math.cos(math.radians(EAsweep[k-1]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))/EI[k-1] + 1/GK[k-1]
                    
                    Czz[i][j] = Czz[i][j] + 0.5 * (A1 + A2)*abs(Y[k]/math.cos(math.radians(EAsweep[k]))-Y[k-1]/math.cos(math.radians(EAsweep[k-1])))
                   
                Ctz[i][j] = - math.sin(math.radians(EAsweep[i])) * Ctz[i][j]

    return Ctt, Ctz, Czz

def Malthopp_quad(Nsect,phi,span):   
    """ Create a matrix required for the Malthopp's quadrature

    References:
    Bisplinghoff, Ashley "Principles of Aeroelasticity"

    """
    
    W    = np.zeros((Nsect,Nsect))
    for i in range(Nsect):
        W[i][i] = math.sin(phi[i])
    W[Nsect-1,Nsect-1] = 0.5 * W[Nsect-1,Nsect-1]
    W = 0.5 * math.pi * span / (2* Nsect) * W

    return W

def row_reshape(mat,Nsect):
    """ Reshape the row matrix into the column 

    """
    
    mat1 = np.zeros((Nsect,1))
    
    for i in range(Nsect):
        mat1[i][0] = mat[i]

    return mat1

def compute_mass_per_length(Nsect,Y_for_AC,mg_data):
    """ Calculate mass per unit length of lumped masses along the wing 

    """

    mg = np.zeros((Nsect,1))
    
    mg[0] = mg_data[0]/(2*Y_for_AC[1])
    for i in range(1,Nsect):
        mg[i] = mg_data[i]/(Y_for_AC[i+1]-Y_for_AC[i])
    mg = mg[::-1]
    
    return mg

def compute_MOI_per_length(Nsect,Y_for_AC,Iyy_data):
    """ Calculate MOI per unit length of lumped masses along the wing 

    """   

    Iyy = np.zeros((Nsect,1))
    
    Iyy[0] = Iyy_data[0]/(2*Y_for_AC[1])
    for i in range(1,Nsect):
        Iyy[i] = Iyy_data[i]/(Y_for_AC[i+1]-Y_for_AC[i])
    Iyy = Iyy[::-1]

    return Iyy
