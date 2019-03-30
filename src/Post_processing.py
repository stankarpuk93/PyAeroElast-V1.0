## Post_processing

# Author  : Stanislav Karpuk
#           (stankarpuk93@gmail.com)
# Created : 01/17/2019
# Modified:

from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import sys
import os

def output(directory, General, Wing, Aerodynamics, Static, Dynamic, Analysis):
    """ Create output documents and plots

    """
    
    sys.path.insert(0, directory)
    os.chdir(directory)
    f = open('Output_summary.dat','w')

    # GENERAL OUTPUT FOR ANY ANALYSIS TYPE
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------

    # Output the geometric results
    f.write('Wing geometric properties: \n')
    f.write('            Span = ' + str(round(Wing.span,2)) + ' in\n')
    f.write('              AR = ' + str(round(Wing.aspect_ratio,2)) + '\n')
    f.write('            Area = ' + str(round(Wing.area, 2)) + ' in\n') 
    f.write('      Root chord = ' + str(round(Wing.chords.chords[0], 2)) + ' in\n')
    f.write('       Tip chord = ' + str(round(Wing.chords.chords[len(Wing.chords.chords)-1], 2)) + ' in\n')
    f.write('  Root incidence = ' + str(round(Wing.incidence[0], 2)) + ' deg\n\n')   

    # Output structural properties of the wing
    try:
        fig, ax1 = plt.subplots()
        ax1.plot(General.structures.EIGJ.Y, General.structures.EIGJ.EI, label = 'EI')
        ax1.plot(General.structures.EIGJ.Y, General.structures.EIGJ.GJ, label = 'GJ')
        ax1.plot(General.structures.EIGJ.Y, General.structures.EIGJ.GK, label = 'GK')
        ax1.set_xlabel('Y, in')
        ax1.set_ylabel('Parameter, psi')
        ax1.set_title('Wing section properties')
        ax1.legend()
        ax1.grid(True)
        plt.savefig('Section_properties.png')
    except:
        pass

    # Plot the wing planform
    plot_3D_wing(Wing,Static,Dynamic)

    
    # ANALYSIS OUTPUT
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    for i in range(len(Analysis)):
        
        # Aerodynamic analysis output
        if Analysis[i] == 1:        
            alpha  = Aerodynamics.alpha
            alpha0 = Aerodynamics.a0
            CLa    = Aerodynamics.CLa
            CL     = Aerodynamics.wing.CL
            CDi    = Aerodynamics.wing.CDi
            ccl    = Aerodynamics.ccl
            cl     = Aerodynamics.cl

            f.write('Aerodynamic analysis output: \n')
            f.write('                     Mach = ' + str(round(General.aerodynamics.Mach, 2)) + '\n')
            f.write('                 Altitude = ' + str(round(General.aerodynamics.altitude, 2)) + ' in\n')   
            f.write('         Lift-curve slope = ' + str(round(CLa,2)) + '\n\n')

            # Plot the lift and drag curves
            fig, ax2 = plt.subplots()
            ax2.plot(alpha,CL, 'r', label = 'CL')
            ax2.grid(True)
            ax2.set_xlabel('AOA, deg')
            ax2.set_ylabel('CL', color = 'r')
            ax2.set_title('Wing Lift Curve')

            ax3 = ax2.twinx()
            ax3.plot(alpha,CDi, 'b',  label = 'CDi')
            ax3.set_ylabel('CDi', color = 'b')
            plt.tight_layout()
            plt.savefig('Lift_and Drag_curves.png')

            # Plot lift distribution
            plot_lift_distribution(alpha,cl,Wing,Static,Analysis[i])
            
        # Wing divergence analysis output
        if Analysis[i] == 2.1:
            f.write('Divergence analysis output: \n')
            f.write('                 Altitude = ' + str(round(General.aerodynamics.altitude, 2)) + ' in\n')
            try:
                f.write('         Divergence speed = ' + str(round(Static.Divergence_speed,2)) + ' fps\n')
                f.write('                            ' + str(round(0.592484*Static.Divergence_speed,2)) + ' KTAS\n')
                f.write('          Divergence Mach = ' + str(round(Static.Divergence_Mach,2)) + '\n\n')
            except:
                f.write('         Divergence speed = ' + Static.Divergence_speed + '\n')
                f.write('          Divergence Mach = ' + Static.Divergence_Mach  + '\n\n')                

        # Modified lift distribution 
        if Analysis[i] == 2.2:
            alpha = Aerodynamics.maneuver_alpha
            cl    = Static.aerodynamics.wing.cl_r
            Y     = Wing.sections.Y_centroid

            f.write('Modified lift distribution due to symmetric maneuver output: \n')
            f.write('                     Mach = ' + str(round(General.aerodynamics.Mach, 2)) + '\n')
            f.write('                 Altitude = ' + str(round(General.aerodynamics.altitude, 2)) + ' in\n')
            f.write('              Load factor = ' + str(round(Static.aerodynamics.load_factor, 2)) + '\n')
            f.write('   Rigid Lift-curve slope = ' + str(np.round(Static.aerodynamics.wing.CLa_r, 2)) + '\n')
            f.write(' Elastic Lift-curve slope = ' + str(np.round(Static.aerodynamics.wing.CLa_el, 2))+ '\n\n')

            # Plot lift distribution
            plot_lift_distribution(alpha,cl,Wing,Static,Analysis[i])
            plot_Elastic_3D_wing(Wing,Static)

        # Wing control reversal analysis output
        if Analysis[i] == 2.3:
            f.write('Control efficeincy analysis output: \n')
            f.write('                     Mach = ' + str(round(General.aerodynamics.Mach , 2)) + '\n')
            f.write('                 Altitude = ' + str(round(General.aerodynamics.altitude, 2)) + ' ft\n')
            f.write('  Control effect-s coef-t = ' + str(round(Static.Aileron_effectiveness[0],2)) + '\n')
            f.write('   Control reversal speed = ' + str(round(Static.Control_reversal_speed[0],2)) + ' fps\n')
            f.write('                            ' + str(round(0.592484*Static.Control_reversal_speed[0],2)) + ' KTAS\n')
            f.write('    Control reversal Mach = ' + str(round(Static.Control_reversal_Mach[0],2)) + '\n')

        if Analysis[i] == 3:
            f.write('Flutter analysis: \n')
            f.write('                 Altitude = ' + str(round(General.aerodynamics.altitude, 2)) + ' ft\n')
            try:
                f.write('            Flutter speed = ' + str(round(Dynamic.Flutter_speed[0],2)/12) + ' fps\n')
                f.write('                          = ' + str(round(Dynamic.Flutter_speed[0],2)*0.0493737) + ' KTAS\n')
                f.write('             Flutter Mach = ' + str(round(Dynamic.Flutter_Mach[0],2)) + '\n\n')
            except:
                f.write('            Flutter speed = ' + Dynamic.Flutter_speed + '\n')
                f.write('             Flutter Mach = ' + Dynamic.Flutter_Mach + '\n\n')                

            plot_flutter_data(Wing,Dynamic)

    f.close()
    
    return

def plot_2D_wing(Wing,Static,ax):
    """ Generates the wing planform plot with an
        AC and EA lines and the aileron (if defined)
    """
    
    #-----------------
    # Unpack inputs
    #-----------------
    Y        = Wing.sections.Y_centroid
    Ytip     = Wing.chords.Y[len(Wing.chords.Y)-1]
    chordtip = Wing.chords.chords[len(Wing.chords.chords)-1]
    chords   = Wing.sections.chords
    LE_sweep = Wing.sections.LEsweep

    # Append the tip parameters
    Y      = np.insert(Y,0,Ytip)
    chords = np.insert(chords,0,chordtip)

    Nsect = len(Y)
    
    # Create the wing 
    X = np.zeros((5,Nsect))
    
    X[0][Nsect-1] = 0
    X[1][Nsect-1] = chords[Nsect-1]
    X[2][Nsect-1] = 0.25 * chords[Nsect-1]

    for i in range(Nsect-2,-1, -1):
        X[0][i] = X[0][i+1] + (Y[i]-Y[i+1])*math.tan(math.radians(LE_sweep[i]))
        X[1][i] = X[0][i] + chords[i]
        X[2][i] = X[0][i] + 0.25 * chords[i]

    # Draw the wing
    ax.plot(Y, X[0][:], 'k')
    ax.plot(Y, X[1][:], 'k')
    ax.plot([Y[0],Y[0]],[X[0][0],X[1][0]], 'k')
    ax.plot([Y[Nsect-1],Y[Nsect-1]],[X[0][Nsect-1],X[1][Nsect-1]], 'k')
    ax.plot(Y, X[2][:], 'r', label = 'Aero center')

    # Draw an elastic axis
    try:
        EAtoAC = Static.structures.distances.to_AC
        EAtoAC = EAtoAC[::-1]
        f      = interpolate.interp1d(Y[1:Nsect], EAtoAC , fill_value='extrapolate')
        EAtoAC = np.insert(EAtoAC,0,f(Ytip))
        
        X[3][Nsect-1] = X[2][Nsect-1] + EAtoAC[Nsect-1]
        for i in range(Nsect-2,-1, -1):
            X[3][i] = X[2][i] + EAtoAC[i]
        ax.plot(Y, X[3][:], 'b', label = 'Elastic axis')
    except:
        pass

    # Draw a cg axis
    try:
        EAtoCG = Static.structures.distances.to_CG
        EAtoCG = EAtoCG[::-1]
        f      = interpolate.interp1d(Y[1:Nsect], EAtoCG , fill_value='extrapolate')
        EAtoCG = np.insert(EAtoCG,0,f(Ytip))

        X[4][Nsect-1] = X[3][Nsect-1] - EAtoCG[Nsect-1]
        for i in range(Nsect-2,-1, -1):
            X[4][i] = X[3][i] - EAtoCG[i]
        ax.plot(Y, X[4][:], 'g', label = 'CG')
    except:
        pass        
        
    # Create the aileron
    try:
        cf_c       = Wing.sections.aileron_ratios
        ail_chords = np.zeros(Nsect)
        ail_X      = np.zeros(1)
        ail_Y      = np.zeros(1)

        cf_c       = np.insert(cf_c,0,cf_c[0])
        ail_chords = cf_c * chords
        for i in range(Nsect):
            if ail_chords[i] > 0:
                ail_X = np.append(ail_X,X[1][i] - ail_chords[i])
                ail_Y = np.append(ail_Y,Y[i])
        ail_Y = np.delete(ail_Y,0)
        ail_X = np.delete(ail_X,0)
        ax.plot(ail_Y,ail_X, 'k')
        ax.plot([ail_Y[len(ail_Y)-1],ail_Y[len(ail_Y)-1]],[ail_X[len(ail_Y)-1],X[1][len(ail_Y)-1]], 'k')
        
    except:
        pass

    ax.set_title('Wing planform')
    ax.set_xlabel('Y, in')
    ax.set_ylabel('X, in')
    ax.axis('equal')
    ax.grid()
    ax.legend(loc='upper right')  

    return ax

def plot_3D_wing(Wing,Static,Dynamic):
    """ Generates the wing planform plot with an
        AC and EA lines and the aileron (if defined)
    """
    
    #-----------------
    # Unpack inputs
    #-----------------
    Y        = Wing.sections.Y_centroid
    Ytip     = Wing.chords.Y[len(Wing.chords.Y)-1]
    chordtip = Wing.chords.chords[len(Wing.chords.chords)-1]
    chords   = Wing.sections.chords
    inc      = Wing.sections.incidence
    inctip   = Wing.incidence[len(Wing.incidence)-1]
    LE_sweep = Wing.sections.LEsweep

    # Append the tip parameters
    Y      = np.insert(Y,0,Ytip)
    chords = np.insert(chords,0,chordtip)
    inc    = np.insert(inc,0,inctip)

    Nsect = len(Y)
    
    # Create the wing 
    X = np.zeros((5,Nsect))
    Z = np.zeros((5,Nsect))
    
    X[0][Nsect-1] = 0
    X[1][Nsect-1] = chords[Nsect-1]*math.cos(math.radians(inc[Nsect-1]))
    X[2][Nsect-1] = 0.25 * chords[Nsect-1]*math.cos(math.radians(inc[Nsect-1]))
    Z[0][Nsect-1] = 0.25 * chords[Nsect-1]*math.sin(math.radians(inc[Nsect-1]))
    Z[1][Nsect-1] = -chords[Nsect-1]*math.sin(math.radians(inc[Nsect-1]))+Z[0][Nsect-1]

    for i in range(Nsect-2,-1, -1):
        X[0][i] = X[0][i+1] + (Y[i]-Y[i+1]) * math.tan(math.radians(LE_sweep[i]))
        X[1][i] = X[0][i] + chords[i] * math.cos(math.radians(inc[i]))
        X[2][i] = X[0][i] + 0.25 * chords[i] * math.cos(math.radians(inc[i]))
        Z[0][i] = 0.25 * chords[i]*math.sin(math.radians(inc[i]))
        Z[1][i] = -chords[i]*math.sin(math.radians(inc[i]))+Z[0][i]

    # Draw the wing
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    ax.plot(Y, X[0][:], Z[0][:], 'k')
    ax.plot(Y, X[1][:], Z[1][:], 'k')
    ax.plot([Y[0],Y[0]],[X[0][0],X[1][0]],[Z[0][0],Z[1][0]], 'k')
    ax.plot([Y[Nsect-1],Y[Nsect-1]],[X[0][Nsect-1],X[1][Nsect-1]],[Z[0][Nsect-1],Z[1][Nsect-1]], 'k')
    ax.plot(Y, X[2][:],Z[2][:], 'r', label = 'Aero center')

    for i in range(1,Nsect-1):
       ax.plot([Y[i],Y[i]],[X[0][i],X[1][i]],[Z[0][i],Z[1][i]], 'grey')

    # Draw an elastic axis
    try:
        try:
            EAtoAC = Static.structures.distances.to_AC
        except:
            EAtoAC = Dynamic.structures.distances.to_AC
        EAtoAC = EAtoAC[::-1]
        f      = interpolate.interp1d(Y[1:Nsect], EAtoAC , fill_value='extrapolate')
        EAtoAC = np.insert(EAtoAC,0,f(Ytip))
        
        X[3][Nsect-1] = (0.25 * chords[Nsect-1] + EAtoAC[Nsect-1])*math.cos(math.radians(inc[Nsect-1])) 
        Z[3][Nsect-1] = -(0.25 * chords[Nsect-1]+EAtoAC[Nsect-1]) \
                        *math.sin(math.radians(inc[Nsect-1])) + Z[0][Nsect-1]
        for i in range(Nsect-2,-1, -1):
            X[3][i] = X[0][i] + (0.25 * chords[i]+ EAtoAC[i]) * math.cos(math.radians(inc[i]))
            Z[3][i] = Z[0][i] - (0.25 * chords[i]+ EAtoAC[i]) * math.sin(math.radians(inc[i]))
        ax.plot(Y, X[3][:], Z[3][:], 'b', label = 'Elastic axis')
    except:
        pass

    # Draw a cg axis
    try:
        EAtoCG = Static.structures.distances.to_CG
        EAtoCG = EAtoCG[::-1]
        f      = interpolate.interp1d(Y[1:Nsect], EAtoCG , fill_value='extrapolate')
        EAtoCG = np.insert(EAtoCG,0,f(Ytip))

        X[4][Nsect-1] = (0.25 * chords[Nsect-1] + EAtoAC[Nsect-1] - EAtoCG[Nsect-1])*math.cos(math.radians(inc[Nsect-1])) 
        Z[4][Nsect-1] = -(0.25 * chords[Nsect-1] + EAtoAC[Nsect-1] - EAtoCG[Nsect-1]) \
                        *math.sin(math.radians(inc[Nsect-1])) + Z[0][Nsect-1]
        for i in range(Nsect-2,-1, -1):
            X[4][i] = X[0][i] + (0.25 * chords[i]+ EAtoAC[i] - EAtoCG[i]) * math.cos(math.radians(inc[i]))
            Z[4][i] = Z[0][i] - (0.25 * chords[i]+ EAtoAC[i] - EAtoCG[i]) * math.sin(math.radians(inc[i]))
        ax.plot(Y, X[4][:], Z[4][:], 'g', label = 'CG')
    except:
        pass        
        
    # Create the aileron
    try:
        cf_c       = Wing.sections.aileron_ratios
        ail_chords = np.zeros(Nsect)
        ail_X      = np.zeros(1)
        ail_Y      = np.zeros(1)
        
        cf_c       = np.insert(cf_c,0,cf_c[0])
        ail_chords = cf_c * chords

        for i in range(Nsect):
            if ail_chords[i] > 0:
                ail_X = np.append(ail_X,X[1][i] - ail_chords[i])
                ail_Y = np.append(ail_Y,Y[i])
        ail_Y = np.delete(ail_Y,0)
        ail_X = np.delete(ail_X,0)
        ax.plot(ail_Y,ail_X, 'k')
        ax.plot([ail_Y[len(ail_Y)-1],ail_Y[len(ail_Y)-1]],[ail_X[len(ail_Y)-1],X[1][len(ail_Y)-1]], 'k')
        
    except:
        pass


    #ax.set_title('Wing planform')
    ax.set_xlim3d(0, Ytip)
    ax.set_ylim3d(0.5*(chords[Nsect-1]-Ytip),0.5*(chords[Nsect-1]+Ytip))
    ax.set_zlim3d(-Ytip/2,Ytip/2)
    ax.grid()
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    #plt.axis('off')
    ax.legend()
    plt.savefig('Wing_planform.png') 

    return ax


def plot_Elastic_3D_wing(Wing,Static):
    """ Generates the wing planform plot with an
        AC and EA lines and the aileron (if defined)
    """
    
    #-----------------
    # Unpack inputs
    #-----------------
    Y        = Wing.sections.Y_centroid
    Ytip     = Wing.chords.Y[len(Wing.chords.Y)-1]
    chordtip = Wing.chords.chords[len(Wing.chords.chords)-1]
    chords   = Wing.sections.chords
    inc      = Wing.sections.incidence
    inctip   = Wing.incidence[len(Wing.incidence)-1]
    LE_sweep = Wing.sections.LEsweep       
    alpha_e  = Static.aerodynamics.wing.elast_alpha
    EAtoAC = Static.structures.distances.to_AC
    
    # Append the tip parameters
    Y      = np.insert(Y,0,Ytip)
    chords = np.insert(chords,0,chordtip)
    inc    = np.insert(inc,0,inctip)
    
    Nsect = len(Y)
    EAtoAC = EAtoAC[::-1]
    f      = interpolate.interp1d(Y[1:Nsect], EAtoAC , fill_value='extrapolate')
    EAtoAC = np.insert(EAtoAC,0,f(Ytip))
    
    # Create the wing 
    X = np.zeros((6,Nsect))
    Z = np.zeros((6,Nsect))
    
    X[0][Nsect-1] = 0
    X[1][Nsect-1] = chords[Nsect-1]*math.cos(math.radians(inc[Nsect-1]))
    X[2][Nsect-1] = 0.25 * chords[Nsect-1]*math.cos(math.radians(inc[Nsect-1]))
    X[3][Nsect-1] = (0.25 * chords[Nsect-1] + EAtoAC[Nsect-1])*math.cos(math.radians(inc[Nsect-1])) 
    X[4][Nsect-1] = X[0][Nsect-1]
    X[5][Nsect-1] = X[1][Nsect-1]   
    Z[0][Nsect-1] = 0.25 * chords[Nsect-1]*math.sin(math.radians(inc[Nsect-1]))
    Z[1][Nsect-1] = -chords[Nsect-1]*math.sin(math.radians(inc[Nsect-1]))+Z[0][Nsect-1]
    Z[3][Nsect-1] = -(0.25 * chords[Nsect-1]+EAtoAC[Nsect-1]) \
                    *math.sin(math.radians(inc[Nsect-1])) + Z[0][Nsect-1]
    Z[4][Nsect-1] = Z[0][Nsect-1]
    Z[5][Nsect-1] = Z[1][Nsect-1]
    
    for i in range(Nsect-2,-1, -1):
        X[0][i] = X[0][i+1] + (Y[i]-Y[i+1]) * math.tan(math.radians(LE_sweep[i]))
        X[1][i] = X[0][i] + chords[i] * math.cos(math.radians(inc[i]))
        X[2][i] = X[0][i] + 0.25 * chords[i] * math.cos(math.radians(inc[i]))
        X[3][i] = X[0][i] + (0.25 * chords[i]+ EAtoAC[i]) * math.cos(math.radians(inc[i]))
        X[4][i] = X[0][i]
        X[5][i] = X[1][i]
        Z[0][i] = 0.25 * chords[i]*math.sin(math.radians(inc[i]))
        Z[1][i] = -chords[i]*math.sin(math.radians(inc[i]))+Z[0][i]
        Z[3][i] = Z[0][i] - (0.25 * chords[i] + EAtoAC[i]) * math.sin(math.radians(inc[i]))
        Z[4][i] = (0.25 * chords[i] + EAtoAC[i]) * math.sin(math.radians(alpha_e[i])) + Z[0][i]
        Z[5][i] = -chords[i]*math.sin(math.radians(inc[i]) +  math.radians(alpha_e[i]))+Z[4][i]       


    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Draw the rigid wing
    ax.plot(Y, X[0][:], Z[0][:], 'k')
    ax.plot(Y, X[1][:], Z[1][:], 'k', label = 'Rigid Wing')
    ax.plot(Y, X[2][:], Z[2][:], 'g', label = 'Aero center')
    ax.plot(Y, X[3][:], Z[3][:], 'b', label = 'Elastic axis')
    for i in range(Nsect):
       ax.plot([Y[i],Y[i]],[X[0][i],X[1][i]],[Z[0][i],Z[1][i]], 'k')
       ax.plot([Y[i],Y[i]],[X[4][i],X[5][i]],[Z[4][i],Z[5][i]], 'r')

    # Draw the elastic wing
    ax.plot(Y, X[4][:], Z[4][:], 'r')
    ax.plot(Y, X[5][:], Z[5][:], 'r', label = 'Elastic Wing')
    
    #ax.set_title('Elastic twist for a static maneuver')
    ax.set_xlim3d(0, Ytip)
    ax.set_ylim3d(0.5*(chords[Nsect-1]-Ytip),0.5*(chords[Nsect-1]+Ytip))
    ax.set_zlim3d(-Ytip/10,Ytip/10)
    ax.grid()
    ax.legend()
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    #plt.axis('off')
    plt.savefig('Wing_planform_Elastic.png') 

    return ax

def plot_lift_distribution(alpha,cl_r,Wing,Static,Analysis):
    """ Generates rigid and elastic lift distribution plots

    """
    
    #-----------------
    # Unpack inputs
    #-----------------
    Y    = Wing.sections.Y_centroid
    Ytip = Wing.chords.Y[len(Wing.chords.Y)-1]

    # Plot lift distribution curves    
    # Append the tip parameters
    Y    = np.insert(Y,0,Ytip)

    fig, ax3 = plt.subplots(2,1)
    if Analysis == 2.2:
        cl    = Static.aerodynamics.wing.cl     
        cl_el = Static.aerodynamics.wing.cl_el

        cl_el = np.insert(cl_el,0,0)
        cl_r  = np.insert(cl_r,0,0)
        cl    = np.insert(cl,0,0)

        ax3[0].plot(Y,cl_r, 'g', label = 'rigid: AoA = ' + str(round(alpha[0],2)) + 'deg')
        ax3[0].plot(Y,cl_el, 'b', label = 'elastic: AoA = ' + str(round(alpha[0],2)) + 'deg')        
        ax3[0].plot(Y,cl, 'r', label = 'total: AoA = ' + str(round(alpha[0],2)) + 'deg')

        ax3[0].legend()
        ax3[0].grid()
        ax3 = plot_2D_wing(Wing,Static,ax3[1])
        
        plt.tight_layout()
        plt.savefig('Elastic_lift_distribution.png')
    
    else:
        N = np.shape(cl_r)
            
        cl_dum   = np.zeros(N[1])
        for i in range(N[0]):
            cl_dum = cl_r[i][:]
            cl_dum = np.insert(cl_dum,0,0)
            ax3[0].plot(Y,cl_dum, c = np.random.rand(3,), label = 'rigid: AoA = ' + str(round(alpha[i],2)) + 'deg')
        ax3[0].set_xlabel('Y, in')
        ax3[0].set_ylabel('cl')
        ax3[0].set_title('Lift distribution')
        ax3[0].legend()
        ax3[0].grid()
        ax3 = plot_2D_wing(Wing,Static,ax3[1])
        
        plt.tight_layout()
        plt.savefig('Rigid_Lift_distribution.png')


    return


def plot_flutter_data(Wing,Dynamic):
    """ Generates flutter output plots

    """
    #-----------------
    # Unpack inputs
    #-----------------
    Y       = Wing.sections.Y_centroid
    M       = Dynamic.Mach
    g       = Dynamic.structural_damping
    phi_h   = Dynamic.bending_modes      
    phi_a   = Dynamic.torsion_modes      
    omega_h = Dynamic.bending_frequency  
    omega_a = Dynamic.torsion_frequency

    N = len(Y)
    
    phih = np.zeros(N)
    phia = np.zeros(N)

    for i in range(N):
        phih[i] = phi_h[i][0]
        phia[i] = phi_a[i][0]
        
    # plot a structural damping diagram
    fig, ax4 = plt.subplots()
    ax4.plot(M[0][:],g[0][:], label = 'root 1')
    ax4.plot(M[1][:],g[1][:], label = 'root 2')
    ax4.set_xlabel('Mach number')
    ax4.set_ylabel('g')
    ax4.set_title('Artificial structural damping diagram')
    ax4.grid()
    ax4.legend()
    plt.savefig('Structural_damping.png')

    # plot natural modes of the wing
    fig, ax5 = plt.subplots()
    ax5.plot(Y,phih, label = 'Bending freq-cy = ' + str(round(omega_h[0],2)) + 'rad/s')
    ax5.plot(Y,phia, label = 'Torsioan freq-cy = ' + str(round(omega_a[0],2)) + 'rad/s')
    ax5.set_xlabel('Y, in')
    ax5.set_ylabel('$\\phi$')
    ax5.set_title('Natural modes and frequencies')    
    ax5.grid()
    ax5.legend()
    plt.savefig('Modes_frequencies')
    
    return
                    

