## PyAeroElast Input

# Author  : Stanislav Karpuk
#           (stankarpuk93@gmail.com)
# Created : 01/17/2019
# Modified:

import os
import math
import numpy as np
from Data_structure  import Data
from Wing_preprocessing import Wing_preprocessing
from Wing_preprocessing import Standard_atmosphere

def PyAeroElast_input():
    
    """
        Description:
            Input-file for an aerodynamic analysis

        Analysis options:       
            1. Aerodynamic analysis 
            2. Static Aeroelasticity
                2.1. Divergence speed 
                2.2. Modified lift distribution for symmetric maneuvers 
                2.3. Aeroelastic control effectiveness
            3. Dynamic Aeroelasticity
    """
    
    # GENERAL INPUT BLOCK
    #-------------------------------------------------------------------

    # Define analysis types
    Analysis = [1]
    
    # GENERAL INPUT BLOCK
    #-------------------------------------------------------------------
    # Define a Wing geometry
    # Requirements: rectangular, const sweep, tapered, linear washout
    Wing                    = Data()
    Wing.chords             = Data()
    Wing.sections           = Data()
    Wing.sweep              = Data()
    Wing.sections           = Data()
    
    Wing.chords.Y           = [0, 300, 500]
    Wing.chords.chords      = [225, 150, 100]
    Wing.sweep.leading_edge = [0, 10]
    Wing.sweep.elastic_axis = [0, 10]
    Wing.incidence          = [3, 3, 3]
    Wing.weight             = 83838
    Wing.sections.Y         = np.linspace(Wing.chords.Y[0], Wing.chords.Y[len(Wing.chords.Y)-1], num=25)

    General                      = Data()
    # Define General Structural Properties
    #-------------------------------------------------------------------
    General.structures           = Data()
    General.structures.EIGJ      = Data()
    General.structures.distances = Data()
   

    # Define General Aerodynamic Properties
    #-------------------------------------------------------------------    
    General.aerodynamics             = Data()
    General.aerodynamics.wing        = Data()
    General.aerodynamics.wing.MAC    = Data()

    General.aerodynamics.altitude    = 27000
    General.aerodynamics.Mach        = 0.7

    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------

    # AERODYNAMIC ANALYSIS INPUT
    #-------------------------------------------------------------------
    Aerodynamics         = Data()
    Aerodynamics.Cla     = 2*math.pi
    Aerodynamics.a0      = [0, 0, 0]
    Aerodynamics.alpha   = [-2, 0, 2, 4, 6]


    # STATIC AEROELASTICITY ANALYSIS INPUT
    #-------------------------------------------------------------------
    Static                       = Data()
    Static.structures            = Data()
    Static.structures.distances  = Data()
    Static.aerodynamics          = Data()
    Static.aerodynamics.wing     = Data()
    Static.aerodynamics.wing.MAC = Data()

    # Structural Data


    # Aerodynamic Data


    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------

    # DYNAMIC AEROELASTICITY ANALYSIS INPUT
    #-------------------------------------------------------------------
    Dynamic                      = Data()
    Dynamic.structures           = Data()
    Dynamic.flutter              = Data()
    Dynamic.flutter.distances    = Data()
    
    # Structural Data (required to calculate the modes)


    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    Wing, General, Static, Aerodynamics = Wing_preprocessing.compute_geometry(Wing, General, Static, Aerodynamics)
    General.aerodynamics = Standard_atmosphere(General.aerodynamics)
    
    General.aerodynamics.airspeed         = General.aerodynamics.Mach * (1.4 * 1716 * General.aerodynamics.temperature)**0.5
    General.aerodynamics.dynamic_pressure = 0.5 * General.aerodynamics.density * General.aerodynamics.airspeed**2

        
    return Analysis, Wing , General, Aerodynamics, Static, Dynamic
