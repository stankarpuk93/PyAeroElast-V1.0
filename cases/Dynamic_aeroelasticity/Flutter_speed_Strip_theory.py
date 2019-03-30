## PyAeroElast Input

# Author  : Stanislav Karpuk
#           (stankarpuk93@gmail.com)
# Created : 01/17/2019
# Modified:

import os
import math
from Data_structure  import Data
from Wing_preprocessing import Wing_preprocessing
from Wing_preprocessing import Standard_atmosphere

def PyAeroElast_input():
    
    """
        Description:
            Input-file for a complete aeroelastic analysis

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
    Analysis = [3]
    
    # Define a Wing geometry
    # Requirements: rectangular, const sweep, tapered, linear washout
    Wing                    = Data()
    Wing.chords             = Data()
    Wing.sweep              = Data()
    Wing.sections           = Data()
    
    Wing.chords.Y           = [0, 500]
    Wing.chords.chords      = [225, 100]
    Wing.sweep.leading_edge = [3.5]
    Wing.sweep.elastic_axis = [0]
    Wing.incidence          = [0, 0]
    Wing.weight             = 83883
    Wing.sections.Y         = [0, 24, 156, 216, 320, 416, 500]
    
    General                 = Data()
    # Define General Structural Properties
    #-------------------------------------------------------------------
    General.structures           = Data()
    General.structures.EIGJ      = Data()
    General.structures.distances = Data()
   
    General.structures.EIGJ.Y   = [0, 10.349, 16.34, 26.144, 40.85, 57.19, 86.057, 122.004,
                                    153.05, 177.56, 200.98, 233.66, 269.063, 301.198, 337.89,
                                    369.826, 403.05, 430.828, 458.606, 480.937, 498.911]

    General.structures.EIGJ.EI  = [67808000000,64767178168,63053995526,60323842995,56396259063,52262212741,
                                    45527646535,38102081698,32488752933,28544699070,25152234279,20989985865,
                                    17168414055,14255569702,11500256067,9533248488,7848860744,6678219219,
                                    5684280119,4988096877,4480033146]

    General.structures.EIGJ.GJ  = [27655600000,27105888700,26869042737,26593086237,26395500801,26404714642,
                                    26757585993,27296384279,27446643992,27184981772,26546926102,24971846293,
                                    22392629759,19390507167,15460587571,11911351861,8466414552,6085501603,
                                    4449146260,3868138277,3989490092]

    General.structures.EIGJ.GK  = [26077000,26100302,26108410,26113169,26100495,26058531,25912669,25602977,
                                    25221192,24844941,24423725,23735199,22856824,21940242,20754878,19602697,
                                    18285061,17090272,15810663,14720470,13803168]


    # Define General Aerodynamic Properties
    #-------------------------------------------------------------------    
    General.aerodynamics             = Data()
    General.aerodynamics.wing        = Data()
    General.aerodynamics.wing.MAC    = Data()

    General.aerodynamics.altitude    = 0
    General.aerodynamics.Mach        = 0.1

    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------

    # AERODYNAMIC ANALYSIS INPUT
    #-------------------------------------------------------------------
    Aerodynamics       = Data()
    Aerodynamics.Cla   = 6.42
    Aerodynamics.a0    = [0, 0]

    # STATIC AEROELASTICITY ANALYSIS INPUT
    #-------------------------------------------------------------------
    Static                       = Data()
    Static.structures            = Data()
    Static.structures.distances  = Data()
    Static.aerodynamics          = Data()
    Static.aerodynamics.wing     = Data()
    Static.aerodynamics.wing.MAC = Data()

    # Structural Data
  

    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------

    # DYNAMIC AEROELASTICITY ANALYSIS INPUT
    #-------------------------------------------------------------------
    Dynamic                      = Data()
    Dynamic.structures           = Data()
    Dynamic.structures.distances = Data()

    Dynamic.structures.mg              = [17400, 6039, 10200, 4200, 3400, 680]
    Dynamic.structures.distances.to_AC = [22.5, 20.25, 17.85, 15.8, 13.3, 11.05]    # positive if ahead of the EA
    Dynamic.structures.distances.to_CG = [0, 7, 2, -2, -2, -4]                      # positive if ahead of the EA
    Dynamic.structures.Iyy             = [8.723*10**6, 8.723*10**6, 5.16*10**6, 3.725*10**6, 2.775*10**6, 0.4*10**6]

    
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    Wing, General, Static, Aerodynamics = Wing_preprocessing.compute_geometry(Wing, General, Static, Aerodynamics)
    General.aerodynamics                = Standard_atmosphere(General.aerodynamics)
    
    General.aerodynamics.airspeed         = General.aerodynamics.Mach * (1.4 * 1716 * General.aerodynamics.temperature)**0.5
    General.aerodynamics.dynamic_pressure = 0.5 * General.aerodynamics.density * General.aerodynamics.airspeed**2
    
    return Analysis, Wing, General, Aerodynamics, Static, Dynamic
