## PyAeroElast

# Run_PyAeroElast.py
# Author  : Stanislav Karpuk
#           (stankarpuk93@gmail.com)
# Created : 03/26/2019
# Modified:

import sys
from PyAeroElast import PyAeroElast

"""
    Description:
        PyAeroElast is a Python set of scripts to estimate static and aeroelastic behavior
        of an airplane conceiver-type wing of medium to high aspect ratios. The code uses aeroelastic
        strip theory and semi-empirical relations to estimate aerodynamic behavior and aeroelastic
        responses of the wing.
        
        The primary use of PyAeroElast is to perform rapid estimations of aeroelastic behavior of the aircraft
        and exploration of the V-n flight envelope to identify potential failure cases.
        Static aeroelastic analysis includes rigid body aerodynamic coefficients and lift distribution,
        divergence speed, lift distribution for symmetric maneuvers due to elastic effects,
        aileron power coefficients and aileron reversal speed, and estimation of flutter speed.    

"""


#ANALYSIS CASE INPUT
#-------------------------------------------------------------------------------------------------------------------------

input_directory  =  r'C:\Users\Станислав\Documents\PyAeroElast\V1.0\cases\Dynamic_aeroelasticity'   # Input file directory
output_directory =  r'C:\Users\Станислав\Documents\PyAeroElast\V1.0\cases\Dynamic_aeroelasticity'   # Output file directory

sys.path.insert(0, input_directory)

import Flutter_speed_Strip_theory as inp                                                             # Input file name


#Execute PyAeroElast
#-------------------------------------------------------------------------------------------------------------------------
PyAeroElast(inp, output_directory)

