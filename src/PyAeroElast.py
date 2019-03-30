## PyAeroElast

# PyAeroElast.py
# Author  : Stanislav Karpuk
#           (stankarpuk93@gmail.com)
# Created : 01/17/2019
# Modified:

import sys
import Aerodynamics as weiss
from Data_structure  import Data
import static_aeroelasticity as stat_ae
from dynamic_aeroelasticity import compute_flutter_speed
import Post_processing as postproc


def PyAeroElast(inp, output_directory):
    """
        Description:
            Executes the PyAeroElast script.
                
            Analysis options:       
                1. Aerodynamic analysis 
                2. Static Aeroelasticity
                    2.1. Divergence speed 
                    2.2. Modified lift distribution for symmetric maneuvers 
                    2.3. Aeroelastic control effectiveness
                3. Dynamic Aeroelasticity
        
        Source:    
                                     
        Assumptions:

    """

    #--------------------------- 
    # Set the analysis up
    #---------------------------------------------------------------------------------------

    Analysis, Wing, General, Aerodynamics, Static, Dynamic = inp.PyAeroElast_input()

    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------


    # RUN THE PyAeroElast
    #--------------------------------------------------------------------------------------

    for i in range(len(Analysis)):
        print('Running case ' + str(i) + ' of ' + str(len(Analysis)) + '\n')

        # Aerodynamics analysis
        if round(Analysis[i]) == 1:       
            print('Running Aerodynamic Analysis...')
            Aerodynamics, Static = weiss.strip_theory(Wing, Aerodynamics, General, Static, Analysis[i])
            print('Aerodynamic Analysis Completed\n')

        # Static Aeroelastic Analysis
        elif round(Analysis[i]) == 2:       
            print('Running Static Aeroelasticity Analysis...')

            # Divergence Speed Estimation
            if Analysis[i] == 2.1:           
                print('Running Divergence Speed Analysis Using Strip Theory...')
                Static  = stat_ae.compute_divergence_speed(General,Aerodynamics,Static,Wing,Analysis[i])
                print('Divergence Speed Analysis Completed\n')

            # Modified Lift Distribution 
            elif Analysis[i] == 2.2:         
                print('Running Modified Lift Distribution Analysis for a symmetric maneuver...')            
                Static = stat_ae.modified_aerodynamics(General,Static,Aerodynamics,Wing,Analysis[i])           
                print('Modified Lift Distribution Analysis Completed\n')

            # Control Reversal Estimation
            elif Analysis[i] == 2.3:         
                print('Running Aeroelastic Control Effectiveness analysis ...')
                Static = stat_ae.elastic_control(General, Aerodynamics, Static, Wing, Analysis[i])           
                print('Control Reversal Analysis Completed\n')

            else:
                print(str(round(Analysis[i],2)) + ' - there is no such index\n')
            print('Static Aeroelasticity Analysis Completed\n')

        # Dynamic Aeroelastic analysis      
        elif round(Analysis[i]) == 3:        
            print('Running Dynamic Aeroelasticity Analysis...')
            Dynamic = compute_flutter_speed(General,Dynamic,Wing)
           
            print('Dynamic Aeroelasticity Analysis Completed\n')
        else:
            print(str(round(Analysis[i],2)) + ' - there is no such index\n')

    print('Calculation Completed')
    print('Outputting the results...')
    postproc.output(output_directory, General, Wing, Aerodynamics, Static, Dynamic, Analysis)

    print('Anslysis is Completed')

    return


