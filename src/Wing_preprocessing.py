## PyAeroElast

# Wing_preprocessing.py
# Author  : Stanislav Karpuk
#           (stankarpuk93@gmail.com)
# Created : 01/19/2019
# Modified:

import numpy as np
import math
from Data_structure  import Data

class Wing_preprocessing:
  
#-----------------------------------------------------------------------------------------------

    # Compute the Wing geometric properties
    #----------------------------------------------

    def compute_geometry(self, General, Static, Aerodynamics):
        """ Computes span, area, aspect ratio, MAC """
        # Unpack
        chords   = self.chords.chords
        Y        = self.chords.Y
        Ypoints  = self.sections.Y
        sweep_LE = self.sweep.leading_edge
        LE_sweep = self.sweep.leading_edge
        EA_sweep = self.sweep.elastic_axis
        inc      = self.incidence
        a0       = Aerodynamics.a0

        # Calculate equivalent wing properties 
        sumr     = 0
        sumt     = 0
        sumsweep = 0
        area     = 0
        span     = 2 * Y[len(Y)-1]

        for i in range(len(chords)-1):
            area     += 0.5 * (chords[i] + chords[i+1]) * (Y[i+1] - Y[i])
            sumr     += chords[i] * 0.5 * (chords[i] + chords[i+1]) * (Y[i+1] - Y[i])
            sumt     += chords[i+1] * 0.5 * (chords[i] + chords[i+1]) * (Y[i+1] - Y[i])
            sumsweep += sweep_LE[i] * 0.5 * (chords[i] + chords[i+1]) * (Y[i+1] - Y[i])
            
        Sw     = span * (sumr + sumt) / (2 * area)
        Cre    = 2 / Sw * sumr
        Cte    = 2 / Sw * sumt
        taper  = Cte / Cre
        sweep  = sumsweep / area
        sweepq = math.atan(math.tan(math.radians(sweep))+0.25/(0.5*span)*(Cte-Cre))

        # calculate wing local quarter-chord sweep
        sweep_C4 = np.zeros(len(Y)-1)
        sweep_C2 = np.zeros(len(Y)-1)
        for i in range(len(sweep_C4)):
            sweep_C4[i] = math.degrees(math.atan(math.tan(math.radians(sweep_LE[i])) + \
                                                 0.25*(chords[i+1]-chords[i])/(Y[i+1] - Y[i])))
            sweep_C2[i] = math.degrees(math.atan(math.tan(math.radians(sweep_LE[i])) + \
                                                  0.5*(chords[i+1]-chords[i])/(Y[i+1] - Y[i])))          

        npoints  = len(Ypoints)

        Ysect        = np.zeros(np.sum(npoints-1))
        EAsweep_sect = np.zeros(np.sum(npoints-1))
        LEsweep_sect = np.zeros(np.sum(npoints-1))
        phi          = np.zeros(np.sum(npoints-1))
        washout      = np.zeros(np.sum(npoints-1))
        chords_sect  = np.zeros(np.sum(npoints-1))
        wing_inc     = np.zeros(np.sum(npoints-1))
        C4_sweep     = np.zeros(np.sum(npoints-1))
        C2_sweep     = np.zeros(np.sum(npoints-1))
        
        # Calculate wing stations and corresponding angles for structural analysis
        phi[0]   = 0.5 * math.pi
        Ysect[0] = 0

        for i in range(1,len(Ysect-1)):
            Ysect[i] = 0.5 * (Ypoints[i+1] + Ypoints[i])
            phi[i]   = math.acos(2 * np.around(Ysect[i], decimals = 4) / span)           

        # Calculate chord length at each spanwise station
        chords_sect  = np.interp(Ysect, Y, chords)
        a0_sect      = np.interp(Ysect, Y, a0)
        wing_inc     = np.interp(Ysect, Y, inc)
        Ysect        = Ysect[::-1]       
        chords_sect  = chords_sect[::-1]      
        wing_inc     = wing_inc[::-1]
        phi          = phi[::-1]
        a0_sect      = a0_sect[::-1]
      
        # Generate a sweep matrix at each section
        for i in range(len(EAsweep_sect)):
            for j in range(len(Y)-1):
              if Ysect[i] >= Y[j] and Ysect[i] <= Y[j+1]:
                EAsweep_sect[i] = EA_sweep[j]
                LEsweep_sect[i] = LE_sweep[j]
                C4_sweep[i]     = sweep_C4[j]
                C2_sweep[i]     = sweep_C2[j]
        
        # Pack
        self.equivalent         = Data()
        self.equivalent.chords  = Data()
        self.equivalent.sweeps  = Data()
        Aerodynamics.sections   = Data()

        self.equivalent.chords.root          = Cre
        self.equivalent.chords.tip           = Cte
        self.span                            = span
        self.equivalent.taper                = taper
        self.equivalent.sweeps.leading_edge  = sweep
        self.equivalent.sweeps.quarter_chord = sweepq
        self.sweep.quarter_chord             = sweep_C4
        self.sweep.half_chord                = sweep_C2
        self.washout                         = washout
        Static.aerodynamics.wing.chords      = chords_sect
        Static.aerodynamics.wing.Y           = Ysect
        self.sections.Y_centroid             = Ysect
        self.sections.incidence              = wing_inc
        self.sections.chords                 = chords_sect
        self.sections.LEsweep                = LEsweep_sect
        Static.aerodynamics.wing.phi         = phi
        General.structures.EAsweep           = EAsweep_sect
        General.structures.C4sweep           = C4_sweep
        General.structures.C2sweep           = C2_sweep
        Aerodynamics.sections.alpha0         = a0_sect
        
        self.area         = 2 * area
        self.aspect_ratio = self.span**2./self.area
        self = Wing_preprocessing.compute_mac(self)
        
        return self, General, Static, Aerodynamics
        

    def compute_mac(self):
        """ Computes mean aerodynamic chord """

        Cr    = self.equivalent.chords.root
        taper = self.equivalent.taper

        self.chords.cmac = 2/3 * Cr * (1 + taper + taper**2)/(1+taper)
        
        return self

    def control_surfaces_calc(self,Static,Aerodynamics):
        cf_data  = self.chords.aileron_ratios
        Cla_data = Aerodynamics.Cla
        Y_data   = self.chords.Y
        Y        = self.sections.Y_centroid 

        N = len(Y)
        cf_c = np.zeros(N)
        Cla  = np.zeros(N)

        for i in range(N):
            for j in range(1,len(Y_data)):
                if Y[i] >= Y_data[j-1] and Y[i] <= Y_data[j]:
                    cf_c[i] = min(cf_data[j-1],cf_data[j])

        #--------------
        # Pack output
        #--------------      
        self.sections.aileron_ratios = cf_c

        return self, cf_c

#--------------------------------------------------------------------------------------------------

def Standard_atmosphere(aerodynamics):

    h = aerodynamics.altitude
    
    aerodynamics.pressure    = 2116 * (1-0.0000068756 * h) ** 5.2561
    aerodynamics.temperature = 518.67 * (1-0.0000068756 * h)
    aerodynamics.density     = aerodynamics.pressure / (1716 * aerodynamics.temperature)

    return aerodynamics
