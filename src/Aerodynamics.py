import math
from math import sqrt, pi, cos, sin, tan
from Data_structure  import Data
from scipy import interpolate
import numpy as np

eps = 1E-10

def strip_theory(Wing, Aerodynamics, General, Static, analysis_type):
  """ Computes lift and aerodynamic influence coefficients based on
      strip theory

  References:
  Bisplinghoff, Ashley "Principles of Aeroelasticity"
  Torenbeek "Synthesis of Subsonic Aircraft Design"

  """

  #---------------------- 
  # Unpack general input
  #----------------------
  Y        = Wing.sections.Y_centroid 
  chords   = Wing.sections.chords  
  Cla      = Aerodynamics.Cla
  M        = General.aerodynamics.Mach
  AR       = Wing.aspect_ratio
  span     = Wing.span
  sweep_eq = Wing.equivalent.sweeps.quarter_chord

  eta      = np.zeros(len(Y))

 
  # Calculate the Wing lift-curve slope ( Diederich app-n, up to M = 0.7, AR = 1.5 )
  # CLa limiter to avoid transonic/supersonic Mach numbers
  if M > 0.7:
    M = 0.7
    
  AA  = Cla * math.cos(math.radians(sweep_eq)) / (1 - M**2)**0.5
  BB  = AR * (1 - M**2)**0.5
  CC  = BB * (1+(Cla * math.cos(math.radians(sweep_eq)) / (math.pi * AR * (1 - M**2)**0.5))**2)**0.5
  DD  = Cla * math.cos(math.radians(sweep_eq)) / math.pi
  CLa = AA * BB / (CC + DD)

  # Generate a matrix of aerodynamic influence coefficients
  A = np.zeros((len(chords),len(chords)))
  A = 1/CLa * np.diag(np.reciprocal(chords))

  # Calculate lift, induced drag and rigid body lift distribution for a given AoA
  #-------------------------------------------------------------------------------
  if analysis_type == 1:
    #--------------- 
    # Unpack input
    #---------------
    alpha = Aerodynamics.alpha
    
    CL  = np.zeros(len(alpha))
    CDi = np.zeros(len(alpha))

    # Normalize spanwise coordinates
    eta = Y / (0.5 * span)
    
    Aerodynamics, Static = Diederich_lift(eta,alpha,CLa,General,Wing,Aerodynamics,Static,analysis_type)

    #------------- 
    # Pack output
    #-------------
    Aerodynamics.alpha = alpha

  # rigid body lift distribution for a given load factor
  #-------------------------------------------------------------------------------
  elif analysis_type == 2.2:
    #--------------- 
    # Unpack input
    #---------------
    q          = General.aerodynamics.dynamic_pressure
    N          = Static.aerodynamics.load_factor
    W          = Wing.weight
    area       = Wing.area
    incedence  = Wing.sections.incidence
    cmac       = Wing.chords.cmac
    sweep      = General.structures.C4sweep 
    alpha0     = Aerodynamics.sections.alpha0

    Ns = len(Y)
    
    inc0  = np.zeros(len(eta))
    alpha = np.zeros(1)
    J     = np.zeros(Ns)
    k1    = np.zeros(Ns)
    C1    = np.zeros(Ns)
    C2    = np.zeros(Ns)
    C3    = np.zeros(Ns)
    La    = np.zeros(Ns)

    # Normalize spanwise coordinates
    for i in range(len(Y)):
      eta[i] = Y[i] / (0.5 * span)
    
    # Calculate free-stream AoA
    for i in range(len(eta)):
      J[i]    = f(sweep[i],eta[i])
      F       = AR/(Cla / (2*math.pi)) * math.cos(math.radians(sweep[i]))
      Cx      = (2*math.pi*AR)/(Cla*math.cos(math.radians(sweep[i]))) 
      C1[i]   = -1.1002*10**(-3)*Cx**2 + 5.8877*10**(-2)*Cx - 8.4084*10**(-3)
      C2[i]   =  2.7595*10**(-3)*Cx**2 - 1.1089*10**(-1)*Cx + 1.0075
      C3[i]   = -1.7613*10**(-3)*Cx**2 + 5.4777*10**(-2)*Cx - 9.1359*10**(-3)
      k1[i]   = (F*(1+4/F**2)**0.5+2)/(F*(1+36/F**2)**0.5+6)    
      La[i]   = C1[i]*chords[i]/cmac + C2[i]*4/math.pi*(1-eta[i]**2)**0.5+C3[i]*J[i]
      inc0[i] = (math.radians(incedence[i])-math.radians(alpha0[i]))*La[i]

    alpha[0] = math.degrees((N*W/(q*area*CLa/144) + np.trapz(inc0, x=eta))/(-np.trapz(La, x=eta)))
    
    Aerodynamics, Static = Diederich_lift(eta,alpha,CLa,General,Wing,Aerodynamics,Static,analysis_type)
    
    #------------- 
    # Pack output
    #-------------
    Aerodynamics.maneuver_alpha = alpha
  

  #------------- 
  # Pack output
  #-------------
  Aerodynamics.CLa           = CLa
  Static.aerodynamics.wing.A = A
  
  return Aerodynamics, Static
    

def f(sweep,eta):
  """ Returns the lift distribution function for
      a particular wing station
  """
  # Input axes for interpolation
  x = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,\
      0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,\
      0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1]

  y = [-45,-30,0,30,45,60]

  xx,yy = np.meshgrid(x,y)

  # read the interpolation data from a text file
  with open('f.dat','r') as textFile:
      z = [line.split() for line in textFile]

  z  = np.asfarray(z,float)
  z1 = interpolate.RectBivariateSpline(x,y,z)
  f  = z1(eta,sweep)
  textFile.close()

  return f
    
def Diederich_lift(eta,alpha,CLa,General,Wing,Aerodynamics,Static,analysis):
  """ Calculates lft distribution and average AoA
      using Diederich's method
  """
  
  #--------------- 
  # Unpack input
  #---------------
  M          = General.aerodynamics.Mach
  sweep      = General.structures.C4sweep
  Y          = Wing.sections.Y_centroid 
  AR         = Wing.aspect_ratio
  cmac       = Wing.chords.cmac
  chords     = Wing.sections.chords  
  incedence  = Wing.sections.incidence
  Cla        = Aerodynamics.Cla
  alpha0     = Aerodynamics.sections.alpha0

  N0 = len(Y)
  N1 = len(alpha)
  N2 = len(eta)
  
  alpha01  = np.zeros(N0)
  ccl_cmac = np.zeros((N1,N0))
  ccl      = np.zeros((N1,N0))
  cl       = np.zeros((N1,N0))
  J        = np.zeros(N0)
  k1       = np.zeros(N0)
  C1       = np.zeros(N0)
  C2       = np.zeros(N0)
  C3       = np.zeros(N0)
  La       = np.zeros(N0)
  Lb       = np.zeros(N0)
  CL       = np.zeros(N1)
  CDi      = np.zeros(N1)

  for i in range(N2):
    J[i]  = f(sweep[i],eta[i])
    F     = AR/(Cla / (2*math.pi)) * math.cos(math.radians(sweep[i]))
    k1[i] = (F*(1+4/F**2)**0.5+2)/(F*(1+36/F**2)**0.5+6)
    
  for j in range(N1):
    for i in range(N2):
      Cx    = (2*math.pi*AR)/(Cla*math.cos(math.radians(sweep[i]))) 
      C1[i] = -1.1002*10**(-3)*Cx**2 + 5.8877*10**(-2)*Cx - 8.4084*10**(-3)
      C2[i] =  2.7595*10**(-3)*Cx**2 - 1.1089*10**(-1)*Cx + 1.0075
      C3[i] = -1.7613*10**(-3)*Cx**2 + 5.4777*10**(-2)*Cx - 9.1359*10**(-3)
      
    for i in range(N2):
      La[i]      = C1[i]*chords[i]/cmac + C2[i]*4/math.pi*(1-eta[i]**2)**0.5+C3[i]*J[i]
      alpha01[i] = (math.radians(alpha[j])+math.radians(incedence[i])-math.radians(alpha0[i]))*La[i]

    alpha01tot = -np.trapz(alpha01, x=eta)
    for i in range(N2):
      Lb[i]          = (1-M**2)**0.5 * k1[i] * CLa * (math.radians(alpha[j])+math.radians(incedence[i])-math.radians(alpha0[i]) - alpha01tot) * La[i]
      ccl_cmac[j][i] = CLa * alpha01tot * La[i] + Lb[i]
      ccl[j][i]      = ccl_cmac[j][i] * cmac
      cl[j][i]       = ccl[j][i] / chords[i]

    # Calculate Lift and Induced Drag coef-s
    CL[j]  = CLa * alpha01tot
    CDi[j] = CL[j]**2 /(math.pi * AR)

  #------------- 
  # Pack output
  #-------------
  if analysis == 1:
    Aerodynamics.wing     = Data()
    Aerodynamics.wing.CL  = CL
    Aerodynamics.wing.CDi = CDi
    Aerodynamics.ccl_cmac = ccl_cmac
    Aerodynamics.ccl      = ccl
    Aerodynamics.cl       = cl
  else:
    Static.aerodynamics.wing.cl       = cl
    Static.aerodynamics.wing.ccl      = ccl
    Static.aerodynamics.wing.ccl_cmac = ccl_cmac

  return Aerodynamics, Static
    

  

  
