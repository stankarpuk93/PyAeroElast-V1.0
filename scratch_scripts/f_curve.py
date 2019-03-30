from scipy import interpolate
import matplotlib as plt
import numpy as np

x = [0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,\
    0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,\
    0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1]

y = [-45,-30,0,30,45,60]

xx,yy = np.meshgrid(x,y)

f = open('new.dat','r')
with open('new.dat','r') as textFile:
    z = [line.split() for line in textFile]

z= np.asfarray(z,float)
print(np.shape(z))
x1 = [5]
#y1 = [0.0, 0.11, 0.24, 0.31, 0.44, 0.5, 0.67, 0.7, 0.8, 0.9]
y1 = [0,0.3,0.6,0.9]
z1 = interpolate.RectBivariateSpline(x,y,z)

test = z1(y1,x1)
print(test)
