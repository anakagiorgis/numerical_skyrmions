from numpy import *
from matplotlib.pyplot import *

m3 = loadtxt("Landau.txt")

m3_i = loadtxt("Initial.txt")

rho = linspace(0, 10, 101)

plot(rho, m3[:,0], 'r')
plot(rho, m3_i[:,0], 'b')
