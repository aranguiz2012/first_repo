import numpy
from matplotlib import pyplot

L = 11
v_max = 136
rho_max = 250

nx = 51
dt = 0.001
dx = L/(nx-1)
nt = 50

x = numpy.linspace(0, L, nx)
rho0 = numpy.ones(nx)*20
rho0[10:20] = 50

inic = rho0.copy()

for i in range(nt):
    rhon = rho0.copy()
    rho0[1:] = rhon[1:]-dt/dx*(v_max-2*v_max/rho_max*rhon[1:])*\
    (rhon[1:]-rhon[:-1])
    

pyplot.plot(x,inic,color='r',lw='2')
pyplot.plot(x,rho0,color='b',lw='2')
pyplot.ylim(0,51)

v1 = v_max*(1-50/rho_max)
print(v1*1000/3600)

rho_max6 = numpy.amax(rho0)
v = v_max*(1-rho_max6/rho_max)
print(v*1000/3600)

rho_avg3 = numpy.average(rho0)
v2 = v_max*(1-rho_avg3/rho_max)
print(v2*1000/3600)