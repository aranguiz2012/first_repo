import numpy
from matplotlib import pyplot
numpy.set_printoptions(suppress = True)

nx = 81
gamma = 1.4
icl = numpy.array([1, 0, 100*1000])
icr = numpy.array([0.125, 0, 10*1000])
nt = 51
dt = 0.0002
dx = 0.25

def initial_condition(nx, icl, icr, gamma):
    
    rho = icr[0]*numpy.ones(nx)
    rho[:(nx-1)/2] = icl[0]
    vel = icr[1]*numpy.ones(nx)
    vel[:(nx-1)/2] = icl[1]
    p = icr[2]*numpy.ones(nx)
    p[:(nx-1)/2] = icl[2]
    e = p/((gamma-1.)*rho)
    et = e + vel**2/2
    u = numpy.vstack((rho,rho*vel,rho*et))
    
    return u
    
u = initial_condition(nx, icl, icr, gamma)

def computeF(gamma, u):
    a = u[1]
    b = u[1]**2/u[0]+(gamma-1)*(u[2]-0.5*u[1]**2/u[0])
    c = (u[2] + (gamma-1)*(u[2]-0.5*u[1]**2/u[0]))*u[1]/u[0]
    f = numpy.vstack((a,b,c))

    return f    

def richtmyer(u, nt, dt, dx):
    u_n = numpy.zeros((nt,u.shape[0],u.shape[1]))
    u_med = numpy.zeros((u.shape[0],u.shape[1]-1))
    u_n[:,:,:] = u.copy()
    #u_med = u.copy()
    for t in range(1,nt):
        F=computeF(gamma,u)
        u_med[:,:]=0.5*(u[:,1:]+u[:,:-1])-dt/(2*dx)*\
        (F[:,1:]-F[:,:-1])
        F_med=computeF(gamma, u_med)
        u_n[t,:,1:-1]=u[:,1:-1]-dt/dx*(F_med[:,1:]-F_med[:,:-1])
        u = u_n[t].copy()
    return u_n
    
u_n = richtmyer(u, nt, dt, dx)
r = 50 #time t=0.01 s

print('Variables at t=0.01 s, x= 2.5 m')
print('density (kg/m3)')
print(u_n[r,0,50])
print('velocity (m/s)')
print(u_n[r,1,50]/u_n[r,0,50])
print('pressure (N/m2)')
print((gamma-1)*(u_n[r,2,50]-0.5*u_n[r,1,50]**2/u_n[r,0,50]))

x = numpy.linspace(-10,10,nx)
pyplot.figure(1)
pyplot.plot(x,(gamma-1)*(u_n[r,2,:]-0.5*u_n[r,1,:]**2/u_n[r,0,:]),lw=3)
pyplot.ylabel('pressure')
pyplot.figure(2)
pyplot.plot(x,u_n[r,1,:]/u_n[r,0,:],lw=3)
pyplot.ylabel('velocity')
pyplot.figure(3)
pyplot.plot(x,u_n[r,0,:],lw=3)
pyplot.ylabel('density')

