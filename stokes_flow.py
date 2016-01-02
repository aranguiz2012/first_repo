from matplotlib import pyplot, cm
import numpy

nx = 41
ny = 41
l = 1.
h = 1.
dx = l/(nx-1)
dy = h/(ny-1)
l1_target = 1e-6
omega_i = numpy.ones((ny,nx))
psi_i = numpy.zeros((ny,nx))

def L1norm(new, old):
    norm = numpy.sum(numpy.abs(new-old))
    return norm
    
def stokes(nx,ny,dx,dy,l,h,psi,omega,l1_target):  
    l1_norm = 1
    iterations = 0
    
    while l1_norm > l1_target:
        psid = psi.copy()
        omegad = omega.copy()
        
        omega[1:-1,1:-1]=0.25*(omegad[2:,1:-1]+omegad[:-2,1:-1]+\
        omegad[1:-1,2:]+omegad[1:-1,:-2])
              
        omega[-1,1:-1] = (8*psi[-2,1:-1]-\
        psi[-3,1:-1])*(-1)/(2*dy**2)-3./dy
        omega[0,1:-1] = (8*psi[1,1:-1]-\
        psi[2,1:-1])*(-1)/(2*dy**2)
        omega[1:-1,-1] = (8*psi[1:-1,-2]-\
        psi[1:-1,-3])*(-1)/(2*dx**2)
        omega[1:-1,0] = (8*psi[1:-1,1]-\
        psi[1:-1,2])*(-1)/(2*dx**2)
        
        psi[1:-1,1:-1]=0.25*(psi[2:,1:-1]+psi[:-2,1:-1]+\
        psi[1:-1,2:]+psi[1:-1,:-2]+omega[1:-1,1:-1]*dy**2)   
        
        l1_normpsi = L1norm(psi,psid)
        l1_normomega = L1norm(omega, omegad)
        l1_norm = max(l1_normpsi,l1_normomega)
        iterations += 1
    print ('Number of GS iterations: {}'.\
        format(iterations))
    return psi, omega
    
psi, omega = stokes(nx,ny,dx,dy,l,h,psi_i.copy(),\
    omega_i.copy(),l1_target)
x = numpy.linspace(0,l,nx)
y = numpy.linspace(0,h,ny)
pyplot.figure(figsize=(8,5))
pyplot.contourf(x,y,psi,10,cmap=cm.viridis)
pyplot.colorbar()
pyplot.figure(figsize=(8,5))
pyplot.contourf(x,y,omega,10,cmap=cm.viridis_r)
pyplot.xlabel('$x$')
pyplot.ylabel('$y$')
pyplot.colorbar()

print(numpy.round(numpy.max(numpy.abs(psi)),4))
print(numpy.round(numpy.max(numpy.abs(omega)),4))
print(numpy.round(psi[32,::8], 4))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
