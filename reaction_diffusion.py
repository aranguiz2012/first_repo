import numpy
from matplotlib import pyplot, cm
numpy.set_printoptions(precision=4)

uvinitial = numpy.load('uvinitial.npz')
ui = uvinitial['U']
vi = uvinitial['V']

n = 192
Du,Dv,F,k = 0.00016, 0.00008, 0.035, 0.065
L = 5.
dh = L/(n-1)
T = 8000
dt = (9./40)*dh**2/max(Du,Dv)
nt = int(T/dt)

print(ui[100,::40])

def ftcs(u,v,nt,Du,Dv,F,k,dh):
    for t in range(nt):
        un = u.copy()
        vn = v.copy()
        
        u[1:-1,1:-1] = un[1:-1,1:-1]+dt*Du/dh**2*\
        (un[1:-1,2:]+un[1:-1,:-2]+un[2:,1:-1]+un[:-2,1:-1]-\
        4*un[1:-1,1:-1])-dt*un[1:-1,1:-1]*vn[1:-1,1:-1]**2+dt*\
        F*(1-un[1:-1,1:-1])
        
        v[1:-1,1:-1] = vn[1:-1,1:-1]+dt*Dv/dh**2*\
        (vn[1:-1,2:]+vn[1:-1,:-2]+vn[2:,1:-1]+vn[:-2,1:-1]-\
        4*vn[1:-1,1:-1])+dt*un[1:-1,1:-1]*vn[1:-1,1:-1]**2-dt*\
        (F+k)*vn[1:-1,1:-1]
        
        u[:,0] = u[:,1]
        u[:,-1] = u[:,-2]
        u[0,:] = u[1,:]
        u[-1,:] = u[-2,:]
        
        v[:,0] = v[:,1]
        v[:,-1] = v[:,-2]
        v[0,:] = v[1,:]
        v[-1,:] = v[-2,:]
        
        
    return u,v


sol = ftcs(ui.copy(),vi.copy(),nt,Du,Dv,F,k,dh)
uans = sol[0]
vans = sol[1]
print(uans[100,::40])

fig = pyplot.figure(figsize=(8,5))
pyplot.subplot(121)
pyplot.imshow(ui,cmap=cm.RdBu)
pyplot.xticks([]),pyplot.yticks([])

pyplot.subplot(122)
pyplot.imshow(vi,cmap=cm.RdBu)
pyplot.xticks([]),pyplot.yticks([])

fig2 = pyplot.figure(figsize=(8,5))
pyplot.imshow(uans,cmap=cm.RdBu)
pyplot.xticks([]),pyplot.yticks([])
