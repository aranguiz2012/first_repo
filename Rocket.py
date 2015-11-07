import numpy

#parametros
m_s = 50
g = 9.81
rho = 1.091
r = 0.50
A = numpy.pi*(r**2)
v_e = 325
c_d = 0.15

#condiciones iniciales
h0 = 0
v0 = 0

def f(u,s):
	
	v = u[1]
	return numpy.array([v,-g + (mpunto[s]*v_e)/(m_s+m_p[s])-\
                        (rho*v*abs(v)*A*c_d)/(2*(m_s+m_p[s]))])

T = 50
dt = 0.1
N = int(T/dt) + 1
t = numpy.linspace(0, T, N)
mpunto = numpy.zeros(t.shape)
for i in range(51):
	mpunto[i] = 20
m_p = numpy.zeros(t.shape)
m_p[0] = 100
for i in range(1,51):
	m_p[i] = m_p[0] - 20*t[i]
u = numpy.zeros((N,2))
u[0] = numpy.array([h0, v0])

for n in range(N-1):
	u[n+1] = u[n] + dt*f(u[n],n)
h = u[:,0]
v = u[:,1]
print(h)
print(numpy.amax(h))
print(h.argmax())