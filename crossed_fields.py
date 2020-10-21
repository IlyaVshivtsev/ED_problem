#! /usr/bin/python3

from electro_magnetic import *
from scipy.integrate import odeint
from matplotlib import pyplot as plt


# initial values
# --------------------------------------------------------------------------------
EH	= 1.0 # commom amplitude
E 	= np.array([0,EH,0])	# electric field
H 	= np.array([0,0,EH])	# magnetic field
field = ElectroMagneticField(E, H, is_hom=True, is_con=True) # homogeneous constant field

x0 	= np.array([0,0,0])	# initial coordinate
p0 	= np.array([0,0,0])	# initial momentum
# --------------------------------------------------------------------------------
	

# firstly, exact solution. Assume x(t=0) = 0 and p2(t=0) = 0
# --------------------------------------------------------------------------------
alpha 	= np.sqrt(1 + np.dot(p0,p0)) - p0[0] # integral of motion
eps 	= 1 + p0[2]*p0[2] # constant parameter

p2_data = np.linspace(p0[1], p0[1]+10, 64) # momentum p2
time_data = p2_data * (eps**2-alpha**2)/(2*EH*alpha**2) + p2_data**3 / (6*EH*alpha**2)

# exact x-,y- and z-component
x1_exact_data = p2_data * (eps**2-alpha**2)/(2*EH*alpha**2) + p2_data**3 / (6*EH*alpha**2)
x2_exact_data = p2_data**2 / (2*EH*alpha)
x3_exact_data = p2_data * p0[2]/(EH*alpha)
# ---------------------------------------------------------------------------------


# numerical solution
# ---------------------------------------------------------------------------------
xp_init	= np.hstack((x0, p0))
xp_data	= odeint(calc_der, xp_init, time_data, args=(field,))	# solves equation

x1_data = xp_data[:,0]	# x-component
x2_data = xp_data[:,1]	# y-component
x3_data = xp_data[:,2]	# z-component
#----------------------------------------------------------------------------------


# output
# ---------------------------------------------------------------------------------
plt.plot(x1_data, x2_data, 'o', x1_exact_data, x2_exact_data)
plt.show()
plt.plot(x2_data, x3_data, 'o', x2_exact_data, x3_exact_data)
plt.show()
plt.plot(x3_data, x1_data, 'o', x3_exact_data, x1_exact_data)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(x1_data, x2_data, x3_data)
plt.show()
# ---------------------------------------------------------------------------------
