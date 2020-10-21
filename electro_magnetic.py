# electro_magnetic.py

import numpy as np
from scipy.integrate import odeint

class ElectroMagneticField():
	'''
	Electro-Magnetic field
	'''
	def __init__(self, E, H, is_hom=True, is_con=True):
		'''
		Constructor
		'''
		self._E = E
		self._H = H
		self.is_hom = is_hom	# is_hom == is homogeneous
		self.is_con = is_con	# is_con == is constant

	def E(self, t, x):
		'''
		Electric field
		'''
		if self.is_hom and self.is_con:
			return self._E
		elif self.is_hom:
			return self._E(t)
		elif self.is_con:
			return self._E(x)
		else:
			return self._E(t,x)

	def H(self, t, x):
		'''
		Magnetic field
		'''
		if self.is_hom and self.is_con:
			return self._H
		elif self.is_hom:
			return self._H(t)
		elif self.is_con:
			return self._H(x)
		else:
			return self._H(t,x)


class Particle():
	'''
	Particle
	'''
	def __init__(self, x, p):
		'''
		Constructor
		'''
		self.x = x # particle coordinate
		self.p = p # particle momentum

	def motion(self, field, time_data):
		'''
		Trajectory of particle and changing its position
		'''
		xp0 = np.hstack((self.x, self.p))
		xp_data = odeint(calc_der, xp0, time_data, args=(field,))

		x_data = xp_data[:,0:3]
		p_data = xp_data[:,3:6]

		self.x = x_data[-1]
		self.p = p_data[-1]

		return x_data, p_data


def calc_der(xp, t, field):
	'''
	Calculates the derivative of the coordinate
	and momentum
	''' 
	x = xp[0:3]	# coordinate
	p = xp[3:6]	# monemtum

	v = p / np.sqrt(1 + np.dot(p,p))	# velocity

	dx = v 	# coordinate derivative
	dp = field.E(t,x) + np.cross(v, field.H(t,x))	# momentun derivative

	return np.hstack((dx,dp))