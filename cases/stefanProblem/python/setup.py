# file taken from H. Scheufler's library https://github.com/DLR-RY/TwoPhaseFlow
from scipy import optimize
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
def fprintf(stream, format_spec, *args):
    stream.write(format_spec % args)



# constants

t0 = 0.03  # [s]
tWall = 378.15  # [K]
tSat = 373.15  # [K]
L = 2260e3  # [J/kg]
cpV = 2030.0  # [J/(kg K)]
rhoV = 0.597  # [kg/m**3]
k = 0.025  # [W/(m K)]
Dv = k/(rhoV*cpV)  # [m**2/s]

t_start = 0.0
t_end = 1
dt = 1e-3



def zita():

    def z(_zita):
        return _zita*math.exp(_zita**2)*math.erf(_zita) - cpV*(tWall - tSat)/(np.sqrt(np.pi)*L)

    risultato = optimize.fsolve(z,1.)
    return risultato[0]

#interface position
t = np.arange(t_start, t_end , dt)
s = np.zeros(len(t))

for i in range(0, len(t)):
    s[i] = 2*zita()*np.sqrt(Dv*t[i])



data = pd.DataFrame()
data['t'] = t
data.set_index('t', inplace=True)
data['s'] = s

plt.figure()
data['s'].plot()
data.to_csv('analytical.csv')
plt.plot(t)


#-------------------------------------------------------------------------------#

x = np.arange(0, 1e-3, 1e-6)
T = np.zeros(len(x))
#Temperature profile
for k in range(0, len(x)):

	T[k] = tWall - ((tWall - tSat)/math.erf(zita()))*math.erf((x[k])/(2.*np.sqrt(Dv*0.03)))
	if T[k] < tSat:
		T[k] = tSat
	else:
		T[k] = T[k]
data = pd.DataFrame()
data['x'] = x
data.set_index('x', inplace=True)
data['T'] = T


plt.figure()
data['T'].plot()
data.to_csv('T0.csv',header=None)

# Write position to file
for id in np.arange(0,t.size):
    if t[id] < t0 + dt/2.0 and t[id] > t0 - dt/2.0:
        f = open("initPos.H", "wr")
        fprintf(f,"xpos %.10f;",s[id])
        break

#t = np.linspace(t_start, t_end + dt, dt)
#plt.plot(t, x(t), label='analytical', color='k')
#plt.show()
