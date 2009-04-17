#!/usr/bin/env python

from __future__ import division
from scipy import *
from scipy.linalg import *
from numpy import *
from pylab import *



##### Magnetic field
#~ A = 1
#~ B = 0
#~ x = 0.5
#~ blist = linspace(0,3)

#~ elev = zeros((blist.size,2))
#~ for k in range(0,blist.size):
    #~ H = array([[A/2, x*blist[k]],[x*blist[k], -A/2]])
    #~ temp = eigvals(H)
    #~ elev[k,0] = temp[0]
    #~ elev[k,1] = temp[1]

#~ plot(blist,elev[:,0],blist,elev[:,1])
#~ xlabel('Magnetic field')
#~ ylabel('Level energy')
#~ show()


d = 1
om = 10
H = array([[d/2, om/2],[om, -d/2]])
temp = eigvals(H)
print diag(temp)
X =     dot(diag(temp).T,dot(H,diag(temp)))
print X/H




#~ omlist = linspace(0,3)

#~ elev = zeros((omlist.size,2))
#~ for k in range(0,omlist.size):
    #~ H = array([[d/2, omlist[k]/2],[omlist[k]/2, -d/2]])
    #~ temp = eigvals(H)
    #~ elev[k,0] = temp[0]
    #~ elev[k,1] = temp[1]

#~ plot(omlist,elev[:,0]-1,omlist,elev[:,1]+1)
#~ xlabel('Rabi-frequency')
#~ ylabel('Level energy')
#~ show()