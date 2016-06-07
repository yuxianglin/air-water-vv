from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put tankVertices and caissonVertices arrays here

#tankVertices
tv=np.array([[0.0, 0.0],
 [1.3890866208485599, 0.0],
 [4.1672598625456798, 0.0],
 [4.6672598625456798, 0.25],
 [4.9172598625456798, 0.25],
 [5.2922598625456798, 0.0],
 [6.6813464833942398, 0.0],
 [9.4595197250913596, 0.0],
 [9.4595197250913596, 0.75],
 [6.6813464833942398, 0.75],
 [1.3890866208485599, 0.75],
 [0.0, 0.75]],
                         )
nt=len(tv)

#caissonVertices
cv=np.array( [ 0.,   0.  ])
nc=len(cv)


xt=[]
yt=[]
xc=[]
yc=[]

for i in range(nt):
    xt.append(tv[i][0])
    yt.append(tv[i][1])
xt.append(tv[0][0])
yt.append(tv[0][1])
"""
for j in range(nc):
    xc.append(cv[j][0])
    yc.append(cv[j][1])
xc.append(cv[0][0])
yc.append(cv[0][1])
"""

# Plot geometry
import matplotlib.pyplot as plt
plt.plot(xt,yt)
#plt.plot(xc,yc)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.suptitle('geometry')
plt.show()
savefig('geometry.png')





