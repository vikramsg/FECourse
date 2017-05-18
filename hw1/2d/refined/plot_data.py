# Need to import the plotting package:
import matplotlib.pyplot as plt
import numpy as np

base   = 'node_u_vs_load_'
kappa  = '3'
newton = '1'
inc    = '2'

filename      = base + kappa + '_' + newton + '_' + inc + '.dat' 

# Read the file. 
f = open(filename, 'r')
lines = f.readlines()
f.close()

# initialize some variable to be lists:
x1 = []
y1 = []

# scan the rows of the file stored in lines, and put the values into some variables:
for line in lines[1:]:
    p = line.split()
    x1.append(float(p[0]))
    y1.append(float(p[1]))

xv1 = np.array(x1)
yv1 = np.array(y1)

# now, plot the data:
plt.plot(xv1, yv1, color = "blue", linewidth =  2.5, \
         linestyle = "dashed", label="Lax u, n = 101")

#plt.xlim([ 0  , 1])
#plt.ylim([-0.6, 1])

plt.title(' Kappa case : '+kappa)
plt.xlabel('load')
plt.ylabel('U')
#plt.legend(loc='upper left', fontsize = 10)
###plt.legend(loc='upper center', bbox_to_anchor = (0.95, 1.12), fontsize = 10)

plt.show()
#filename      = 'Velocities_Re_' + Re 
#plt.savefig(filename + '.png')
#
