# Need to import the plotting package:
import matplotlib.pyplot as plt
import numpy as np

base   = 'normR_vs_iteration_'
kappa  = '2'
newton = '1'
inc    = '1'

filename1     = base + kappa + '_' + newton + '_' + inc + '.dat' 

# Read the file. 
f1 = open(filename1,'r')
lines1 = f1.readlines()
f1.close()

newton = '1'; inc    = '2'
filename2     = base + kappa + '_' + newton + '_' + inc + '.dat' 

# Read the file. 
f2 = open(filename2,'r')
lines2 = f2.readlines()
f2.close()

newton = '2'; inc    = '1'
filename3     = base + kappa + '_' + newton + '_' + inc + '.dat' 

# Read the file. 
f3 = open(filename3,'r')
lines3 = f3.readlines()
f3.close()

newton = '2'; inc    = '2'
filename4     = base + kappa + '_' + newton + '_' + inc + '.dat' 

# Read the file. 
f4 = open(filename4,'r')
lines4 = f4.readlines()
f4.close()

# initialize some variable to be lists:
x1 = []; y1 = []
x2 = []; y2 = []
x3 = []; y3 = []
x4 = []; y4 = []

# scan the rows of the file stored in lines, and put the values into some variables:
for line in lines1[2:]:
    p = line.split()
    x1.append(float(p[0])); y1.append(float(p[1]))

for line in lines2[2:]:
    p = line.split()
    x2.append(float(p[0])); y2.append(float(p[1]))

for line in lines3[2:]:
    p = line.split()
    x3.append(float(p[0])); y3.append(float(p[1]))

for line in lines4[2:]:
    p = line.split()
    x4.append(float(p[0])); y4.append(float(p[1]))

xv1 = np.array(x1); yv1 = np.array(y1)
xv2 = np.array(x2); yv2 = np.array(y2)
xv3 = np.array(x3); yv3 = np.array(y3)
xv4 = np.array(x4); yv4 = np.array(y4)

# now, plot the data:
plt.plot(xv1, yv1, color = "blue", linewidth =  2.5, \
         linestyle = 'dashed', marker = 'h', markersize =10, label="FN")
plt.plot(xv2, yv2, color = "blue", linewidth =  2.5, \
         linestyle = 'dotted', marker = 's', label="FNIL")
plt.plot(xv3, yv3, color = "blue", linewidth =  2.5, \
         linestyle = 'dashed', marker = 'o', label="MN")
plt.plot(xv4, yv4, color = "blue", linewidth =  2.5, \
         linestyle = 'dashdot', marker = 'D', label="MNIL")

plt.yscale('log')
plt.xscale('log')

#plt.xlim([ 0  , 1])
#plt.ylim([-0.6, 1])

plt.title(' Kappa case : '+kappa)
plt.xlabel('iterations')
plt.ylabel('norm(R)')
#plt.legend(loc='upper left', fontsize = 10)
plt.legend(loc='upper center', bbox_to_anchor = (0.80, 0.30), fontsize = 12)

#plt.show()
filename      = 'normR_kapp_' + kappa 
plt.savefig(filename + '.png')

