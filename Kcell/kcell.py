#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys


#fname="kcell_dry.dat";
fname="kcell.dat";

if len(sys.argv)>1:
    fname=sys.argv[1]

fp=open(fname,"r");

fp.readline();
time=float(fp.readline().strip());
fp.readline(); # computational domain
Xa=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline(); # computational domain
Xb=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline(); # computational domain
Ndiv=list(map(int,fp.readline().lstrip().split(" ")));

fp.readline();	# Imaging area
K=list(map(float,fp.readlines()));	# Imaging area
K=np.transpose(np.reshape(K,Ndiv))

print(np.shape(K));

fig=plt.figure();
ax=fig.add_subplot(111)
#im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=2,cmap="gray");
im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=2)
plt.colorbar(im);
ax.set_title(fname);
plt.show();


fp.close();
