#!/home/kazushi/Enthought/Canopy_64bit/User/bin/python

import numpy as np
import matplotlib.pyplot as plt

def ramp(t1,t2,v1,v2,time):
	amp=[];
	T=t2-t1;
	V=v2-v1;
	for t in time:
		if t < t1:
			amp.append(v1);
			continue;
		if t > t2:
			amp.append(v2);
			continue;

		amp.append( (t-t1)/T*V+v1);
	amp=np.array(amp);

	return amp
#----------------------------------

	
PI=np.pi;

a0=60.0;
b0=60.0;
exx_end=-0.0;
eyy_end=-0.0;
exy_end=0.0/180.*PI;
eyx_end=45.0/180.*PI;

t1=0.0; t2=30.0; dt=0.01;
Nt=int((t2-t1)/dt)+1;

time=np.linspace(t1,t2,Nt)

ta=0.0; tb=30.0;
exx=ramp(ta,tb,0.0,exx_end,time);
eyy=ramp(ta,tb,0.0,eyy_end,time);

exy=ramp(ta,tb,0.0,exy_end,time);
eyx=ramp(ta,tb,0.0,eyx_end,time);

a=a0*(1.0+exx);
#a=a0/np.cos(exy)
#b=b0*(1.0+eyy);
b=b0/np.cos(eyx)

fp=open("cell_const.dat","w")
fp.write("# nwt\n");
fp.write(str(Nt)+"\n");
fp.write("# W[0]=a, W[1]=b, W[2]=alph, W[3]=beta \n");
for k in range(len(a)):
	data=str(a[k])+", "+str(b[k])+", "+str(exy[k])+", "+str(eyx[k])+"\n"	
	fp.write(data)

fp.close()
	

ismp=range(0,Nt,100)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_xlim([-a0*0.1, a0*1.2])
ax.set_ylim([-b0*0.1, b0*1.2])
ax.grid(True)


x=[0.0, a[0], a[0], 0.0, 0.0]
y=[0.0,  0.0, b[0],b[0], 0.0]
grph,=ax.plot(x,y,"-o",lw=3,ms=6)
for k in ismp:
	print(k)
	alp=exy[k]
	bet=eyx[k]
	av=np.array([np.cos(alp),np.sin(alp)])*a[k];
	bv=np.array([np.sin(bet),np.cos(bet)])*b[k];

	cv=av+bv;

	x=[0, av[0],cv[0], bv[0],0]
	y=[0, av[1],cv[1], bv[1],0]

	grph.set_data(x,y);
	plt.pause(0.2)

	
"""
time=np.linspace(0,100,200);
amp=ramp(10,30,5,15,time)

ax.plot(time,amp)
ax.grid(True)
plt.show()
"""
