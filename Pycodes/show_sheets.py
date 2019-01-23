#!/home/kazushi/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt


fp=open("gen_sheet.out","r");

fp.readline();

x=[]; y=[]; th=[];
Wd=[];
for row in fp:
	dat=list(map(float,row.strip().split(",")));

	x.append(dat[0]);
	y.append(dat[1]);
	Wd.append(dat[2]);
	th.append(dat[3]);


x=np.array(x);
y=np.array(y);
Wd=np.array(Wd);
th=np.array(th)/180.*np.pi;
nst=len(x);


fig=plt.figure(figsize=(8,4))
ax=fig.add_subplot(121)
bx=fig.add_subplot(122)
ax.grid(True);
bx.grid(True);
for k in range(nst):
	xc=x[k];
	yc=x[k];

	et=np.array([np.cos(th[k]),np.sin(th[k])]);
	xc=np.array([x[k],y[k]]);

	x1=xc-0.5*Wd[k]*et;
	x2=xc+0.5*Wd[k]*et;

	xx=[x1[0],x2[0]]
	yy=[x1[1],x2[1]]


	ax.plot(xx,yy,'-');

ax.set_xlabel("x[nm]")
ax.set_ylabel("y[nm]")
bx.set_xlabel("length [nm]")
bx.set_ylabel("count")
bx.hist(Wd,bins=15)
fig.savefig("model.png",bbox_inches="tight")
plt.show()
