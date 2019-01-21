#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class SHEET:
	def __init__(self,npnt):
		self.x=[];
		self.y=[];
		self.npnt=npnt;
		self.irev=[];
		self.jrev=[];
	def plot(self,ax,clr):
		#xx=self.x+Wd[0]*i;
		#yy=self.y+Wd[0]*j;
		indx= abs(np.diff(self.irev)) >0   
		indx=np.append(indx,True);
		id=np.where(indx == True)
		id=np.append([0],id);
		print(id) 
		print(np.shape(id));
		#indx = indx.tolist();
		ax.plot(self.x,self.y,"-"+clr,lw=2);

class PTCS:
	def __init__(self,fname):
		fp=open(fname,"r");
		fp.readline();
		tt=float(fp.readline());
		fp.readline();
		dat=fp.readline();
		dat=dat.split(" ");
		rho_d=float(dat[0]);	
		poro=float(dat[1]);
		print(str(tt)+"[ps], "+str(rho_d)+"[g/cm3]")
		self.time=tt;

		fp.readline();
		dat=fp.readline();
		dat=dat.strip().split(" ");
		Xc=list(map(float,dat));	

		fp.readline();
		dat=fp.readline();
		dat=dat.strip().split(" ");
		Wd=list(map(float,dat));	

		fp.readline();
		Np=int(fp.readline());

		fp.readline();

		x=[];
		y=[];
		irev=[];
		jrev=[];
		sigp=[];
		sigm=[];
		vx=[];
		vy=[];
		#for row in fp:
		for k in range(Np):
			dat=fp.readline();
			dat=dat.split(" ");
			irev.append(int(dat[0]));	
			jrev.append(int(dat[1]));	
			x.append(float(dat[2]));	
			y.append(float(dat[3]));	
			vx.append(float(dat[4]));	
			vy.append(float(dat[5]));	
			sigp.append(float(dat[6]));	
			sigm.append(float(dat[7]));	
		
		self.x=x;
		self.y=y; #print("x=",x); input("pause")
		self.irev=irev;
		self.jrev=jrev; self.sigp=sigp; self.sigm=sigm;
		self.Np=Np;
		self.vx=vx;
		self.vy=vy;
		fp.close();

	def plot(self,ax,nps,Movie=False):
		if Movie == False:
			ax.cla()
	
		clrs=["r","b","g","c","y","m","k"];
		nclrs=len(clrs);
		n1=0;
		st=0;


		plts=[];
		for n in nps:
			n2=n1+n;
			irev=self.irev[n1:n2];
			jrev=self.jrev[n1:n2];
			x=self.x[n1:n2];
			y=self.y[n1:n2];
			itmp=np.abs(np.diff(irev));
			jtmp=np.abs(np.diff(jrev));
			itmp=np.append(itmp,1);
			jtmp=np.append(jtmp,1);

			tmp=itmp+jtmp;
			indx,=np.where(tmp>0)
			indx+=1;

			m1=0;
			for m2 in indx:
				plt2,=ax.plot(x[m1:m2],y[m1:m2],"-",color="skyblue",ms=2,lw=6);
				plt,=ax.plot(x[m1:m2],y[m1:m2],"-"+clrs[st%nclrs],ms=2,lw=1);
				plts.append(plt);
				m1=m2;
			n1=n2;
			st+=1;

		return plts;

#------------------------ MAIN ROUTINE -------------------------------

if __name__=="__main__":

	nfile1=0;
	nfile2=180;

	args=sys.argv
	narg=len(args)
	if  narg > 1:
		nfile1=int(args[1]);
	if  narg > 2:
		nfile2=int(args[2]);

	fp=open("ptc_nums.dat");
	nums=fp.readlines();
	nums=list(map(int,nums));

	fig=plt.figure();
	ax=fig.add_subplot(111);
	ax.set_aspect(1.0)
	ax.grid(True);

	fnum=np.arange(nfile1,nfile2+1,1);
	fnum=fnum.astype(int)
	nfile=len(fnum);
	Sigp=np.array([]);
	Sigm=np.array([]);
	Vx=np.array([]);
	Vy=np.array([]);
	for k in fnum:
		fname="x"+str(k)+".dat";
		ptc=PTCS(fname);
		Sigp=np.hstack( (Sigp,ptc.sigp) );
		Sigm=np.hstack( (Sigm,ptc.sigm) );
		Vx=np.hstack( (Vx,ptc.vx) );
		Vy=np.hstack( (Vy,ptc.vy) );
		print("sigb=",0.5*(np.mean(ptc.sigp)+np.mean(ptc.sigm)));

	#plt.savefig("x"+str(k)+".png");
	Np=ptc.Np;
	Sigp=np.reshape(Sigp,[nfile,Np])
	Sigm=np.reshape(Sigm,[nfile,Np])
	Vx=np.reshape(Vx,[nfile,Np])
	Vy=np.reshape(Vy,[nfile,Np])
	fig=plt.figure()
	ax=fig.add_subplot(211)
	bx=fig.add_subplot(212)
	im1=ax.imshow(Sigm,aspect="auto",origin="lower",cmap="jet");
	im2=bx.imshow(Sigp,aspect="auto",origin="lower",cmap="jet");
	#im2=bx.imshow(np.sqrt(Vx*Vx+Vy*Vy),aspect="auto",origin="lower",cmap="jet");
	#cbar1=fig.colorbar(im1)	
	#cbar2=fig.colorbar(im2)	
	plt.show()
