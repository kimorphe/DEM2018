#!/home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt

class KCELL:
    def __init__(self):
        self.nshow=0;
    def load(self,fname):
        fp=open(fname,"r")
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
        self.K=np.transpose(np.reshape(K,Ndiv))
        self.Ndiv=np.array(Ndiv);
        self.time=time;
        self.Xa=np.array(Xa);
        self.Xb=np.array(Xb);
        self.fname=fname;
        self.Wd=self.Xb-self.Xa;
        self.dx=self.Wd/self.Ndiv;
        fp.close();

    def show(self,ax,tlt=""):
        #im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=2,cmap="gray");
        if self.nshow>0:
            plt.cla()
        Xa=self.Xa;
        Xb=self.Xb;
        #im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",vmin=-30,vmax=70,cmap="gnuplot")
        ax.set_xlabel("x [nm]");
        im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=10000,cmap="jet")
        #im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",cmap="jet")
        ax.set_ylabel("y [nm]");
        if self.nshow==0:
            plt.colorbar(im);
            self.Xa0=Xa;
            self.Xb0=Xb;
        ax.set_xlim((self.Xa0[0],self.Xb0[0]))
        ax.set_ylim((self.Xa0[1],self.Xb0[1]))
        ax.set_title(self.fname+", t="+str(self.time)+"[ps]")
        self.nshow+=1;

    def showK(self,ax,tlt=""):
        Ka=self.Ka;
        Kb=self.Kb;
        im=ax.imshow(np.abs(self.K),origin="lower",extent=[Ka[0],Kb[0],Ka[1],Kb[1]],interpolation="none",cmap="gray");

    def export(self,fig,fname):
        fig.savefig(fname,bbox_inches="tight")

    def bin_filt(self,val):
        indx=np.argwhere(self.K[:]<2);
        self.K[indx[:,0],indx[:,1]]=0
        #indx=(self.K<=1);
        #self.K[indx]=0

    def FFT2D(self):
        K=np.fft.fft2(self.K)
        K=np.fft.fftshift(K)
        #self.K=np.abs(K)
        self.K=K

        dk=np.zeros(2)
        #dk=2.*np.pi/Wd;
        dk=1./self.Wd;
        self.Ka=np.zeros(2)
        self.Kb=np.zeros(2)
        self.Kb=dk*self.Ndiv;

        self.Ka=-0.5*self.Kb;
        self.Kb= 0.5*self.Kb;
        self.dk=dk;

    

if __name__=="__main__":

    fig=plt.figure();
    ax=fig.add_subplot(111)

    nums=range(240,251,10);

    lmb=0.15418;
    kin=1./lmb;
    th2=np.linspace(0,np.pi*0.15,300)

    K=KCELL();
    for k in nums:
        fname="k"+str(k)+".dat"
        K.load(fname);
        val=2
        K.bin_filt(val)
        """
        K.show(ax)
        fnimg=fname.replace(".dat",".png");
        print(fname+" --->"+fnimg)
        K.export(fig,fnimg)
        """
        K.FFT2D()
        #K.show(ax)
        K.showK(ax)
        fnimg=fname.replace(".dat","k.png");
        print(fname+" --->"+fnimg)
        K.export(fig,fnimg)

    F=np.zeros(len(th2))
    th0s=np.linspace(0.0,np.pi*0.25,30)
    for th0 in th0s:
        xi0=np.array([np.cos(th0),np.sin(th0)])
        iad=0
        for th in th2:
            xi1=np.array([np.cos(th+th0),np.sin(th+th0)])
            xi=xi1-xi0;
            kv=kin*xi;
            ii=int((kv[0]-K.Ka[0])/K.dk[0])
            jj=int((kv[1]-K.Ka[1])/K.dk[1])
            f=np.abs(K.K[ii][jj])
            F[iad]+=f;
            iad+=1

    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(th2/np.pi*180.,F/len(th0s))
    ax2.grid(True)
    ax2.set_ylim([0,10000])
    plt.show();
