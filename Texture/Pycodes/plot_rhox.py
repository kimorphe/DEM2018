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

        rho_count,bins=np.histogram(K,bins=60,density=True);

        nbin=len(rho_count)
        indx=np.arange(nbin)
        rho_bins=0.5*(bins[indx]+bins[indx+1])

        self.count=rho_count;
        self.bins=rho_bins;

        self.K=np.transpose(np.reshape(K,Ndiv))
        self.Ndiv=Ndiv;
        self.time=time;
        self.Xa=Xa;
        self.Xb=Xb;
        self.fname=fname;
        fp.close();

    def show(self,ax,tlt=""):
        #im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=2,cmap="gray");
        if self.nshow>0:
            plt.cla()
        Xa=self.Xa;
        Xb=self.Xb;
        im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",cmap="jet",vmin=0,vmax=1.0)
        #im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",cmap="gray");
        ax.set_xlabel("x [nm]");
        ax.set_ylabel("y [nm]");
        if self.nshow==0:
            plt.colorbar(im);
            self.Xa0=Xa;
            self.Xb0=Xb;
        ax.set_xlim((self.Xa0[0],self.Xb0[0]))
        ax.set_ylim((self.Xa0[1],self.Xb0[1]))
        ax.set_title(self.fname)
        ax.set_xlabel("x [nm]",fontsize=12)
        ax.set_ylabel("y [nm]",fontsize=12)
        ax.tick_params(labelsize=12)
        self.nshow+=1;
    def export(self,fig,fname):
        fig.savefig(fname,bbox_inches="tight")


if __name__=="__main__":

    fig=plt.figure();
    ax=fig.add_subplot(111)
    nums=range(130,251,10);
    K=KCELL();
    for k in nums:
        fname="x"+str(k)+".rho"
        K.load(fname);
        K.show(ax)
        fnimg=fname.replace(".rho","_rho.png");
        print(fname+" --->"+fnimg)
        K.export(fig,fnimg)
    
    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    bx.grid(True)
    for k in nums:
        fname="x"+str(k)+".rho"
        K.load(fname);
        bx.plot(K.bins,K.count)
    bx.set_ylim([0,10])
    bx.set_xlabel(r"number density of particles [nm$^{-2}$]",fontsize=14)
    bx.set_ylabel(r"probability density",fontsize=14)
    bx.tick_params(labelsize=14)
    fig2.savefig("rhox.png",bbox_inches="tight")

    plt.show();
