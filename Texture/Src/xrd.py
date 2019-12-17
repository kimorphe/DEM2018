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
        #im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",vmin=-30,vmax=70,cmap="gnuplot")
        ax.set_xlabel("x [nm]");
        im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",vmin=0,vmax=500,cmap="jet")
        ax.set_ylabel("y [nm]");
        if self.nshow==0:
            plt.colorbar(im);
            self.Xa0=Xa;
            self.Xb0=Xb;
        ax.set_xlim((self.Xa0[0],self.Xb0[0]))
        ax.set_ylim((self.Xa0[1],self.Xb0[1]))
        ax.set_title(self.fname+", t="+str(self.time)+"[ps]")
        self.nshow+=1;
    def export(self,fig,fname):
        fig.savefig(fname,bbox_inches="tight")

    def bin_filt(self):
        indx=self.K>=0;
        self.K[indx]=1
        indx=self.K<0;
        self.K[indx]=0

    def FFT2D(self):
        K=np.fft.fft2(self.K)
        K=np.fft.fftshift(K)
        self.K=np.abs(K)

    

if __name__=="__main__":

    fig=plt.figure();
    ax=fig.add_subplot(111)

    nums=range(240,251,10);

    K=KCELL();
    for k in nums:
        fname="kn"+str(k)+".dat"
        K.load(fname);
        """
        K.show(ax)
        fnimg=fname.replace(".dat",".png");
        print(fname+" --->"+fnimg)
        K.export(fig,fnimg)
        """

        K.bin_filt()
        K.FFT2D()
        K.show(ax)
        fnimg=fname.replace(".dat","k.png");
        print(fname+" --->"+fnimg)
        K.export(fig,fnimg)

    
    #plt.show();
