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
        iprd=list(map(int,fp.readline().lstrip().split(" ")));
        fp.readline();	# Imaging area
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
        im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=2)
        ax.set_xlabel("x [nm]");
        ax.set_ylabel("y [nm]");
        #if self.nshow==0:
        if self.nshow>=0:
            #plt.colorbar(im);
            self.Xa0=Xa;
            self.Xb0=Xb;
        ax.set_xlim((self.Xa0[0],self.Xb0[0]))
        ax.set_ylim((self.Xa0[1],self.Xb0[1]))
        ax.set_title(self.fname+", t="+str(self.time)+"[ps]")
        self.nshow+=1;
    def export(self,fig,fname):
        fig.savefig(fname,bbox_inches="tight")
    

if __name__=="__main__":

    fig=plt.figure();
    ax=fig.add_subplot(111)

    nums=range(0,251,10);

    K=KCELL();
    for k in nums:
        fname="k"+str(k)+".dat"
        K.load(fname);
        K.show(ax)
        fnimg=fname.replace(".dat",".png");
        print(fname+" --->"+fnimg)
        K.export(fig,fnimg)
    
    #plt.show();
