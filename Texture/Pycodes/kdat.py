import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class Kdat:
    def __init__(self):
        self.nshow=0

    def load(self,fname):
        fp=open(fname,"r")
        fp.readline();
        time=float(fp.readline().strip());
        fp.readline(); # computational domain
        Ka=list(map(float,fp.readline().lstrip().split(" ")));
        fp.readline(); # computational domain
        Kb=list(map(float,fp.readline().lstrip().split(" ")));
        fp.readline(); # computational domain
        Ndiv=list(map(int,fp.readline().lstrip().split(" ")));
        self.Ndiv=Ndiv;
        self.time=time;
        self.Ka=Ka;
        self.Kb=Kb;
        self.fname=fname;
        fp.close();

        self.dat=pd.read_csv(fname,skiprows=9);
        Kr=np.array(self.dat['Real'])
        Ki=np.array(self.dat['Imag'])
        Z=Kr+1j*Ki;

        self.Z=np.transpose(np.reshape(Z,Ndiv));
        #self.Z=np.fft.fftshift(Z);


    def show(self,ax):
        Ka=self.Ka;
        Kb=self.Kb;
        #kb=np.mean(np.abs(self.Z[:,:]))
        #ks=np.std(np.abs(self.Z[:,:]))
        kmax=np.max(np.abs(self.Z[:,:]))
        #print(kb,ks,kmax)
        im=ax.imshow(np.abs(self.Z),origin="lower",extent=[Ka[0],Kb[0],Ka[1],Kb[1]],interpolation="bilinear",cmap="jet",vmin=0,vmax=kmax*0.05);
        #im=ax.imshow(np.abs(self.Z),origin="lower",extent=[Ka[0],Kb[0],Ka[1],Kb[1]],interpolation="bilinear",cmap="jet")
        #plt.colorbar(im)
        if self.nshow==0:
            self.Ka0=Ka;
            self.Kb0=Kb;
            ax.set_xlabel("$k_x$ [nm$^{-1}$]",fontsize=12)
            ax.set_ylabel("$k_y$ [nm$^{-1}$]",fontsize=12)
        ax.set_xlim([self.Ka0[0],self.Kb0[0]])
        ax.set_ylim([self.Ka0[1],self.Kb0[1]])
        ax.set_title(self.fname)
        ax.tick_params(labelsize=12)

        self.nshow+=1



if __name__=="__main__":


    Dat=Kdat();
    fig=plt.figure()
    ax=fig.add_subplot(111)

    nums=range(0,251,10)
    for k in nums:
        fname="k"+str(k)+".fft"
        print(fname)
        Dat.load(fname);
        Dat.show(ax)
        fname_out="fft"+str(k)+".png"
        fig.savefig(fname_out,bbox_inches="tight")


    #plt.show()
