import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class Kdat:
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
        print(self.dat.columns)
        Kr=np.array(self.dat['Real'])
        Ki=np.array(self.dat['Imag'])
        Z=Kr+1j*Ki;

        self.Z=np.transpose(np.reshape(Z,Ndiv));
        #self.Z=np.fft.fftshift(Z);


    def show(self,ax):
        Ka=self.Ka;
        Kb=self.Kb;
        im=ax.imshow(np.abs(self.Z),origin="lower",extent=[Ka[0],Kb[0],Ka[1],Kb[1]],interpolation="bilinear",cmap="jet",vmin=0,vmax=0.005);
        #im=ax.imshow(np.abs(self.Z),origin="lower",extent=[Ka[0],Kb[0],Ka[1],Kb[1]],interpolation="bilinear",cmap="jet")
        plt.colorbar(im)



if __name__=="__main__":

    fname="kn240.dat"

    Dat=Kdat();
    Dat.load(fname);
    fig=plt.figure()
    ax=fig.add_subplot(111)

    Dat.show(ax)
    plt.show()
