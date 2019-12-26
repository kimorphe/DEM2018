import numpy as np
import matplotlib.pyplot as plt

class XRD:
    def load(self,fname):
        fp=open(fname,"r")
        th2=[];
        Ith2=[];
        for row in fp:
            dat=list(map(float,row.strip().split(" ")))
            th2.append(dat[0]);
            Ith2.append(dat[1]);

        
        self.th2=np.array(th2);
        self.Ith2=np.array(Ith2);

    def show(self,ax,i_start=0):
        ax.plot(self.th2[i_start:-1],self.Ith2[i_start:-1],linewidth=1.5)


if __name__=="__main__":

    #----------- XRD Intensity Plot ---------------------
    fig1=plt.figure()
    ax1=fig1.add_subplot(111)
    ax1.grid(True)
    th2_max=15.0;
    Imax=15.0;
    ax1.set_xlim([0,th2_max])
    ax1.set_ylim([0,Imax])
    ax1.set_xlabel(r'2$\theta$[deg]',fontsize=14);
    ax1.set_ylabel("Intensity ",fontsize=14);
    ax1.tick_params(labelsize=14)

    Ith=XRD();
    nums=range(130,251,20);
    for num in nums:
        fname="k"+str(num)+".xrd"
        Ith.load(fname);
        Ith.show(ax1,i_start=1)

    #----------- Radial Distribution Function ------------
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.grid(True)
    r_max=15.0;
    ax2.set_xlim([0,r_max])
    ax2.set_xlabel(r'radial distance [nm]',fontsize=14);
    ax2.set_ylabel("radial ditribution function f(r) ",fontsize=14);
    ax2.tick_params(labelsize=14)
    fr=XRD();
    nums=range(130,251,20);
    for num in nums:
        fname="x"+str(num)+".rad"
        fr.load(fname);
        fr.show(ax2)
    #------------ Normal Vector Histogram---------------
    fig3=plt.figure()
    #ax3=fig3.add_subplot(111,projection="polar")
    ax3=fig3.add_subplot(111)
    ax3.grid(True)
    r_max=15.0;
    ax3.set_xlabel(r'normal direction [deg]',fontsize=14);
    ax3.set_ylabel("density",fontsize=14);
    ax3.tick_params(labelsize=14)

    nr=XRD();
    nums=range(130,251,20);
    for num in nums:
        fname="x"+str(num)+".nml"
        nr.load(fname);
        nx=nr.th2
        ny=nr.Ith2

        thp=np.arccos( nx)
        thm=np.arccos(-nx)
        indxp=np.argwhere(ny<0);
        thp[indxp]=2.*np.pi-thp[indxp];
        indxm=np.argwhere(-ny<0);
        thm[indxm]=2.*np.pi-thm[indxm];
        th=np.hstack([thp,thm])

        th_hist,th_bins=np.histogram(th,bins=120,density=True)
        th_hist=np.array(th_hist)
        th_bins=np.array(th_bins)
        n=len(th_bins)-1
        indx1=np.arange(n)
        indx2=np.arange(n)+1
        th_h=(th_bins[indx1]+th_bins[indx2])*0.5
        #ax3.hist(th,bins=60)
        #ax3.bar(th_h,th_hist,width=0.1)
        ax3.plot(th_h/np.pi*180.,th_hist)
        ax3.set_xlim([0,180.])
        ax3.set_ylim([0,0.4])

    #---------------------------------------------------

    fig1.savefig("xrd.png",bbox_inches="tight")
    fig2.savefig("rad.png",bbox_inches="tight")
    fig3.savefig("nml.png",bbox_inches="tight")
    plt.show()

