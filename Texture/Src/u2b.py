#!/home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt

fname="Sr1/u2b.out"

class DISP:
    def load(self,fname):
        fp=open(fname,"r")

        u2=[];
        ux=[];
        uy=[];
        for row in fp:
            dat=row.strip().split(" ")
            u2.append(float(dat[0]))
            ux.append(float(dat[1]))
            uy.append(float(dat[2]))
        fp.close()

        self.u2=np.array(u2)/2
        self.ux=np.array(ux)
        self.uy=np.array(uy)
        self.nt=len(u2);

    def set_xaxis(self,Lev=8):
        self.Lev=8
        self.dt=2**(-2*(Lev+1))
        self.tt=np.arange(self.nt)*self.dt;

    def fit(self,a1=0.1,a2=1.1, na=41):
        aexp=np.linspace(a1,a2,na)
        Ks=[];
        Res=[];
        tt=self.tt
        for a in aexp:
            K=np.sum(self.u2*(tt**a))/np.sum(tt**(2*a))*0.5
            u2fit=2*K*(tt**a)
            err=(self.u2-u2fit)
            res=np.sum(err*err)/np.sum(self.u2*self.u2)
            Res.append(np.sqrt(res))
            Ks.append(K);
        Res=np.array(Res)
        Ks=np.array(Ks)
        indx=np.argmin(np.abs(Res))
        self.K=Ks[indx];
        self.a=aexp[indx]
    def lin_fit(self):
        #n2=int(self.nt*3/5)
        n2=int(self.nt*1/2)
        x=self.tt[n2:-1];
        y=self.u2[n2:-1];
        X=np.mean(x)
        Y=np.mean(y)
        XX=np.mean(x*x)
        XY=np.mean(x*y)
        A=np.array([[XX,X],[X,1]])
        B=np.array([XY,Y])

        sol=np.linalg.solve(A,B)

        self.alph=sol[0];
        self.beta=sol[1];
                

    def u2td(self):
        return(self.a*self.K*(self.tt**self.a)/(self.tt+0.0))
    def u2t(self):
        return(2*self.K*(self.tt**self.a))
    def u2t_lin(self):
        return(self.alph*self.tt+self.beta)


if __name__=="__main__":


    fig=plt.figure()
    ax=fig.add_subplot(121)
    bx=fig.add_subplot(122)
    ax.grid(True)
    bx.grid(True)


    nums=range(1,10,1)

    K=[];
    a=[];
    D=[];
    for k in nums:
        fname="Sr"+str(k)+"/u2b.out"
        print(fname)
        u2b=DISP()
        u2b.load(fname)
        u2b.set_xaxis(Lev=8)
        u2b.fit()
        u2b.lin_fit()

        ax.plot(u2b.tt,u2b.u2)
        ax.plot(u2b.tt,u2b.u2t())
        ax.plot(u2b.tt,u2b.u2t_lin(),"--")
        print("K=",u2b.K)
        print("a=",u2b.a)
        print("alph=",u2b.alph)
        print("beta=",u2b.beta)

        bx.plot(u2b.tt,u2b.u2td())
        bx.set_ylim([0,1])

        K.append(u2b.K);
        a.append(u2b.a);
        D.append(u2b.alph*0.5);

    fig2=plt.figure()
    cx=fig2.add_subplot(121)
    dx=fig2.add_subplot(122)
    sr=np.array(nums)*0.1+0.1;
    cx.plot(sr,D,"-o")
    cx.plot(sr,K,"-s")
    cx.grid(True)
    dx.plot(sr,a,"-s")
    dx.grid(True)


plt.show()


