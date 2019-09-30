# -*- coding: utf-8 -*-
"""
Created on Mon Feb 07 15:42:59 2011

@author: Sasha
"""

import mpmath as mm
import numpy as np
import json
   
class barWithDispersion:
    def __init__(self, r0=5e-3, rho=7800.,
                 E=200e9, nu=0.29):
        self.r0=r0
        self.rho=rho
        self.E=E
        self.nu=nu
        self.l=E*nu/(1+nu)/(1-2*nu)
        self.mu=E/2./(1+nu)
        self.cb=np.sqrt(E/rho)
        self.ct=np.sqrt(self.mu/rho)
        self.cl=np.sqrt((self.l+2*self.mu)/rho)
        self.cr=(0.862+1.14*nu)/(1+nu)*self.ct

    def f(self, ksi, omega):
        alpha2=self.rho*omega**2/(self.l+2*self.mu)-ksi**2
        beta2=self.rho*omega**2/self.mu-ksi**2
        alpha=mm.sqrt(alpha2)
        beta=mm.sqrt(beta2)
        rez=2.*alpha/self.r0*(beta2+ksi**2)*mm.besselj(1,alpha*self.r0)*mm.besselj(1,beta*self.r0)
        rez-=(beta2-ksi**2)**2*mm.besselj(0,alpha*self.r0)*mm.besselj(1,beta*self.r0)+4*ksi**2*alpha*beta*mm.besselj(1,alpha*self.r0)*mm.besselj(0,beta*self.r0)
        return abs(rez)
    
    def solve(self, omegamax=5e6, domega0=lambda x: 10e3):
        omega=[0]#[0,1000.,2000.]
        ksi=[0]#[0,omega[1]/cb, omega[2]/cb]
        domega=domega0(0)
        while omega[-1]<=omegamax and domega>1e-6*domega0(omega[-1]):
            xx = omega[-1]+domega
            q=lambda x: self.f(x, xx)/xx**2
            try:
                if len(omega)>=3: 
                    y0=(xx-omega[-2])/(omega[-1]-omega[-2])*(ksi[-1]-ksi[-2])+ksi[-2]
                else:
                    y0=1.0
                y=mm.findroot(q, y0, solver = 'mnewton')
                omega.append(omega[-1]+domega)
                ksi.append(y.real)
                print("converged:", omega[-1])
                domega=domega0(omega[-1])
            except ValueError:
                    print("not converged:", omega[-1]+domega, domega)
                    domega*=0.5
        
        self.ksi_omega=lambda x: np.interp(x, omega, ksi)
        self.omega=np.array(omega)
        self.ksi=np.array(ksi)
        self.omegamax=max(self.omega)
    
    @property
    def cp_c0(self):
        rez=[]
        for i in range(len(self.omega)):
            if self.ksi[i]==0:
                rez.append(1)
            else:
                rez.append(self.omega[i]/self.ksi[i]/self.cb)
        return rez
    
    @property
    def r_l(self):
        return self.r0*self.ksi/2/np.pi
     
    def Ksi(self, o):
        if o<=self.omegamax:
            return self.ksi_omega(o)
        else:
            return self.ksi[-1]+(o-self.omega[-1])/self.cr
#            return self.ksi[-1]+(omega-omega[-1])/self.cr

    def save(self, fname):
        rez={}
        rez['r0']=self.r0
        rez['rho']=self.rho
        rez['E']=self.E
        rez['nu']=self.nu
        rez['omegamax']=self.omegamax
        rez['omega']=list(map(float, self.omega))
        rez['ksi']=list(map(float, self.ksi))
        json.dump(rez, open(fname, 'w'))
            
    def load(self, fname):
        rez=json.load(open(fname, 'r'))
        self.r0=rez['r0']
        self.rho=rez['rho']
        self.E=rez['E']
        self.nu=rez['nu']
        self.omegamax=rez['omegamax']
        self.omega=np.array(rez['omega'])
        self.ksi=np.array(rez['ksi'])
        self.l=self.E*self.nu/(1+self.nu)/(1-2*self.nu)
        self.mu=self.E/2./(1+self.nu)
        self.cb=np.sqrt(self.E/self.rho)
        self.ct=np.sqrt(self.mu/self.rho)
        self.cl=np.sqrt((self.l+2*self.mu)/self.rho)
        self.cr=(0.862+1.14*self.nu)/(1+self.nu)*self.ct        
        self.ksi_omega=lambda x: np.interp(x, self.omega, self.ksi)

    def disp_shift_wave(self, t, y, dz, backshift=False):
        fy=np.fft.rfft(y)
        n=len(fy)
        dz=-dz
        omega=2*np.pi*np.arange(n)/max(t)
        fy2=[]
        for i in range(n):
            fy2.append(np.exp(complex(0,self.Ksi(omega[i])*dz))*fy[i])
        rez = np.fft.irfft(fy2).tolist()
        if backshift:    
            dt=dz/self.cb
            rez = self.simple_shift_wave(t,rez,dt)
        rez = rez+[rez[-1]]*(len(y)-len(rez))
        return rez
    
    def simple_shift_wave(self, t, y, tshift):
        dt = t[-1]-t[0]
        tshift = tshift-np.floor(tshift/dt)*dt
        dt=t[1]-t[0]
        n=int(np.round(tshift/dt))
        if isinstance(y,np.ndarray): y=y.tolist()
        return y[-n:]+y[:-n]

if __name__=='__main__':
    import matplotlib.pylab as plt
    c=barWithDispersion(r0=30e-3)
    c.solve()
    c.save('bar2.txt')
#    c1=barWithDispersion()
#    c1.load('bar1.txt')
    plt.plot(c.omega, c.ksi)
    plt.show()

