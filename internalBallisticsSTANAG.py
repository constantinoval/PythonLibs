# coding: cp1251
# v 2.0 - added output of Pproj and Pbreech
# v 2.1 - fixed bug in saving data for post shot period
# v 2.2 - added saving of mean density of propelant gases
# v 2.3 - implicitSwitchTime fix
# v 3.0 - added water in barrel
# v 3.1 - added air resistant Pair output.
# v 3.2 - result output fix
# v 3.3 - Stop criteria fix
# v 3.4 os.path.join fix
# 01.12.2018 - STANAG_solver object
# 14.02.2019 - fixed combustion model for tube and multichannel propellant grains

from __future__ import division, print_function, generators, with_statement
import json
import os
import sys
from math import pi
from functools import partial
from scipy.optimize import fsolve
from numpy import interp
import argparse

class configReader:
    def __init__(self, configFile):
        self.configFile=configFile
        self.conf=json.loads(open(configFile).read())
        self.check()
        if self.errorCode:
            print('Error in task configuration.\n'+self.errorMessage)
            exit()
    def check(self):
        self.errorCode=0
        self.errorMessage=""
    def getFEpath(self):
        return self.conf.get('FEpath', '')
    def getBarrelArea(self):
        return self.conf["barrelArea"]
    def getBulletMass(self):
        return self.conf["bulletMass"]
    def getBulletIy(self):
        return self.conf["bulletIy"]
    def getBarrelLength(self):
        return self.conf["barrelLength"]
    def getChamberVolume(self):
        return self.conf["chamberVolume"]
    def getRiffleStep(self):
        return self.conf["riffleStep"]
    def getBarrelT0(self):
        return self.conf.get("barrelT0",293)
    def getPowderMass(self):
        return self.conf["powderMass"]
    def getPowderForce(self):
        return self.conf["powderForce"]
    def getPowderBirnT(self):
        return self.conf["powderBirnT"]
    def getSpecificCovolum(self):
        return self.conf.get("specificCovolum",0)                      
    def getGama(self):
        return self.conf.get("gama",0)               
    def getPowderDens(self):
        return self.conf["powderDens"]
    def getW(self):
        return self.conf["w"]            
    def getBeta(self):
        return self.conf["beta"]              
    def getMinPressure(self):
        return self.conf.get("minPressure", 1e5)
    def getTimeStep(self):
        return self.conf["timeStep"]
    def getPowderGrainType(self):
        return self.conf["powderForm"]["type"]
    def getPowderGrainParams(self):
        return self.conf["powderForm"]["params"]
    def getTotalTime(self):
        return self.conf.get('totalTime', 0)
    def getFormLSpressures(self):
        return self.conf.get('formLSpressures', 0)
    def getWaterStart(self):
        return self.conf.get('waterStart', 0)
    def getWaterEnd(self):
        return self.conf.get('waterEnd', self.getBarrelLength())
    def getWaterDens(self):
        return self.conf.get('waterDens', 0)
    def getRezPoints(self):
        rp = self.conf.get('rezPoints', 0)
        if rp==0:
           return 0
        if type(rp)==int:
            rp = [rp]
        if len(rp)==1:
            rp.append(0)
        return rp
    def getAtmPressure(self):
        return self.conf.get('atmPressure', 1e5)

def tubeInitialGeom(args, h=0):
    D=args["D"]-2*h
    L=args["L"]
    D0=args["D0"]+2*h
    S0=pi*D*L+pi*D0*L+0.5*pi*(D**2-D0**2)
    V0=pi*(D**2/4.-D0**2/4.)*L
    return S0, V0   
# def tubeFormFunction(z, args):
#     D0=args["D0"]
#     L=args["L"]
#     return (1-0.57*z)**0.5#D0-L+(D0**2-4*D0*L*z+2*D0*L+L**2)**0.5/2./D0
def tubeCheckConsistency(args):
    err=0
    errorCode="Powder grain parameters error: "
    if args["D0"]<=0:
        err=1
        errorCode+="D0<=0. "
    if args["D"]<=0:
        err=1
        errorCode+="D<=0. "
    if args["D0"]>=args["D"]:
        err=1
        errorCode+="D<=D0. "
    return err, errorCode
    
def sphereInitialGeom(args):
    D=args["D"]
    S0=pi*D**2
    V0=4./3.*pi*(D/2.)**3
    return S0, V0
def sphereFormFunction(z, args):
    return (1-z)**0.6666666666666666    
def sphereCheckConsistency(args):
    err=0
    errorCode="Powder grain parameters error: "
    if args["D"]<=0:
        err=1
        errorCode+="D<=0. "
    return err, errorCode

def cylInitialGeom(args):
    D=args["D"]
    L=args["L"]
    S0=pi*L*D
    V0=0.25*pi*D**2*L
    return S0, V0
def cylFormFunction(z, args):
    return (1-z)**0.5  
def cylCheckConsistency(args):
    err=0
    errorCode="Powder grain parameters error: "
    if args["D"]<=0:
        err=1
        errorCode+="D<=0. "
    if args["L"]<=0:
        err=1
        errorCode+="L<=0. "
    return err, errorCode
        

def multiChannelGeom(args, h=0):
    D=args["D"]-2*h
    L=args["L"]
    D0=args["D0"]+2*h
    N=args["N"]
    S=pi*L*(D+N*D0)+0.5*pi*(D**2-N*D0**2)
    V=0.25*pi*L*(D**2-N*D0**2)
    return S, V
#def multiChannelFormFunction(z, args):
#    S0, V0 = multiChannelGeom(args)
#    fi=lambda u: multiChannelGeom(args,u)[0]/S0
#    zz_z=lambda u: 1-multiChannelGeom(args,u)[1]/V0-z
#   uu=fsolve(zz_z, 0)[0]#mpmath.findroot(zz_z,0)
#    return fi(uu)
def multiChannelCheckConsistency(args):
    err=0
    errorCode="Powder grain parameters error: "
    NN=[20,19,17,16,15,14,13,12,11,10,9,8,7,5,4,3,2,1]
    d_d0=[5.122, 4.863, 4.792, 4.615, 4.521, 4.328, 4.236, 4.029, 3.923,
          3.813, 3.613, 3.304, 3, 2.701, 2.414, 2.154, 2, 1]
    D=args["D"]
    D0=args["D0"]
    if D<=0:
        err=1
        errorCode+="D<=0. "
    if D0<=0:
        err=1
        errorCode+="D0<=0. "
    if D<=D0:
        err=1
        errorCode+="D<=D0. "
    crit=D/D0
    N=args["N"]
    for nnn, ddd in zip(NN, d_d0):
        if crit>ddd:
            break
    if N>nnn:
        err=1
        errorCode+="too many holes. Max number=%d. " % (nnn,)
    return err, errorCode

def diskGeom(args, u=0):
    D=args["D"]
    w=args["w"]
    S=pi*(D-u)*(0.5*(D-u)+(w-u))
    V=0.25*pi*(D-u)**2*(w-u)
    return S, V
def diskFormFunction(z, args):
    S0, V0 = diskGeom(args)
    fi=lambda u: diskGeom(args,u)[0]/S0
    zz_z=lambda u: 1-diskGeom(args,u)[1]/V0-z
    uu=fsolve(zz_z, 0)[0]
    return fi(uu)
def diskCheckConsistency(args):
    err=0
    errorCode="Powder grain parameters error: "
    if args["D"]<=0:
        err=1
        errorCode+="D<=0. "
    if args["w"]<=0:
        err=1
        errorCode+="w<=0. "
    return err, errorCode
    
def solidPrismGeom(args, u=0):
    D=args["D"]
    w=args["w"]
    L=args["L"]
    S=2*((L-u)*(D-u)+(L-u)*(w-u)+(D-u)*(w-u))
    V=(L-u)*(D-u)*(w-u)
    return S, V
def solidPrismFormFunction(zz, args):
    S0, V0 = solidPrismGeom(args)
    fi=lambda u: solidPrismGeom(args,u)[0]/S0
    zz_z=lambda u: 1-solidPrismGeom(args,u)[1]/V0-zz
    uu=fsolve(zz_z, z[-1])[0]
    return fi(uu)
def solidPrismCheckConsistency(args):
    err=0
    errorCode="Powder grain parameters error: "
    if args["D"]<=0:
        err=1
        errorCode+="D<=0. "
    if args["w"]<=0:
        err=1
        errorCode+="w<=0. "
    if args["L"]<=0:
        err=1
        errorCode+="L<=0. "
    return err, errorCode

powderGrains={"tube": [tubeInitialGeom, None, tubeCheckConsistency],
              "sphere": [sphereInitialGeom, sphereFormFunction, sphereCheckConsistency],
              "cyl": [cylInitialGeom, cylFormFunction, cylCheckConsistency],
              "multiChannel": [multiChannelGeom, None, multiChannelCheckConsistency],
              "disk": [diskGeom, diskFormFunction, diskCheckConsistency],
              "solidPrism": [solidPrismGeom, solidPrismFormFunction, solidPrismCheckConsistency]}

class STANAG_solver:
    def __init__(self, cfg_file):
        self.cfg = configReader(cfg_file)
    def solve(self, dt=None):
        #parameters of weapon
        if dt==None:
            dt=self.cfg.getTimeStep()
        A=self.cfg.getBarrelArea()
        m=self.cfg.getBulletMass()
        Iy=self.cfg.getBulletIy()
        xm=self.cfg.getBarrelLength()
        Vchamber=self.cfg.getChamberVolume()
        htw=self.cfg.getRiffleStep()
        Tw=self.cfg.getBarrelT0()
        #parametets of Powder
        C0=self.cfg.getPowderMass()
        F=self.cfg.getPowderForce()
        T0=self.cfg.getPowderBirnT()
        b=self.cfg.getSpecificCovolum()
        gama=self.cfg.getGama()
        rhoprop=self.cfg.getPowderDens()
        w=self.cfg.getW()
        beta=self.cfg.getBeta()

        pmin=self.cfg.getMinPressure()
        grainType=self.cfg.getPowderGrainType()
        grainPars=self.cfg.getPowderGrainParams()
        if powderGrains[grainType][2]:
            err, errorCode=powderGrains[grainType][2](grainPars)
            if err:
                print(errorCode)
                exit()
            
        S0, V0 = powderGrains[grainType][0](grainPars)
        fi=partial(powderGrains[grainType][1], args=grainPars) if powderGrains[grainType][1] else None
        viv=w*S0/V0
        #------------
        x0=Vchamber/A
        patm=self.cfg.getAtmPressure()
        p=[patm]
        pproj=[(p[-1]+C0/3./m*patm)/(1+C0/3./m)]
        pbreech=[pproj[-1]*(1+C0/2./m)-C0/2./m*patm]
        z=[0.0]
        T=[T0]
        xp=[x0]
        vp=[0]
        t=[0]
        omegap=[0]
        rhog=[0]
        Eair=0
        Pair=[patm]
        waterStart=self.cfg.getWaterStart()
        waterLength=self.cfg.getWaterEnd()-self.cfg.getWaterStart()
        waterDens=self.cfg.getWaterDens()
        waterStatus=0
        hhh = 0
        while (xp[-1]-x0)<=xm:
            addedMass=0
            x1=xp[-1]-x0
            if waterDens and x1>=waterStart:
                x2=min(x1+waterLength, xm)
                addedMass=waterDens*(x2-x1)*A
                if waterStatus==0:
                    vp[-1]*=m/(m+addedMass)
                waterStatus=1
            M=m+addedMass
            if z[-1]>=1:
                zz=1
            else:
                if not fi:
                    hhh+=w*p[-1]**beta*dt
                    S, V = powderGrains[grainType][0](grainPars, h=hhh)
                    zz=z[-1]+w*p[-1]**beta*dt*S/V0
                else:
                    zz=z[-1]+viv*fi(z[-1])*p[-1]**beta*dt
            z.append(zz)
            pair=patm*(1+0.84*(vp[-1]/340.)**2)
            Ep=0.5*M*vp[-1]**2
            Erot=0.5*Iy*omegap[-1]**2
            Egas=1./6.*C0*vp[-1]**2
            Et = 0
            if len(xp)>=2: Eair+=pair*A*(xp[-1]-xp[-2])
            TT=T0*(1-(Ep+Egas+Eair+Erot+Et)/C0/z[-1]/F*(gama-1))
            T.append(TT)
            pp=patm+C0*zz*F*T[-1]/T0/(Vchamber+A*(xp[-1]-x0)-C0*(1-z[-1])/rhoprop-C0*z[-1]*b)
            rhog.append(C0*z[-1]/(Vchamber+A*(xp[-1]-x0)-C0*(1-z[-1])/rhoprop-C0*z[-1]*b))
            p.append(pp)
            ppproj=(p[-1]+C0/3./M*pair)/(1+C0/3./M)
            ppbreech=ppproj*(1+C0/2.0/M)-C0/2.0/M*pair
            ap=(ppproj-pair)*A/M
            vvp=vp[-1]+ap*dt
            xxp=xp[-1]+vvp*dt
            pproj.append(ppproj)
            pbreech.append(ppbreech)
            xp.append(xxp)
            vp.append(vvp)
            omegap.append(vvp/htw)
            t.append(t[-1]+dt)
            Pair.append(pair)

        for i, pp in enumerate(p):
            if pp>=pmin:
                break
        p=p[i:]
        Pair=Pair[i:]
        xp=xp[i:]
        vp=vp[i:]
        omegap=omegap[i:]
        z=z[i:]
        T=T[i:]
        pproj=pproj[i:]
        pbreech=pbreech[i:]
        rhog=rhog[i:]
        t=t[i:]
        t0=t[0]
        t=[tt-t0 for tt in t]
        tswitch=t[-1]
        totTime=self.cfg.getTotalTime()
        if totTime-t[-1]-t0>0:
            tm=t[-1]
            pm=p[-1]
            Tm=T[-1]
            tmax=totTime-t[-1]
            rhom=C0*z[-1]/(Vchamber+A*(xp[-1]-x0)-C0*(1-z[-1])/rhoprop-C0*z[-1]*b)
            v_m=vp[-1]/xp[-1]
            rho=lambda t_: rhom/(1+v_m*t_)
            Pafter=lambda t_: pm*((1.0/rhom-b)/(1/rho(t_)-b))**gama
            Tafter=lambda t_: Tm*((1.0/rhom-b)/(1/rho(t_)-b))**(gama-1)
            tt=dt
            while tt<=tmax:
                t.append(tm+tt)
                rhog.append(rhom)
                xp.append(xp[-1]+vp[-1]*dt)
                vp.append(vp[-1])
                omegap.append(omegap[-1])
                p.append(Pafter(tt))
                pproj.append((p[-1]+C0/3./m*patm)/(1+C0/3./m))
                pbreech.append(pproj[-1]*(1+C0/2./m)-C0/2./m*patm)
                Pair.append(Pair[-1])
                T.append(Tafter(tt))
                z.append(z[-1])
                tt+=dt
        # end solve
        self.t = t
        self.xp = xp
        self.vp = vp
        self.omegap = omegap
        self.p = p
        self.T = T
        self.z = z
        self.pproj = pproj
        self.pbreech = pbreech
        self.rhog = rhog
        self.tswitch = tswitch
        self.totTime = totTime
        self.x0 = x0

    def remap(self, rezPoints=None):
        """Перенос результата на сетку, определяемую параметром rezPoints
        
        Keyword Arguments:
            rezPoints {[N1, N2]} -- если None, то сетка берется из конфигурационного файла.
                                    N1 - число точек на основной участок (активный выстрел),
                                    N2 - число точек на период последействия
        """

        if rezPoints==None:
            rezPoints=self.cfg.getRezPoints()
        if rezPoints==0: return
        newT=map(lambda x: x*self.tswitch/rezPoints[0], range(rezPoints[0]+1))
        if rezPoints[1] and self.cfg.getTotalTime():
            newT+=map(lambda x: x*(self.totTime-self.tswitch)/rezPoints[1]+self.tswitch, range(1, rezPoints[1]+1))
        self.xp=interp(newT, self.t, self.xp)
        self.vp=interp(newT, self.t, self.vp)
        self.omegap=interp(newT, self.t, self.omegap)
        self.p=interp(newT, self.t, self.p)
        self.T=interp(newT, self.t, self.T)
        self.z=interp(newT, self.t, self.z)
        self.pproj=interp(newT, self.t, self.pproj)
        self.pbreech=interp(newT, self.t, self.pbreech)
        self.rhog=interp(newT, self.t, self.rhog)
        self.t=newT
    def save_k_file(self, task_dir):
        task_dir = os.path.abspath(task_dir)
        if not os.path.exists(task_dir):
            os.mkdir(task_dir)
        kfilePath=os.path.join(task_dir,"pressures.k")
        try:
            mesherData=json.loads(open(os.path.join(task_dir,"mesherData.dat")).read())
            y0=mesherData["bulletY0"]
        except:
            y0=0
            mesherData={}
        with open(kfilePath, 'w') as fout:
            fout.write('*keyword\n')
            fout.write('*DEFINE_FUNCTION_TABULATED\n       100pressure_vs_time\npress\n')
            fout.write('$#                a1                  o1\n')
            for ttt, ppp in zip(self.t, self.p):
                fout.write('%20.10e%20.10e\n' % (ttt, ppp))
            fout.write('*DEFINE_FUNCTION_TABULATED\n       101position_vs_time\nposition\n')
            fout.write('$#                a1                  o1\n')
            for ttt, lll in zip(self.t, self.xp):
                fout.write('%20.10e%20.10e\n' % (ttt, lll+y0-self.x0))
            fout.write('''*DEFINE_FUNCTION
$#     fid
    102                                                                      
float brlpres(float t, float x, float y, float z,
float x0, float y0, float z0)
{
float sc;
sc=1.0;
if(y>position(t)) sc=0.0;
return sc*press(t);
}
*end
''')
        mesherData["bulletInitialVelocity"]=self.vp[0]
        mesherData["implicitSwitchTime"]=self.tswitch
        f=open(os.path.join(task_dir, "mesherData.dat"),'w')
        json.dump(mesherData, f)
        f.close()

    def save(self, fname):
        fname = os.path.abspath(fname)
        if not os.path.exists(os.path.split(fname)[0]):
            os.mkdir(os.path.split(fname)[0])
        with open(fname, 'w') as fout:
            fs="%s\t"*10+'\n'
            fout.write(fs % ("time", "displecement", "velocity", "omega",
                    "pressure", "temperature", "Z", "Pproj", "Pbreech", 'rhogas'))
            fs="%e\t"*10+'\n'
            for tt, xx, vv, om, pp, tmp, zz, ppr, pbr, rhom in zip(self.t,self.xp,self.vp,self.omegap,self.p,self.T,self.z,self.pproj,self.pbreech, self.rhog):
                fout.write(fs % (tt, xx-self.x0, vv, om, pp, tmp, zz, ppr, pbr, rhom))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', help='Path to configuration file.')
    parser.add_argument('-od', '--output_dir', dest= 'out_dir', help = 'Path to save results.')
    args = parser.parse_args()
    if not os.path.exists(args.cfg):
        print('Cant open configuration file.')
        sys.exit()
    solver = STANAG_solver(args.cfg)
    wd = args.out_dir if args.out_dir else solver.cfg.getFEpath()
    solver.solve()
    solver.remap()
    solver.save(os.path.join(wd, 'STANAG.rez'))
    if solver.cfg.getFormLSpressures():
        solver.save_k_file(wd)