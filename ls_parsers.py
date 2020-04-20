# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:59:15 2015

@author: Konstantinov
change log
v2.0 - Change log started. Added rcforc and ncforc parsers for contact analisys.
v2.1 - Simplified parser for big ncforc files.
v2.2 - Parse Ncforce reading file by parts, not all at time (memory error issue).
v2.3 - parseNodout modification to rotational degree of freedom.
v2.4 - parseRcforc moments output.
v2.5 - spcforc parser.
v2.6 - parseNodout fix. SplitByNseq procedure.
v2.7 - parseNodout bug fixes + speed optimization + progress.
       parseRcforc progress.
"""
from pyparsing import *

def tofloat(x):
    rez=0
    try:
        rez=float(x)
    except:
        pass
    return rez
    
def toint(x):
    rez=0
    try:
        rez=int(x)
    except:
        pass
    return rez

def splitByNseq(s, n):
    while sum(n)+n[-1]<=len(s):
        n.append(n[-1])
    pos=[sum(n[:i+1]) for i in range(len(n))]
    pos.insert(0,0)
    rez=[]
    for i in range(len(pos)-1):
        rez.append(s[pos[i]:pos[i+1]].strip())
    return rez

def progress(cur, total, scaleLength=20, text='Progress='):
    prog=int(100.*cur/total)
    rez='\r{0}{1:4}% '.format(text, prog)
    numD=int(prog*scaleLength/100.)
    rez+='['+'#'*numD+'-'*(scaleLength-numD)+']'
    rez+=' Total={0:10}, current={1:10}'.format(total,cur)
    return rez

floatArrayWord=Word(nums+'+-.Ee').setParseAction(lambda x: [tofloat(xx) for xx in x])
floatWord=Word(nums+'+-.Ee').setParseAction(lambda x: tofloat(x[0]))
intWord=Word(nums).setParseAction(lambda x: toint(x[0]))

def parseElout(fname='elout'):
    '''
    Procedure to parse LS-DYNA elout file.
    Return t    - list of times
           stresses - dic of lists for all elements (element numbers - keys)
           stresses[elnum] - list of 7 lists:
                              0 - sxx
                              1 - syy
                              2 - szz
                              3 - sxy
                              4 - syz
                              5 - szx
                              6 - plastic strain
           if strains in elout, then also return strains={}
           strains[elnum]=[exx=[], eyy=[], ...]
    '''
    header=Suppress('e l e m e n t   s t r e s s   c a l c u l a t i o n s   f o r   t i m e  s t e p')
    header+=Suppress(Word(nums)+'( at time ')+floatWord+Suppress(')')
    header+=Suppress('''element  materl(global)
 ipt-shl  stress       sig-xx      sig-yy      sig-zz      sig-xy      sig-yz      sig-zx       plastic
           state                                                                                 strain ''')
    elNum=intWord+Suppress(restOfLine)
    data=Suppress(SkipTo('stic'))+Suppress('stic')+Group(floatArrayWord*7)
    strainHeader=Suppress('strains (global)      eps-xx      eps-yy      eps-zz      eps-xy      eps-yz      eps-zx')
    strainElNum=intWord+Suppress(restOfLine)
    strainLine1=Suppress('lower ipt')+floatArrayWord*6
    strainLine2=Suppress('upper ipt')+floatArrayWord*6
    strainPattern=strainHeader+OneOrMore(Group(strainElNum+Group(strainLine1+strainLine2)))
    pattern=header+Group(OneOrMore(Group(elNum+data)))+ZeroOrMore(Group(strainPattern))
    rez=pattern.searchString(open(fname).read())
    t=[]
    stresses={}
    strains={}
    for s in rez[0][1]:
        stresses[s[0]]=[[],[],[],[],[],[],[]]
    for r in rez:
        t.append(r[0])
        for s in r[1]:
            for j in range(7):
                stresses[s[0]][j].append(s[1][j])
    if len(rez[0])==3:
        for e in rez[0][2]:
            strains[e[0]]=[[],[],[],[],[],[]]
        for r in rez:
            for e in r[2]:
                for j in range(6):
                    strains[e[0]][j].append(0.5*(e[1][j]+e[1][j+6]))
    return t, stresses, strains
    
def parseNodfor(fname='nodfor'):
    '''
    Procedure to parse LS-DYNA nodfor.
    Return t=[] - times
           forces[ndgroup][ndnum] - [fx=[], fy=[], fz=[]]
           ndgroup - nodal group
           ndnum   - node number or 'total' for total force for this group.
    '''
    timeHeader=Suppress('n o d a l   f o r c e   g r o u p    o u t p u t  t=')+floatWord
    groupHeader=Suppress('nodal group output number')+intWord
    indNodes=Group(Suppress('nd#')+intWord+Group(Suppress('xforce=')+floatWord+Suppress('yforce=')\
            +floatWord+Suppress('zforce=')+floatWord+Suppress(restOfLine)))
    totalForces=Group(Suppress('xtotal=')+floatWord+Suppress('ytotal=')+floatWord+Suppress('ztotal=')+floatWord+Suppress(restOfLine))
    groupForce=Group(groupHeader+Group(OneOrMore(indNodes)+totalForces))
    pattern=timeHeader+Group(OneOrMore(groupForce))
    rez=pattern.searchString(open(fname).read())
    t=[]
    forces={}
    for gr in rez[0][1]:
        forces[gr[0]]={}
        for nd in gr[1][:-1]:
            forces[gr[0]][nd[0]]=[[],[],[]]
        forces[gr[0]]['total']=[[],[],[]]
    for r in rez:
        t.append(r[0])
        for gr in r[1]:
            for nd in gr[1][:-1]:
                for j in range(3):
                    forces[gr[0]][nd[0]][j].append(nd[1][j])
            for j in range(3):
                forces[gr[0]]['total'][j].append(gr[1][-1][j])
            
    return t, forces
    
def parseMatsum2(fname='matsum'):
    '''
    ls-dyna matsum parser.
    t=[] - times
    vels={}
    vels[mat]=[vx=[], vy=[], vz=[]]
    mat - material
    '''
    timeHeader=Suppress('time =')+floatWord
    matHeader=Suppress('mat.#=')+intWord+Suppress(restOfLine)
    momLine=Suppress('x-mom='+restOfLine)
    velLine=Suppress('x-rbv=')+floatWord+Suppress('y-rbv=')+floatWord+Suppress('z-rbv=')+floatWord
    lastLine=Suppress('hgeng='+restOfLine)
    dataBlock=Group(matHeader+momLine+Group(velLine)+lastLine)
    pattern=timeHeader+Group(OneOrMore(dataBlock))
    rez=pattern.searchString(open(fname).read())
    t=[]
    vels={}
    for d in rez[0][1]:
        vels[d[0]]=[[],[],[]]
    for r in rez:
        t.append(r[0])
        for d in r[1]:
            for j in range(3):
                vels[d[0]][j].append(d[1][j])
    return t, vels

def parseDbfsi(fname='dbfsi'):
    '''
    Parse dbfsi file.
    return t=[] - times
    f={}
    f[surfaceNum]=[p=[], fx=[], fy=[], fx=[]]
    '''
    timeLine=Suppress('time=')+floatWord
    dataLine=Group(intWord+Group(floatArrayWord*4)+Suppress(restOfLine))
    skippedLine=Suppress(6*Word(nums+'+-.eE'))
    pattern=timeLine+Group(OneOrMore(dataLine+skippedLine))
    rez=pattern.searchString(open(fname).read())
    t=[]
    forces={}
    for d in rez[0][1]:
        forces[d[0]]=[[],[],[],[]]
    for r in rez:
        t.append(r[0])
        for d in r[1]:
            for j in range(4):
                forces[d[0]][j].append(d[1][j])
    return t, forces

def parseNodout2(fname='nodout'):
    '''
    Parse nodout file.
    Return t=[] - times.
    nodouts={}
    nodouts[node]=[Ux=[], Uy=[], Uz=[], Vx=[],..., ax=[],..., x=[],...z=[]]
    '''
    timeLine=Suppress('at time')+floatWord+Suppress(restOfLine)
    displLine=Suppress('nodal point  x-disp     y-disp      z-disp'+restOfLine)
    rotLine=Suppress('nodal point  x-rot      y-rot       z-rot'+restOfLine)
    dataLine1=Group(intWord+Group(floatArrayWord*12))
    dataLine2=Group(intWord+Group(floatArrayWord*9))
    pattern=timeLine+skippedLine+Group(OneOrMore(dataLine))
    rez=pattern.searchString(open(fname).read())
    t=[]
    nodouts={}
    for d in rez[0][1]:
        nodouts[d[0]]=[[] for i in range(12)]
    for r in rez:
            t.append(r[0])
            for d in r[1]:
                for j in range(12):
                    nodouts[d[0]][j].append(d[1][j])
    return t, nodouts

def parseNodout(fname='nodout'):
    '''
    Parse nodout file.
    Return t=[] - times.
    nodouts={}
    nodouts[node]=[Ux=[], Uy=[], Uz=[], Vx=[],..., ax=[],..., x=[],...z=[], rotX=[], ...]
    '''
    print('Reading data from ', fname)
    data=open(fname).readlines()
#    while "\n" in data:
#        data.remove("\n")
    print('Done...')
    t=[]
    rez={}
    idxs=[]
    for i, l in enumerate(data):
        if "n o d a l   p r i n t   o u t   f o r   t i m e  s t e p" in l:
            idxs.append(i)
    idxs.append(len(data)+1)
    lastStep=''
    currentStep=''
    def parseOneRecord(data1):
        ll=data1[0].split()
        t1=float(ll[-2])
        currentStep=ll[-6]
        rez1={}
        for l in data1[2:]:
            if l=='\n': continue
            ss=splitByNseq(l, [9,13,12])
            nn=toint(ss[0])
            rez1[nn]=[tofloat(sss) for sss in ss[1:]]
        return t1, rez1, currentStep
    total=len(idxs)-1
    print('Parsind data...')
    for i in range(total):
        tt, rr, currentStep = parseOneRecord(data[idxs[i]:idxs[i+1]])
        rot=True
        if currentStep!=lastStep:
            t.append(tt)
            rot=False
        if not rez:
            for k in list(rr.keys()):
                rez[k]=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        j=12 if rot else 0
        try:
            for k, v in list(rr.items()):
                for ii, vv in enumerate(v):
                    rez[k][ii+j].append(vv)
        except:
            pass
        lastStep=currentStep
        print(progress(i,total-1), end='')
    print('Done...')
    return t, rez

def parseD3HSP(fname='d3hsp'):
    """
    Parse d3hsp file for mass properties of parts.\n
    Return dict rez[partNumber]=massPropertiesData.\n
    massPropertiesData.keys:\n
         "mass"                 - part mass\n
         "cm"                   - mass center\n
         "inertiaTensor"        - inertia tensor\n
         "principalInertias"    - principal inertias [I11, I22, I33]\n
         "principalDirections"  - principal directions\n
    """
    def parseOnePart(data1):
        rez={}
        partN=int(data1[0].split("#")[-1])
        rez["mass"]=float(data1[1].split("=")[-1])
        rez["cm"]=[float(data1[i].split("=")[-1]) for i in (2, 3, 4)]
        rez["inertiaTensor"]=[[float(s) for s in data1[i][10:].split()] for i in (7,8,9)]
        rez["principalInertias"]=[float(data1[i].split("=")[-1]) for i in (12, 13, 14)]
        rez["principalDirections"]=[[float(s) for s in data1[i][10:].split()] for i in (17,18,19)]
        return partN, rez
    data=open(fname).readlines()
    rez={}
    for i, l in enumerate(data):
        if "m a s s   p r o p e r t i e s   o f   p a r t" in l:
            n, r = parseOnePart(data[i:i+20])
            rez[n]=r
    return rez

def parseMatsum(fname='matsum'):
    '''
    ls-dyna matsum parser.
    t=[] - times
    vels={}
    vels[mat]=[vx=[], vy=[], vz=[]]
    mat - material
    '''
    data=open(fname).readlines()
    while "\n" in data:
        data.remove("\n")
    t=[]
    rez={}
    idxs=[]
    for i, l in enumerate(data):
        if "time =" in l:
            idxs.append(i)
    idxs.append(len(data)+1)
    def parseOneRecord(data1):
        t1=float(data1[0].split("=")[-1])
        rez1={}
        for i, l in enumerate(data1):
            if "mat.#=" in l:
                nn=toint(l.split()[1])
                rez1[nn]=[tofloat(ss) for ss in data1[i+2].split()[1::2]]
        return t1, rez1
    for i in range(len(idxs)-1):
        tt, rr = parseOneRecord(data[idxs[i]:idxs[i+1]])
        t.append(tt)
        if not rez:
            for k in list(rr.keys()):
                rez[k]=[[],[],[]]
        for k, v in list(rr.items()):
            for i, vv in enumerate(v):
                rez[k][i].append(vv)
    return t, rez
    
def parseTrhist(fname='trhist'):
    '''
    ls-dyna trhist parser.
    return:
        times - array of times
        rez - array of particles data, i.e. rez[0] - data
        for particle 1 and so on.
        rez[i]=[x,y,z,vx,vy,vz,sx,sy,sz,sxy,syz,szx]
    trhist header:
       x       y       z      vx      vy      vz
      sx      sy      sz     sxy     syz     szx
     efp     rho    rvol  active
    '''
    f=open(fname, 'r')
    data=f.readlines()
    nParticles=int(data[1].split()[0])
    times=[]
    rez=[]
    for i in range(nParticles):
        rez.append([[],[],[],[],[],[],[],[],[],[],[],[]])
    for i in range(5, len(data), 1+3*nParticles):
        times.append(tofloat(data[i]))
        for j in range(nParticles):
            vrow=data[i+1+j*3].split()
            srow=data[i+2+j*3].split()
            for k in range(6):
                rez[j][k].append(tofloat(vrow[k]))
                rez[j][k+6].append(tofloat(srow[k]))
    f.close()
    return times, rez

def parseRcforc(fname='rcforc'):
    '''
    ls-dyna rcforc parser.
    return:
        t - array of times.
        rez={'m1': [fx, fy, fz], 's1': [fx, fy, fz], ...}
        m1 - master surface 1
        s1 - slave surface 1
        m2 - master surface 2
        ...
        example fx for master surface 1: rez['m1'][0] 
    '''
    print('Reading data from ', fname)
    f=open(fname, 'r')
    data=f.readlines()
    f.close()
    print('Done...')
    t=[]
    lastTime=None
    def parceOne(l):
        ll=l.split()
        tp=ll[0]
        n=ll[1]
        tm=ll[3]
        f=list(map(tofloat, [ll[5], ll[7], ll[9], ll[13], ll[15], ll[17]]))
        return tp[0]+n, tm, f
    i=0    
    while not '{END LEGEND}' in data[i]:
        i+=1
    i+=2
    rez={}
    print('Parsing data...')
    total=len(data[i:])
    for ii, d in enumerate(data[i:]):
        tp, tm, f = parceOne(d)
        if tm!=lastTime:
            lastTime=tm
            t.append(tofloat(tm))
        if tp not in rez:
            rez[tp]=[[],[],[],[],[],[]]
        for j in range(6):
            rez[tp][j].append(f[j])
        print(progress(ii, total), end='')
    print('Done...')
    return t, rez

def parseNcforc2(fname='ncforc'):
    '''
    ls-dyna rcforc parser.
    return:
        t - array of times.
        rez={'m-1-10': [fx, fy, fz, p, x, y, z], 'm-1-11': [fx, fy, fz, p, x, y, z], ...}
        m-1-10 - master surface 1 node 10
        ...
        example p for slave surface 2 node 5: rez['s-2-5'][3] 
    '''
    f=open(fname, 'r')
    data=f.readlines()
    f.close()
    t=[]
    lastTime=None
    def parceOne(lns):
        tm=lns[0].split('t=')[1].split(')')[0]
        n=lns[0].split()[4]
        tp=lns[0].split()[5][0]
        rez={}
        for i in range(4, len(lns),2):
            if lns[i]=='\n':
                break
            ll1=lns[i].split()
            ll2=lns[i+1].split()
            rez['-'.join([tp, n, ll1[0]])]=list(map(tofloat, ll1[1:]+ll2))
        return tm, rez
    idxs=[]
    for i, l in enumerate(data):
        if 'forces' in l:
            idxs.append(i)
    idxs.append(len(data))
    rez={}
    for i in range(len(idxs)-1):
        tm, d = parceOne(data[idxs[i]:idxs[i+1]])
        if tm!=lastTime:
            lastTime=tm
            t.append(float(tm))
        for k, f in d.items():
            if k not in rez:
                rez[k]=[[],[],[],[],[],[],[]]
            for j in range(7):
                rez[k][j].append(f[j])
    return t, rez

def parseBigNcforc(fname='ncforc', idx=-1):
    '''
   ls-dyna rcforc parser.
    return:
        t - array of times.
        rez={'m-1-10': [fx, fy, fz, p, x, y, z], 'm-1-11': [fx, fy, fz, p, x, y, z], ...}
        m-1-10 - master surface 1 node 10
        ...
        example p for slave surface 2 node 5: rez['s-2-5'][3]
    '''
    f=open(fname, 'r')
    t=[]
    rez={}
    lastTime=None
    dataStart=False
    for l in f:
        if 'forces' in l:
            dataStart=True
            tm=l.split('t=')[1].split(')')[0]
            n=l.split()[4]
            tp=l.split()[5][0]
            if lastTime!=tm:
               lastTime=tm
               t.append(float(tm))
        if l[0:10].lstrip().isdigit() and dataStart:
            ll=l.split()
            k='-'.join([tp, n, ll[0]])
            if k not in rez:
                rez[k]=[]
            rez[k].append(tofloat(ll[idx]))
    f.close()
    return t, rez

def parseNcforc(fname='ncforc'):
    '''
    ls-dyna rcforc parser.
    return:
        t - array of times.
        rez={'m-1-10': [fx, fy, fz, p, x, y, z], 'm-1-11': [fx, fy, fz, p, x, y, z], ...}
        m-1-10 - master surface 1 node 10
        ...
        example p for slave surface 2 node 5: rez['s-2-5'][3] 
    '''
    t=[]
    lastTime=None
    dataBlock=[]
    
    def parceOne(lns):
        if not lns:
            return None, None
        if not 'forces' in lns[0]:
            return None, None
        tm=lns[0].split('t=')[1].split(')')[0]
        n=lns[0].split()[4]
        tp=lns[0].split()[5][0]
        rez={}
        for i in range(4, len(lns),2):
            if lns[i]=='\n':
                break
            ll1=lns[i].split()
            ll2=lns[i+1].split()
            rez['-'.join([tp, n, ll1[0]])]=list(map(tofloat, ll1[1:]+ll2))
        return tm, rez

    rez={}
    fl=open(fname, 'r')
    for l in fl:
        if 'forces' in l:
            tm, d = parceOne(dataBlock)
            dataBlock=[l]
            if not tm:
                continue
            if tm!=lastTime:
                lastTime=tm
                t.append(float(tm))
            for k, f in d.items():
                if k not in rez:
                    rez[k]=[[],[],[],[],[],[],[]]
                for j in range(7):
                    rez[k][j].append(f[j])
            continue
        dataBlock.append(l)
    fl.close()    
    return t, rez

def parseSpcforc(fname='spcforc'):
    '''
    ls-dyna spcforc parser
    return t, rez
    t - times
    rez={}
    rez[0]=[fx, fy, fz] - rezultant forces.
    rez[node]=[fx, fy, fz, mx, my, mz] - forces at node.
    '''
    t=[]
    rez={}
    rez[0]=[[],[],[]]
    def parseOne(lns):
        t=tofloat(lns[0].split()[-1])
        f={}
        m={}
        fr=map(tofloat, lns[-1].split()[-3:])
        for l in lns[1:-1]:
            ll=l.split()
            if 'forces' in l:
                f[int(ll[1])]=map(tofloat, ll[-5:-2])
            if 'moments' in l:
                m[int(ll[1])]=map(tofloat, ll[-5:-2])
        return t, f, m, fr
    data=open(fname, 'r').readlines()
    idxs=[]
    for i, l in enumerate(data):
        if 'output at time' in l:
            idxs.append(i)
    idxs.append(len(data))
    for i in range(len(idxs)-1):
        tt, ff, mm, ffr = parseOne(data[idxs[i]:idxs[i+1]])
        t.append(tt)
        for i in range(3):
            rez[0][i].append(ffr[i])
        for nd in ff.keys():
            if not rez.has_key(nd):
                rez[nd]=[[],[],[],[],[],[]]
            for i in range(3):
                rez[nd][i].append(ff[nd][i])
                if mm.has_key(nd):
                    rez[nd][i+3].append(mm[nd][i])
    return t, rez


if __name__=="__main__":
#    import time
#    t0=time.time()
#    t, n =parseMatsum("../task/matsum")
#    t1=time.time()
#    print(t1-t0)
#    t0=time.time()
#    t, n =parseMatsum2("../task/matsum")
#    t1=time.time()
#    print t1-t0
    t, d=parseNcforc()
