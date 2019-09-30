# -- coding: cp1251
import numpy as np
import pyodbc
from scipy.optimize import minimize_scalar, fmin, differential_evolution
#matplotlib.use('PS')
#%%
def runningMeanFast(x, N):
    rez = np.roll(np.convolve(x, np.ones((N,))/N)[(N-1):], N//2)
    rez[:N//2] = rez[N//2]
    return rez

def toReal(diag, exp_type):
    sign=[-1,1][exp_type=='t']
    diag['der']=np.array(diag['det'])/(1+sign*np.array(diag['et']))
    diag['sr']=np.array(diag['st'])*(1+sign*np.array(diag['et']))
    diag['er']=sign*np.log(1+sign*np.array(diag['et']))

def meanDE(data):
    maxe=data['et'].max()
    N1=(data['de']*(1-data['et']/maxe)).argmax()
    N2=(data['st']*(1+data['et']/maxe)).argmax()
#    plt.plot(data[3]*(1-data[1]/maxe))
#    plt.show()
    if N2<N1:
        N1,N2=N2,N1
    N=N2-N1
    return data['de'][N1+N//5:N2-N//5].mean()#, N1+N//3, N2-N//3

def sync(p):
    def shiftTr(p, n=0, rez=0):
        n=int(n)
        i1=p[0]-p[1]
        if n==0:
            i2=p[2]
        if n<0:
            i2=np.array(p[2][-n:].tolist()+[p[2][-1]]*(-n))
        if n>0:
            i2=np.array([p[2][0]]*n+p[2][:-n].tolist())
        if rez==0:
            return ((i1-i2)**2).sum()
        else:
            return i2

    def shiftRef(p, dy=0):
        i1=p[0]-p[1]-dy
        i2=p[2]
        return ((i1-i2)**2).sum()
    pp=p.copy()
    N=len(p[0])
    amp=max(p[0])
    r=minimize_scalar(lambda n: shiftTr(p,n,0), bounds=(-N//3,N//3), method='bounded')
    pp[2]=shiftTr(p,int(r.x),1)
    r=minimize_scalar(lambda dy: shiftRef(pp,dy), bounds=(-amp/10,amp/10), method='bounded')
    pp[1]+=r.x
    return pp

#def syncPulses(p):
#    tmax=p[0][-1]
#    emax=max(p[1][0])
#    ref=lambda t: np.interp(t, p[0], p[1][1])
#    tr=lambda t: np.interp(t, p[0], p[1][2])
#    def residual(params):
#        dtref, dttr, dyref = params
#        return ((p[1][0]-(ref(p[0]-dtref)+dyref)-tr(p[0]-dttr))**2).sum()
#    r=fmin(residual,  (tmax/10, tmax/10, emax/10), ftol=1e-6, xtol=1e-8)
#    return p[0], [p[1][0], ref(p[0]-r[0])+r[2], tr(p[0]-r[1])]

def syncPulses(p):
    tmax=p[0][-1]
    emax=max(p[1][0])
    ref=lambda t: np.interp(t, p[0], p[1][1])
    tr=lambda t: np.interp(t, p[0], p[1][2])
    def residual(params):
        dtref, dttr, dyref = params
        return ((p[1][0]-(ref(p[0]-dtref)+dyref)-tr(p[0]-dttr))**2).sum()
    r=differential_evolution(residual, bounds=((-tmax/3, tmax/3), (-tmax/3, tmax/3), (-emax/10, emax/10)))
    return p[0], [p[1][0], ref(p[0]-r.x[0])+r.x[2], tr(p[0]-r.x[1])]

def correctE(e,s, Ecor=100000., window=800):
    k=[]
    a=[]
    N=len(e)
    window=max(2, N//10)
    step=max(1,window//100)
    for i in range(0,N//3,step):
        r=np.polyfit(e[i:i+window],s[i:i+window],1)
        k.append(r[0])
        a.append(r[1])
    NN=np.argmax(k)
    kk=k[NN]
    aa=a[NN]
    return e-s*(Ecor-kk)/Ecor/kk+aa/kk

def integrate(y, dx):
    rez=[0]
    for i in range(1,len(y)):
        rez.append(rez[-1]+0.5*(y[i]+y[i-1])*dx)
    return np.array(rez)      

def calcDiagram(dt, ein, eref, etr, cfg, isCorrE=False, Ecor=100000.):
    n=len(ein)
    et=[]
    st=[]
    det=[]
    for i in range(n):
        V1=cfg['c1']*(ein[i]+eref[i])
        V2=cfg['c2']*etr[i]
        det.append((V1-V2)/cfg['Lsp'])
        F=cfg['E2']*cfg['S2']*etr[i]
        st.append(F/cfg['Ssp'])
        et.append(integrate(det[i], dt))
    for i in range(1,n):
        et[i]+=et[i-1][-1]
    et=np.array(et).flatten()
    det=np.array(det).flatten()
    st=np.array(st).flatten()
    if isCorrE:
        et = correctE(et, st, Ecor, len(ein)//10)
    return et,st,det

def calcDiagram2(db, exp_code):
    experiment = db.getExperimentData(exp_code)
    b1 = db.getBarData(experiment.bars[0])
    b2 = db.getBarData(experiment.bars[1])
    t = experiment.pulses['t']
    return t, calcDiagram(dt=(t[1]-t[0]),
                          ein=[experiment.pulses['pulses'][0]],
                          eref=[experiment.pulses['pulses'][1]],
                          etr=[experiment.pulses['pulses'][2]],
                          cfg={'E2': b2.E, 'c1': b1.c, 'c2': b2.c, 'S2': b2.S,
                               'Ssp': experiment.d0**2/4.*3.14, 'Lsp': experiment.l0*1e-3}
                          )

def tofloat(s):
    if s==None:
        return 0
    if type(s)!=str:
        return float(s)
    i=0
    for c in s:
        if not c.isdigit() and c!='.':
            break
        i+=1
    if not i:
        return 0.0
    return float(s[:i])

def unpackTable(tbl):
    if not tbl:
        return [],[]
    tmp=tbl.split()
    N=int(tmp[0])
    dt=float(tmp[1])
    t=np.arange(N)*dt
    NN=(len(tmp)-2)//N
    cols=[]
    for i in range(NN):
        cols.append(np.array(list(map(float, tmp[2+N*i:2+N*(i+1)]))))
    return t, cols

def packTable(t, cols):
    rez='{0:d}\t\n{1}\t\n'.format(len(t), t[1]-t[0])
    for c in cols:
        for num in c:
            rez+='{}\t'.format(num)
        rez+='\t\n'
    return rez

class odbc:
    def __init__(self, dbFile):
        self.conn=pyodbc.connect(r"DRIVER={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=%s;" % dbFile)
        self.cursor=self.conn.cursor()
        
    def getInfo(self, table, fieldsCond='', fieldsCondValues='', getFields='*', cond=' and '):
        if type(getFields) not in [list, tuple]:
            getFields=[getFields]
        getFields=list(getFields)
        for i in range(len(getFields)):
            if getFields[i]!='*':
                getFields[i]='['+getFields[i]+']'
        getFields=' , '.join(getFields)
        if fieldsCond and fieldsCondValues:
            if type(fieldsCond)==str: fieldsCond=[fieldsCond]
            if type(fieldsCondValues) not in [list, tuple]: fieldsCondValues=[fieldsCondValues]
            if len(fieldsCond)!=len(fieldsCondValues):
                return []
            fieldsCond=list(fieldsCond)
            N=len(fieldsCond)
            for i in range(N):
                fieldsCond[i]='['+fieldsCond[i]+']'
            fieldsCond=list(fieldsCond)
            for i, f in enumerate(fieldsCond):
                fieldsCond[i]=f+'=?'
            fieldsCond=cond.join(fieldsCond)
        s='select {0} from {1} '.format(getFields, table)
        if fieldsCond and fieldsCondValues:
            s+='where {0}'.format(fieldsCond)
        if fieldsCondValues:
            data=self.cursor.execute(s, *fieldsCondValues).fetchall()
        else:
            data=self.cursor.execute(s).fetchall()
        allrez=[]
        for d in data:
            rez={}
            for i, c in enumerate(d.cursor_description):
                rez[c[0]]=d[i]
            allrez.append(rez)
        return allrez        

    def putInfo(self, table, putFields='', putFieldsValues='', fieldsCond='', fieldsCondValues='', cond=' and '):
        if type(putFields)==str:
            putFields=[putFields]
        if type(putFieldsValues) not in [list, tuple]:
            putFieldsValues=[putFieldsValues]
        if len(putFields)!=len(putFieldsValues):
            return 0
        putS=[]
        for f in putFields:
            putS.append('[{}] = ?'.format(f))
        putS=' , '.join(putS)
        if fieldsCond and fieldsCondValues:
            if type(fieldsCond)==str: fieldsCond=[fieldsCond]
            if type(fieldsCondValues) not in [list, tuple]: fieldsCondValues=[fieldsCondValues]
            if len(fieldsCond)!=len(fieldsCondValues):
                return 0
            fieldsCond=list(fieldsCond)
            N=len(fieldsCond)
            for i in range(N):
                fieldsCond[i]='['+fieldsCond[i]+']'
            for i, f in enumerate(fieldsCond):
                fieldsCond[i]=f+'=?'
            fieldsCond=cond.join(fieldsCond)
        s='update {0} set {1} '.format(table, putS)
        if fieldsCond and fieldsCondValues:
            s+='where {0}'.format(fieldsCond)
        if fieldsCondValues:
            vals=list(putFieldsValues)+list(fieldsCondValues)
            self.cursor.execute(s, *vals)
        else:
            self.cursor.execute(s, *putFieldsValues)
        self.cursor.commit()
        return 1

    def insertInfo(self, table, putFields='', putFieldsValues='', commit=True):
        if type(putFields)==str:
            putFields=[putFields]
        for i in range(len(putFields)):
            putFields[i]='['+putFields[i]+']'
        if type(putFieldsValues) not in [list, tuple]:
            putFieldsValues=[putFieldsValues]
        if len(putFields)!=len(putFieldsValues):
            return 0
        putS=' , '.join(putFields)
        s='insert into {0}({1}) values ({2})'.format(table, putS, ' , '.join(['?']*len(putFields)))
        self.cursor.execute(s, *putFieldsValues)
        if commit:
            self.cursor.commit()
        return 1

    def close(self):
        self.conn.close()

class bar(object):
    def __init__(self, odbcRez):
        self.E = tofloat(odbcRez['���������������(���)'])
        self.code = str(odbcRez['������������']) 
        self.mat = str(odbcRez['��������'])
        self.d = tofloat(odbcRez['�������(��)'])
        self.d0 = tofloat(odbcRez['�����������������'])
        self.c = tofloat(odbcRez['�������������(����)'])
        self.l = tofloat(odbcRez['�����(��)'])
        self.S = np.pi*self.d**2/4.
    def __repr__(self):
        rez=''
        rez+=f'��� �������: {self.code}\n'
        rez+=f'�������� �������: {self.mat}\n'
        rez+=f'E = {self.E} ���\n'
        rez+=f'c = {self.c} �/c\n'
        rez+=f'l = {self.l} ��\n'
        rez+=f'd = {self.d} ��\n'
        rez+=f'd0 = {self.d0} ��\n'
        return rez
class striker(object):
    def __init__(self, odbcRez):
        self.code = str(odbcRez['�����������']) 
        self.mat = str(odbcRez['����������������'])
        self.d = tofloat(odbcRez['���������������(��)'])
        self.l = tofloat(odbcRez['�������������(��)'])
        self.S = np.pi*self.d**2/4.
    def __repr__(self):
        rez=''
        rez+=f'��� ��������: {self.code}\n'
        rez+=f'�������� ��������: {self.mat}\n'
        rez+=f'l = {self.l} ��\n'
        rez+=f'd = {self.d} ��\n'
        return rez
class experimentalData(object):
    def __init__(self, odbcRez):
        self.data = str(odbcRez['����']).split()[0]
        self.code = str(odbcRez['����������'])
        self.striker = str(odbcRez['�������'])
        self.expType = str(odbcRez['���������������'])
        self.T = tofloat(odbcRez['�����������'])
        self.P = tofloat(odbcRez['�����������'])
        self.V = tofloat(odbcRez['����������������'])
        self.d0 = tofloat(odbcRez['�������'])
        self.l0 = tofloat(odbcRez['�����'])
        self.l = tofloat(odbcRez['���������������'])
        self.d = tofloat(odbcRez['�����'])
        self.note = str(odbcRez['����������'])
        self.osc = {}
        self.osc['t'], self.osc['rays'] = unpackTable(odbcRez['�������������'])
        self.pulses = {}
        self.pulses['t'], self.pulses['pulses'] = unpackTable(odbcRez['��������������������'])
        self.tarir = [
            tofloat(odbcRez['��������������������������']),
            tofloat(odbcRez['��������������������������']),
            tofloat(odbcRez['��������������������������2(������)'])
        ]
        self.datPosition = [
            tofloat(odbcRez['�������������������(��)']),
            tofloat(odbcRez['�������������������(��)'])
        ]
        self.bars = [
            str(odbcRez['�������������������']),
            str(odbcRez['���������������']),
            str(odbcRez['���������������2(������)'])
        ]
    def __repr__(self):
        rez = ''
        rez+= f'��� ������������: {self.code}\n'
        return rez

class expODBC(odbc):
    def __init__(self, dbFile):
        super().__init__(dbFile)
    def getExpTypes(self):
        return self.getInfo(getFields=('���������������', '���������������'), table='���������������')
    def getMaterials(self):
        return self.getInfo(getFields=('��������', '������������'), table='�������������������')
    def getNumbers(self, expType, materialCode):
        materialCode = tofloat(materialCode)
        return self.getInfo(getFields='������������', table='�����������',
                            fieldsCond=('���������������', '������������'),
                            fieldsCondValues=(expType, materialCode))   
    def getExperimentData(self, sampleCode):
        rez = self.getInfo(table='�����������', fieldsCond='����������',
                            fieldsCondValues=sampleCode)
        if rez:
            return experimentalData(rez[0])
        else:
            print(f'� ���� �� ������ ����������� � ����� {sampleCode}')
            return None
    def getStrickerData(self, strickerCode):
        rez = self.getInfo(table='�������', fieldsCond='�����������',
                            fieldsCondValues=strickerCode)       
        if rez:
            return striker(rez[0])
        else:
            print(f'� ���� �� ������ ������� � ����� {strickerCode}')
            return None
            
    def getBarData(self, barCode):
        rez = self.getInfo(table='��������������', fieldsCond='������������',
                            fieldsCondValues=barCode)
        if rez:
            return bar(rez[0])
        else:
            print(f'� ���� �� ������ �������� � ����� {barCode}')
            return None
    def putOsc(self, sampleCode, data):
        self.putInfo(table='�����������', putFields='�������������', putFieldsValues=data,
                     fieldsCond='����������', fieldsCondValues=sampleCode)
    def putPulses(self, sampleCode, data):
        self.putInfo(table='�����������', putFields='��������������������', putFieldsValues=data,
                     fieldsCond='����������', fieldsCondValues=sampleCode)

if __name__=='__main__':
    pass
    dbFile=r"d:\experiments\���������������������������������.accdb"
    db=expODBC(dbFile)
    d=db.getExperimentData('c620-62')
    print(d['�������������'][0:50])
    #db.insertInfo(table='�����������', putFields=('�����������', '������������'), putFieldsValues=(25, '013'))
    

#%%
