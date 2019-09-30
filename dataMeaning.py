import matplotlib.pylab as plt
import pandas
import numpy as np
import scipy as sp
import scipy.stats
import json


def mean_confidence_interval(data, confidence=0.95):
    N=len(data[0])
    c=len(data)
    m=[]
    h=[]
    for i in range(N):
        a=[]
        for j in range(c):
            a.append(data[j][i])
        a = np.array(a)
        mm, se = np.mean(a), scipy.stats.sem(a)
        hh = se * sp.stats.t._ppf((1+confidence)/2., c-1)
        m.append(mm)
        h.append(hh)
    return m, h

def readDataFromXLS(xls, expCode, expNumbers):
    print(f'Reading data from {xls}')
    data = {}
    for N in expNumbers:
        print(f'\tReading experiment {expCode}-{N:02d}')
        d=pandas.read_excel(xls, sheet_name=expCode + '-{:02d}'.format(N))
        exp={}
        for k in ['время(мкс)', 'деформация', 'напряжение(МПа)', 'скорость деформации(1/c)']:
            exp[k]=list(d[k])
        data[expCode+'-{:02d}'.format(N)]=exp
    return data

def meanCurvesT(data, confidence=0.95, plot=True, numPoints=200, epslim=None, stressLevel=0.05, picPrefix=None):
    expKeys=list(data.keys())
    tmax = data[expKeys[0]]['время(мкс)'][-1]
    tmin = data[expKeys[0]]['время(мкс)'][0]
    stress0 = stressLevel * max(data[expKeys[0]]['напряжение(МПа)'])
    for j, s in enumerate(data[expKeys[0]]['напряжение(МПа)']):
        if s >= stress0:
            t0 = data[expKeys[0]]['время(мкс)'][j]
            break
    for exp in data:
        for j, s in enumerate(data[exp]['напряжение(МПа)']):
            if s >= stress0:
                tt = data[exp]['время(мкс)'][j]
                break
        data[exp]['время(мкс)'] = list(map(lambda x: x-tt+t0, data[exp]['время(мкс)']))
        tmin = max(tmin, data[exp]['время(мкс)'][0])
        tmax = min(tmax, data[exp]['время(мкс)'][- 1])
    tt = np.linspace(tmin, tmax, numPoints)
    ee = []
    ss = []
    dee = []
    for d in data:
        ee.append(np.interp(tt, data[d]['время(мкс)'], data[d]['деформация']))
        ss.append(np.interp(tt, data[d]['время(мкс)'], data[d]['напряжение(МПа)']))
        dee.append(np.interp(tt, data[d]['время(мкс)'], data[d]['скорость деформации(1/c)']))
    e, he = mean_confidence_interval(ee, confidence)
    s, hs = mean_confidence_interval(ss, confidence)
    de, hde = mean_confidence_interval(dee, confidence)
    if plot:
        c = []
        for i, d in enumerate(data):
            l, = plt.plot(data[d]['деформация'], data[d]['напряжение(МПа)'], label=d)
            c.append(l.get_color())
        plt.grid()
        plt.legend(bbox_to_anchor=(1.5, 1))
        plt.xlabel('деформация')
        plt.ylabel('напряжение, МПа')
        plt.twinx()
        plt.ylabel('скорость деформации, 1/c')

        for i, d in enumerate(data):
            l, = plt.plot(data[d]['деформация'], data[d]['скорость деформации(1/c)'], '--',
                          color=c[i], label=d)
        if epslim: plt.xlim(0, epslim)
        if picPrefix:
            plt.gcf().savefig(picPrefix+'-all.png')
        plt.figure()
        for d in data:
            plt.plot(data[d]['деформация'], data[d]['напряжение(МПа)'])
        plt.errorbar(e, s, yerr=hs, xerr=he, color='k', errorevery=3)
        plt.grid()
        plt.xlabel('деформация')
        plt.ylabel('напряжение, МПа')
        if picPrefix:
            plt.gcf().savefig(picPrefix+'-es.png')

        plt.figure()
        for d in data:
            plt.plot(data[d]['время(мкс)'], data[d]['деформация'])
        plt.grid()
        plt.errorbar(tt, e, yerr=he, color='k', errorevery=3)
        plt.xlabel('время, мкс')
        plt.ylabel('деформация')
        if picPrefix:
            plt.gcf().savefig(picPrefix+'-te.png')

        plt.figure()
        for d in data:
            plt.plot(data[d]['время(мкс)'], data[d]['напряжение(МПа)'])
        plt.grid()
        plt.errorbar(tt, s, yerr=hs, color='k', errorevery=3)
        plt.xlabel('время, мкс')
        plt.ylabel('напряжение, МПа')
        if picPrefix:
            plt.gcf().savefig(picPrefix+'-ts.png')

        plt.figure()
        for d in data:
            plt.plot(data[d]['время(мкс)'], data[d]['скорость деформации(1/c)'])
        plt.grid()
        plt.errorbar(tt, de, yerr=hde, color='k', errorevery=3)
        plt.xlabel('время, мкс')
        plt.ylabel('скорость деформации, 1/c')
        if picPrefix:
            plt.gcf().savefig(picPrefix+'-tde.png')

        plt.show()
    return {'et': list(e), 'st': list(s), 'det': list(de),
            'he': list(he), 'hs': list(hs), 'hde': list(hde),
            't': list(tt)}

if __name__=='__main__':
    NN = [6, 2, 10, 4, 3, 9, 5, 1, 8, 7]
    data=readDataFromXLS('с633.xls', 'c633', NN)
    r = meanCurvesT(data, stressLevel=0.25, picPrefix='c633-450-20')
    #json.dump(r, open('c633-450-20.json', 'w'))