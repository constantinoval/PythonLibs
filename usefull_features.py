import numpy as np


def group_by(data, by, tol=1e-6):
    tmp = data.copy()
    tmp.sort(order=by)
    ii = [0] + (np.where(np.diff(tmp[by]) >= tol)
                [0]+1).tolist() + [len(tmp[by])]
    groups = []
    for i in range(len(ii)-1):
        groups.append(tmp[ii[i]:ii[i+1]])
    return groups


def runningMeanFast(x, N):
    rez = np.roll(np.convolve(x, np.ones((N,))/N)[(N-1):], N//2)
    rez[:N//2] = rez[N//2]
    return rez

def splitByN(l, n):
    if n>=len(l):
        return [l]
    nn=len(l)//n
    rez=[]
    for i in range(nn):
        rez.append(l[i*n:(i+1)*n])
    if len(l)%n:
        rez.append(l[nn*n:])
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

class ConsoleProgress(object):
    def __init__(self, total, scaleLength=50, text='Progress=',
                 symbol='-#'):
        self.current = 0
        self.total = total
        self.chekRedraw = lambda i: not (i % (total//scaleLength))
        self.scaleLength = scaleLength
        self.text = text
        self.symbol = symbol

    def redraw(self):
        if self.chekRedraw(self.current):
            print(progress(self.current, self.total, self.scaleLength, self.text, self.symbol), end='')#, flush=True)
        self.current+=1

    def finalize(self):
        print(progress(self.total, self.total, self.scaleLength, self.text, self.symbol))


def progress(cur, total, scaleLength=20, text='Progress=', symbol='-#'):
    if total!=0:
        prog=int(100.*cur/total)
    else:
        prog=100
    rez='\r{0}{1:4}% '.format(text, prog)
    numD=int(prog*scaleLength/100.)
    rez+=symbol[1]*numD+symbol[0]*(scaleLength-numD)
    rez+=' Total={0:10}, current={1:10}'.format(total,cur)
    return rez