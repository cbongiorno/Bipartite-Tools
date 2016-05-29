#!/usr/bin/python

import scipy.stats as st
from collections import Counter,defaultdict
import sys
import pandas as pd

def awi(A,B):
    
    mb = defaultdict(list)
    for i in xrange(len(A)):
        mb[A[i]].append(i)
    mb = [mb[e] for e in mb if len(mb[e])>1]
    
    mr = defaultdict(list)
    for i in xrange(len(B)):
        mr[B[i]].append(i)
    mr = [mr[e] for e in mr if len(mr[e])>1]
    Bp = sum(map(lambda x:len(x)*(len(x)-1)/2,mr))

    N = len(A)*(len(A)-1)/2
    
    n = 0
    m = 0
    for x in mb:
        for i in xrange(len(x)):
            for j in range(i+1,len(x)):
                n+=1
                if B[x[i]]==B[x[j]]:
                    m+=1
                    
    Ap = float(n)
    if (Ap-(Bp*Ap)/N)==0:
        return None,1.
    if n!=0:
        return (m-(Bp*Ap)/N)/(Ap-(Bp*Ap)/N),st.hypergeom.sf(m-1,N,Ap,Bp)
    else:
        return None,1.

if __name__=='__main__':

	x = sys.argv[1]

	x = pd.read_csv(x,sep='\t',header=None)
	A = awi(x[0],x[1])
	
	print "AWI: ",A[0]
	print "pvalue: ",A[1]
