#!/usr/bin/python

import scipy.stats as st
from collections import Counter,defaultdict
import sys
import pandas as pd
from itertools import combinations

def AssortativeMixing(g,c):
	'It wants as input an (g) igraph weighted objected, and an list/array with the membership of the nodes. It returns the Assortativity Coefficent'
	
    cl = defaultdict(list)
    for i,x in enumerate(c):
        cl[x].append(i)
    cl = [np.array(cl[i]) for i in cl if cl[i]>1]

    s = g.strength(weights=g.es["weight"])

    kk= sum(sum(s[a]*s[b] for a,b in combinations(xc, 2) ) for xc in cl)

    m = sum(s)

    return g.modularity(c,weights=g.es["weight"])/(1.-(2.*kk)/m**2)

def awi(A,B):
	'It wants as input two list/array  (A,B) with the membership of the nodes. It will evaluate the AWI  A in B. That is an estimator of the inclusion of A in B'
    
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
