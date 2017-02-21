#!/usr/bin/python

import igraph as ig
from Validate import SVN
import scipy.stats as st
import metric as mtr
import random as rd
import numpy as np
import pandas as pd
import sys, getopt

def ADD_NOISE(Net,N,prew,Pb=None):
	A = set(range(N))
	for x in xrange(len(Net)):
		n = st.binom.rvs(len(Net[x]),1.-prew)
		s = rd.sample(Net[x],n)
		sel = A - set(s)

		if Pb==None: Net[x] = s + rd.sample(sel,len(Net[x])-n)
		else:
			p = Pb[list(sel)]
			p = p/sum(p)
			
			Net[x] = s+list(np.random.choice(list(sel),replace=False,size=len(Net[x])-n,p=p))

			
	NOD = list(np.repeat(0,len(Net)))+list(np.repeat(1,N))
	EDG = [(i,len(Net)+p) for i in range(len(Net)) for p in Net[i]]
		
	gb = ig.Graph.Bipartite(NOD,EDG)
	gb.vs["name"] = map(str,range(gb.vcount()))

	return gb 

def CREATE_NET(SA,SB,pc,pr,degree_dist=None,Pb=None):
	
	q = len(SA)
	
	if degree_dist==None:
		'Create Community'
		N = []
		Net = []
		for i in xrange(q):
			if N==[]: M = 0
			else: M = max(N)+1
			n = map(lambda x: x+M,xrange(SB[i]))
			k = st.binom.rvs(len(n),pc[i],size=SA[i])
			Net.extend([rd.sample(n,k[i]) for i in xrange(len(k))])
			N.extend(n)
	else:
		Nb = sum(SB)
		P = degree_dist
		'Degree Sequance'
		N = []
		Net = []
		for i in range(q):
			if N==[]: M = 0
			else: M = max(N)+1
			n = map(lambda x: x+M,range(SB[i]))
			k2 = P[sum(SA[:i]):sum(SA[:i+1])]
			k = [a  if a<=SB[i] else len(n) for a in k2]
			r = np.array(k2)-np.array(k)
			CMP = list(set(range(Nb)) -set(n))
			s = [rd.sample(n,k[i]) for i in xrange(len(k))]
			l = [rd.sample(CMP,r[i]) for i in xrange(len(r))]
			ad = [a+b for a,b in zip(s,l)]
			Net.extend(ad)
			N.extend(n)

	'Reference Membership'
	O =[i for i,c in enumerate(SA) for h in xrange(c)]
	
	gb = ADD_NOISE(Net,len(N),pr,Pb)
	
	return gb,O

def main(argv):
	inputfile = ''
	outputfile = ''
	noise = 0.
	sim = 1

	try:
	  opts, args = getopt.getopt(argv,"hi:o:n",["ifile=","ofile=","noise="])
	except getopt.GetoptError:
		print 'Benchmark.py -i <inputfile> -o <outputfile> --noise <noise>' 
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Benchmark.py -i <inputfile> -o <outputfile> --noise <noise>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt in ("-n", "--noise"):
			n = float(arg)

	return inputfile,outputfile,n


def write_network(g,file_w):
	nm = g.vs["name"]
	to_ck = [e.tuple for e in g.es]
	
	X = ["%s %s"%(nm[a],nm[b]) for a,b in to_ck]
	
	with open(file_w,"w") as fw:
		fw.write("\n".join(X)+"\n")
	return

if __name__=='__main__':
	
	inputfile,outputfile,pr = main(sys.argv[1:])
	X = pd.read_table(inputfile)
	
	SA = list(X['SA'].astype(int))
	SB = list(X['SB'].astype(int))
	pc = list(X['pc'].astype(float))
	
	gb,O = CREATE_NET(SA,SB,pc,pr)
	write_network(gb,outputfile)

	
