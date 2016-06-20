#!/usr/bin/python

import igraph as ig
from Validate import SVN,CONTRACT_COMMUNITY
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

	
def COMPARE(SA,SB,pc,pr,nsim):
	Res = []
	for isim in xrange(nsim):

		gb,O = CREATE_NET(SA,SB,pc,pr)

		g,gbonf,gfdr = SVN(gb,False)
			
		'classic community detection'
		Cfull = g.community_multilevel(weights=g.es["weight"]).membership
		Cbonf = gbonf.community_multilevel(weights=gbonf.es["weight"]).membership
		Cfdr = gfdr.community_multilevel(weights=gfdr.es["weight"]).membership

		'Community detection on contracted node'
		Cfull_bonf = CONTRACT_COMMUNITY(g,Cbonf)
		Cfull_fdr = CONTRACT_COMMUNITY(g,Cfdr)

		'Estimate ARI'
		AriFull = ig.compare_communities(Cfull,O,method='ari')
		AriBonf = ig.compare_communities(Cbonf,O,method='ari')
		AriFdr = ig.compare_communities(Cfdr,O,method='ari')
		AriFullBonf = ig.compare_communities(Cfull_bonf,O,method='ari')
		AriFullFdr = ig.compare_communities(Cfull_fdr,O,method='ari')

		'Estimate AWI'
		AwiFull,pAwiFull = mtr.awi(Cfull,O)
		AwiBonf,pAwiBonf = mtr.awi(Cbonf,O)
		AwiFdr,pAwiFdr = mtr.awi(Cfdr,O)
		AwiFullBonf,pAwiFullBonf = mtr.awi(Cfull_bonf,O)
		AwiFullFdr,pAwiFullFdr  = mtr.awi(Cfull_fdr,O)

		'Estimate Modularity'
		ModRef = g.modularity(O,weights=g.es["weight"])
		ModFull = g.modularity(Cfull,weights=g.es["weight"])
		ModFullBonf = g.modularity(Cfull_bonf,weights=g.es["weight"])
		ModFullFdr = g.modularity(Cfull_fdr,weights=g.es["weight"])

		'Aggregate solution'
		X = [AriFull,AriBonf,AriFdr,AriFullBonf,AriFullFdr]
		X.extend([AwiFull,AwiBonf,AwiFdr,AwiFullBonf,AwiFullFdr])
		X.extend([ModRef,ModFull,ModFullBonf,ModFullFdr])
		X.extend([pAwiFull,pAwiBonf,pAwiFdr,pAwiFullBonf,pAwiFullFdr])
		
		Res.append(X)
	
	T = ['ARI Full','ARI Bonf','ARI FDR','ARI Full+Bonf','ARI Full+FDR']
	T.extend(['AWI Full','AWI Bonf','AWI FDR','AWI Full+Bonf','AWI Full+FDR'])
	T.extend(['Modularity Reference','Modularity Full','Modularity Full+Bonf','Modularity Full+FDR'])
	T.extend(['pvalue AWI Full','pvalue AWI Bonf','pvalue AWI FDR','pvalue AWI Full+Bonf','pvalue AWI Full+FDR'])
	
	Res = pd.DataFrame(Res,columns=T)
	
	return Res
	

def main(argv):
	inputfile = ''
	outputfile = ''
	noise = 0.
	sim = 1

	try:
	  opts, args = getopt.getopt(argv,"hi:o:n:r",["ifile=","ofile=","noise=","run="])
	except getopt.GetoptError:
		print 'Benchmark.py -i <inputfile> -o <outputfile> --noise <noise> --run <run>' 
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Benchmark.py -i <inputfile> -o <outputfile> --nrun <run> --ncpu <ncpu> --hier <bool>' 
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt in ("-n", "--noise"):
			n = float(arg)
		elif opt in ("-r", "--run"):
			sim = int(arg)

	return inputfile,outputfile,noise,sim


if __name__=='__main__':
	
	inputfile,outputfile,pr,sim = main(sys.argv[1:])
	X = np.array(pd.read_table(inputfile,header=None))
	
	SA = list(X[0].astype(int))
	SB = list(X[1].astype(int))
	pc = list(X[2].astype(float))
	
	x = COMPARE(SA,SB,pc,pr,sim)

	x.to_csv(outputfile,index=False)
	
