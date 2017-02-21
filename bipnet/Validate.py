#!/usr/bin/python

import sys, getopt
import random as rd
import scipy.stats as st
import igraph as ig
from collections import Counter,defaultdict
import numpy as np
import pandas as pd
    
def Pvalue(gb,g,to_ck,which=False):
	I = g.vs["Tid"]
	
	if which==False:  
		Nb =  sum(gb.vs["type"])
	else:
		Nb = gb.vcount()-sum(gb.vs["type"])

	survivor = st.hypergeom.sf 
	D = set(tuple([g[(a,b)]]+sorted([gb.vs[I[a]].degree(),gb.vs[I[b]].degree()])) for a,b in to_ck)
	PV = {(a,b,c): survivor(a-1,Nb,b,c) for a,b,c in D}
	PB = [(PV[tuple([g[(a,b)]]+sorted([gb.vs[I[a]].degree(),gb.vs[I[b]].degree()]))],(a,b),g[(a,b)]) for a,b in to_ck ]
	return PB

def SVN(gb,which=False,alpha=0.01):
		
	gb.vs["Tid"] = range(gb.vcount())	
	g = gb.bipartite_projection(multiplicity=True,which=which)
	
	to_ck = [e.tuple for e in g.es]
	PB = Pvalue(gb,g,to_ck,which)
	g.es["pvalue"] = zip(*PB)[0]
	
	'Bonferroni'
	alpha = 0.01
	bnf = alpha/(g.vcount()*(g.vcount()-1))/2
	PVf = filter(lambda x:x[0]<bnf,PB)

	'FDR'
	PB = sorted(PB)
	p,ed,w = zip(*PB)
	sfdr = bnf*np.arange(1,len(PB)+1)
	s = np.where(p<sfdr)[0]
	if len(s)>0:
		s = s[-1]
		PBf = PB[:s]
	else:
		PBf = []
	

	if len(PVf)>0:
		ED = list(zip(*PVf)[1])
		gbonf = ig.Graph(g.vcount(),edges=ED)
		gbonf.vs["name"] = g.vs["name"]
		gbonf.es["weight"]=list(zip(*PVf)[2])			
	else:
		gbonf = ig.Graph(g.vcount())
		gbonf.vs["name"] = g.vs["name"]
		gbonf.es["weight"] = 1.0
	   
	if len(PBf)>0:
		ED = list(zip(*PBf)[1])
		gfdr = ig.Graph(g.vcount(),edges=ED)
		gfdr.vs["name"] = g.vs["name"]
		gfdr.es["weight"]=list(zip(*PBf)[2])

	else:
		gfdr = ig.Graph(g.vcount())
		gfdr.vs["name"] = g.vs["name"]
		gfdr.es["weight"] = 1.0
		
	del g.vs["Tid"]
	del gb.vs["Tid"]

	return g,gbonf,gfdr


def main(argv):
	inputfile = ''
	outputfile = ''
	side = 1

	try:
		opts, args = getopt.getopt(argv,"hi:o:s:t",["ifile=","ofile=","side=","thr="])
	except getopt.GetoptError:
		print 'Validate.py -i <inputfile> -o <outputfile> --side <bool> --thr <float>' 
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Validate.py -i <inputfile> -o <outputfile> --side <bool> --thr <float>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt in ("-s", "--side"):
			side = int(arg)
		elif opt in ("-t", "--thr"):
			alpha = float(arg)



	return inputfile,outputfile,side,alpha

def get_BipartiteFromFile(file_r):
	
	ED = pd.read_table(file_r,sep=' ',header=None)
	
	A = list(set(ED[0]))
	B = list(set(ED[1]))
	
	iA = dict(zip(A,range(len(A))))
	iB = dict(zip(B,range(len(A),len(A)+len(B))))
	
	ed = [(iA[a],iB[b]) for a,b in zip(ED[0],ED[1])]
	nod = [0]*len(A)+[1]*len(B)
	
	gb = ig.Graph.Bipartite(nod,ed)
	gb.vs["name"] = A+B
	
	return gb

def write_validation(g,file_w):
	nm = g.vs["name"]
	to_ck = [e.tuple for e in g.es]
	
	X = ["%s %s %d"%(nm[a],nm[b],w) for (a,b),w in zip(to_ck,g.es["weight"])]
	
	with open(file_w,"w") as fw:
		fw.write("\n".join(X)+"\n")
	return

if __name__=='__main__':
	
	file_r,file_w,side,alpha = main(sys.argv[1:])

	gb = get_BipartiteFromFile(file_r)
		
	g,gbonf,gfdr = SVN(gb,side,alpha)
	
	fw_full = file_w.split('.')[0]+'_full.net'
	fw_bonf = file_w.split('.')[0]+'_bonf.net'
	fw_fdr = file_w.split('.')[0]+'_fdr.net'
	
	write_validation(g,fw_full)
	write_validation(gbonf,fw_bonf)
	write_validation(gfdr,fw_fdr)
	





