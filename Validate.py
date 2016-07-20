#!/usr/bin/python

import sys, getopt
import random as rd
import scipy.stats as st
import igraph as ig
from collections import Counter,defaultdict
import numpy as np
import pandas as pd
from copy import deepcopy


def CONTRACT_COMMUNITY(g,CB):
    h = deepcopy(g)
    h.contract_vertices(CB, combine_attrs=list)
    h.simplify(loops=False,combine_edges=dict(weight=sum))
    CN = h.community_multilevel(weights=h.es["weight"])
    NN = dict([(a,e) for e,b in zip(CN.membership,h.vs["name"]) for a in b])
    R =[NN[b] for b in g.vs["name"]]
    return R
    
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
	
	bnf = alpha/(g.vcount()*(g.vcount()-1))/2
	PVf = filter(lambda x:x[0]<bnf,PB)
	
	pv = np.array(zip(*PB)[0])
	pv.sort()

	bnf = 2*alpha/(g.vcount()*(g.vcount()-1))
	ind_fdr =np.where(pv < np.arange(1,len(pv)+1)*bnf)[0][-1]

	PBf = PB[:ind_fdr]
	
	#PBf,PB = filter(lambda (i,x):x[0]<(i+1)*bnf,enumerate(sorted(PB))) ,PVf

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
		#PBf = list(zip(*PBf)[1])
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
	  opts, args = getopt.getopt(argv,"hi:o:s",["ifile=","ofile=","side="])
	except getopt.GetoptError:
	  print 'Validate.py -i <inputfile> -o <outputfile> --side <bool>' 
	  sys.exit(2)
	for opt, arg in opts:
	  if opt == '-h':
		 print 'Validate.py -i <inputfile> -o <outputfile> --side <bool>' 
		 sys.exit()
	  elif opt in ("-i", "--ifile"):
		 inputfile = arg
	  elif opt in ("-o", "--ofile"):
		 outputfile = arg
	  elif opt in ("-s", "--side"):
		 side = int(arg)


	return inputfile,outputfile,side

def _get_BipartiteFromFile(file_r):
	
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

def _write_validation(g,file_w):
	nm = g.vs["name"]
	to_ck = [e.tuple for e in g.es]
	
	X = ["%s\t%s\t%d"%(nm[a],nm[b],w) for (a,b),w in zip(to_ck,g.es["weight"])]
	
	with open(file_w,"w") as fw:
		fw.write("\n".join(X)+"\n")
	return

if __name__=='__main__':
	
	file_r,file_w,side = main(sys.argv[1:])
	
	gb = _get_BipartiteFromFile(file_r)
		
	g,gbonf,gfdr = SVN(gb,side)
	
	fw_full = file_w.split('.')[0]+'_full.net'
	fw_bonf = file_w.split('.')[0]+'_bonf.net'
	fw_fdr = file_w.split('.')[0]+'_fdr.net'
	
	_write_validation(g,fw_full)
	_write_validation(gbonf,fw_bonf)
	_write_validation(gfdr,fw_fdr)
	
	'classic community detection'
	Cfull = g.community_multilevel(weights=g.es["weight"]).membership
	Cbonf = gbonf.community_multilevel(weights=gbonf.es["weight"]).membership
	Cfdr = gfdr.community_multilevel(weights=gfdr.es["weight"]).membership

	'Community detection on contracted node'
	Cfull_bonf = CONTRACT_COMMUNITY(g,Cbonf)
	Cfull_fdr = CONTRACT_COMMUNITY(g,Cfdr)
	
	'write community'
	T = ['Full','Bonf','FDR','Full+Bonf','Full+FDR']
	X = np.array([Cfull,Cbonf,Cfdr,Cfull_bonf,Cfull_fdr]).transpose()
	

	X = pd.DataFrame(X,columns=T)
	X.to_csv(file_w.split('.')[0]+'_CommunityMembership.net',index=False,sep='\t')
	





