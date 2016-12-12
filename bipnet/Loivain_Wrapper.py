#!/usr/bin/python

import igraph as ig
import numpy as np
from matplotlib.pylab import flatten
from collections import OrderedDict,Counter
import os
import pandas as pd


def LoivainModified(g,n,seed,part=None,DIRs=None):
	"It wants mandatory a igraph object (g), an integer number (n) that must be different if you computer in paraller different network, an integer number (seed) the it is the seed used to initialize the random function, Optionally you can provide a seed partition (part) the must be different from the trivail one i.e. each node in a separate partition, the last argument is the path of the directory of the Community_latest C code. It returns the membership of the nodes"

	

	
	if not os.path.exists(DIRs+'convert'):
		print "Error please cheak the Directory %s"%DIRs
		print "Maybe you need to compile the Loivain code"
		exit(0)
	
	to_ck = [e.tuple for e in g.es]

	Od= OrderedDict.fromkeys(list(flatten(to_ck)))

	Od = dict(zip(Od.keys(),range(len(Od))))

	x = [(Od[a],Od[b],w) for (a,b),w in zip(to_ck,g.es["weight"])]

	x = pd.DataFrame(x)

	x.to_csv('/tmp/graph_%d.txt'%n,sep=' ',index=False,header=False)
	
	
	os.system(DIRs+'convert  -i /tmp/graph_%d.txt -o /tmp/graph_%d.bin -w /tmp/graph_%d.weights'%(n,n,n))
	#print DIRs+'convert  -i /tmp/graph_%d.txt -o /tmp/graph_%d.bin -w /tmp/graph_%d.weights'%(n,n,n)
	
	if part!=None:
		Part = np.zeros(g.vcount()).astype(int)
		for i,p in enumerate(part):
			Part[Od[i]] = p

		Part = zip(range(g.vcount()),Part)

		pd.DataFrame(Part).to_csv('/tmp/Part_%d.part'%n,sep=' ',header=False,index=False)


	if part!=None:
		os.system(DIRs+'community %d /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights -p /tmp/Part_%d.part > /tmp/graph_%d.tree'%(seed,n,n,n,n))
		#print DIRs+'community %d /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights -p /tmp/Part_%d.part > /tmp/graph_%d.tree'%(seed,n,n,n,n)
	else:
		os.system(DIRs+'community %d /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights  > /tmp/graph_%d.tree'%(seed,n,n,n))
		#print DIRs+'community %d /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights  > /tmp/graph_%d.tree'%(seed,n,n,n)

	lv = int(os.popen(DIRs+'hierarchy /tmp/graph_%d.tree'%n).read().split('\n')[0].split(': ')[1])-1
	#print DIRs+'hierarchy /tmp/graph_%d.tree'%n

	os.system(DIRs+'hierarchy /tmp/graph_%d.tree -l %d > /tmp/graph_node2comm_level2_%d'%(n,lv,n))
	#print DIRs+'hierarchy /tmp/graph_%d.tree -l %d > /tmp/graph_node2comm_level2_%d'%(n,lv,n)
	
	cf = pd.read_csv('/tmp/graph_node2comm_level2_%d'%n,header=None,sep=' ')[1]

	Odrv = dict(map(list,map(reversed,Od.items())))

	c = np.array([-1]*g.vcount())
	if cf[0]!='of':
		for i,a in enumerate(cf):
			c[Odrv[i]] = a

	iso = np.where(c==-1)[0]
	c[iso] = range(max(c)+1,max(c)+len(iso)+1)

	c = c.astype(int)

	
	os.system('rm /tmp/graph_%d.bin /tmp/graph_%d.txt /tmp/graph_%d.tree /tmp/graph_%d.weights /tmp/graph_node2comm_level2_%d'%(n,n,n,n,n))
	if part!=None:
		os.system('rm  /tmp/Part_%d.part'%n)
	
	return c


def LoivainOriginal(g,n=0,part=None):
	"It wants mandatory a igraph object (g), an integer number (n) that must be different if you computer in paraller different network (optional), ,Optionally you can provide a seed partition (part) the must be different from the trivail one i.e. each node in a separate partition. It returns the membership of the nodes"


	to_ck = [e.tuple for e in g.es]

	Od= OrderedDict.fromkeys(list(flatten(to_ck)))

	Od = dict(zip(Od.keys(),range(len(Od))))

	x = [(Od[a],Od[b],w) for (a,b),w in zip(to_ck,g.es["weight"])]

	x = pd.DataFrame(x)

	x.to_csv('/tmp/graph_%d.txt'%n,sep=' ',index=False,header=False)
	
	
	os.system('convertX  -i /tmp/graph_%d.txt -o /tmp/graph_%d.bin -w /tmp/graph_%d.weights'%(n,n,n))
	#print DIRs+'convert  -i /tmp/graph_%d.txt -o /tmp/graph_%d.bin -w /tmp/graph_%d.weights'%(n,n,n)
	
	if part!=None:
		Part = np.zeros(g.vcount()).astype(int)
		for i,p in enumerate(part):
			Part[Od[i]] = p

		Part = zip(range(g.vcount()),Part)

		pd.DataFrame(Part).to_csv('/tmp/Part_%d.part'%n,sep=' ',header=False,index=False)


	if part!=None:
		os.system('loivainX /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights -p /tmp/Part_%d.part > /tmp/graph_%d.tree'%(n,n,n,n))
		#print DIRs+'community %d /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights -p /tmp/Part_%d.part > /tmp/graph_%d.tree'%(seed,n,n,n,n)
	else:
		os.system('louvainX /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights  > /tmp/graph_%d.tree'%(n,n,n))
		#print DIRs+'louvain /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights  > /tmp/graph_%d.tree'%(n,n,n)

	lv = int(os.popen('hierarchyX /tmp/graph_%d.tree'%n).read().split('\n')[0].split(': ')[1])-1
	#print DIRs+'hierarchy /tmp/graph_%d.tree'%n
	#print os.popen(DIRs+'hierarchy /tmp/graph_%d.tree'%n).read()
	#print lv
	#print DIRs+'hierarchy /tmp/graph_%d.tree'%n

	os.system('hierarchyX /tmp/graph_%d.tree -l %d > /tmp/graph_node2comm_level2_%d'%(n,lv,n))
	#print DIRs+'hierarchy /tmp/graph_%d.tree -l %d > /tmp/graph_node2comm_level2_%d'%(n,lv,n)
	
	cf = pd.read_csv('/tmp/graph_node2comm_level2_%d'%n,header=None,sep=' ')[1]

	Odrv = dict(map(list,map(reversed,Od.items())))

	c = np.array([-1]*g.vcount())
	if cf[0]!='of':
		for i,a in enumerate(cf):
			c[Odrv[i]] = a

	iso = np.where(c==-1)[0]
	c[iso] = range(max(c)+1,max(c)+len(iso)+1)

	c = c.astype(int)

	
	os.system('rm /tmp/graph_%d.bin /tmp/graph_%d.txt /tmp/graph_%d.tree /tmp/graph_%d.weights /tmp/graph_node2comm_level2_%d'%(n,n,n,n,n))
	if part!=None:
		os.system('rm  /tmp/Part_%d.part'%n)
	
	
	return ig.VertexClustering(g,c)

