#!/usr/bin/python

import numpy as np
from matplotlib.pylab import flatten
from collections import OrderedDict,Counter
import os
import igraph as ig
import pandas as pd

def Louvain(g,n=0,part=None):
	"It wants mandatory a igraph object (g), an integer number (n) that must be different if you computer in paraller different network (optional), ,Optionally you can provide a seed partition (part) the must be different from the trivail one i.e. each node in a separate partition. It returns the membership of the nodes"


	to_ck = [e.tuple for e in g.es]

	Od= OrderedDict.fromkeys(list(flatten(to_ck)))

	Od = dict(zip(Od.keys(),range(len(Od))))

	x = [(Od[a],Od[b],w) for (a,b),w in zip(to_ck,g.es["weight"])]

	x = pd.DataFrame(x)

	x.to_csv('/tmp/graph_%d.txt'%n,sep=' ',index=False,header=False)
	
	
	os.system('convertX  -i /tmp/graph_%d.txt -o /tmp/graph_%d.bin -w /tmp/graph_%d.weights'%(n,n,n))
	
	if part!=None:
		Part = np.zeros(g.vcount()).astype(int)
		for i,p in enumerate(part):
			Part[Od[i]] = p

		Part = zip(range(g.vcount()),Part)

		pd.DataFrame(Part).to_csv('/tmp/Part_%d.part'%n,sep=' ',header=False,index=False)


	if part!=None:
		os.system('loivainX /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights -p /tmp/Part_%d.part > /tmp/graph_%d.tree'%(n,n,n,n))
	else:
		os.system('louvainX /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights  > /tmp/graph_%d.tree'%(n,n,n))

	lv = int(os.popen('hierarchyX /tmp/graph_%d.tree'%n).read().split('\n')[0].split(': ')[1])-1

	os.system('hierarchyX /tmp/graph_%d.tree -l %d > /tmp/graph_node2comm_level2_%d'%(n,lv,n))
	
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

