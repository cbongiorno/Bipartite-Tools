import igraph as ig
import numpy as np
from matplotlib.pylab import flatten
from collections import OrderedDict,Counter
import os


def TestVar(part,vm):
    p2 = part.copy()
    np.random.shuffle(p2)
    return ig.compare_communities(part,p2,method='vi')>vm

def Lypou(graph,n,seed_rand,seed=None,vm=None,var=True,cmpOut=False):
    'Devi corregere il seed della random che manca'
    
    if seed!=None:
        seed_var = seed.copy()

    diff_seed = 0
    s = Counter(seed)
 
    while var:
        sx = np.where(np.array([s[a]  for a in seed])>1)[0]
        a = np.random.randint(0,len(seed))
        b = np.random.choice(sx)
        if seed[a]==seed[b]: continue
        seed_var[a],seed_var[b] = seed_var[b],seed_var[a]
        diff_seed = ig.compare_communities(seed_var,seed,method='vi')
        if  diff_seed > vm: break
    
    if seed!=None:
        comm = LovivenOriginal(graph,n,seed_rand,seed)
        comm_var = LovivenOriginal(graph,n,seed_rand,seed_var)
    else:
        comm = LovivenOriginal(graph,n,seed_rand)
        comm_var = LovivenOriginal(graph,n,seed_rand)
     
    diff_out = ig.compare_communities(comm,comm_var,method='vi')
    if cmpOut==False:
        return (diff_seed,diff_out)
    if cmpOut==True:
        return (diff_seed,diff_out,seed,seed_var,comm,comm_var)
        
        

def LoivenOriginal(g,n,seed,part=None):
	DIRs = 'Community_latest/'
	to_ck = [e.tuple for e in g.es]

	Od= OrderedDict.fromkeys(list(flatten(to_ck)))

	Od = dict(zip(Od.keys(),range(len(Od))))

	x = [(Od[a],Od[b],w) for (a,b),w in zip(to_ck,g.es["weight"])]

	x = pd.DataFrame(x)

	x.to_csv('/tmp/graph_%d.txt'%n,sep=' ',index=False,header=False)


	os.system(DIRs+'convert  -i /tmp/graph_%d.txt -o /tmp/graph_%d.bin -w /tmp/graph_%d.weights'%(n,n,n))

	if part!=None:
		Part = np.zeros(g.vcount()).astype(int)
		for i,p in enumerate(part):
			Part[Od[i]] = p

		Part = zip(range(g.vcount()),Part)

		pd.DataFrame(Part).to_csv('/tmp/Part_%d.part'%n,sep=' ',header=False,index=False)


	if part!=None:
		os.system(DIRs+'community %d /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights -p /tmp/Part_%d.part > /tmp/graph_%d.tree'%(seed,n,n,n,n))
	else:
		os.system(DIRs+'community %d /tmp/graph_%d.bin -l -1 -w /tmp/graph_%d.weights  > /tmp/graph_%d.tree'%(seed,n,n,n))

	print DIRs+'hierarchy /tmp/graph_%d.tree'%n
	lv = int(os.popen(DIRs+'hierarchy /tmp/graph_%d.tree'%n).read().split('\n')[0].split(': ')[1])-1


	os.system(DIRs+'hierarchy /tmp/graph_%d.tree -l %d > /tmp/graph_node2comm_level2_%d'%(n,lv,n))

	cf = pd.read_csv('/tmp/graph_node2comm_level2_%d'%n,header=None,sep=' ')[1]

	Odrv = dict(map(list,map(reversed,Od.items())))

	c = np.array([-1]*g.vcount())
	if cf[0]!='of':
		for i,a in enumerate(cf):
			c[Odrv[i]] = a

	iso = np.where(c==-1)[0]
	c[iso] = range(max(c)+1,max(c)+len(iso)+1)

	c = c.astype(int)


	os.system('rm /tmp/graph_%d.bin /tmp/graph_%d.tree /tmp/graph_%d.weights /tmp/graph_node2comm_level2_%d'%(n,n,n,n))
	if part!=None:
		os.system('rm  /tmp/Part_%d.part'%n)

	return c
