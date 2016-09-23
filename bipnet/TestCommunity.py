#!/usr/bin/python

import sys
from multiprocessing import Pool
import pandas as pd

import Loivain_Wrapper as loi
import Validate as val
import metric as mtr
from functools import partial
import numpy as np
import igraph as ig


def Mod(g,com):
    return g.modularity(com,weights=g.es["weight"])

def Estimate_Partition(g,ncpu,Ntest=1000):
	
	func = partial(loi.LoivainOriginal,g)
	#cn = [(i,np.random.randint(0,1e9)) for i in xrange(Ntest)]
	p = Pool(ncpu)
	c = p.map(func,range(Ntest))
	p.close()
	
	return c

def Stability_and_BestP(g,c,ncpu):
	func = partial(Mod,g)

	p = Pool(ncpu)
	Modx =p.map(func,c)
	p.close()
	
	c = np.array(c)
	
	sel = np.percentile(Modx,99)
	Modx = np.array(Modx)
	SEL = c[Modx>=sel]
	
	ari = [ig.compare_communities(SEL[i],SEL[j],method='ari') for i in xrange(len(SEL)) for j in xrange(i)]
	
		
	if len(ari)==0:
		ari_g = np.nan,np.nan,np.nan
		#best_g = c[Modx==max(Modx)][0]
		
	else:	
		ari_g = np.mean(ari),min(ari),max(ari) 
	
	best_g = c[Modx==max(Modx)][0]
	
	
	
	return best_g,ari_g

def TestNet(g,ncpu,Ntest=1000):
	
	c = Estimate_Partition(g,ncpu,Ntest)
	best_g,ari_g = Stability_and_BestP(g,c,ncpu)
	
	return best_g,ari_g
	

def write_stats(file_w, selfAri, Ass, Awi):
	
	with open(file_w+'.stats','w') as fw:
		
		fw.write('------- SELF ARI -----------\n')
		fw.write('FULL: %lf\t(%lf, %lf)\n'%selfAri['full'])
		fw.write('BONF: %lf\t(%lf, %lf)\n'%selfAri['bonf'])
		fw.write('FDR: %lf\t(%lf, %lf)\n'%selfAri['fdr'])
		
		fw.write('\n------- Assortativity Coeff -----------\n')
		fw.write('FULL: %lf\n'%Ass['full'])
		fw.write('BONF: %lf\n'%Ass['bonf'])
		fw.write('FDR: %lf\n'%Ass['fdr'])
		
		fw.write('\n------- AWI -----------\n')
		fw.write('BONF->FULL: %lf\n'%Awi['full'])
		fw.write('FDR->FULL: %lf\n'%Awi['bonf'])
		fw.write('BONF->FDR: %lf\n'%Awi['fdr'])
		
	return
		
def write_membership(file_r,g,best_g,best_gbonf,best_gfdr):
	
	cl = ['Name','FULL','BONF','FDR']
	x = pd.DataFrame(zip(g.vs["name"],best_g,best_gbonf,best_gfdr))
	x.to_csv(file_r+'_membership.dat',index=False,sep='\t')
	
	return
		
def Test_Bipartite(gb,ncpu=1,side=0,alpha=0.01):
	
	'Statitical Validation'
	g,gbonf,gfdr = val.SVN(gb,which=bool(side),alpha=alpha)
	
	proj = {'full': g,'bonf': gbonf,'fdr': gfdr}
	
	'Estimate BestPartition and Stability'
	best_g, ari_g = TestNet(g,ncpu)
	best_gbonf, ari_gbonf = TestNet(gbonf,ncpu)
	best_gfdr, ari_gfdr = TestNet(gfdr,ncpu)
	
	best = {'full': best_g, 'bonf': best_gbonf, 'fdr': best_gfdr }
	
	selfAri = {'full':ari_g,'bonf': ari_gbonf,'fdr': ari_gfdr }
	
	'Estimate the Assortativity coefficient'
	ass_g = mtr.AssortativeCoefficient(g,best_g)
	ass_gbonf = mtr.AssortativeCoefficient(gbonf,best_gbonf)
	ass_gfdr = mtr.AssortativeCoefficient(gfdr,best_gfdr)
	
	Ass = {'full': ass_g,'bonf': ass_gbonf,'fdr': ass_gfdr }
	
	'Estimate di AWI of the bestpartitions'
	awi_gbonf_g = mtr.awi(best_gbonf,best_g)
	awi_gfdr_g = mtr.awi(best_gfdr,best_g)
	awi_gbonf_gfdr = mtr.awi(best_gbonf,best_gfdr)
	
	Awi = {'bonf_full': awi_gbonf_g,'fdr_full': awi_gfdr_g,'bonf_fdr': awi_gbonf_gfdr }

	return {'selfARI': selfAri,'Ass.Coef': Ass,'AWI':Awi,'BestPartition': best,'projection':proj}

if __name__=='__main__':
	
	file_r = sys.argv[1]
	side = sys.arv[2]
	alpha = sys.arv[3]
	
	ncpu = sys.argv[4]
			
	gb = val.get_BipartiteFromFile(file_r)
		
	x = Test_Bipartite(gb,side,alpha)
	selfAri,Ass,Awi,best,proj = x['selfARI'],x['Ass.Coef'],x['AWI'],x['BestPartition'],x['projection']
	
	g,gbonf,gfdr = proj['full'],proj['bonf'],proj['fdr']
	best_g,best_gbonf,best_gfdr = best['full'],best['bonf'],best['fdr']
	
	'Write the output'
	write_stats(file_r, selfAri, Ass, Awi)
	write_membership(file_r,g,best_g,best_gbonf,best_gfdr)
	
	val.write_validation(g,file_r+'_full_proj.dat')
	val.write_validation(gbonf,file_r+'_bonf_proj.dat')
	val.write_validation(gfdr,file_r+'_fdr_proj.dat')

