import bipnet

SA,SB,q,pc,pr = 50,50,50,0.8,0.7

gb,O = bipnet.CREATE_NET([SA]*q,[SB]*q,[pc]*q,pr)

# Statistical Validation
g,gbonf,gfdr = bipnet.SVN(gb,alpha=0.01)

# community detection
c = bipnet.Louvain(g)
cf = bipnet.Louvain(gfdr)

print "AWI - Full projection"
print bipnet.awi(c.membership, O)
print "-----------"
print "AWI - FDR projection"
print bipnet.awi(cf.membership, O)
