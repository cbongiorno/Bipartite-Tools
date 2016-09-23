#!/usr/bin/python

import sys
import os

DIRs = os.getcwd()
DIRs = DIRs + "/bipnet/gen-louvain"


os.chdir(DIRs)
os.system('make clean')
os.system('make')

os.system('cp louvain /usr/local/bin/louvainX')
os.system('cp convert /usr/local/bin/convertX')
os.system('cp hierarchy /usr/local/bin/hierarchyX')

print "Now run:",'python setup.py install'



