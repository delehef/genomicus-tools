#! /usr/bin/env python

__doc__ = """
	Cree la base Genomicus (table SpeciesTree)
"""

import sys

import utils.myTools
from utils.myTools import file
import utils.myPhylTree


arguments = utils.myTools.checkArgs( [("phylTree.conf",file)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

left = 0

# Utilise left comme gauche
# Met a jour left comme la prochaine valeur libre
def do(node):
	global left
	
	thisleft = left
	left += 1
	
	if node in phylTree.items:
		for (anc,_) in phylTree.items[node]:
			do(anc)
	
	print(utils.myFile.myTSV.MySQLFileWriter([phylTree.indNames[node], thisleft, left, phylTree.indNames[phylTree.parent[node][0]] if node in phylTree.parent else None]))
	
	left += 1

do(phylTree.root)

