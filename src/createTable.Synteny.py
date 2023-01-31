#! /usr/bin/env python

__doc__ = """
	Cree la base Genomicus (table Synteny)
	Necessite la table Gene pour positionner les blocs sur les chromosomes
	Lit les fichiers de diagonales + genes ancestraux
"""

import sys
import collections

import utils.myTools
from utils.myTools import file
import utils.myGenomes
import utils.myPhylTree



def storeAncDiags(anc):

	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])
	genome = utils.myGenomes.Genome(arguments["diags"] % phylTree.fileName[anc], ancGenes=ancGenes)
	global synt_id

	# Parcours des chromosomes
	print("Inserting ancestral diags", anc, "...", end=' ', file=sys.stderr)
	for (chrom,l) in genome.lstGenes.items():
		if len(l) == 1:
			continue
		trans = [dicPos[phylTree.indNames[anc]][gene.names[0]] for gene in l]
		transS = sorted(trans)
		# Tout le monde est sur un seul chromosome
		assert transS[0][0] == transS[-1][0]
		assert transS[0][1] + len(trans) - 1 == transS[-1][1]
		if trans != transS:
			assert list(reversed(trans)) == transS
		print(utils.myFile.myTSV.MySQLFileWriter((synt_id, phylTree.indNames[anc], trans[0][0], transS[0][1], transS[-1][1])))
		synt_id += 1
	print("OK", file=sys.stderr)


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("GeneTable",file)], \
	[("ancGenesFiles",str,""), ("diags",str,"")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Lecture de la table des genes
print("Loading table 'Gene' ...", end=' ', file=sys.stderr)
dicPos = collections.defaultdict(dict)
for line in utils.myFile.myTSV.readTabular(arguments["GeneTable"], [str]*11):
	dicPos[int(line[1])][line[2]] = (line[3],int(line[4]))
print("OK", file=sys.stderr)

synt_id = 0
for anc in phylTree.listAncestr:
	storeAncDiags(anc)



