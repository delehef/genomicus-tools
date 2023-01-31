#! /usr/bin/env python

__doc__ = """
	Cree la base Genomicus (tables Tree & Orthologs)
	Necessite les genomes modernes pour rajouter les genes qui ne sont pas dans des arbres
	Necessite les genes ancestraux pour creer le dictionnaire gene_name <-> gene_id
"""

import sys
import itertools
import collections

import utils.myFile
import utils.myTools
from utils.myTools import file
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree



# Lit chaque arbre de proteines et le stocke dans la base
# Definit la liste des genes de chaque espece (ancestrale/moderne)
###################################################################
def storeProteinTrees():

	def doProt(node, parent_name, parent_dup, parent_id, parent_distance, last_orth_id):

		global gene_id, nextLeft, orth_id

		dup = tree.info[node]["Duplication"]

		if last_orth_id is None:
			last_orth_id = orth_id
			orth_id += 1

		# Les ancetres qui s'intercalent entre le precedent et le courant
		newAnc = tree.info[node]['taxon_name']
		if parent_name != None:
			links = phylTree.dicLinks[parent_name][newAnc]
			if not parent_dup:
				links = links[1:]
		else:
			links = [newAnc]

		# L'id du gene courant
		curr_id = gene_id
		gene_id += len(links)
		newparent_id = gene_id-1

		# La borne gauche evolue
		left0 = nextLeft
		nextLeft += len(links)

		ids = []
		allids = []
		if node in tree.data:

			names = []
			# Appels recursifs sur les fils
			for (i,(fils,d)) in enumerate(tree.data[node]):
				# Appel
				if dup >= 2:
					neworth = orth_id
					allids.append(neworth)
					orth_id += 1
				else:
					neworth = last_orth_id
				x = doProt(fils, newAnc, dup>=2, newparent_id, d, neworth)
				names.extend( x[0] )
				ids.append( x[1] )
				allids.extend( x[1] )

		else:
			# On garde le gene en memoire pour plus tard
			allGenes[newAnc].add(tree.info[node]['gene_name'])
			names = [tree.info[node]['gene_name']]

		# Il faut enregistrer les genes des ancetres intermediaires
		for (i,anc) in enumerate(links):

			if anc == newAnc:
				tmpdist = parent_distance
				tmpdup = dup
			else:
				tmpdist = 0
				tmpdup = 0

			# Il y a un gene uniquement si ce n'est pas une duplication
			if tmpdup < 2:
				#print >> sys.stderr,anc
				if anc in phylTree.listSpecies:
					assert len(names) == 1
					print(utils.myFile.myTSV.MySQLFileWriter([anc, names[0], curr_id, root_id, last_orth_id]), file=outputFiles["dicGeneID"])
				else:
					print(utils.myFile.myTSV.MySQLFileWriter([anc, dicAncGeneNames[anc][frozenset(names)], curr_id, root_id, last_orth_id]), file=outputFiles["dicGeneID"])
			print(utils.myFile.myTSV.MySQLFileWriter((curr_id, left0+i, nextLeft+len(links)-2*i-1, parent_id, tmpdist, tmpdup, phylTree.indNames[anc], root_id)), file=outputFiles["Tree"])

			(curr_id, parent_id, parent_name) = (curr_id+1, curr_id, anc)
			nextLeft += 1
		assert parent_id == newparent_id

		if dup < 2:
			# Si pas de duplication: les fils sont orthologues entre eux
			for (l1,l2) in itertools.combinations(ids, 2):
				for (x,y) in itertools.product(l1, l2):
					print(utils.myFile.myTSV.MySQLFileWriter([x,y]), file=outputFiles["Orthologs"])
					print(utils.myFile.myTSV.MySQLFileWriter([y,x]), file=outputFiles["Orthologs"])
		else:
			# Si duplication, il faut creer les liens peres/fils
			for x in allids:
				print(utils.myFile.myTSV.MySQLFileWriter([x,last_orth_id]), file=outputFiles["Orthologs"])
				print(utils.myFile.myTSV.MySQLFileWriter([last_orth_id,x]), file=outputFiles["Orthologs"])

		return (names,allids)

	# Parcours du fichier d'arbres
	for tree in utils.myProteinTree.loadTree(arguments["proteinTree"]):
		root_id = gene_id
		doProt(tree.root, None, None, None, None, None)
	print(file=sys.stderr)



# Enregistre un genome moderne
################################
def storeModernGenome(esp):
	print("Inserting modern genome", esp, "...", end=' ', file=sys.stderr)
	global gene_id, nextLeft

	# Lecture du genome
	f = utils.myFile.openFile(arguments["genesFiles"] % phylTree.fileName[esp], "r")
	for l in f:
		name = l[:-1].split("\t")[4]
		if name not in allGenes[esp]:
			curr_id = gene_id
			gene_id += 1
			print(utils.myFile.myTSV.MySQLFileWriter((curr_id, nextLeft, nextLeft+1, None, None, 0, phylTree.indNames[esp], curr_id)), file=outputFiles["Tree"])
			print(utils.myFile.myTSV.MySQLFileWriter((esp, name, curr_id, curr_id, None)), file=outputFiles["dicGeneID"])
			nextLeft += 2
	f.close()
	print("OK", file=sys.stderr)



# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("proteinTree",file)], \
	[("genesFiles",str,""), ("ancGenesFiles",str,""), ("outputFile",str,"dump/%s.txt")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

outputFiles = {}
for x in ["Tree", "Orthologs", "dicGeneID"]:
	outputFiles[x] = utils.myFile.openFile(arguments["outputFile"] % x, 'w')

dicAncGeneNames = {}
for anc in phylTree.listAncestr:
	dic = {}
	genome = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])
	for gene in genome:
		dic[frozenset(gene.names[1:])] = gene.names[0]
	dicAncGeneNames[anc] = dic

nextLeft = 0
gene_id = 0
orth_id = 0
allGenes = collections.defaultdict(set)
storeProteinTrees()

for x in range(orth_id):
	print(utils.myFile.myTSV.MySQLFileWriter([x,x]), file=outputFiles["Orthologs"])

for esp in phylTree.listSpecies:
	storeModernGenome(esp)

for f in outputFiles.values():
	f.close()

