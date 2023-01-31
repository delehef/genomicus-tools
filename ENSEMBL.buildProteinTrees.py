#! /usr/bin/env python

__doc__ = """
	Corrige les arbres d'Ensembl en fonction du seuil minimal de duplication_score et de l'arbre des especes desire
		1: score par defaut (0 -> 1)
		2: coef multiplicateur d'un score reference 6X_species / all_species
		3: duplication_confidence_score calcule sur uniquement 6X_species
"""

import sys
import itertools
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree
from utils.myTools import file
sys.setrecursionlimit(20000) 

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("ensemblTree",file)], \
	[("cutoff",str,"-1"), ("defaultFamName",str,"FAM%08d"), ("scoreMethod",int,[1,2,3]), ("newNodeID",int,1000000000), ("recurs",bool,False)], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Limites automatiques de score de duplication
if arguments["scoreMethod"] in [1, 3]:
	def calc(anc, val):
		return val
elif arguments["scoreMethod"] == 2:
	def calc(anc, val):
		nesp = len(phylTree.species[anc])
		n2X = len(phylTree.lstEsp2X.intersection(phylTree.species[anc]))
		# La moitie des especes non 2X a vu la duplication (au minimum 1 espece)
		return round(max(1., val*(nesp-n2X)) / nesp, 3) - 2e-3

minDuplicationScore = {}
try:
	# Une limite pour tout le monde
	val = float(arguments["cutoff"])
	for anc in phylTree.listAncestr:
		minDuplicationScore[anc] = calc(anc, val)
except ValueError:
	f = utils.myFile.openFile(arguments["cutoff"], "r")
	for l in f:
		t = l.split()
		anc = phylTree.officialName[t[0]]
		minDuplicationScore[anc] = calc(anc, float(t[1]))
	f.close()
print("minDuplicationScore:", minDuplicationScore, file=sys.stderr)

# Les scores dans l'abre pour les especes modernes valent toujours 1, on doit toujours les accepter
for esp in phylTree.listSpecies:
	minDuplicationScore[esp] = 0


utils.myProteinTree.nextNodeID = arguments["newNodeID"]

@utils.myTools.memoize
def goodSpecies(anc):
	return phylTree.species[anc].difference(phylTree.lstEsp2X)

def alwaysTrue(tree, rnode):
	return True

def hasLowScore(tree, rnode):

	@utils.myTools.memoize
	def getSpeciesSets(node):
		if node in tree.data:
			return set().union(*(getSpeciesSets(x) for (x,_) in tree.data[node]))
		else:
			print(tree.info[node]["taxon_name"], file=sys.stderr)
			assert tree.info[node]["taxon_name"] in phylTree.listSpecies
			return set([tree.info[node]["taxon_name"]])

	if rnode not in tree.data:
		return False

	speciessets = [getSpeciesSets(x) for (x,_) in tree.data[rnode]]
	inters = set()
	for (s1,s2) in itertools.combinations(speciessets, 2):
		inters.update(s1.intersection(s2))
	all = set().union(*speciessets)
	anc = tree.info[rnode]["taxon_name"]

	if arguments["scoreMethod"] == 3:
		inters.intersection_update(goodSpecies(anc))
		all.intersection_update(goodSpecies(anc))
	#print >> sys.stderr,rnode
	return ((len(inters) == 0) and (minDuplicationScore[anc] == 0)) or (len(inters) < (minDuplicationScore[anc] * len(all)))


nbEdit = {"dubious": 0, "toolow": 0, "good": 0}

for (nb,tree) in enumerate(utils.myProteinTree.loadTree(arguments["ensemblTree"])):

	assert max(tree.info) < arguments["newNodeID"]

	# On trie les bonnes duplications des mauvaises
	################################################
	for (node,inf) in tree.info.items():
		print(node,inf, file=sys.stderr)
		if inf['Duplication'] != 0:
		
			if 'dubious_duplication' in inf:
				# On considere que les duplications 'dubious' ne sont pas valables pour une duplication
				assert inf['Duplication'] == 1
				del inf['dubious_duplication']
				nbEdit["dubious"] += 1

			if hasLowScore(tree, node):
				inf['Duplication'] = 1
				nbEdit["toolow"] += 1
			else:
				# Il faut la passer a 2 si le score est suffisant
				# Attention: pour les arbres d'Ensembl dont la racine est une duplication, celle-ci vaut 1 (parce qu'elle n'a pas d'outgroup)
				if inf['Duplication'] == 1:
					inf['Duplication'] = 3
				else:
					assert inf['Duplication'] in [2,3]
				nbEdit["good"] += 1

	tree.flattenTree(phylTree, True)
	#print >> sys.stderr,node,inf
	tree.rebuildTree(phylTree, hasLowScore if arguments["recurs"] else alwaysTrue)
	
	if "tree_name" not in tree.info[tree.root]:
		tree.info[tree.root]["tree_name"] = arguments["defaultFamName"] % nb

	tree.printTree(sys.stdout)

print(nbEdit, file=sys.stderr)

