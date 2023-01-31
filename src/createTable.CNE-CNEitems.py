#! /usr/bin/env python

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import re
import sys
import bisect
import operator
import itertools
import collections

import utils.myTools
from utils.myTools import file
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("KX_HDMC_blocs.list",file), ("KX_HDMC_especes.list",file)], \
	[("ancGenesFiles",str,""), ("genesFiles",str,""), ("syntenyCutoff",int,5), ("output",str,""), ("cneTarget",str,""), ("chromID", str, ""), ("geneID", str, "")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

dicScientificName = {
	'Cow': 'Bos taurus', \
	'Mouse': 'Mus musculus', \
	'hg19': 'Homo sapiens', \
	'Dog': 'Canis familiaris', \
	'Chicken': 'Gallus gallus', \
	'Zebrafish': 'Danio rerio' \
}

# Load chromosome id
print("Load chrom id...", file=sys.stdout)
f = utils.myFile.openFile(arguments["chromID"], "r")
dicChrom = collections.defaultdict(lambda: collections.defaultdict())
for line in f:
	t = line.split()
	# species id - chromname - chromid
	dicChrom[t[0]][t[1]] = t[2]

# Load gene id
print("Load gene id ...", file=sys.stdout)
f = utils.myFile.myTSV.readTabular(arguments["geneID"], [int,str,str,str,str,str,str,str,str,str,str])
dicGeneID = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
for line in f:
	# species id - chrom id - gene name - gene id
	# take only gene name outside the braces
	if "(" in line[2]:
		matchobj = re.match('(.*)\((.*)\)', line[2])
		if matchobj:
			dicGeneID[line[1]][line[3]][matchobj.group(1).strip()] = line[0]
	else:
		dicGeneID[line[1]][line[3]][line[2].strip()] = line[0]
print("OK", file=sys.stdout)

# Chargement des genomes des especes modernes si pas encore fait
dicGenomes = {}
posCNE1 = {}
posCNE2 = {}
def loadSpecies(species_name):
	if species_name in dicGenomes:
		return

	# Chargement du genome
	lst = collections.defaultdict(list)
	for gene in utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[species_name]):
		lst[str(gene.chromosome)].append( (gene.beginning,gene.end,gene.names[0]) )

	dicGenomes[species_name] = lst
	posCNE1[species_name] = collections.defaultdict(list)
	posCNE2[species_name] = collections.defaultdict(list)



# Chargement des CNE
CNE = {}
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_blocs.list"], [str,int,int,float,str]):
	# Les infos generales du CNE
	CNE[line[0]] = (line[1:],[])

levels = collections.defaultdict(set)
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_especes.list"], [str,str,int,int,str,str,int,str]):
	# Les infos de chaque CNE pour chaque espece
	species_name = dicScientificName[line[0]]
	loadSpecies(species_name)

	# Il faut que le chromosome de l'alignement multiple soit dans le genome de reference
	if line[1] in dicGenomes[species_name]:
		# On recherche tout de suite l'index de position sur le chromosome
		i = bisect.bisect_left(dicGenomes[species_name][line[1]], (line[2],line[3]))
	else:
		i = 0
		#print >> sys.stderr, "no chromosome *%s* for *%s*" % (line[1],species_name)
	CNE[line[5]][1].append( (species_name,line[1],i,line) )
	levels[line[6]-1].add(species_name)

# Chargement des genes ancestraux
ancNames = [(int(x),phylTree.lastCommonAncestor(list(levels[x]))) for x in levels]
ancNames.sort()
ancGenes = [utils.myGenomes.Genome(arguments["ancGenesFiles"] % x.replace(" ", ".")).dicGenes for (_,x) in ancNames]

# Filtre des CNE sur la syntenie
goodCNE = []
for (cne_id,(_,cne_items)) in CNE.items():
	found = []
	ancGenesID = []
	# Parcours de chaque espece
	for (species_name,chrom,pos,line) in cne_items:
		lst = set()
		refGenes = ancGenes[int(line[6])-1]
		# On recupere les genes voisins
		for g in dicGenomes[species_name][chrom][max(0,pos-arguments["syntenyCutoff"]):pos+arguments["syntenyCutoff"]+1]:
			if g[-1] in refGenes:
				lst.add(refGenes[g[-1]])
		ancGenesID.append( (species_name,lst) )
	goodSpecies = set()
	# Comparaison des especes
	for ((e1,l1),(e2,l2)) in itertools.combinations(ancGenesID, 2):
		# Si deux especes ont un gene syntenique, on les garde
		if len(l1.intersection(l2)) > 0:
			goodSpecies.add(e1)
			goodSpecies.add(e2)
	goodCNE.append( (cne_id,goodSpecies) )

####### special: for only human
# Load target genes
print("Load cne target ...", file=sys.stdout)
f = utils.myFile.myTSV.readTabular(arguments["cneTarget"], [str,str,str])
cne_target = {}
human_id = str(phylTree.indNames["Homo sapiens"])
for t in f:
	if t[0] not in CNE:
		continue
	if t[1] not in dicChrom[human_id]:
		continue
	
	chrom_id = dicChrom[human_id][t[1]]
	targets = t[2].split()
	tg_ids = ""
	if targets != "-":
		for tg in targets:
			if tg in dicGeneID[human_id][str(chrom_id)]:
				gene_id = dicGeneID[human_id][str(chrom_id)][tg]
				tg_ids += str(gene_id) + " "
	# cne name - list(targets names)
	if tg_ids != "":
		cne_target[t[0]] = tg_ids.strip()
	else:
		cne_target[t[0]] = None
print("OK", file=sys.stdout)

dataCNE = {}
f = utils.myFile.openFile(arguments["output"] % "CNE", "w")

# Parcours des CNE restants
for (realID,(cne_id,cne_spec)) in enumerate(goodCNE):

	# Parcours des especes
	for (species_name,chrom,pos,line) in CNE[cne_id][1]:

		level = line[6]
		syntenic = int(species_name in cne_spec)
		strand = int(line[4] + '1')
		intronic = -1

		# Si on a une info sur le chromosome
		#if line[1] in dicGenomes[species_name]:
		#	i = bisect.bisect_left(dicGenomes[species_name][line[1]], (line[2],line[3]))
		#	if (i >= 1) and (dicGenomes[species_name][line[1]][i-1][1] > line[2]):
		#		intronic = i-1
		#else:
		#	i = 0
		if pos >= 1 and (dicGenomes[species_name][line[1]][pos-1][1] > line[2]):
			intronic = pos-1
		
		spec_id = phylTree.indNames[species_name]
		# On enregistre la position
		if line[1].strip() in dicChrom[str(spec_id)]:
			chrom_id = dicChrom[str(spec_id)][line[1]]
		else:
			chrom_id = line[1]

		if line[5] in cne_target:
			target = cne_target[line[5]]
		else:
			target = None
		dataCNE[(realID,species_name)] = [realID,spec_id, chrom_id,None,None,None,None,None,None,  line[2],line[3],strand,intronic,syntenic,level,line[7],line[5], target]
		posCNE1[species_name][line[1]].append( (pos,line[2],realID) )
		if intronic >= 0:
			posCNE2[species_name][line[1]].append( (pos-1,line[2],realID) )
		else:
			posCNE2[species_name][line[1]].append( (pos,line[2],realID) )

	# On imprime les infos generales
	print(utils.myFile.myTSV.printLine((realID,) + CNE[cne_id][0] + (level,)), file=f)
f.close()


def writePositions(pos, indexP, index, indexR, reverse):

	# Parcours des CNE par espece
	for (species_name,x) in pos.items():
		# Puis par chromosome
		for (chrom,lst) in x.items():
			lst.sort(reverse=reverse)
			# Puis par espace intergenique
			for (pos0,y) in itertools.groupby(lst, operator.itemgetter(0)):
				# On positionne les CNE dans le bon ordre et a intervalle regulier
				l = list(y)
				for (i,y) in enumerate(l):
					dataCNE[(y[2],species_name)][indexP] = pos0 - 1 + int(reverse)
					dataCNE[(y[2],species_name)][index] = i
					dataCNE[(y[2],species_name)][indexR] = len(l)

writePositions(posCNE1, 3, 4, 5, False)
writePositions(posCNE2, 6, 7, 8, True)

f = utils.myFile.openFile(arguments["output"] % "CNE_items", "w")
for x in dataCNE.values():
	print(utils.myFile.myTSV.MySQLFileWriter(x), file=f)
f.close()

