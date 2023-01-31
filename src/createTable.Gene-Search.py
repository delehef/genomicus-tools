#! /usr/bin/env python

__doc__ = """
	Cree la base Genomicus (table Gene)
	Necessite la liste d'association Gene <-> gene_id definie par les arbres
	Lit les fichiers de genomes (modernes / diags+genes ancestraux)
"""

import sys
import collections

import utils.myFile
import utils.myTools
from utils.myTools import file
import utils.myGenomes
import utils.myPhylTree
import multiprocessing
#from joblib import Parallel, delayed

def chromSize():
	f = utils.myFile.openFile(arguments["chromSize"], "r")
	chroms = collections.defaultdict(lambda: collections.defaultdict())
	currSpec = -1
	for line in f:
		if "####" in line:
			if line[4:].strip() in phylTree.officialName:
				currSpec = phylTree.indNames[phylTree.officialName[line[4:].strip()]]
			elif line[4:5].upper() + line[5:] in phylTree.officialName:
				currSpec = phylTree.indNames[phylTree.officialName[line[4:5].upper() + line[5:]]]
		else:
			chrom = line.split()
			chroms[currSpec][chrom[0].strip()] = chrom[1].strip()
	return chroms

# Enregistre un genome ancestral
#################################
def storeAncGenes(anc):
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])
	genome = utils.myGenomes.Genome(arguments["diags"] % phylTree.fileName[anc], ancGenes=ancGenes)

	spec_id = phylTree.indNames[anc]
	# Parcours des chromosomes
	print("Inserting ancestral genome", anc, "...", end=' ', file=sys.stderr)
	
	#################
	# Sorting chroms
	#################
	sortedChroms = sorted(list(genome.lstGenes.items()), key=lambda x: len(x[1]), reverse=True)
	count = 1
	for (x,y) in sortedChroms:
		if x != None:
			genome.lstGenes[count] = y
			count += 1

	count = 0
	for (chrom,l) in genome.lstGenes.items():
		if len(l) == 1:
			gene = l[0]
			(gene_id,root_id,orth_id) = dicGeneID[anc].pop(gene.names[0])
			#res = [gene_id, phylTree.indNames[anc], gene.names[0], None, chrom, 0, None, None, root_id, orth_id, None]
			res = [gene_id, spec_id, gene.names[0], None, chrom, 0, None, None, root_id, orth_id, None]
			print(utils.myFile.myTSV.MySQLFileWriter(res), file=fG)
			print(utils.myFile.myTSV.MySQLFileWriter([gene.names[0], gene_id, spec_id]), file=fS)
		else:
			countG = 0
			for (i,gene) in enumerate(l):
				(gene_id,root_id,orth_id) = dicGeneID[anc].pop(gene.names[0])
				res = [gene_id, spec_id, gene.names[0], count, i, gene.strand, gene.beginning, gene.end, root_id, orth_id, None]
				print(utils.myFile.myTSV.MySQLFileWriter(res), file=fG)
				print(utils.myFile.myTSV.MySQLFileWriter([gene.names[0], gene_id, spec_id]), file=fS)
				countG += 1
			
			resC = [spec_id, arguments["ancBlockName"] % chrom, count, None, countG]
			print(utils.myFile.myTSV.MySQLFileWriter(resC), file=fC)
			count += 1

	print("OK", file=sys.stderr)


# Enregistre un genome moderne
################################
def storeModernGenome(esp):
	spec_id = phylTree.indNames[esp]
	# Lecture du genome
	genome = utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[esp])
	
	print("Loading names & descr ", esp, "...", end=' ', file=sys.stderr)
	dicNames = {}
	for t in utils.myFile.myTSV.readTabular(arguments["namesFiles"] % phylTree.fileName[esp], (str,str)):
		if len(t[1]) > 0:
			#print >> sys.stderr, t[0]
			assert t[0] not in dicNames
			dicNames[t[0]] = t[1]
	dicDescr = {}
	for t in utils.myFile.myTSV.readTabular(arguments["descrFiles"] % phylTree.fileName[esp], (str,str)):
		if len(t[1]) > 0:
			assert t[0] not in dicDescr
			dicDescr[t[0]] = t[1]
	print("OK", file=sys.stderr)

	# Cherche un autre nom dans la table dicNames, et l'ecrit dans Search
	def newName(name, gene_id):
		print(utils.myFile.myTSV.MySQLFileWriter([name, gene_id, spec_id]), file=fS)
		if name in dicNames:
			print(utils.myFile.myTSV.MySQLFileWriter([dicNames[name], gene_id, spec_id]), file=fS)
			return "%s (%s)" % (dicNames[name],name)
		else:
			return name

	# Parcours des chromosomes
	print("Inserting modern genome", esp, "...", end=' ', file=sys.stderr)
	
	count = 0
	for (chrom,l) in genome.lstGenes.items():
		countG = 0
		for (i,gene) in enumerate(l):
			if gene.names[0] in dicGeneID[esp]:
				(gene_id,root_id,orth_id) = dicGeneID[esp].pop(gene.names[0])
				#print >> sys.stderr, gene_id,root_id,orth_id,gene.names[0]
				#res = [gene_id, spec_id, newName(gene.names[0], gene_id), count, i, gene.strand, gene.beginning, gene.end, root_id, orth_id, dicDescr.get(gene.names[0])]
				res = [gene_id, spec_id, newName(gene.names[0], gene_id), count, countG, gene.strand, gene.beginning, gene.end, root_id, orth_id, dicDescr.get(gene.names[0])]
				print(utils.myFile.myTSV.MySQLFileWriter(res), file=fG)
				countG += 1
		
		if str(chrom) in chromSize[spec_id]:
			size = chromSize[spec_id][str(chrom)]
		else:
			size = None
		resC = [spec_id, chrom, count, size, countG]
		print(utils.myFile.myTSV.MySQLFileWriter(resC), file=fC)
		count += 1

	if arguments["strict"]:
		assert len(dicGeneID[esp]) == 0, list(dicGeneID[esp].keys())
	elif len(dicGeneID[esp])  != 0:
		for (i,(name,(gene_id,root_id,orth_id))) in enumerate(dicGeneID[esp].items()):
			res = [gene_id, spec_id, newName(name, gene_id), None, i, 0, None, None, root_id, orth_id, dicDescr.get(name)]
			print(utils.myFile.myTSV.MySQLFileWriter(res), file=fG)
		print("OK (%d genes in trees, but not in genome)" % len(dicGeneID[esp]), file=sys.stderr)
	else:
		print("OK", file=sys.stderr)


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("dicGeneID",file)], \
	[("genesFiles",str,""), ("ancGenesFiles",str,""), ("diags",str,""), ("namesFiles",str,""), ("descrFiles",str,""), ("output",str,"dump/%s.txt"), ("ancBlockName",str,"block_%s"), ("chromSize",str,""),("strict",bool,False)], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

print("Loading gene_id dictionary ...", end=' ', file=sys.stderr)
dicGeneID = collections.defaultdict(dict)
for (esp,name,gene_id,root_id,orth_id) in utils.myFile.myTSV.readTabular(arguments["dicGeneID"], [str,str,int,int,str]):
	dicGeneID[esp][name] = (gene_id,root_id,orth_id)
print("OK", file=sys.stderr)

print("Loading chromosome size ...", end=' ', file=sys.stderr)
chromSize = chromSize()
print("OK.", file=sys.stderr)

fG = utils.myFile.openFile(arguments["output"] % "Gene", "w")
fS = utils.myFile.openFile(arguments["output"] % "Search", "w")
fC = utils.myFile.openFile(arguments["output"] % "Chromosome", "w")

for anc in phylTree.listAncestr:
	storeAncGenes(anc)
#	assert len(dicGeneID[anc]) == 0, dicGeneID[anc]

for esp in phylTree.listSpecies:
	storeModernGenome(esp)
	#assert len(dicGeneID[esp]) == 0, dicGeneID[esp]

fG.close()
fS.close()
fC.close()

