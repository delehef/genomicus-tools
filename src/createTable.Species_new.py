#! /usr/bin/env python

__doc__ = """
	Cree la base Genomicus (table Species)
"""
import re
import sys
import utils.myFile
from utils.myTools import file
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("genome_db.txt",file)], [], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Infos supplementaires d'assemblage et de genebuild
dicEsp = {}
for t in utils.myFile.myTSV.readTabular(arguments["genome_db.txt"], (str,str,str,str,str,str,str,str,str,str,str,str,str)):
	dicEsp[t[2]] = (t[1],t[3],t[6],t[5])

print(dicEsp, file=sys.stderr)

# Les groupes pour le menu
allgrp = []
todo = set(phylTree.indNames)


# Ajoute un nouveau groupe avec toutes les especes descedantes de node, qu'on n'a pas encore vues
def addModernSpecies(node, name):
	grp = todo.intersection(phylTree.species[node])
	todo.difference_update(grp)
	allgrp.append( (node,name,grp) )

addModernSpecies("Primates", "Primates")
addModernSpecies("Euarchontoglires", "Rodents etc.")
addModernSpecies("Laurasiatheria", "Laurasiatherias")
addModernSpecies("Xenarthra", "Xenarthras etc.")
addModernSpecies("Afrotheria", "Afrotherias etc.")
addModernSpecies("Mammalia", "Marsupials &amp; Monotremes")
addModernSpecies("Sauria", "Birds &amp; Reptiles")
addModernSpecies("Tetrapoda", "Tetrapodes")
#addModernSpecies("Actinopterygii", "Fish")
# Avant ou apres l'insertion d'amphioxus
if "Deuterostomia" in phylTree.items:
	addModernSpecies("Deuterostomia", "Other chordates")
else:
	addModernSpecies("Chordata", "Other chordates")
addModernSpecies(phylTree.root, "Other eukaryotes")


# Ajoute un nouveau groupe avec tous les ancetres descendants de node, qui n'ont pas encore ete vus
def addAncestors(node, name):
	grp = todo.intersection(phylTree.allDescendants[node])
	todo.difference_update(grp)
	allgrp.append( (node,name,grp) )

addAncestors("Primates", "Ancestors in Primates")
addAncestors("Glires", "Ancestors in Rodents etc.")
addAncestors("Laurasiatheria", "Ancestors in Laurasiatherias")
addAncestors("Mammalia", "Ancestors in mammals")
#addAncestors("Actinopterygii", "Ancestors in fish")
addAncestors("Sauria", "Ancestors in Birds &amp; Reptiles")
addAncestors("Vertebrata", "Ancestors in vertebrates")
addAncestors(phylTree.root, "Other ancestors")


# Impression finale des groupes
maxid = len(phylTree.allNames)
for (i,(node,name,lst)) in enumerate(allgrp):

	# Entete
	res = [maxid+i, None, None, None, name, None,None, None, 100*i, 0]
	print(utils.myFile.myTSV.MySQLFileWriter(res))

	# Contenu
	for (j,esp) in enumerate(lst):
		print(j, esp, file=sys.stderr)
		vers = "/".join(dicEsp[esp][1:]) if esp in dicEsp else ""
		print(vers, file=sys.stderr)
		othernames = [x for x in phylTree.commonNames[esp] if (x != esp) and (type(x) != int) and not(re.search("_",x))]
		print(othernames, file=sys.stderr)
		if len(othernames) == 0:
			scientific_name = None
			common_name = esp
		else:
			scientific_name = esp
			common_name = othernames[0]
		res = [phylTree.indNames[esp], scientific_name, phylTree.ages[esp], vers, common_name, None, None, None if esp in phylTree.items else esp.replace(' ', '_'), 100*i+j+1, int(esp in phylTree.lstEsp2X), int(dicEsp[esp.lower().replace(' ', '_')][-1]) if esp.lower().replace(' ', '_') in dicEsp else 1]
		print(utils.myFile.myTSV.MySQLFileWriter(res))


