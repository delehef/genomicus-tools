#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright  2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
        Convert gene trees from NHX to my protTree forest
"""

import sys

import utils.myProteinTree as myProteinTree
import utils.myTools as myTools
from utils.myTools import file

arguments = myTools.checkArgs( [("tree",file)], [], __doc__)

for tree in myProteinTree.loadTree(arguments["tree"]):
    myProteinTree.ProteinTree(tree.data, tree.info, tree.root).printTree(sys.stdout)
