
import sys
import collections

import utils.myTools

##########################################################
# Gere les arbres de proteines sous la forme (data,info) #
##########################################################

# Indique les noms des ancetres depuis le dernier lu
######################################################
def getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, isDuplication):

    # Les noeuds ou ecrire les familles
    if lastWrittenAnc == None:
        if previousAnc == None:
            toWrite = [newAnc]
        else:
            toWrite = phylTree.dicLinks[previousAnc][newAnc]
    else:
        toWrite = phylTree.dicLinks[lastWrittenAnc][newAnc][1:]

    root = False
    if lastWrittenAnc == None:
        if isDuplication:
            toWrite = []
        else:
            root = True
            assert len(toWrite) > 0
            toWrite = toWrite[-1:]
    else:
        if isDuplication:
            toWrite = toWrite[:-1]

    if len(toWrite) == 0:
        newLastWritten = lastWrittenAnc
    else:
        newLastWritten = toWrite[-1]

    return (toWrite,newLastWritten,root)


#
# Gere les arbres de proteines sous la forme (data,info)
##########################################################
class ProteinTree:

    def __init__(self, data=None, info=None, root=None):
        self.data = {} if data is None else data
        self.info = {} if info is None else info
        self.root = root


    def doBackup(self):
        self.backRoot = self.root
        self.backData = dict((node,values[:]) for (node,values) in self.data.items())
        self.backInfo = dict((node,values.copy()) for (node,values) in self.info.items())



    # Imprime l'arbre sous forme indentee avec les infos
    ######################################################
    def printTree(self, f, node=None):

        def rec(n, node):
            indent = "\t" * n
            # L'id du noeud
            print("%sid\t%d" % (indent, node), file=f)
            # Les infos
            print("%sinfo\t{%s}" % (indent, ", ".join(repr(key) + ": " + repr(value) for (key, value) in sorted(self.info[node].items()))), file=f)
            # Les enfants
            for (g,d) in self.data.get(node,[]):
                print("%s\tlen\t%f" % (indent, d), file=f)
                rec(n+1, g)

        rec(0, self.root if node is None else node)
        try:
            f.flush()
        except AttributeError:
            pass

    def printNewick(self, f, root=None, withDist=True, withTags=False, withAncSpeciesNames=False):
        NHX = {"Duplication": "D", "Bootstrap": "B", "taxon_name": "S", "duplication_confidence_score": "SIS", "dubious_duplication": "DD"}
        def rec(node):
            if node in self.data:
                return "(" + ",".join(
                        rec(g)
                        + ((":%f" % l) if withDist else "")
                        + ("[&&NHX:" + ":".join(("%s=%s" % (NHX[tag],self.info[g][tag])).replace(" ", ".") for tag in NHX if tag in self.info[g]) + "]" if withTags else "")
                        for (g,l) in self.data[node]
                ) + ")" + (self.info[node]["taxon_name"].replace(' ', '.') if withAncSpeciesNames and ("taxon_name" in self.info[node]) else '')
            else:
                return self.info[node]['gene_name'].split("/")[0]

        if root is None:
            root = self.root
        print(rec(root) + ("[&&NHX:" + ":".join(("%s=%s" % (NHX[tag],self.info[root][tag])).replace(" ", ".") for tag in NHX if tag in self.info[root]) + "]" if withTags else "") + ";", file=f)
        try:
            f.flush()
        except AttributeError:
            pass


    # Imprime l'arbre au format Newick
    ####################################
    def printNewickTree(self, f, node=None):
        genes = []
        def rec(node):
            if node not in self.data:
                genes.append(self.info[node]['gene_name'])
                return self.info[node]['gene_name']
            else:
                return "(" + ",".join([rec(x) + ":" + str(l) for (x,l)  in self.data[node]]) + ") " + info[node]['family_name']

        tr = rec(self.root if node is None else node)
        print(" ".join(genes), file=f)
        print(tr, ";", file=f)


    ######################################################################################
    # Compacte un arbre en supprimant les noeuds intermdiaires qui n'ont qu'un seul fils #
    ######################################################################################
    def compactTree(self, phylTree, node=None):

        def do(node):
            # Fin du process sur une feuille
            if node not in self.data:
                return False

            flag = False
            # Appels recursifs
            for (gg,_) in self.data[node]:
                flag |= do(gg)

            if len(self.data[node]) > 1:
                return flag

            # Edition du noeud courant
            (g,l) = self.data[node][0]
            if g in self.data:
                self.data[node] = [(gg,ll+l) for (gg,ll) in self.data[g]]
                del self.data[g]
            else:
                del self.data[node]
            self.info[node] = self.info[g]
            del self.info[g]
            return True

        return do(self.root if node is None else node)


    ###############################################################################################
    # Renomme tous les noeuds pour qu'ils correspondent aux ancetres communs de leurs descendants #
    ###############################################################################################
    def renameTree(self, phylTree, node=None):

        def do(node):
            # Fin du process sur une feuille
            if node not in self.data:
                return False

            flag = False
            # Appels recursifs
            for (g,_) in self.data[node]:
                flag |= do(g)

            # Renommage du noeud courant
            newName = phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )
            flag |= (self.info[node]['taxon_name'] != newName)
            self.info[node]['taxon_name'] = newName

            return flag

        return do(self.root if node is None else node)


    ###################################################################################################################
    # Aplatit un noeud et des descendants directs si ils representent le meme taxon et qu'il n'y a pas de duplication #
    # Par exemple:                                                                                                    #
    #  ((Eutheria1,Eutheria2)XA,(Eutheria3,Eutheria4)XB)XC se transforme en (Eutheria1,Eutheria2,Eutheria3,Eutheria4) #
    #    uniquement si XA, XB et XC sont des speciations                                                              #
    ###################################################################################################################
    def flattenTree(self, phylTree, rec, node=None):

        def do(node):

            # Fin du process sur une feuille
            if node not in self.data:
                return False
            assert len(self.data[node]) > 0
           # print >> sys.stderr, 'nooode', node, self.info[node], self.data[node]


            flag = False
            # Appels recursifs
            if rec:
                for (g,_) in self.data[node]:
                       # print >> sys.stderr, 'g', g, self.info[g]
                    flag |= do(g)


           # print >> sys.stderr, 'essai', node, self.info[node]['taxon_name']
            self.info[node]['taxon_name'] = phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )
            #print >> sys.stderr, 'new', node, self.info[node]['taxon_name']
            #print >> sys.stderr, "ca marche"
            # Si c'est une vraie duplication, on n'a plus rien a faire
            #print >> sys.stderr, 'dup', node, self.info[node]['Duplication']
            if self.info[node]['Duplication'] >= 2:
                   # print >> sys.stderr, "oooookkkkkk"
                return flag

            newData = []
            taxonName = self.info[node]['taxon_name']
            for (g,d) in self.data[node]:
                inf = self.info[g]
                # 2x le meme taxon et pas de duplication

                if (inf['taxon_name'] == taxonName) and (inf['Duplication'] < 2) and (g in self.data):

                    newData.extend([(g2,d+d2) for (g2,d2) in self.data[g]])
                    #print >> sys.stderr, 'NNNO', node, newData
                    del self.data[g]

                    self.info[node].update(self.info[g])
                    del self.info[g]
                    flag = True
                else:
                    newData.append( (g,d) )
                       # print >> sys.stderr, 'ELSE', node, newData

            assert len(newData) == len(set(g for (g,_) in newData)), newData
            assert len(newData) > 0, (node,self.data[node],newData)

            self.info[node]['Duplication'] = 0
            self.data[node] = newData
            #print >> sys.stderr, 'end', node, newData
            if len(self.data[node]) > 0:
                assert self.info[node]['taxon_name'] == phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )

            return flag

        return do(self.root if node is None else node)


    ###############################################################################
    # Redonne la topologie attendue a l'arbre (par rapport a l'arbre des especes) #
    #   Rassemble les noeuds equivalents sous le meme fils                        #
    ###############################################################################
    def rebuildTree(self, phylTree, hasLowScore, node=None):
           # print >> sys.stderr, 'REBUILDTREE'
        def do(node):
            #print >> sys.stderr, 'node1', node
            # Fin du process sur une feuille
            if node not in self.data:
                return False

            flag = False

            # On ne change les fils que si ce n'est pas une vraie duplication
            if self.info[node]['Duplication'] < 2:
                #print >> sys.stderr, 'pas true dup'
                # Le noeud sera une speciation sauf cas special ci-dessous
                self.info[node]['Duplication'] = 0

                # On redefinit les fils pour -notre- arbre phylogenetique en triant les enfants en paquets
                fils = collections.defaultdict(list)
                #print >> sys.stderr, 'fils', fils
                anc = self.info[node]['taxon_name']
                #print >> sys.stderr, 'anc', anc
                lfils = phylTree.items.get(anc, [])
                #print >> sys.stderr, 'lfils', lfils
                for (g,d) in self.data[node]:
                    #print >> sys.stderr, 'g', g
                    gname = self.info[g]['taxon_name']
                    for (a,_) in lfils:
                        #print >> sys.stderr, 'a', a
                        if phylTree.isChildOf(gname, a):
                            #print >> sys.stderr, 'ouiiiii', a
                            fils[a].append( (g,d) )
                            break
                    else:
                        #print >> sys.stderr,'else'
                        # Ne peut apparaitre que si g est le meme ancetre que node et que g est une duplication,
                        # ce qui a normalement implique Duplication >= 2, sauf pour les nouveaux noeuds crees
                        assert (gname == anc), "ERREUR: name!=anc [%s / %s / %s]" % (node, anc, gname)
                        # Le noeud courant sera donc un noeud de duplication
                        #print >> sys.stderr, 'noooooooooooooo', anc
                        self.info[node]['Duplication'] = 3
                        fils[anc].append( (g,d) )

                # len(fils):
                #  1 -> uniquement anc
                #  2 ou 3 parmi F1/F2/anc
                assert (len(fils) != 1) or (anc in fils), "ERREUR: 1=anc [%s / %s]" % (node, fils)
                assert (len(fils) <= (1+len(lfils))), "ERREUR: len>(1+nbFils) [%s / %s]" % (node, fils)

                todo = []

                if len(fils) > 1:
                    if anc in fils:
                        lst1 = fils.pop(anc)
                        lst2 = []
                        for tmp in fils.values():
                            lst2.extend(tmp)
                        items = [(anc,lst1), (anc,lst2)]
                    else:
                        items = list(fils.items())

                    newData = set()
                    for (anc,lst) in items:
                        if len(lst) == 1:
                            newData.add( lst[0] )
                        elif len(lst) > 1:
                            for (g,l) in self.data[node]:
                                if (g in self.data) and (self.data[g] == lst):
                                    newData.add( (g,l) )
                                    break
                                if g in self.data:
                                    assert sorted(self.data[g]) != sorted(lst)
                            else:
                                global nextNodeID
                                #print >> sys.stderr, 'nextnode', nextNodeID
                                nextNodeID += 1
                                length = min([d for (_,d) in lst]) / 2
                                self.data[nextNodeID] = [(g,d-length) for (g,d) in lst]
                                anc = phylTree.lastCommonAncestor([self.info[g]['taxon_name'] for (g,_) in lst])
                                self.info[nextNodeID] = {'taxon_name':anc}
                                self.info[nextNodeID]["Duplication"] = 1 if hasLowScore(self, nextNodeID) else 3
                                todo.append(nextNodeID)
                                newData.add( (nextNodeID,length) )
                                #print >> sys.stderr, 'FLATTENTREE\n'
                                self.flattenTree(phylTree, False,  nextNodeID)
                                flag = True
                    assert len(newData) == len(set(g for (g,_) in newData)), newData
                    self.data[node] = [x for x in self.data[node] if x in newData] + list(newData.difference(self.data[node]))
                    for x in todo:
                        if hasLowScore(self, x):
                            self.info[x]["Duplication"] = 0
                            self.flattenTree(phylTree, False, x)

            # Appels recursifs
            for (g,_) in self.data[node]:
                #print >> sys.stderr,'gggggggggggdebut', g
                flag |= do(g)



            return flag

        return do(self.root if node is None else node)




    ###################################################################################################################
    # Lorsq'une duplication a ete suprimee parce qu'elle ne correspond pas aux criteres voulu, on cree une duplication#
    # chez l'ancetre suivant pour tester apres si la duplication suprimee est en faite soutenue chez l'ancetre suivant#
    # et ainsi on peut 'deplacer' des duplications                                                                    #
    #                                                                                                                 #
    ###################################################################################################################
    def CreateDuplication(self, phylTree, rec, node=None):

        def do(node):
            #print >> sys.stderr, 'nodeBIS', node
            # Fin du process sur une feuille
            if node not in self.data:
                return False
            assert len(self.data[node]) > 0

            flag = False
            # Appels recursifs
            if rec:
                for (g,_) in self.data[node]:
                    #print >> sys.stderr, 'g', g, self.info[g]
                    flag |= do(g)

            self.info[node]['taxon_name'] = phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )
            # Si c'est une vraie duplication, on n'a plus rien a faire
            #print >> sys.stderr, 'dup', node, self.info[node]['Duplication']
            if self.info[node]['Duplication'] >= 2:
                   # print >> sys.stderr, "oooookkkkkk"
                return flag
            else:
                return flag


        return do(self.root if node is None else node)






nextNodeID = -1


@utils.myTools.deprecated
def printTree(ft, data, info, root):
    ProteinTree(data, info, root).printTree(ft)



# Renvoie le suffixe associe a un numero de duplication ('a' -> 'z', 'aa' -> 'az', 'ba' ...)
##############################################################################################
def getDupSuffix(n, upper):
    base = 64 if upper else 96
    assert 1 <= n
    s = "."
    while n > 26:
        s = s + chr(base + (n-1)%26)
        n =  1 + (n-1)/26
    return s + chr(base + n)


# Charge l'arbre depuis un fichier
###################################
def loadTree(name):
    #, 'name', name
    import utils.myTools
    ns = utils.myTools.Namespace()

    # Lit la prochaine ligne du fichier (et bufferise la prochaine)
    def nextLine():

        old = ns.curr
           # print >> sys.stderr, 'old', old
        try:
            l = ""
           # print >> sys.stderr, 'llll', l
            while (l == "") or l.startswith("#"):
                # On enleve le \n final et on coupe suivant les \t
                l = f.readline().replace('\n', '')
            l = l.split('\t')
            # On stocke le triplet (indentation,key,value)
            ns.curr = (len(l)-2, l[-2], l[-1])
           # #, 'ns', ns.curr
        except StopIteration:
            ns.curr = None
        return old

    # La procedure d'analyse des lignes du fichier
    def recLoad(tree, indent):

        # l'ID du point
        currID = int(nextLine()[2])
           # print >> sys.stderr, 'currID', currID

        # Les infos associees
        tree.info[currID] = eval(nextLine()[2])
           # print >> sys.stderr, 'treeinfo', tree.info[currID]

        # Des fils ?
        child = []
           # print >> sys.stderr, '1111111111111', ns.curr
        while (ns.curr != None) and (ns.curr[0] == indent+1):
           # print >> sys.stderr, '2222222222'
            length = float(nextLine()[2])
            child.append( (recLoad(tree, indent+1), length) )
           # print >> sys.stderr, 'child1', child, currID
        if len(child) > 0:
            tree.data[currID] = child
          #  print >> sys.stderr, 'child2', child, currID
           # print >> sys.stderr, 'treedata', tree.data[currID]
        #print >> sys.stderr, 'treedata', tree.data
        return currID

    import utils.myFile

    print("Chargement du fichier d'arbres %s ..." % name, end=' ', file=sys.stderr)
    f = utils.myFile.openFile(name, "r") if isinstance(name, str) else name
    ns.curr = None
   # print >> sys.stderr, 'DDDDDDDDDDDDDDDDDDDEB'
    nextLine()
  #  print >> sys.stderr, 'EEEEEEEEEEEEEEEEEEEND'
    n = (0,0,0)
    while True:
        tree = ProteinTree()
      #  print >> sys.stderr, 'tree', tree
        tree.root = recLoad(tree, 0)
           # print >> sys.stderr, 'treeroot', tree.root
        yield tree
        n = (n[0]+1, n[1]+len(tree.data), n[2]+len(tree.info)-len(tree.data))
        #print >> sys.stderr, 'n', n
        if ns.curr == None:
            break
    print("%d racines, %d branches, %d noeuds OK" % n, file=sys.stderr)

    f.close()
