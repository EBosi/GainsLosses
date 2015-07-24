from IPython import embed

def getName(clade):
	""" by default BioPhy clades are unnamed if not leaves """
	if clade.name != None: return clade.name
	return '_'.join(i.name for i in clade.get_terminals())
	
def getDirectChildren(clade):
	clade=clade.get_nonterminals()[0]
	dChildren=[i for i in getChildren(clade) if i in clade] # should work also for tree
	return dChildren

def isLeaf(t):
	return t.is_terminal()
	
def getChildren(t):
	return t.get_nonterminals()[1:] + t.get_terminals()

def spawnLeaf(clade,container,isPresent=False):
	name=getName(clade)
	EEs=EvolutionaryEvents(name)
	EEs.startLeaf(isPresent)
	container[clade]=EEs
	return EEs

def spawnParent(clade,children,container):
	name=getName(clade)
	EEs=EvolutionaryEvents(name)
	container[clade]=EEs
	children=[container[i] for i in children]
	EEs.startParent(children)
	return EEs
		
#########################################
###		Algorithm specific objects/fxns #
#########################################
	
def assumptionAi(parent, children):
	""" under assumption Ai, there are two alternative scenarios:
	i) the gene has been lost in the parent or ii) the gene has not been lost
	in the parent."""
	# compute events according to rules in Mirkin et al. 2003
	scenarioI = {'Gi': {e for child in children for e in child.eventPatterns['Gn']} , 'Li':{e for child in children for e in child.eventPatterns['Ln']}.union({parent.name})}
	scenarioII = {'Gi': {e for child in children for e in child.eventPatterns['Gi']} , 'Li':{e for child in children for e in child.eventPatterns['Li']}}
	# choose most parsimonious event and set parent's events
	mostParsimonious = min(scenarioI,scenarioII,key=lambda x: sum(map(len,x.values())))
	parent.eventPatterns.update(mostParsimonious)
	
def assumptionAn(parent, children, g=1):
	""" under assumption An, there are two alternative scenarios:
	i) the gene has been gained in the parent or ii) the gene has not been gained
	in the parent.
	The g parameter measures the likelihood of gains over losses [default=1 : equiprobable] """
	# compute events according to rules in Mirkin et al. 2003
	scenarioI = {'Gn': {e for child in children for e in child.eventPatterns['Gi']}.union({parent.name}) , 'Ln':{e for child in children for e in child.eventPatterns['Li']}}
	scenarioII = {'Gn': {e for child in children for e in child.eventPatterns['Gn']} , 'Ln':{e for child in children for e in child.eventPatterns['Ln']}}
	# choose most parsimonious event and set parent's events
	scenarios=[scenarioI, scenarioII]
	numberOfEvents = [sum(map(len,scenarioI.values())) -1 + g, sum(map(len,scenarioII.values()))]
	mostParsimonious = scenarios[scenarios.index(min(scenarioI,scenarioII, key=lambda x: numberOfEvents[scenarios.index(x)]))]
	parent.eventPatterns.update(mostParsimonious)

class EvolutionaryEvents(object):
	""" class accounting for evolutionary events at a given node """
	def __init__(self, name):
		self.name = name
		self.eventPatterns = {'Gi':set(),'Li':set(),
								'Gn':set(),'Ln':set()}
	def startLeaf(self,isPresent):
		if isPresent: self.eventPatterns['Gn'] = self.eventPatterns['Gn'].union(set([self.name]))
		else: self.eventPatterns['Li'] = self.eventPatterns['Li'].union(set([self.name]))

	def startParent(self,children,g=1):
		self.children = children
		assumptionAi(self,self.children)
		assumptionAn(self,self.children,g)

	def backTracking(self,g=1):
		d={k:len(v) for k,v in self.eventPatterns.iteritems()}
		d['Gi'], d['Gn'] = g*d['Gi'], g*d['Gn']
		mostParsimoniousAssumption = min(['Gi','Li'],['Gn','Ln'],key=lambda x: (sum([d[i] for i in x]),sum([d[i] for i in x if i.startswith('G')])))
		solution = {k: self.eventPatterns[k] for k in mostParsimoniousAssumption}
		return solution

	def hasChildren(self):
		return hasattr(self,'children')

	def propagate(self,container,g=1):
		""" return clades with gains and losses """
		clades={'gains':set(),'losses':set()}
		solution=self.backTracking(g)
		namesToClades={getName(k):k for k in container.keys()}
		if 'Gi' in solution.keys(): clades['gains'].add(namesToClades[self.name])
		for k,v in solution.iteritems():
			if k.startswith('G'): k_='gains'
			else: k_='losses'
			for v_ in v: clades[k_].add(namesToClades[v_])
		return clades

def walk(t):
	""" 
	general algorithm to walk trees. needs a global variable named 'visited' 
		WORKING!
		"""
	if t in visited: return
	if isLeaf(t):
		visited.add(t)
		print t
		return t
	for c in getDirectChildren(t):
		walk(c)
	visited.add(t)
	print t
	return t

def tryWalk():
	from Bio import Phylo
	visited=set()
	t=Phylo.read('prova.nwk','newick')
	walk(t) 

def phyleticWalk(t,presDict):
	""" 
	specific algorithm to walk trees while recovering phyletic patterns.
	Still needs a global dict named *visited* and a set of nodes with the protein
	"""
	if t in visited: return
	if isLeaf(t):
		hasIt=t.name in presDict
		node=spawnLeaf(t,visited,hasIt)
		visited[t]=node
		return node
	children=getDirectChildren(t)
	for c in children:
		phyleticWalk(c,presDict)
	node=spawnParent(t,children,visited)
	visited[t]=node
	return node

###########################################

def main(presence_list,tree_file,fmt='newick',protein='protein',g=1):
	""" presence list is the list/set of taxa in which the protein is present """
	from Bio.Phylo import read
	try: tree_file.trace
	except: t=read(tree_file,fmt)
	else: t=tree_file
	visited={}
	global visited
	clade=t.get_nonterminals()[0]
	phyleticWalk(clade,presence_list)
	return visited,clade
	
def test_main():
	""" Figure 5 of mirkin et al. 2003 """
	presence_list={'y','o','b','l','w','d','c','g','h','e','t','i','v'}
	t='prova.nwk'
	return main(presence_list,t)

def test_main2():
	""" Figure 6 of mirkin et al. 2003  """
	presence_list={'y','o','b','l','w','d','c','g','h','e','t','i','v'}
	t='prova2.nwk'
	return main(presence_list,t)

###########################################

def extendedMain(pres_matrix,tree_file,fmt='newick',g=1,transpose=False,sep='\t'):
	""" matrix should be rows -> proteins; columns -> strains; else use transpose.
		First row/column should be headers.
	"""
	import numpy as np
	from Bio.Phylo import read
	pres_matrix=[l.strip().split(sep) for l in open(pres_matrix)]
	pres_matrix,colnames,rownames=np.array([i[1:] for i in pres_matrix[1:]]),np.array(pres_matrix[0]),np.array([i[0] for i in pres_matrix[1:]])
	if transpose: pres_matrix,colnames,rownames=pres_matrix.T,rownames,colnames
	proteins,strains=rownames,colnames
	prot_dict={}
	tree=read(tree_file,fmt)
	# compute parsimony events for each protein
	for i,prot in enumerate(proteins):
		pres_abs=np.array(map(int,pres_matrix[i]))
		pres_list=set(strains[pres_abs==1])
		container,root=main(pres_list,tree)
		prot_dict[prot]=container[root].propagate(container)
	# get luca genes
	root=tree.get_nonterminals()[0]
	LUCA_genes=[prot for k,v in prot_dict.iteritems() if root in v['gains']]
	# get a count of gain/losses for each clade
	clades=tree.get_terminals()+tree.get_nonterminals()
	counts={}
	for c in clades:
		counts[c]={'gains':set(),'losses':set()}
		for p in prot_dict:
			if c in prot_dict[p]['gains']: counts[c]['gains'].add(p)
			if c in prot_dict[p]['losses']: counts[c]['losses'].add(p)
	return LUCA_genes,counts
		
	
# test_main()

############################
# toy example: toyTree.csv #
############################





















