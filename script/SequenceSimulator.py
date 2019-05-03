import pyvolve
import os
import ete3
import numpy
import AuxiliarFunctions as af
import numpy as np



class SequenceSimulator():
	
	def __init__(self, parameters,gene):
		self.parameters = parameters
		size = self.parameters["SIZE"]
		numgene= self.parameters["GENE_NUMBER"]
		if (gene==True) :
			self.size_seq = size*numgene
		else : 
			self.size_seq= self.parameters["HT_SIZE"]
		
		self.model = self.get_codon_model()
	
	def run(self, tree_file, sequences_folder):
		with open(tree_file) as f:
			line = f.readline().strip()
			if "(" not in line or line == ";":
				return None
			else:
				my_tree = ete3.Tree(line, format=1)
		tree = pyvolve.read_tree(tree=my_tree.write(format=5), scale_tree = self.parameters["SCALING"])
		name_mapping = self.get_mapping_internal_names(tree, my_tree)
		partition = pyvolve.Partition(models=self.model, size=self.size_seq)
		evolver = pyvolve.Evolver(tree=tree, partitions=partition)
		fasta_file = tree_file.split("/")[-1].split(".")[0] + "_complete.fasta"
		evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None, write_anc=True)

		# Correct the names
		self.correct_names(os.path.join(sequences_folder, fasta_file), name_mapping)

	def get_codon_model(self):
		codon_params = {}
		for param in ["ALPHA", "BETA", "KAPPA"]:
			codon_params[param.lower()] = float(self.parameters[param])
		codon_fitness_file = open("Input/codon_fitness.txt","r")
		codon_fitness =[]
		for line in codon_fitness_file :
			codon_fitness.append(float(line))
		#codon_fitness = np.random.normal(size = 61) 
		#f = pyvolve.ReadFrequencies("codon", file = "FBgn0034744.fasta")
		#frequencies = f.compute_frequencies()
		return pyvolve.Model("MutSel", {"fitness": codon_fitness})
		#return pyvolve.Model(self.parameters['CODON_MODEL'], codon_params, neutral_scaling=True)
	
	def get_mapping_internal_names(self, pytree, ettree):

		pyvolvemap = dict()

		def traverse(root):
			if root:
				if len(root.children) == 0:
					return None
				else:
					traverse(root.children[0])
					traverse(root.children[1])
					pyvolvemap[root.children[0].name + "+" + root.children[1].name] = root.name
		traverse(pytree)
		good_mapping = dict()
		etroot = ettree.get_tree_root().name
		good_mapping["myroot"] = etroot
		
		for n in ettree.traverse(strategy="postorder"):
			if not n.is_leaf():
				c1, c2 = n.get_children()
				n1 = c1.name + "+" + c2.name
				n2 = c2.name + "+" + c1.name
				if n1 in pyvolvemap:
					good_mapping[pyvolvemap[n1]] = n.name
					n.name = pyvolvemap[n1]
				if n2 in pyvolvemap:
					good_mapping[pyvolvemap[n2]] = n.name
					n.name = pyvolvemap[n2]
					
		return good_mapping
		
		
	def correct_names(self,fasta_file, good_mapping):
		entries = list()
		for n,v in af.fasta_reader(fasta_file):
			if "root" in n:
				entries.append((">"+good_mapping["myroot"], v))
			elif "internal" in n:
				entries.append((">"+good_mapping[n[1:]], v))
			else:
				entries.append((n, v))
				
		af.fasta_writer(fasta_file, entries)











