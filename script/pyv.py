import argparse
import os
import AuxiliarFunctions as af
from SequenceSimulator import SequenceSimulator

from collections import OrderedDict
from Bio import SeqIO


def gene(parameters_file, sequences_folder,tree):
	gene=True
	print("\nPreparing simulator of sequences")
	ss = SequenceSimulator(parameters, gene)
	ss.run(tree, sequences_folder)
	surviving_nodes=af.write_pruned_sequences(tree, sequences_folder)
	sequences=cut_seq(parameters,sequences_folder,tree)
	print_seq(parameters,sequences,surviving_nodes)

def HT(parameters,HT_folder,HT_tree):
	gene=False
	print("\nPreparing simulator of HT")
	ss = SequenceSimulator(parameters,gene)
	ss.run(HT_tree, HT_folder)
	
	af.write_pruned_sequences(HT_tree, HT_folder)
	

def cut_seq(parameters,sequences_folder,tree):
	dico={}
	size = parameters["SIZE"]
	numgene= parameters["GENE_NUMBER"]
	file_fasta=tree.split("/")[-1].split(".")[0]+"_pruned.fasta"
	for seq_record in SeqIO.parse(os.path.join(sequences_folder, file_fasta), "fasta"):
		sequence = str(seq_record.seq).upper()
		seq=([sequence[i:i+size] for i in range(0, len(sequence), size)])
		for i in range(1, numgene+1, 1):
			seq_id=str(i)+"_"+seq_record.id
			dico[seq_id]= seq[(i-1)]
	dico=OrderedDict(sorted(dico.items(), key=lambda t: t[0]))
	return (dico)

def print_seq(parameters,sequences,surviving_nodes) :
	numgene= parameters["GENE_NUMBER"]
	out_gene= parameters["OUTPUT_GENE"]
	af.check_folder(out_gene)
	i=1  
	nb=0    
	filename = open("%s/group%s.fasta" % (out_gene,i), "w")
	for key, v in sequences.items():
		nb+=1
		key= str(key)
		key=key.split("_")[1]
		filename.write(">%s\n%s\n" % (key,v))
		if ((nb==len(surviving_nodes))and(i!=numgene)) :
			nb=0
			print("Simulating sequence for gene family %s" % (i))
			filename.close()
			i+=1
			filename = open("%s/group%s.fasta" % (out_gene,i), "w")
	print("Simulating sequence for gene family %s" % (i))
	filename.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("params",  type=str, help="Parameters file")
	args = parser.parse_args()
	
	parameters_file =  args.params
	
	parameters = af.prepare_sequence_parameters(af.read_parameters(parameters_file))
	
	sequences_folder = parameters["OUTPUT_SEQ"]
	tree=parameters["TREE"]
	HT_folder = parameters["OUTPUT_HT"]
	HT_tree=parameters["HT_TREE"]
	option=parameters["OPTION"]
	if option==1 :
		af.check_folder(sequences_folder)
		af.check_folder(HT_folder)
		gene(parameters, sequences_folder, tree)
		HT(parameters,HT_folder,HT_tree)
	elif option==2 :
		af.check_folder(sequences_folder)
		gene(parameters, sequences_folder, tree)
	else : 
		af.check_folder(HT_folder)
		HT(parameters,HT_folder,HT_tree)
		 
	


