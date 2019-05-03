import ete3
import numpy
import copy
import sys
import os

def write_pruned_sequences(tree_file, fasta_folder):

    with open(tree_file) as f:
        line = f.readline().strip()
        if "(" not in line or line == ";":
            return None
        else:
            my_tree = ete3.Tree(line, format=1)

    surviving_nodes = {x.name for x in my_tree.get_leaves()}
    file_name = tree_file.split("/")[-1].split(".")[0]
    entries = fasta_reader(fasta_folder + "/" + file_name + "_complete.fasta")
    clean_entries = list()
    for h, seq in entries:
        if h[1:] in surviving_nodes:
            clean_entries.append((h, seq))

    fasta_writer(fasta_folder + "/" + file_name + "_pruned.fasta", clean_entries)
    return(surviving_nodes)


def prepare_sequence_parameters(parameters):

    for parameter, value in parameters.items():

        if parameter == "VERBOSE" or parameter =="SIZE" or parameter =="GENE_NUMBER" or parameter =="HT_SIZE" or parameter =="OPTION":
            parameters[parameter] = int(value)

        if parameter == "SCALING":
            parameters[parameter] = float(value)
        

    return parameters


def read_parameters(parameters_file):

    parameters = dict()

    with open(parameters_file) as f:

        for line in f:

            if line[0] == "#" or line == "\n":
                continue

            if "\t" in line:
                parameter, value = line.strip().split("\t")
                parameters[parameter] = value
            elif " " in line:
                parameter, value = line.strip().split(" ")
                parameters[parameter] = value

    return parameters


def fasta_reader(fasta_file):

    with open(fasta_file) as f:

        seq = ""
        for line in f:
            if ">" == line[0]:
                if seq != "":
                    yield header, seq
                    header = line.strip()
                    seq = ""
                else:
                    header = line.strip()
                    seq = ""
            else:
                seq += line.strip()

        yield header, seq
        

def fasta_writer(outfile, entries):

    x = 80
    with open(outfile, "w") as f:
        for h, seq in entries:
            f.write(h + "\n")
            lines = [seq[i: i + x] for i in range(0, len(seq), x)]
            for line in lines:
                f.write(line +"\n")

def check_folder(folder):
	if not os.path.isdir(folder):
		os.mkdir(folder)
	else:
		os.system("rm -fr " + os.path.join(folder))
		os.mkdir(folder)
