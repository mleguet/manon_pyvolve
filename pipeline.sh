##!/bin/sh

rm -fr gene
rm -fr HT
rm -fr Output
rm -fr gene_ali/*
rm -fr HT_ali/*
rm -fr group*

python3 script/pyv.py ./Input/SequenceParameters.tsv

ls -l gene | awk '{system("java -jar script/macse_v2.03.jar -prog alignSequences -seq gene/"$9"");}'

rm gene/*_AA*
mv gene/*_NT* gene_ali/  

rename -v "s/_NT//g" gene_ali/*


sed -i s/'!'/N/g gene_ali/*
cp gene_ali/* ./

#seaview gene_ali/group1.fasta

java -jar script/macse_v2.03.jar -prog alignSequences -seq HT/*_pruned.fasta

rm HT/*_AA*
mv HT/*_NT* HT_ali/  


rename -v "s/_pruned_NT//g" HT_ali/*


sed -i s/'!'/N/g HT_ali/*
cp HT_ali/* ./

