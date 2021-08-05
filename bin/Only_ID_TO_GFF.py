#!/usr/bin/env python
#coding=utf-8
'''
Assemble_ncRNA_mergeGff1.py 根据genetype 来整理的注释，Assemble_ncRNA_mergeGff.py 把gene_type都是rRNA
'''

import sys

gff=open(sys.argv[1],'r')
gff_new=open(sys.argv[2],'w')
ncRNA_num=1
for line in gff:
	lines=line.strip().split("\t")
	col_head2="\t".join(lines[0:2])
	col_4_to_8="\t".join(lines[3:8])
	id_name=lines[8]
	if lines[2] == "CDS":
		gff_new.write(col_head2 + "\tgene\t"+ col_4_to_8 \
		+"\tID="+ id_name+";Name="+id_name+";gbkey=Gene;gene_biotype="+"protein_coding"+"\n")
		gff_new.write(col_head2 + "\tCDS\t"+ col_4_to_8 \
		+"\tID="+ id_name+".cds;Parent="+id_name+";Name="+id_name+".cds;gbkey=CDS"+"\n")
	if lines[2] =="rRNA":
		ID_name=id_name+".rRNA"
		ID_name1=ID_name+".1"
		gff_new.write(col_head2 + "\tgene\t"+ col_4_to_8 +"\tID="+ id_name+";Name="+id_name+";gbkey=Gene;gene_biotype="+lines[2]+"\n")
		gff_new.write(col_head2 + "\t"+lines[2]+"\t"+ col_4_to_8 +"\tID="+ID_name +";Parent=" +id_name+";gbkey="+lines[2]+"\n")
		gff_new.write(col_head2 + "\texon\t"+ col_4_to_8 +"\tID="+ID_name1 +";Parent=" +ID_name+";gbkey="+lines[2]+"\n")
	if lines[2] =="tRNA":
		ID_name=id_name+".tRNA"
		ID_name1=ID_name+".1"
		gff_new.write(col_head2 + "\tgene\t"+ col_4_to_8 +"\tID="+ id_name+";Name="+id_name+";gbkey=Gene;gene_biotype="+lines[2]+"\n")
		gff_new.write(col_head2 + "\t"+lines[2]+"\t"+ col_4_to_8 +"\tID="+ID_name +";Parent=" +id_name+";gbkey="+lines[2]+"\n")
		gff_new.write(col_head2 + "\texon\t"+ col_4_to_8 +"\tID="+ID_name1 +";Parent=" +ID_name+";gbkey="+lines[2]+"\n")
	if lines[2] =="sRNA":
		ID_name=id_name+".sRNA"
		ID_name1=ID_name+".1"
		gff_new.write(col_head2 + "\tgene\t"+ col_4_to_8 +"\tID="+ id_name+";Name="+id_name+";gbkey=Gene;gene_biotype="+lines[2]+"\n")
		gff_new.write(col_head2 + "\t"+lines[2]+"\t"+ col_4_to_8 +"\tID="+ID_name +";Parent=" +id_name+";gbkey="+lines[2]+"\n")
		gff_new.write(col_head2 + "\texon\t"+ col_4_to_8 +"\tID="+ID_name1 +";Parent=" +ID_name+";gbkey="+lines[2]+"\n")
	if lines[2] =="repeat":
		ID_name=id_name+".repeat"
		ID_name1=ID_name+".1"
		gff_new.write(col_head2 + "\tgene\t"+ col_4_to_8 +"\tID="+ id_name+";Name="+id_name+";gbkey=Gene;gene_biotype="+lines[2]+"\n")
		gff_new.write(col_head2 + "\t"+lines[2]+"\t"+ col_4_to_8 +"\tID="+ID_name +";Parent=" +id_name+";gbkey="+lines[2]+"\n")
		gff_new.write(col_head2 + "\texon\t"+ col_4_to_8 +"\tID="+ID_name1 +";Parent=" +ID_name+";gbkey="+lines[2]+"\n")
	if lines[2] =="ncRNA":
		ID_name=id_name+".ncRNA"
		ID_name1=ID_name+".1"
		gff_new.write(col_head2 + "\tgene\t"+ col_4_to_8 +"\tID="+ id_name+";Name="+id_name+";gbkey=Gene;gene_biotype="+lines[2]+"\n")
		gff_new.write(col_head2 + "\t"+lines[2]+"\t"+ col_4_to_8 +"\tID="+ID_name +";Parent=" +id_name+";gbkey="+lines[2]+"\n")
		gff_new.write(col_head2 + "\texon\t"+ col_4_to_8 +"\tID="+ID_name1 +";Parent=" +ID_name+";gbkey="+lines[2]+"\n")
gff.close()
gff_new.close()
