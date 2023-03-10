#!/usr/bin/python3

import sys
import warnings
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
import re
from math import log
import tempfile

print(sys.argv[1:])


#this is to take a DNA string and return it as a list of codons within the string
def Split_Into_Codons(DNAseq):
    codon_list = []
    for i in range(0, len(DNAseq)-2, 3):
        codon = DNAseq[i] + DNAseq[i+1] + DNAseq[i+2]
        codon_list.append(codon)
    return codon_list



#this takes a genbank file and returns the name of the organism
#"genus species"
#this is literally so I can see it at the top of the file and be
#less confused
def OrganismName(fileR):

    with open(fileR, "r") as fh:
        for line in fh:
            if re.search(r"^DEFINITION", line):
                c = 0
                name = ""
                for let in line:
                    if let == " ":
                         c += 1
                    if c == 4:
                        break
                    if c > 1:
                        if let == " ":
                            name += "_"
                        else:
                            name += let
                break
    return name[1:]


#takes genbank file to read from
def averageproteinlength(fileR):

    AllGeneIDs = {}

    Organism = OrganismName(fileR)
    for record in SeqIO.parse(fileR, "gb"):
        try:
            Acc = record.annotations["accessions"][0]
            print("On new gb file: ", Acc)
            full_chrom_seq = str(record.seq)       
        except: #some sort of Bio.SeqError
            print("cannot get full chrom seq")
            #now we need to pull from the fasta files for the unplaced genes
            ID = record.id
            if Organism == "Homo_sapiens":
                filepath = "/scratch/Human/Unplaced.fna"
            elif Organism == "Drosophila_melanogaster":
                filepath = "/scratch/d.melanogaster/drosophilaUnplaced.fna"
            else:
                raise FileNotFoundError("need a fasta file path for "+ Organism)
            
            for newrecord in SeqIO.parse(filepath, "fasta"):
                if newrecord.id == ID:
                    full_chrom_seq = str(newrecord.seq)
                
        for i in range(len(record.features)):
            P = record.features[i]
            if P.type == "CDS":
                strand = P.strand                            #will be +1 or -1
                GeneStart = P.location.nofuzzy_start + 1     #1 indexed
                GeneEnd = P.location.nofuzzy_end  
                intron = False
                exonPos = {}

                #now get CDS information
                try:
                    for j in P.qualifiers["db_xref"]:
                        find = re.search(r"GeneID:\d+", j)
                        if find != None:
                            GeneID = find.group()
                            break
                except KeyError:
                    GeneID = "LocusTag:" + P.qualifiers["locus_tag"][0]
                    warnings.warn(Acc + "  i: " + str(i) + "  no gene ID")
                try:
                    ProID = P.qualifiers["protein_id"][0]
                    aaSeq = P.qualifiers["translation"][0]
                except KeyError:
                    continue  #if there is no proteinID, this will not be checked for LCRs

                if GeneID in AllGeneIDs: #if same Gene, then check aaSeq length.
                    a = AllGeneIDs[GeneID][-1]  #this indexes the aaSeq
                    b = aaSeq
                    if len(a) >= len(b): #if previous aaSeq(withLCR)is longer, keep it and continue, otherwise check if it has LCR
                        continue     

                if P.location_operator == "join": #otherwise add in new information for that geneID
                    intron = True                 #or get info for that unique gene

                if intron == True:
                    for exon in P.location.parts:
                        ex = re.search(r"(\d+):[><]*(\d+)", str(exon)) #< or < cuz some exons are partial sequences
                        exstrand = re.search(r"\([+-]\)", str(exon))
                        exonPos[(int(ex.group(1))+1, int(ex.group(2)))] = exstrand.group() #1 indexed

                if (strand != 1 and strand != -1) and intron == False:    #exons on different strands
                    raise ValueError("There is no strand type and intron == False")
                if P.location_operator != "join" and P.location_operator != None:
                    raise ValueError("intron is neither 'join nor None")

                #get geneSeq
                if intron == True:    #concatenate the exons together into GeneSeq variable
                    GeneSeq = ""
                    for tup, s in exonPos.items():
                        forward = full_chrom_seq[tup[0]-1:tup[1]]
                        if s == "(+)":
                            GeneSeq += forward
                        elif s == "(-)":
                            GeneSeq += str(Seq(forward).reverse_complement())
                    #if there are no introns 
                else:
                    GeneSeq = full_chrom_seq[GeneStart-1: GeneEnd]
                    if strand == -1:
                        GeneSeq = str(Seq(GeneSeq).reverse_complement())   
                GeneSeq = GeneSeq[:-3].lower()  #remove stop codon and ensure all in lower case

                #Ensure DNAseq is 3X length of aaSeq
                if len(aaSeq)*3 != len(GeneSeq):
                    continue  #these will not be checked for LCRs

                #check for stop codons in ORF, these sequences will not be included
                codon_list = Split_Into_Codons(GeneSeq)
                if "tag" in codon_list or "tga" in codon_list or "taa" in codon_list:
                    continue  #do not run seg on it

                AllGeneIDs[GeneID] = [ProID, aaSeq]


    totalproteins = 0
    sumproteinlength = 0
    for LIST in AllGeneIDs.values():
        sumproteinlength += len(LIST[1])
        totalproteins += 1

    averageprotlength = round(sumproteinlength/totalproteins, 3)
    print("average protein length: ", averageprotlength)

if len(sys.argv[1]) == 0:
    print("takes file containing all genbank files in an organism; will give average protein length")
else:
    averageproteinlength(sys.argv[1])

