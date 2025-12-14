#!/usr/bin/python3

from Bio import SeqIO
import re
from Bio.Seq import Seq
import sys

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


#this is to take a DNA string and return it as a list of codons within the string
def Split_Into_Codons(DNAseq):

    DNAseq = DNAseq.lower()
    codon_list = []
    for i in range(0, len(DNAseq)-2, 3):
        codon = DNAseq[i] + DNAseq[i+1] + DNAseq[i+2]
        codon_list.append(codon)
    return codon_list


#takes the fullchromosome sequence as a string and returns a dictionary of the number of a,c,g, and ts
#int it
#if ambig_nts are encountered, a fraction of a nt will be added
def count_nts(DNAseq):

    DNAseq = DNAseq.lower()
    ambig_nts = {"a":["a"], "c":["c"], "g":["g"], "t":["t"], "n":["a","c","g","t"], "m":["a","c"],\
                "r":["a","g"], "w":["a","t"], "s":["c","g"], "y":["c","t"], "k":["g","t"],\
                "v":["a","c","g"], "h":["a","c","t"], "d":["a","g","t"], "b":["c","g","t"]}
    count = {}
    count["a"] = 0; count["c"] = 0; count["g"] = 0; count["t"] = 0
    for nt1 in DNAseq:
        if nt1 in count:
            count[nt1] += 1
        else:
            list_nts = ambig_nts[nt1]
            for nt2 in list_nts:
                count[nt2] += 1/len(list_nts)
    return count

def count_aas(fullaastring):
   
    fullaastring = fullaastring.lower()
    ambig_aas = {"x":["a","r","n","d","c","q","e","g","h","i","l","k","m","f","p","s","t","w","y","v"],\
                 "b":["d","n"], "z":["e","q"], "j":["i","l"], "u":["c"], "o":["k"]}
    count = {}
    count["a"] = 0; count["c"] = 0; count["d"] = 0; count["e"] = 0
    count["f"] = 0; count["g"] = 0; count["h"] = 0; count["i"] = 0
    count["k"] = 0; count["l"] = 0; count["m"] = 0; count["n"] = 0
    count["p"] = 0; count["q"] = 0; count["r"] = 0; count["s"] = 0
    count["t"] = 0; count["v"] = 0; count["w"] = 0; count["y"] = 0
    for aa in fullaastring:
        if aa in count:
            count[aa] += 1
        else:
            ambig_aa_list = ambig_aas[aa]
            for aa2 in ambig_aa_list:
                count[aa2] += 1/len(ambig_aa_list)
    return count
   

def count_codons(codon_list):

    ambig_nts = {"a":["a"], "c":["c"], "g":["g"], "t":["t"], "n":["a","c","g","t"], "m":["a","c"],\
                "r":["a","g"], "w":["a","t"], "s":["c","g"], "y":["c","t"], "k":["g","t"],\
                "v":["a","c","g"], "h":["a","c","t"], "d":["a","g","t"], "b":["c","g","t"]}
    count = {}
    count["ttt"] = 0; count["ttc"] = 0; count["tta"] = 0; count["ttg"] = 0
    count["tct"] = 0; count["tcc"] = 0; count["tca"] = 0; count["tcg"] = 0
    count["tat"] = 0; count["tac"] = 0; count["tgt"] = 0; count["tgc"] = 0
    count["tgg"] = 0; count["ctt"] = 0; count["ctc"] = 0; count["cta"] = 0
    count["ctg"] = 0; count["cct"] = 0; count["ccc"] = 0; count["cca"] = 0
    count["ccg"] = 0; count["cat"] = 0; count["cac"] = 0; count["caa"] = 0
    count["cag"] = 0; count["cgt"] = 0; count["cgc"] = 0; count["cga"] = 0
    count["cgg"] = 0; count["att"] = 0; count["atc"] = 0; count["ata"] = 0
    count["atg"] = 0; count["act"] = 0; count["acc"] = 0; count["aca"] = 0
    count["acg"] = 0; count["aat"] = 0; count["aac"] = 0; count["aaa"] = 0
    count["aag"] = 0; count["agt"] = 0; count["agc"] = 0; count["aga"] = 0
    count["agg"] = 0; count["gtt"] = 0; count["gtc"] = 0; count["gta"] = 0
    count["gtg"] = 0; count["gct"] = 0; count["gcc"] = 0; count["gca"] = 0
    count["gcg"] = 0; count["gat"] = 0; count["gac"] = 0; count["gaa"] = 0
    count["gag"] = 0; count["ggt"] = 0; count["ggc"] = 0; count["gga"] = 0
    count["ggg"] = 0

    for codon in codon_list:
        if codon in count:
            count[codon] += 1
        elif codon == "taa" or codon == "tag" or codon == "tga":
            pass
        else:                                   #if there are ambiguous nt in codon
            possible_cods = []
            possib_nt_1 = ambig_nts[codon[0]]   #list of all possible nts for 1st pos in codon
            possib_nt_2 = ambig_nts[codon[1]]   #list of all possible nts for 2nd pos in codon
            possib_nt_3 = ambig_nts[codon[2]]   #list of all possible nts for 3rd pos in codon
            for char1 in possib_nt_1:
                for char2 in possib_nt_2:
                    for char3 in possib_nt_3:
                        possible_cods.append(char1 + char2 + char3)
            if "taa" in possible_cods:
                possible_cods.remove("taa")
            if "tag" in possible_cods:
                possible_cods.remove("tag")
            if "tga" in possible_cods:
                possible_cods.remove("tga")
            #now add the fraction for proportion of each codon observed
            for new_cod in possible_cods:
                count[new_cod] += 1/len(possible_cods)
    return count


def determine_proportions(fileR):

    print("started running")
    AllGeneIDs = {}  #only take one of each geneID
    AllChromSeqs = ""

    organism = OrganismName(fileR)
    for record in SeqIO.parse(fileR, "gb"):
        try:
            #Acc = record.annotations["accessions"][0]
            full_chrom_seq = str(record.seq)
        except:
            ID = record.id
            if organism == "Homo_sapiens":
                filepath = "/scratch/Human/Unplaced.fna"
            elif organism == "Drosophila_melanogaster":
                filepath = "/scratch/d.melanogaster/drosophilaUnplaced.fna"
            else:
                raise FileNotFoundError("need a fasta file path for "+organism)

            for newrecord in SeqIO.parse(filepath, "fasta"):
                if newrecord.id == ID:
                    full_chrom_seq = str(newrecord.seq)
        AllChromSeqs += full_chrom_seq  #make one giant sequence (count nts later)

        for i in range(len(record.features)):          #now get coding seq info
            P = record.features[i]
            if P.type == "CDS":
                try:              #get GeneID
                    for j in P.qualifiers["db_xref"]:
                        find = re.search(r"GeneID:\d+", j)
                        if find != None:
                            GeneID = find.group()
                            break
                except KeyError:
                    raise KeyError("There is no GeneID. Accession = "+Acc+"feature index ="+str(i))
                try:
                    aaSeq = P.qualifiers["translation"][0]
                except KeyError:
                    pass #some CDS do not have a ProID

                if GeneID in AllGeneIDs: #ensure we are only counting for 1 splice variant
                    a = AllGeneIDs[GeneID][-1] #indexes old aaSeq
                    b = aaSeq
                    if len(a) >= len(b): #if previous aaSeq is longer, keep it, otherwise count!
                        continue

                strand = P.strand
                GeneStart = P.location.nofuzzy_start + 1     #1 indexed
                GeneEnd = P.location.nofuzzy_end  
                intron = False
                exonPos = {}

                if P.location_operator == "join": #otherwise add in new information for that geneID
                    intron = True                 #or get info for that unique gene
                if intron == True:
                    for exon in P.location.parts:
                        ex = re.search(r"(\d+):[><]*(\d+)", str(exon)) #< or < cuz some exons are partial sequences
                        exstrand = re.search(r"\([+-]\)", str(exon))
                        exonPos[(int(ex.group(1))+1, int(ex.group(2)))] = exstrand.group() #1 indexed
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
                GeneSeq = GeneSeq[:-3]  #remove stop codon

                #Ensure DNAseq is 3X length of aaSeq
                if len(aaSeq)*3 == len(GeneSeq):     #just becuase I don't trust that I would be counting the codons correctly otherwise
                    AllGeneIDs[GeneID] = [GeneSeq, aaSeq]
                    print("adding seqs to dictionary")
   
    allcds = ""
    allaas = ""
    for seq_list in AllGeneIDs.values():  #put geneseq and aaseq into 1 big seq
        allcds += seq_list[0] #GeneSeqs (all multiples of 3)
        allaas += seq_list[1] #aaseq

    print("counting codons, aas, nt in cds")
    allnt_count = count_nts(AllChromSeqs) #dict of a,c,g,t for all sequences            #F1
    ntcds_count = count_nts(allcds)             #dict of nt frequencies in CDS          #F2
    aa_count =  count_aas(allaas)               #dict of aa frequencies in proteome     #F3
    codon_list = Split_Into_Codons(allcds)
    codon_count = count_codons(codon_list)      #dict of all codon frequencies          #F4

    print("writing files")
    
    #now put into nice table
    NOTE = "#"+organism+"\n"
    F1 = organism+"_allnt_count"
    F2 = organism+"_ntcds_count"
    F3 = organism+"_aa_count"
    F4 = organism+"_codon_count"
    #make the headers for 3 tables
    with open(F1,"w") as fh:
        fh.write(NOTE+"#count of all nts in genome\n.\ta\tc\tg\tt\n")
    with open(F2, "w") as fh:
        fh.write(NOTE+"#count of nts in all CDS (same seqs used as CombinedUsingSysCall)\n.\
    \ta\tc\tg\tt\n")
    with open(F3, "w") as fh:
        fh.write(NOTE+"#count of AAs (same seqs used as CombinedUsingSysCall)\n.\
    \ta\tc\td\te\tf\tg\th\ti\tk\tl\tm\tn\tp\tq\tr\ts\tt\tv\tw\ty\n")
    with open(F4, "w") as fh:
        fh.write(NOTE+"#count of codons in all coding sequences\n.\
    \tttt\tttc\ttta\tttg\ttct\ttcc\ttca\ttcg\ttat\ttac\ttgt\ttgc\
    \ttgg\tctt\tctc\tcta\tctg\tcct\tccc\tcca\tccg\tcat\tcac\tcaa\
    \tcag\tcgt\tcgc\tcga\tcgg\tatt\tatc\tata\tatg\tact\tacc\taca\
    \tacg\taat\taac\taaa\taag\tagt\tagc\taga\tagg\tgtt\tgtc\tgta\
    \tgtg\tgct\tgcc\tgca\tgcg\tgat\tgac\tgaa\tgag\tggt\tggc\tgga\tggg\n")

    to_F1 = "freq\t"  #all nt table   (F1)
    total_allnt = 0
    for num_val in allnt_count.values():
        total_allnt += num_val
        to_F1 += str(num_val)+"\t"
    to_F1 += "\nprop\t"
    for num_val in allnt_count.values():
        to_F1 += str(round(num_val/total_allnt,5))+"\t"

    to_F2 = "freq\t"   #nt in cds table   (F2)
    total_cdsnt = 0
    for num_val in ntcds_count.values():
        total_cdsnt += num_val
        to_F2 += str(num_val)+"\t"
    to_F2 += "\nprop\t"
    for num_val in ntcds_count.values():
        to_F2 += str(round(num_val/total_cdsnt,5))+"\t"

    to_F3 = "freq\t"   #aa freq table   (F3)
    total_aa = 0
    for num_val in aa_count.values():
        total_aa += num_val
        to_F3 += str(num_val)+"\t"
    to_F3 += "\nprop\t"
    for num_val in aa_count.values():
        to_F3 += str(round(num_val/total_aa,5))+"\t"

    to_F4 = "freq\t"   #aa freq table   (F3)
    total_cod = 0
    for num_val in codon_count.values():
        total_cod += num_val
        to_F4 += str(num_val)+"\t"
    to_F4 += "\nprop\t"
    for num_val in codon_count.values():
        to_F4 += str(round(num_val/total_cod, 5))+"\t"


    with open(F1, "a") as fh:
        fh.write(to_F1)
    with open(F2, "a") as fh:
        fh.write(to_F2)
    with open(F3, "a") as fh:
        fh.write(to_F3)
    with open(F4, "a") as fh:
        fh.write(to_F4)


determine_proportions(sys.argv[1])   #takes gb file
