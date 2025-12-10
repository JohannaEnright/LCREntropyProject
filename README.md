This project investigates the evolution of low complexity regions (LCRs) in proteins by comparing the entropy, as calculated by Shannon's Entropy equation, of the LCR segment of a protein to its encoding DNA sequence and vice versa. Simulations are used to compare the observed entropy correlations to what we would expect under various evolutionary mechanisms.
The project looks at LCRs in the entire proteome and genome of five model organisms: Homo sapiens, Saccharomyces cerevisiae, Drosophila melanogaster, Caenorhabditis elegans, and Arabidopsis thaliana. 
Results from the analysis are published in the Journal of Molecular Biology and Evolution and can be found at DOI: 10.1093/molbev/msad084 (Enright, Dickson, & Golding, 2023). 

The main script components include:

1. CombinedUsingSysCall_Codon3.py
   This script was used to determine the entropies of LCRs from real biological sequences and requires the following inputs:
   i) a genbank file for a particular organism, which can be multiple genbank files for different chromosomes, concatenated together
   ii) the direction of comparison, which can be "pro_dna", "pro_cod", "dna_pro", "dna_cod", "cod_dna", or "cod_pro".
       For example, pro_dna scans for LCRs at the protein level and compares the sequence entropy to that of its underlying DNA sequence,
       while dna_cod scans for LCRs at the DNA level and compares the entropy to that of the codon entropy it encodes.
   iii) segA parameter W ,which is the window length as an integer
   iv) segA parameter K1, which is the "trigger complexity" as a float.
       Trigger complexity is the maximum entropy threshold within the sliding window that triggers segA to stop and expand the window in search of the full LCR.
   v) segA parameter K2, which is the "extension complexity" as a float.
      Upon finding an LCR of the given window length below the trigger complexity threshold, segA will scan outwards from this window until it determines the full LCR with an entropy below the extension complexity threshold.

   This script outputs a plain text files parsed on whitespace, name as inputfilename_Output.
   It contains the information for protein ID, genbank accession number, the DNA strand the LCR was found on, if there were introns found in the gene, if the LCR spanned an exon juncion, the start position of the gene, the end position of the gene, the starting amino acid position of the protin LCR, the end amino acid position of the protein LCR, the entropy of the LCR from the sequence type you were searching for, the starting LCR position on the chromosome, the end LCR position on the chromosome, and the entropy of the corresponding sequence type which was being compared. The key information which was used for further analysis was entropy y ~ entropy x.
   
   Other plain text files are also outputed and include more information on the LCRs, information on proteins which did not contain LCRs (for counting/checking purposes), and information if the LCRs/genes were more atypical. For example, if the LCR was found in a gene which was located on both strands of the DNA, or if the LCR spanned an exon junction.
   
2. segA
   This script is a modified version of the original Seg algorithm created by Wootton and Federhen (1993). It has been modified to extend to longer LCR sequences than were previously capable of being handled by Seg. As well, it has an added alphabet to identify codon LCRs and deal with ambiguous letters. segA is called automatically in the CombinedUsingSysCall_Codon3.py scipt.
   
3. GiveCombSysCallArgs.sh
   This is a bash file which is used to feed CombinedUsingSysCall_Codon3.py a combination of set parameters for segA. This was important to ensure that we were observing the same trends given slightly different LCR identifying parameters, considering the somewhat arbitrary nature of these parameters. In total, it will feed CombinedUsingSysCall_Codon3.py 27 unique parameter combinations as shown in the supplementary section of the publication (Enright, Dickson, & Golding, 2023).
   
4. averageproteinlength.py
   This script requires a genbank file as input, and will return the average protein length for that file.
   Again, it also accepts multiple concatenated genbank files, hence it was used to calculate the average protein length for each organism.
   This information was used in the generaterandomproteome script.

6. codonclass_categorizeLCRs.py
   This script requires the _Seq file outputed by CombinedUsingSysCall_Codon3.py which contains the sequences of all of the detected LCRs for a given organism, and a file output location.
   It categorizes each LCR into one of three categories:
     i) containing three unique nucleotides, ii) containing two unique nucleotides, and iii) containing one unique nucleotide and outputs this information as counts of LCRs in each category.
   
8. countcodontypes_allCDS.py
   This script takes a genbank (or multiple) file as input. It outputs a text file with information for the total codons of each codon class summed for all sequences in a genome.
   This information is used for the generaterandomproteome script
    
9. generaterandomproteome
   This script is to simulate 100 000 sequences outputed in fasta file format. 
   The input values are the simulation type ('mutate' or 'random'), species ('sc', 'dm', 'hs', 'ce', 'at'), fasta file of sequences to mutate (only if you choose the 'mutate' option),
   write file (name of file to write to).
   The species input is necessary to generate sequences which have the average length, the correct codon class proportions, and the same proportion of LCRs as the corresponding species.
   The 'random' simulation will randomly generate 100 000 random sequences given the various parameters attributed with the species.
   The 'mutate' simulation will add 1000 mutations to the sequences generated from the 'random' simulation.
   Only the 'mutate' simulation was used for the final publication. ******** DOUBLE CHECK ****
   
10. improvedcodonclassLCRsimulation.py
    
  
11. incslipsim_speciesspecific.py
    
12. nt_pro_cod_genomebias
    
13. periodicityANDEntropy


