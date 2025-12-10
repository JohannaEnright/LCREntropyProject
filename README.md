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
   This script is a modified version of the original Seg algorithm created by Wootton and Federhen (1993). It has been modified to extend to longer LCR sequences than were previously capable of being handled by Seg. As well, it has an added alphabet to identify codon LCRs and deal with ambiguous letters. segA is called automatically in the CombinedUsingSysCall_Codon3.py script.
   
3. GiveCombSysCallArgs.sh
   This is a bash file which is used to feed CombinedUsingSysCall_Codon3.py a combination of set parameters for segA. This was important to ensure that we were observing the same trends given slightly different LCR identifying parameters, considering the somewhat arbitrary nature of these parameters. In total, it will feed CombinedUsingSysCall_Codon3.py 27 unique parameter combinations as shown in the supplementary section of the publication (Enright, Dickson, & Golding, 2023).
   
4. averageproteinlength.py
   This script requires a genbank file as input, and will return the average protein length for that file.
   Again, it also accepts multiple concatenated genbank files, hence it was used to calculate the average protein length for each organism.
   This information was used in the generaterandomproteome script.

5. codonclass_categorizeLCRs.py
   This script requires the _Seq file outputed by CombinedUsingSysCall_Codon3.py which contains the sequences of all of the detected LCRs for a given organism, and a file output location.
   It categorizes each LCR into one of three categories:
     i) containing three unique nucleotides, ii) containing two unique nucleotides, and iii) containing one unique nucleotide and outputs this information as counts of LCRs in each category.
   This information is used for the calculation of the preference coefficient, see Enright, Dickson, & Golding (2023).
   
6. countcodontypes_allCDS.py
   This script takes a genbank (or multiple) file as input. It outputs a text file with information for the total codons of each codon class summed for all sequences in a genome.
   This information is used for the calculation of the preference coefficient, see Enright, Dickson, & Golding (2023).
    
7. generaterandomproteome.py
   This script is to simulate 100 000 sequences outputed in fasta file format. 
   The input values are the simulation type ('mutate' or 'random'), species ('sc', 'dm', 'hs', 'ce', 'at'), fasta file of sequences to mutate (only if you choose the 'mutate' option),
   write file (name of file to write to).
   The species input is necessary to generate sequences which have the average length, the correct codon class proportions, and the same proportion of LCRs as the corresponding species.
   The 'random' simulation will randomly generate 100 000 random sequences given the various parameters attributed with the species.
   The 'mutate' simulation will add 1000 mutations to the sequences generated from the 'random' simulation.
   Only the 'mutate' simulation was used in the final publication to add 1000 synonomous mutations to create the 'Slip+CC+Syn' model.   
   
8. improvedcodonclassLCRsimulation.py
    It will simulate 100 000 LCR sequences using an increasing DNA polyermase slippage model that also takes into account the codonclass proportions (class 1,2 or 3).
    Codon class proportion biases are based off of precalculated preference coefficients for each species, see Enright, Dickson, & Golding (2023).
    This model is reffered to as the 'Slip+CC' model.
    The input requirements are a file of secies specific codon bias proportions Outputed by the nt_cod_aa_bias.py script and species type ('hs', 'sc', 'dm', 'at', 'ce').
    It generates sequences using the average protein length, codon proportion, and preference coefficient for biological relevance.
  
9. incslipsim_speciesspecific.py
    This simulation generates 100 000 sequences protein and corresponding DNA sequences, outputed as a genbank file. It creates LCRs simulating the DNA polyermase slippage mechanism,
    where the probability of resampling the same codon increases with increasing number of repeats of that codon. Sequences are generated such that they include the same proportion of LCRs, the same protein length, and same codon bias as are biologically found in the organism.
    The script requires the inputs:
    i) simulation type ('random', 'biased' or 'mutate') where 'random' randomly generates sequences not using an increasing slippage model,
     'biased' generates sequences using an increasing slippage model (referred to as the 'Slip model' in the publication),
      and 'mutate' adds 1000 random mutations to the third position of codons from the generated sequences (reffered to as 'Slip + Syn model' in the publication),
      using biologially relevant weights (5:1) for transition/tranversion mutations.
    ii) species ('hs', 'sc', 'dm', 'at', 'ce')
    iii) file of codon proportions as generated from nt_pro_cod_genomebias.py

10. nt_pro_cod_genomebias
    This script takes a genbank file (or multiple concatenated) as input and outputs 4 files: each containing the nucleotide frequency, the nucleotide in coding sequence frequency, the codon frequency, and the amino acid frequency.
    The codon frequency output was used for seqeunce simulations such as improvedcodonclassLCRsimulation.py and incslipsim_speciesspecific.py.
    
11. periodicityANDEntropy
    This script requires a _Seq file outputed by CombinedUsingSysCall_Codon3.py containing the sequences of the detected biological LCRs and the sequence type, either 'protein' or 'dna'.
    Optional inputs include:
       i) wholeseq="True" or "False, which outputs either the entropies of only the repeat segments ("False") or the entropies of the entire LCR sequence ("True").
       ii) The lengths that each repeat has to be to be considered a repeat. For example, for mono-amino acid repeats, the releat size needed to be 4 (rep1p="4"), and the length of tri-nucleotide repeats needed to be 4 (rep3d="4").
    It outputs 2 files:
       i) one contains the DNA, codon, and protein entropy information for LCR sequences which did not contain periodic repeats.
       ii) the second contains this information for sequences which do contain periodic reapeats, the index positions where the repeat begins and ends, as well as the repeat type, for example mono-, di-, or tri- amino acid repeats.
    This information was used to compare the correlation of entropies between LCRs which contained periodic repeats and those which did not.


