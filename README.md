# MLPA-Probe-Generating
This project includes two programs written in Python for generating Multiplex Ligation-dependent Probe Amplification (MLPA) probes. MLPA.py generates hybridising sequences (HS) for both halves. It can take multiple sequences of multiple exons as an intake. Criteria for the probe generation are primarly taken from the MLPA probe design protocol from MRC-Holland. LPO_RPO.py takes the hybridising sequences from MLPA.py as an input and attaches primer binding site and stuffer sequence to them, to make the probes fit into the refernce probemix P300-100R by MRC Holland.

The sequences of the exons have to be copied into the text-file input.txt in FASTA format.
The input must be in this format for the code to work properly.
This an example of an possible input:

>Exon17
ATGTCTTCATGACCCAGTTCTCTGCCCTGCAGACAGCTCGATCTGTTCGAACAAGACGGTTGGCAGCTGCAGAGGAAAATATTGAAGTGGCTCGGGCAGCCCGCCTAGCCCAGATCTTCAAAGAAATTTGTGATGGTATCATCTCTTATAAAG
>Exon18
ttttgcagATTCTTCCCGGCAAGCACTGGCAGCTCCACTTTTGAACCTTCCCCCAAAGAAAAAgtaggttc
>Exon19
GAATGCTGATTATTATGAGAAGATCTCTGATCCCCTAGATCTTATCACCATAGAGAAGCAGATCCTCACTGGTTACTATAAGACAGTGGAAGCTTTTGATGCTGACATGCTCAAAGTCTTTCGGAATGCTGAG

MLPA.py can than be started. The generated probes can be seen in the text-file blast_result_final.txt.
The output of MLPA.py would look like this:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exon17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>LHS(Exon17, LHS-Länge: 24, GC-Gehalt: 50%, Tm: 68°C)	CAGCTCGATCTGTTCGAACAAGAC
>RHS(Exon17, RHS-Länge: 24, GC-Gehalt: 50%, Tm: 71°C)	GGTTGGCAGCTGCAGAGGAAAATA

>LHS(Exon17, LHS-Länge: 26, GC-Gehalt: 50%, Tm: 71°C)	CAGACAGCTCGATCTGTTCGAACAAG
>RHS(Exon17, RHS-Länge: 26, GC-Gehalt: 50%, Tm: 74°C)	ACGGTTGGCAGCTGCAGAGGAAAATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exon18 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>LHS(Exon18, LHS-Länge: 20, GC-Gehalt: 50%, Tm: 65°C)	CTGGCAGCTCCACTTTTGAA
>RHS(Exon18, RHS-Länge: 23, GC-Gehalt: 48%, Tm: 66°C)	CCTTCCCCCAAAGAAAAAGTAGG

>LHS(Exon18, LHS-Länge: 22, GC-Gehalt: 50%, Tm: 69°C)	CACTGGCAGCTCCACTTTTGAA
>RHS(Exon18, RHS-Länge: 23, GC-Gehalt: 48%, Tm: 66°C)	CCTTCCCCCAAAGAAAAAGTAGG

One probe consists of the LHS and RHS (Left/Right Hybridising Sequence). The hybridising sequence and additional information about the sequence is shown. For each Exon multiple HS are generated and sorted by the GC-content closest to 50%.

Of each Exon one pair of LHS and RHS can be picked and copied to the text-file auswahl.txt as an input. LPO_RPO.py reads auswahl.txt as an input and attaches primer binding and stuffer sequence to the probes. The lenght of the probes is set to fit in the refernce probemix P300-100R by MRC Holland (p300_geordnet.txt) during the MLPA reaction.

A possible input for auswahl.txt:

>LHS(Exon17, LHS-Länge: 24, GC-Gehalt: 50%, Tm: 68°C)	CAGCTCGATCTGTTCGAACAAGAC
>RHS(Exon17, RHS-Länge: 24, GC-Gehalt: 50%, Tm: 71°C)	GGTTGGCAGCTGCAGAGGAAAATA
>LHS(Exon18, LHS-Länge: 20, GC-Gehalt: 50%, Tm: 65°C)	CTGGCAGCTCCACTTTTGAA
>RHS(Exon18, RHS-Länge: 23, GC-Gehalt: 48%, Tm: 66°C)	CCTTCCCCCAAAGAAAAAGTAGG
>LHS(Exon19, LHS-Länge: 24, GC-Gehalt: 46%, Tm: 61°C)	CTCTGATCCCCTAGATCTTATCAC
>RHS(Exon19, RHS-Länge: 24, GC-Gehalt: 50%, Tm: 65°C)	CATAGAGAAGCAGATCCTCACTGG


The text-file MLPA_SONDEN.txt shows the final MLPA probes:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exon18 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>LPO:
GGGTTCCCTAAGGGTTGGACGCTACTACTACTGGCAGCTCCACTTTTGAA
-------------------+++++++++++********************

>RPO:
CCTTCCCCCAAAGAAAAAGTAGGCTACTCTAGATTGGATCTTGCTGGCAC
***********************++++-----------------------
Gesamtlänge: 100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exon17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>LPO:
GGGTTCCCTAAGGGTTGGACGCTACTACTATTAGCAGCTCGATCTGTTCGAACAAGAC
-------------------+++++++++++++++************************

>RPO:
GGTTGGCAGCTGCAGAGGAAAATAACTAAATCTACTCTAGATTGGATCTTGCTGGCAC
************************+++++++++++-----------------------
Gesamtlänge: 116

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exon19 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>LPO:
GGGTTCCCTAAGGGTTGGACGCTACTACTATTAGTAGAATTGATGCTCTGATCCCCTAGATCTTATCAC
-------------------++++++++++++++++++++++++++************************

>RPO:
CATAGAGAAGCAGATCCTCACTGGCTAATGGTCAAACTAAATCTACTCTAGATTGGATCTTGCTGGCAC
************************++++++++++++++++++++++-----------------------
Gesamtlänge: 138

* marks the hybridising sequence
+ marks the stuffer sequence
- marks the primer binding sequence
