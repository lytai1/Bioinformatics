from Utilities import Utilities
from GuideRNA import GuideRNA
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import csv

def main():
    file = open("data/huntingtin_gene.fasta") 
    target = Utilities.readTargetSeq(file)
    target.setMutation(51,141)
    file = open("data/cds_huntingtin.txt")
    Utilities.readCDS(target, file)

    # file = open("data/Wu_result.txt", "r")
    # gRNAs = Utilities.readgRNASeq(target, file)
    # # Utilities.gRNAFasta(gRNAs, "Wu")
    # Utilities.openBlastRecords("Wu-Alignment.xml", gRNAs)
    # Utilities.writeGRNAAnalysis("WU", gRNAs)
    
    file = open("data/ChopChopresult.txt", "r")
    gRNAs = Utilities.readgRNASeq(target, file)
    # Utilities.gRNAFasta(gRNAs, "ChopChop")
    Utilities.openBlastRecords("ZGUG078K01N-Alignment.xml", gRNAs)

    Utilities.writeGRNAAnalysis("ChopChop_1181_1246", gRNAs[1181:])
    
    # file = open("data/gpp_result.txt", "r")
    # gRNAs = Utilities.readgRNASeq(target, file)
    # # Utilities.gRNAFasta(gRNAs,"gpp")
    # Utilities.openBlastRecords("gpp-Alignment.xml", gRNAs)
    # Utilities.writeGRNAAnalysis("gpp", gRNAs)
    
    # gRNAs = target.searchGRNA()
    # # Utilities.gRNAFasta(gRNAs,"gRNA")
    # Utilities.openBlastRecords("gRNA-Alignment.xml", gRNAs)
    # Utilities.writeGRNAAnalysis("gRNA", gRNAs)
   
main()
