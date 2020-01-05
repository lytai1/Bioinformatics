from TargetGene import TargetGene
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 
from GuideRNA import GuideRNA
import csv
from Bio.Blast import NCBIXML
from Bio import Entrez

class Utilities:

    @staticmethod
    def readTargetSeq(file):
        records = parse(file, "fasta")
        for record in records:
            return TargetGene(record)

    @staticmethod
    def readgRNASeq(target, file):
        gRNAtxt = file.readlines()
        gRNA = []
        for line in gRNAtxt:
            gRNA.append(GuideRNA(target, Seq(line.strip().upper())))
        return gRNA

    @staticmethod
    def readCDS(target, file):
        cdsTxt = file.readlines()
        cds = ""
        for line in cdsTxt:
            cds += line.strip()

        start = cds.index("(")
        end = cds.index(")")

        cds = cds[start+1: end]
        cdsList = cds.split(",")
        for cds in cdsList:
            items = cds.split("..")
            target.addCDS(int(items[0]),int(items[1]))

    @staticmethod    
    def writeGRNAAnalysis(fileName, gRNAs):
        with open('result/'+ fileName +'.csv', 'w', newline='') as csvfile:
            fieldnames = ['sequence', 'length', 'sense', 'start', 'end', 'distance', 'percent cds', 'repeat bases',
             'No of offtarget', "off targets", "no of match nucleaotides", 'gc content', 'score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for gRNA in gRNAs:
                score = gRNA.score()
                print('sequence'+ gRNA.seq)
                print(gRNA.otResult)
                writer.writerow({'sequence': gRNA.seq, 'length': len(gRNA.seq),'sense': gRNA.sense, 'start': gRNA.start, 'end': gRNA.end,
                                'distance': gRNA.distance(), 'percent cds': gRNA.percentCDS(), 'repeat bases': gRNA.repeativeBases(),
                                 'No of offtarget':gRNA.otResult[0], "off targets":gRNA.otResult[1], "no of match nucleaotides": gRNA.otResult[2],
                                 'gc content': gRNA.getGCContent(), 'score': score})

    @staticmethod
    def gRNAFasta(gRNAs, fileName):
        file = open("gRNA/"+fileName + ".fasta", 'w')
        for i ,gRNA in enumerate(gRNAs):
            file.write(">" + str(i)+"\n")
            file.write(str(gRNA.seq) + "\n")

    @staticmethod
    def openBlastRecords(fileName, gRNAs):
        with open("blastResult/" + fileName, "r") as f:
            for i, record in enumerate(NCBIXML.parse(f)):
                # print("query:" + record.query)
                gRNAs[i].saveBlastRecord(record)
                # print(record.alignments[0].accession)
                # print(record.alignments[0].hsps[0].sbjct, record.alignments[0].hsps[0].sbjct_start, record.alignments[0].hsps[0].sbjct_end)

    
