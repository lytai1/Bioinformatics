from Bio.Seq import Seq 
from GuideRNA import GuideRNA

class TargetGene:
    def __init__(self, gene):
        super().__init__()
        self.gene = gene
        self.cds = []
        self.mutation = ()
        self.gRNAs = set()
    
    def addCDS(self,start, end):
        self.cds.append((start,end))

    def setMutation(self, start, end):
        self.mutation = (start, end)

    def searchGRNA(self):
        gRNAs = []
        #search for sense gRNA
        index = self.gene.seq.find(Seq("GG"))
        while index >0:
            if index>20:
                gRNA = GuideRNA(self, "", index - 20, index)
                if gRNA.scoreForSearching() > 0:
                    gRNAs.append(gRNA)
            index = self.gene.seq.find("GG", index+2)

        # search for antisense gRNA
        index = self.gene.seq.find(Seq("CC"))
        while (index >0):
            if((index + 20) <len(self.gene.seq)):
                gRNA = GuideRNA(self, "", index + 22, index+2)
                if gRNA.scoreForSearching() > 0:
                    gRNAs.append(gRNA)
            index = self.gene.seq.find(Seq("CC"), index +2)

        return gRNAs

    def addgRNA(self, GuideRNAs):
        gRNAs.update(GuideRNAs)
    
