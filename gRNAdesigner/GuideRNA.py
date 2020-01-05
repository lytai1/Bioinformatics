from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import timeit
from Bio import Entrez

class GuideRNA:
    def __init__(self, gene, seq = None, start = -1, end = -1):
        super().__init__()
        self.record = None
        self.otResult = [-1,0,0] #off target result
        self.target = gene #initialize target gene

        if start >= 0  and end >= 0: #if start and end are given
            self.start = start
            self.end = end
            if start<end:
                self.seq = self.target.gene.seq[start:end]
            else: self.seq = self.target.gene.seq[end:start].reverse_complement()
        elif seq != None: # if sequence is given
            self.seq = seq
            self.start = gene.gene.seq.find(seq)
            self.end = self.start + len(seq)
            if self.start <0:
                seq = seq.reverse_complement()
                self.end = gene.gene.seq.find(seq)
                if self.end <0:
                    raise LookupError("cannot locate sequence " + self.seq)
                self.start = self.end + len(seq)
        else: #error 
            raise ValueError("sequence is null or start, end does not input correctly")

        #set sense
        if self.start < self.end:
            self.sense = True
        else:
            self.sense = False

    def saveBlastRecord(self, record):
        self.record = record

    def getGCContent(self):
        # print(GC(self.seq))
        return GC(self.seq)

    def repeativeBases(self):
        if self.seq.count("AAAAA") > 0:
            return True
        elif self.seq.count("TTTT") > 0:
            return True
        elif self.seq.count("CCCCC") > 0:
            return True
        elif self.seq.count("GGGG") > 0:
            return True
        else: return False

    def percentCDS(self):
        if self.sense:
            for (start, end) in self.target.cds:
                if self.start >= start and self.start<=end:
                    return min(end-self.start, len(self.seq))/len(self.seq)
                elif self.end >= start and self.end<=end:
                    return min(self.end - start, len(self.seq))/len(self.seq)
            return 0
        else:
            for (start, end) in self.target.cds:
                if self.end >=start and self.end<=end:
                    return min(end-self.end, len(self.seq))/len(self.seq)
                elif self.start >= start and self.start<=end:
                    return min(self.start - start, len(self.seq))/len(self.seq)
            return 0
        
    def distance(self):
        location = -1
        for cds in self.target.cds:
            if self.start >=cds[0] and self.start<=cds[1]:
                location = self.cdsLocation(self.start)
                return min(abs(self.target.mutation[0] - location),abs(self.target.mutation[1] - location))
            elif self.end >= cds[0] and self.end<=cds[1]:
                location = self.cdsLocation(self.end)
                return min(abs(self.target.mutation[0] - location),abs(self.target.mutation[1] - location))
        raise ValueError("gRNA not in cds")
       
        
    
    def cdsLocation(self, location):
        position = 0
        i = 0
        while i<len(self.target.cds) and location>self.target.cds[i][0]:
            if location <= self.target.cds[i][0] and location>=self.target.cds[i][1]:
                position += location - self.target.cds[i][0]
                break
            else:
                position += self.target.cds[i][1] - self.target.cds[i][0]
            i +=1

        return position 
    
    def alignWithExon(self, aln):
        accession = aln.accession
        start = aln.hsps[0].sbjct_start
        end = aln.hsps[0].sbjct_end

        #read from genebank for gene record
        email = "kk0kathleen@gmail.com"
        handles = Entrez.efetch(db="nucleotide",
                                email=email,
                                id=accession,
                                rettype="gb",
                                retmode="txt")
        data = handles.read()
        handles.close()
        print(accession)
        print("start: " + str(start))
        print("end: " + str(end))
        for i,line in enumerate(data.split()):
            if "exon" in line:
                line = data.split()[i+1]
                if line[0] in "1234567890" and ".." in line:
                    exon = line.split("..")
                    exon[0] = int(exon[0])
                    exon[1] = int(exon[1])
                    if end>start:
                        if start > exon[0] and start < exon[1]:
                            intersect = min(exon[1] - start, end-start)
                            return [True, intersect]
                        elif end > exon[0] and end <exon[1]:
                            intersect = min(end - exon[0], end - start)
                            return [True, intersect]
                    else:
                        if end > exon[0] and end < exon[1]:
                            intersect = min(exon[1] - end, start-end)
                            return [True, intersect]
                        elif start > exon[0] and start <exon[1]:
                            intersect = min(start - exon[0], start - end)
                            return [True, intersect]
        return [False, 0]

    def offTarget(self):
        totalOffTargets = 0
        matches = []
        matchesNo = []
        for algn in self.record.alignments:
            result = self.alignWithExon(algn)
            if result[0]:
                totalOffTargets += 1
                matches.append(algn.accession)
                matchesNo.append(result[1])

        self.otResult = [totalOffTargets, matches, matchesNo]
            
    def scoreForSearching(self):
        if self.repeativeBases():
            return 0
        gc = self.getGCContent()
        if gc<30 or gc>70:
            return 0
        percentCDS = self.percentCDS()
        if percentCDS<=0:
            return 0
        score = 10000
        return score

    def score(self):
        print(self.seq)
        if self.repeativeBases():
            return 0
        gc = self.getGCContent()
        if gc<30 or gc>70:
            return 0
        percentCDS = self.percentCDS()
        if percentCDS<=0:
            return 0
        self.offTarget()
        match = self.otResult[0]
        if match == 0:
            score = 100000*percentCDS/self.distance()
        else:
            score = 10000*percentCDS/self.distance()/match
        return score
        