import gzip
import pysam
import numpy as np

from .seq import Seq

class Locus:
    def __init__(self, name, start, end, strand, offset=1):
        self.name = name
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.offset = offset
        
    def length(self):
        return self.end - self.start + self.offset
        
    def region_string(self):
        return '{}:{}-{}'.format(self.name, self.start, self.end)

        
class Gene:
    def __init__(self, gene_id, gene_name, gene_type, locus):
        self.id = gene_id
        self.name = gene_name
        self.type = gene_type
        self.pos = locus

        
class Exon:
    def __init__(self, exon_id, exon_number, locus):
        self.id = exon_id
        self.number = exon_number
        self.pos = locus
        self.cds_pos = None
        self.cds_frame = -1
        
    def set_cds(self, locus, frame):
        self.cds_pos = locus
        self.cds_frame = frame


class Transcript:
    def __init__(self, gene, transcript_id, transcript_name, transcript_type, locus):
        self.gene = gene
        self.id = transcript_id
        self.name = transcript_name
        self.type = transcript_type
        self.pos = locus
        self.exons = dict()
        self.is_coding = False
        
    def get_exon_order(self):
        """ Return exon order based on strand """
        if self.pos.strand == "-":
            return range(len(self.exons), 0, -1)
        else:
            return range(1, len(self.exons)+1)

    def get_mrna_seq(self, genome):
        """ Return mRNA sequence """
        seq = ''
        fa = pysam.FastaFile(genome)
        for i in self.get_exon_order():
            seq += fa.fetch(region=self.exons[i].pos.region_string())
        fa.close()
        mrna = Seq(seq)
        if self.pos.strand == "-":
            mrna.reverse_complement()
        fa.close()
        return mrna.seq

    
    def get_cds_seq(self, genome):
        """ Return coding sequence """
        seq = ''
        if not self.is_coding:
            return seq
        fa = pysam.FastaFile(genome)
        for i in self.get_exon_order():
            if self.exons[i].cds_pos == None:
                continue
            seq += fa.fetch(region=self.exons[i].cds_pos.region_string())
        fa.close()
        cds = Seq(seq)
        if self.pos.strand == "-":
            cds.reverse_complement()
        fa.close()
        return cds.seq        

    
    def raw_cds_with_downstream(self, genome):
        """
        Return CDS and the downstream UTR genomic
        sequences and positions as arrays, and
        exon boundary string as a string
        """
        seq = []
        pos = []
        exon_str = ''
        if not self.is_coding:
            return seq, pos, exon_str

        # find an exon where the CDS ends
        exon_cds_end = -1
        exon_numbers = list(self.get_exon_order())
        if self.pos.strand == '+':
            exon_numbers = exon_numbers[::-1]
        for i in exon_numbers:
            if self.exons[i].cds_pos == None:
                continue
            exon_cds_end = i
            break
        if exon_cds_end == -1:
            return seq, pos, exon_str
        
        fa = pysam.FastaFile(genome)
        for i in self.get_exon_order():            
            if i > exon_cds_end:
                region = self.exons[i].pos
            elif i == exon_cds_end:
                if self.pos.strand == '+':
                    region = Locus(self.pos.name,
                                   self.exons[i].cds_pos.start,
                                   self.exons[i].pos.end,
                                   self.pos.strand)
                else:
                    region = Locus(self.pos.name,
                                   self.exons[i].pos.start,
                                   self.exons[i].cds_pos.end,
                                   self.pos.strand)
            elif self.exons[i].cds_pos != None:
                region = self.exons[i].cds_pos
            else:
                continue

            seq.append(list(fa.fetch(region=region.region_string())))
            pos.append(region)
            exon_str += '|{}|'.format(' '*(region.length()-2))

        fa.close()

        # remove the first exon boundary character
        if self.pos.strand == '+':
            exon_str = ' ' + exon_str.lstrip('|')
        else:
            exon_str = exon_str.rstrip('|') + ' '
            
        return seq, pos, exon_str

    
class Annotation:
    """
    Parse GTF file and create transcript data
    """
    def __init__(self, gtf_path, transcript_list=[]):
        self.transcripts = dict()
        
        try:
            if gtf_path.endswith('.gz'):
                gtf = gzip.open(gtf_path, 'rt')
            else:
                gtf = open(gtf_path, 'r')
        except IOError:
            raise Exception('cannot read: {}'.format(gtf_path))

        skip_transcript = False
        transcript_id = None
        gene_read = 0
        for row in gtf:
            # skip comment lines
            if row.startswith('#'):
                continue
            row = row.rstrip('\n').split('\t')
            
            # Locus(chrom, start, end, strand)
            pos = Locus(row[0], row[3], row[4], row[6])
            feature = row[2]
            frame = row[7]

            attr = dict()
            for a in row[8].replace('"', '').replace('_biotype', '_type').split(';')[:-1]:
                v = a.strip().split(' ')
                if v[0] != 'tag':
                    attr.setdefault(v[0], v[1])
                else:
                    attr.setdefault(v[0], []).append(v[1])

            if feature == 'gene':
                if 'gene_id' not in attr:
                    raise Exception('"gene_id" not found: {}'.format('\t'.join(row)))
                if 'gene_name' not in attr:
                    attr['gene_name'] = attr['gene_id']
                g = Gene(attr['gene_id'], attr['gene_name'], attr['gene_type'], pos)
                gene_read += 1
                
            elif feature == 'transcript':
                if 'transcript_id' not in attr:
                    raise Exception('"transcript_id" not found: {}'.format('\t'.join(row)))
                transcript_id = attr['transcript_id']
                if 'transcript_name' not in attr:
                    attr['transcript_name'] = transcript_id
                t = Transcript(g, transcript_id, attr['transcript_name'], attr['transcript_type'], pos)
                # include all transcripts or specified transcripts
                if transcript_list != [] and transcript_id not in transcript_list:
                    skip_transcript = True
                else:
                    skip_transcript = False
                    self.transcripts.setdefault(transcript_id, t)

            elif feature == 'exon':
                if skip_transcript:
                    continue
                if 'exon_id' not in attr:
                    attr['exon_id'] = str(len(t.exons)+1)
                if 'exon_number' not in attr:
                    attr['exon_number'] = len(t.exons)+1
                e = Exon(attr['exon_id'], attr['exon_number'], pos)
                self.transcripts[transcript_id].exons[int(attr['exon_number'])] = e

            elif feature == 'CDS':
                if skip_transcript:
                    continue
                if 'exon_number' not in attr:
                    attr['exon_number'] = len(t.exons)
                try:
                    frame = int(frame)
                except ValueError:
                    raise Exception('CDS frame not defined: {}'.format('\t'.join(row)))
                self.transcripts[transcript_id].exons[int(attr['exon_number'])].set_cds(pos, frame)
                self.transcripts[transcript_id].is_coding = True

