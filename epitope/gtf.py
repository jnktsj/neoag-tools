from .seq import Seq
import gzip
from typing import Any, Dict, List, Literal, Optional, Tuple
import pysam
from dataclasses import dataclass, field

"""
GTF organizes the transcript data hierarchically, with the organization levels being:
1. Gene
2. Transcript
3. Exon

The order by which the rows introduce new elements is a DFS on this tree.
Every element is associated with a span of the genome, corresponding to the coding sequence for an exon, 
or the span containing all child nodes for an object higher in the hierarchy.
"""


@dataclass
class Locus:
    """
    A substring of one of the DNA strands.
    """

    contig: str
    """The name of the locus."""
    start: int
    """The start position."""
    end: int
    """The end position."""
    strand: Literal["+", "-"]
    """Whether the locus is on the + or - strand."""
    offset: int = field(default=1)

    def length(self):
        return self.end - self.start + self.offset

    def region_string(self) -> str:
        return "{}:{}-{}".format(self.contig, self.start, self.end)


@dataclass
class Gene:
    id: str
    name: str
    type: str
    pos: Locus


@dataclass
class Exon:
    id: str
    exon_number: int
    pos: Locus
    cds_pos: Optional[Locus] = field(default=None)
    cds_frame: int = field(default=-1)

    def set_cds(self, locus: Locus, frame: int):
        self.cds_pos = locus
        self.cds_frame = frame


@dataclass
class Transcript:
    gene: Gene
    id: str
    name: str
    type: str
    pos: Locus
    exons: Dict[int, Exon] = field(default_factory=dict)
    is_coding: bool = field(default=False)

    def get_exon_order(self) -> range:
        """Return exon order based on strand"""
        if self.pos.strand == "-":
            return range(len(self.exons), 0, -1)
        else:
            return range(1, len(self.exons) + 1)

    def get_mrna_seq(self, genome):
        """Return mRNA sequence"""
        seq = ""
        fa = pysam.FastaFile(genome)
        for i in self.get_exon_order():
            seq += fa.fetch(region=self.exons[i].pos.region_string())
        fa.close()
        mrna = Seq(seq)
        if self.pos.strand == "-":
            mrna = mrna.reverse_complement()
        fa.close()
        return mrna.seq

    def get_cds_seq(self, genome: str) -> str:
        """
        Reads the coding sequence of the given transcript from the given genome file.
        """
        seq = ""
        if not self.is_coding:
            return seq
        fa = pysam.FastaFile(genome)
        for i in self.get_exon_order():
            cds_pos = self.exons[i].cds_pos
            if cds_pos is None:
                continue
            seq += fa.fetch(region=cds_pos.region_string())
        fa.close()
        cds = Seq(seq)
        if self.pos.strand == "-":
            cds = cds.reverse_complement()
        fa.close()
        return cds.seq

    def raw_cds_with_downstream(
        self, genome: str
    ) -> Tuple[List[List[str]], List[Locus], str]:
        """
        Args:
            genome (str): The path to the reference genome.

        Returns:
            Tuple:
             * The CDS and downstream UTR, represented as a list of the sequences, each of which is a list of nucleotides.
             * The list of loci of the CDS and downstream UTR.
             * The exon boundary string as a string.
        """
        seq: List[List[str]] = []
        pos: List[Locus] = []
        exon_str = ""
        if not self.is_coding:
            return seq, pos, exon_str

        # find an exon where the CDS ends
        exon_cds_end = -1
        exon_numbers = list(self.get_exon_order())
        if self.pos.strand == "+":
            exon_numbers = exon_numbers[::-1]
        for i in exon_numbers:
            if self.exons[i].cds_pos is None:
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
                cds_pos = self.exons[i].cds_pos
                assert cds_pos is not None
                if self.pos.strand == "+":
                    region = Locus(
                        self.pos.contig,
                        cds_pos.start,
                        self.exons[i].pos.end,
                        self.pos.strand,
                    )
                else:
                    region = Locus(
                        self.pos.contig,
                        self.exons[i].pos.start,
                        cds_pos.end,
                        self.pos.strand,
                    )
            elif self.exons[i].cds_pos is not None:
                region = self.exons[i].cds_pos
            else:
                continue
            assert region is not None
            seq.append(list(fa.fetch(region=region.region_string())))  # type: ignore
            pos.append(region)
            exon_str += "|{}|".format(" " * (region.length() - 2))

        fa.close()

        # remove the first exon boundary character
        if self.pos.strand == "+":
            exon_str = " " + exon_str.lstrip("|")
        else:
            exon_str = exon_str.rstrip("|") + " "

        return seq, pos, exon_str


class Annotation:
    """
    Parse GTF file and create transcript data
    """

    transcripts: Dict[str, Transcript]

    def __init__(self, gtf_path, transcript_list=[]):
        self.transcripts = dict()

        try:
            if gtf_path.endswith(".gz"):
                gtf = gzip.open(gtf_path, "rt")
            else:
                gtf = open(gtf_path, "r")
        except IOError:
            raise Exception("cannot read: {}".format(gtf_path))

        skip_transcript = False
        transcript_id = None
        gene_read = 0
        g = None
        t = None
        for row in gtf:
            # skip comment lines
            if row.startswith("#"):
                continue
            row = row.rstrip("\n").split("\t")

            # Locus(chrom, start, end, strand)
            assert row[6] in ("-", "+"), f"Illegal strand: {row[6]}"
            pos = Locus(row[0], int(row[3]), int(row[4]), row[6])
            feature = row[2]
            frame = row[7]

            attr = dict()
            for a in (
                row[8]
                .replace('"', "")
                .replace("_biotype", "_type")
                .split(";")[:-1]
            ):
                v = a.strip().split(" ")
                if v[0] != "tag":
                    attr.setdefault(v[0], v[1])
                else:
                    attr.setdefault(v[0], []).append(v[1])

            if feature == "gene":
                if "gene_id" not in attr:
                    raise Exception(
                        '"gene_id" not found: {}'.format("\t".join(row))
                    )
                if "gene_name" not in attr:
                    attr["gene_name"] = attr["gene_id"]
                g = Gene(
                    attr["gene_id"], attr["gene_name"], attr["gene_type"], pos
                )
                gene_read += 1

            elif feature == "transcript":
                if "transcript_id" not in attr:
                    raise Exception(
                        '"transcript_id" not found: {}'.format("\t".join(row))
                    )
                transcript_id = attr["transcript_id"]
                if "transcript_name" not in attr:
                    attr["transcript_name"] = transcript_id
                assert g is not None, "We have not yet set encountered a gene!"
                t = Transcript(
                    g,
                    transcript_id,
                    attr["transcript_name"],
                    attr["transcript_type"],
                    pos,
                )
                # include all transcripts or specified transcripts
                if (
                    transcript_list != []
                    and transcript_id not in transcript_list
                ):
                    skip_transcript = True
                else:
                    skip_transcript = False
                    self.transcripts.setdefault(transcript_id, t)

            elif feature == "exon":
                assert t is not None
                if skip_transcript:
                    continue
                if "exon_id" not in attr:
                    attr["exon_id"] = str(len(t.exons) + 1)
                if "exon_number" not in attr:
                    attr["exon_number"] = len(t.exons) + 1
                e = Exon(attr["exon_id"], attr["exon_number"], pos)
                self.transcripts[transcript_id].exons[
                    int(attr["exon_number"])
                ] = e

            elif feature == "CDS":
                if skip_transcript:
                    continue
                if "exon_number" not in attr:
                    attr["exon_number"] = len(t.exons)
                try:
                    frame = int(frame)
                except ValueError:
                    raise Exception(
                        "CDS frame not defined: {}".format("\t".join(row))
                    )
                self.transcripts[transcript_id].exons[
                    int(attr["exon_number"])
                ].set_cds(pos, frame)
                self.transcripts[transcript_id].is_coding = True
