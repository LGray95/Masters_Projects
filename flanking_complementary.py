import pysam
import os.path
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from io import StringIO

# Creates fasta file with both upstream and downstream sequences
def upstream_and_downstream_seq(coords, out_path):
    chromosome = coords.split(":")[0]
    start = str(int(coords.split(":")[1].split("-")[0]))
    downstream = str(int(coords.split(":")[1].split("-")[0])-500)
    end = str(int(coords.split(":")[1].split("-")[1]))
    upstream = str(int(coords.split(":")[1].split("-")[1])+500)

    #using the samtools faidx function to take the appropriate sequence from a reference genome
    downstream_fa = Seq(pysam.faidx("/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/ICM/GRCh38_genome.fa", chromosome+":"+downstream+"-"+start), generic_dna)

    upstream_fa = Seq(pysam.faidx("/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/ICM/GRCh38_genome.fa", chromosome+":"+end+"-"+upstream), generic_dna)

    # Selecting only the sequence and converting to uppercase
    downstream_seq = downstream_fa[(len(downstream_fa.split('\n')[0])):-1].upper()
    # Selecting only the sequence, converting to uppercase, reversing and then getting the complementary sequence
    reverse_compliment_upstream_seq = upstream_fa[(len(upstream_fa.split('\n')[0])):-1].upper().reverse_complement()

    # Making sequence records with ID header and sequence
    downstream_seq = SeqRecord(downstream_seq, id="downstream_sequence")
    reverse_compliment_upstream_seq = SeqRecord(reverse_compliment_upstream_seq, id="upstream_sequence")

    # Writing sequences to fasta file
    downstream_outfile = open(os.path.join(out_path, "downstream.fa"), "w")
    downstream_outfile.write(">"+str(downstream_seq.id) + "\n" + str(downstream_seq.seq))

    upstream_outfile = open(os.path.join(out_path, "upstream.fa"), "w")
    upstream_outfile.write(">"+str(reverse_compliment_upstream_seq.id) + "\n" + str(reverse_compliment_upstream_seq.seq))

infile = open("/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/ICM/Fibroblast/enriched_IDs.txt").read().splitlines()
for line in infile[0:-1]:
    upstream_and_downstream_seq(line, '/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/ICM/Fibroblast/tmp/')

    def complementary_pairing(inpath, outpath):
        output = NcbiblastnCommandline(query=inpath + 'downstream.fa', subject=inpath + 'upstream.fa', word_size=7, outfmt=5,
                                       task="blastn")()[0]
        blast_result_record = NCBIXML.read(StringIO(output))


        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                # print('****Alignment****')
                # print('length:', alignment.length)
                # print('e value:', hsp.expect)
                # print(hsp.query)
                # print(hsp.match)
                # print(Seq(hsp.sbjct).reverse_complement())

                # Setting match variable
                match = hsp.query

                # Calculating GC content of each match
                def get_GC_content(sequence):
                    gc_count = sequence.count("G") + sequence.count("C")
                    gc_fraction = float(gc_count) / len(sequence)
                    ll = round((gc_fraction * 100),2)
                    return ll

                GC_content = get_GC_content(match)

                # Set conditional argument for match length if >= 20 nt
                if len(match) >= 20:

                    # Parsing in FASTA files and selecting the sequence
                    # Removing '-' from string and replacing each match with lowercase
                    for SeqRecord in SeqIO.parse(inpath + 'downstream.fa', 'fasta'):
                        downstream_inseq = SeqRecord.seq
                        downstream_match = str(match.replace("-", ""))
                        downstream_outseq = str(downstream_inseq).replace(downstream_match, downstream_match.lower())

                    for SeqRecord in SeqIO.parse(inpath + 'upstream.fa', 'fasta'):
                        upstream_inseq = str(SeqRecord.seq)
                        upstream_outseq = str(Seq(upstream_inseq).reverse_complement())
                        upstream_match = str(Seq(hsp.sbjct).reverse_complement()).replace("-", "")
                        masked_upstream = str(upstream_outseq).replace(upstream_match, upstream_match.lower())
                    filename = line.replace(":", "_")
                    outfile = open(outpath + filename + ".fa", 'w')
                    outfile.write("e value: " + str(hsp.expect) + "\n"
                                  + "Match length: " + str(len(match)) + '\n'
                                  + "GC_content: " + str(GC_content) + "\n"
                                  + "Downstream match: " + "\n"
                                  + str(hsp.query).replace("-", "") + "\n"
                                  + "Upstream match: " + "\n"
                                  + str(Seq(hsp.sbjct).reverse_complement()).replace("-", "") + "\n"
                                  + ">Downstream_seq" + "\n" + str(downstream_outseq) + "\n"
                                  + ">Upstream_seq" + "\n" + str(masked_upstream))
    complementary_pairing('/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/ICM/Fibroblast/tmp/', '/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/ICM/Fibroblast/enriched_fasta/')
