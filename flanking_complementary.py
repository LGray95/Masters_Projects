import argparse
from split_coordinates import split_coords
import os
import pysam
import os.path
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from io import StringIO


def argparser():
    parser = argparse.ArgumentParser(description='Ask for path')
    parser.add_argument("-c", "--coords", help = "circRNA coordinates")
    parser.add_argument("-g", "--genome", help="path to genome file")
    parser.add_argument("-d", "--directory", help="Path to working directory")
    args = parser.parse_args()
    return args


# Calculating GC content of each match
def get_GC_content(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    gc_fraction = float(gc_count) / len(sequence)
    gc_content = round((gc_fraction * 100),2)
    return gc_content


# Creates fasta file with both upstream and downstream sequences
def upstream_and_downstream_seq(args):
    chromosome = split_coords(args.coords)[0]
    start = str(split_coords(args.coords)[1])
    downstream = str(int(start)-1000)
    end = str(split_coords(args.coords.replace('"', ""))[2])
    upstream = str(int(end)+1000)

    #using the samtools faidx function to take the appropriate sequence from a reference genome
    downstream_fa = Seq(pysam.faidx(args.genome, chromosome+":"+downstream+"-"+start), generic_dna)

    upstream_fa = Seq(pysam.faidx(args.genome, chromosome+":"+end+"-"+upstream), generic_dna)

    # Selecting only the sequence and converting to uppercase
    downstream_seq = downstream_fa[(len(downstream_fa.split('\n')[0])):-1].upper()
    # Selecting only the sequence, converting to uppercase, reversing and then getting the complementary sequence
    reverse_compliment_upstream_seq = upstream_fa[(len(upstream_fa.split('\n')[0])):-1].upper().reverse_complement()

    # Making sequence records with ID header and sequence
    downstream_seq = SeqRecord(downstream_seq, id="downstream_sequence")
    reverse_compliment_upstream_seq = SeqRecord(reverse_compliment_upstream_seq, id="upstream_sequence")

    if os.path.isdir(args.directory+"tmp/") == False:
        os.mkdir(args.directory+"tmp/")

    # Writing sequences to fasta file
    downstream_outfile = open(os.path.join(args.directory+"tmp/", "downstream.fa"), "w")
    downstream_outfile.write(">"+str(downstream_seq.id) + "\n" + str(downstream_seq.seq))

    upstream_outfile = open(os.path.join(args.directory+"tmp/", "upstream.fa"), "w")
    upstream_outfile.write(">"+str(reverse_compliment_upstream_seq.id) + "\n" + str(reverse_compliment_upstream_seq.seq))


def complementary_pairing(args):
    output = NcbiblastnCommandline(query=args.directory+"tmp/" + 'downstream.fa',
                                   subject=args.directory+"tmp/" + 'upstream.fa',
                                   word_size=7, outfmt=5, task="blastn")()[0]
    blast_result_record = NCBIXML.read(StringIO(output))

    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:

            # Setting match variable
            match = hsp.query

            gc_content = get_GC_content(match)

            # Set conditional argument for match length if >= 20 nt
            if len(match) >= 20:

                # Parsing in FASTA files and selecting the sequence
                # Removing '-' from string and replacing each match with lowercase
                for SeqRecord in SeqIO.parse(args.directory+"tmp/" + 'downstream.fa', 'fasta'):
                    downstream_inseq = SeqRecord.seq
                    downstream_match = str(match.replace("-", ""))
                    downstream_outseq = str(downstream_inseq).replace(downstream_match, downstream_match.lower())

                for SeqRecord in SeqIO.parse(args.directory+"tmp/" + 'upstream.fa', 'fasta'):
                    upstream_inseq = str(SeqRecord.seq)
                    upstream_outseq = str(Seq(upstream_inseq).reverse_complement())
                    upstream_match = str(Seq(hsp.sbjct).reverse_complement()).replace("-", "")
                    masked_upstream = str(upstream_outseq).replace(upstream_match, upstream_match.lower())
                filename = args.coords.replace(":", "_")
                outfile = open(args.directory + filename + ".fa", 'w')
                outfile.write("e value: " + str(hsp.expect) + "\n"
                              + "Match length: " + str(len(match)) + '\n'
                              + "GC_content: " + str(gc_content) + "\n"
                              + "Downstream match: " + "\n"
                              + str(hsp.query).replace("-", "") + "\n"
                              + "Upstream match: " + "\n"
                              + str(Seq(hsp.sbjct).reverse_complement()).replace("-", "") + "\n"
                              + ">Downstream_seq" + "\n" + str(downstream_outseq) + "\n"
                              + ">Upstream_seq" + "\n" + str(masked_upstream))


def main():

    args = argparser()

    upstream_and_downstream_seq(args)

    complementary_pairing(args)


if __name__ == '__main__':
    main()