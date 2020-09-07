import argparse
import pysam
from split_coordinates import split_coords
from remove_duplicates import remove_duplicates
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from Bio.SeqUtils import MeltingTemp as mt

# Set up command line interface
def argparser():
    parser = argparse.ArgumentParser(description='Ask for files')
    parser.add_argument("-G", help = "Path to Reference Genome")
    parser.add_argument("-A", help = "Path to Annotation file")
    parser.add_argument("-ID", help = "circRNA coordinates chr:start-end")
    parser.add_argument("-s", help = "circRNA strand")
    #parser.add_argument("-o", help = "path to output file")
    args = parser.parse_args()
    
    return args


# indexing gtf for first and last exons
def index_gtf(args):
    line_number1 = []
    line_number2 = []
    gtf_file = open(args.A).read().splitlines()
    for line in gtf_file:
        col = line.split("\t")
        if col[0] == split_coords(args.ID)[0] and col[2] == "exon" and col[3] == str(split_coords(args.ID)[1]):
            line_number1.append(gtf_file.index(line))
        elif col[0] == split_coords(args.ID)[0] and col[2] == "exon" and col[4] == str(split_coords(args.ID)[2]):
            line_number2.append(gtf_file.index(line))
    
    return line_number1, line_number2, gtf_file


# Making list of exons composing circRNA
def exon_lists(line_number1, line_number2, gtf_file):
    exon_list = []
    for line in gtf_file[line_number1[0]:line_number2[0]+1]:
        col = line.split("\t")
        if col[2] == "exon":
            exon_list.append(col[0] + ":" + col[3] + "-" + col[4])
    exon_list = remove_duplicates(exon_list)
    
    return exon_list

# Pulling sequence from reference and joining together
# Check user input for forward or reverse strand
def exonic_circRNA(exon_list, args):
    circRNA_seq = ""
    for exon in exon_list:
        exon_seq = Seq(pysam.faidx(args.G, exon), generic_dna)
        if args.s == "-":
            exon_seq = str(exon_seq[(len(exon_seq.split('\n')[0])):-1].upper().reverse_complement())
            circRNA_seq += exon_seq
        elif args.s == "+":
            exon_seq = str(exon_seq[(len(exon_seq.split('\n')[0])):-1].upper())
            circRNA_seq += exon_seq
    # Making a one line sequence
    circRNA_seq = circRNA_seq.replace("\n", "")
    
    return circRNA_seq

# Selecting the primers based on size
def select_primers(circRNA_seq):
    for i in range(7,14):
        print(i)
        five_end = circRNA_seq[0:i]
        three_end = circRNA_seq[(i-20):]
        forward_primer = Seq(circRNA_seq[(i-150):(i-150)+20])

        Rev_Comp = Seq(three_end + five_end)
        Rev_Primer = Rev_Comp.reverse_complement()

        print(forward_primer, mt.Tm_Wallace(forward_primer))
        print(Rev_Primer, mt.Tm_Wallace(Rev_Primer))
        print(Rev_Comp)
        difference = mt.Tm_Wallace(forward_primer) - mt.Tm_Wallace(Rev_Primer)
        print(abs(difference))
        
 
def main():
    # args define
    args = argparser()
    
    print("Designing Primers for circRNA: " + args.ID + "" + args.s)
    
    # gather information
    line_number1 = index_gtf(args)[0]
    line_number2 = index_gtf(args)[1]
    gtf_file = index_gtf(args)[2]
    
    # generate exon list
    exon_list = exon_list(line_number1, line_number2, gtf_file)
    
    # generate exonic circRNA
    circRNA_seq = exonic_circRNA(exon_list, args)
    
    select_primers(circRNA_seq)

if __name__ == '__main__':
    main()
    
    
    

    
