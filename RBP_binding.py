import argparse
from split_coordinates import split_coords
from remove_duplicates import remove_duplicates


# Set arguments
def argparser():
    parser = argparse.ArgumentParser(description="Ask for paths")
    parser.add_argument("-c", "--coords", help="circRNA coordinates")
    parser.add_argument("-i", "--infile", help="ENCORI infile")
    parser.add_argument("-A", "--annotation", help="path to intron file")
    parser.add_argument("-E", "--exon", help="Path to exon file")
    parser.add_argument("-o", "--outfile", help="Out Directory")
    args = parser.parse_args()

    return args


# Set coordinates of the downstream intron
def down_intron(args):
    annotation = open(args.annotation, "r").read().splitlines()
    for line in annotation:
        col = line.split("\t")
        if int(col[2]) == split_coords(args.coords)[1] - 1:
            downstream_intron = (int(col[1]))

            return downstream_intron


# Set coordinates of upstream intron
def up_intron(args):
    annotation = open(args.annotation, "r").read().splitlines()
    for line in annotation:
        col = line.split("\t")
        if int(col[1]) == split_coords(args.coords)[2]:
            upstream_intron = (int(col[2]))

            return upstream_intron


# indexing gtf for first and last exons
def exon_coords(args):
    line_number1 = []
    line_number2 = []
    exon_file = open(args.exon).read().splitlines()
    for line in exon_file:
        col = line.split("\t")
        if col[0] == split_coords(args.coords)[0] and col[2] == "exon" and col[3] == str(split_coords(args.coords)[1]):
            line_number1.append(exon_file.index(line))
        if col[0] == split_coords(args.coords)[0] and col[2] == "exon" and col[4] == str(split_coords(args.coords)[2]):
            line_number2.append(exon_file.index(line))


    # Making list of exons composing circRNA
    exon_list = []
    for line in exon_file[line_number1[0]:line_number2[0]+1]:
        col = line.split("\t")
        if col[2] == "exon":
            exon_list.append(col[0] + ":" + col[3] + "-" + col[4])

    exon_list = remove_duplicates(exon_list)

    return exon_list



# Search within ENCORI output to find binding sites of RBP
def ENCORI(args, downstream_intron, upstream_intron):
    outlist = []
    infile = open(args.infile, "r").read().splitlines()
    for line in infile:
        if "#" not in line:
            col = line.split("\t")
            if col[8] == split_coords(args.coords)[0] \
                    and int(col[9]) >= downstream_intron + 1 \
                    and int(col[10]) <= upstream_intron:

                outlist.append(line)

    return outlist


# Determine where proteins are binding to; Downstream intron, within circRNA, BSJ, exon, or Upstream intron
def location_of_binding(args, outlist, exon_list):
    finaloutlist = []
    for line in outlist:
        col = line.split("\t")
        if int(col[9]) < split_coords(args.coords)[1]:
            line = line + "\t" + "Downstream"
        elif int(col[9]) > split_coords(args.coords)[2]:
            line = line + "\t" + "Upstream"
        elif int(col[9]) > split_coords(args.coords)[1] and int(col[10]) < split_coords(args.coords)[2]:
            line = line + "\t" + "circRNA"
        elif int(col[9]) == split_coords(args.coords)[1] or int(col[10]) == split_coords(args.coords)[2]:
            line = line + "\t" + "BSJ"
        elif int(col[10]) > split_coords(args.coords)[1] or int(col[9]) < split_coords(args.coords)[2]:
            line = line + "\t" + "BSJ"
        finaloutlist.append(line)
        for exon in exon_list:
            if int(col[9]) >= split_coords(exon)[1] and int(col[10]) <= split_coords(exon)[2]:
                line = line + "\t" + "Exon Binding"
            finaloutlist.append(line)

    return finaloutlist


# Count how many are binding within introns
def intron_binding(args, finaloutlist):
    yes_count = 0
    for element in finaloutlist:
        col = element.split("\t")
        if int(col[9]) < split_coords(args.coords)[1] \
                or int(col[9]) > split_coords(args.coords)[2]:
            yes_count += 1

    return yes_count


# Write output to new file
def write_out(args, finaloutlist, yes_count):
    outfile = open(args.outfile, "w")
    for element in finaloutlist:
        outfile.write(element + "\n")
    outfile.write("Number of intronic Binding sites" + "\t" + str(yes_count) + "\n")
    outfile.write("Number of circRNA Binding sites" + "\t" + str(len(finaloutlist)-yes_count) + "\n")


# Set main arguments
def main():

    args = argparser()

    downstream_intron = down_intron(args)

    upstream_intron = up_intron(args)

    exon_list = exon_coords(args)

    outlist = ENCORI(args, downstream_intron, upstream_intron)

    finaloutlist = location_of_binding(args, outlist, exon_list)

    yes_count = intron_binding(args, finaloutlist)

    write_out(args, finaloutlist, yes_count)

if __name__ == '__main__':
    main()




