import argparse
from split_coordinates import split_coords


# Set arguments
def argparser():
    parser = argparse.ArgumentParser(description="Ask for paths")
    parser.add_argument("-c", "--coords", help="circRNA coordinates")
    parser.add_argument("-i", "--infile", help="ENCORI infile")
    parser.add_argument("-A", "--annotation", help="path to genome file")
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


# Determine where proteins are binding to; Downstream intron, within circRNA, BSJ or Upstream intron
def location_of_binding(args, outlist):
    finaloutlist = []
    for line in outlist:
        col = line.split("\t")
        if int(col[9]) < split_coords(args.coords)[1]:
            line = line + "\t" + "Downstream"
        elif int(col[9]) > split_coords(args.coords)[2]:
            line = line + "\t" + "Upstream"
        elif int(col[9]) > split_coords(args.coords)[1] and int(col[10]) < split_coords(args.coords)[2]:
            line = line + "\t" + "circRNA"
        else:
            line = line + "\t" + "BSJ"
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

    outlist = ENCORI(args, downstream_intron, upstream_intron)

    finaloutlist = location_of_binding(args, outlist)

    yes_count = intron_binding(args, finaloutlist)

    write_out(args, finaloutlist, yes_count)

if __name__ == '__main__':
    main()




