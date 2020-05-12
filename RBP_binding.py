import argparse
from split_coordinates import split_coords

def argparser():
    parser = argparse.ArgumentParser(description="Ask for paths")
    parser.add_argument("-c", "--coords", help="circRNA coordinates")
    parser.add_argument("-i", "--infile", help="ENCORI infile")
    parser.add_argument("-A", "--annotation", help="path to genome file")
    parser.add_argument("-o", "--outfile", help="Out Directory")
    args = parser.parse_args()

    return args


def down_intron(args):
    annotation = open(args.annotation, "r").read().splitlines()
    for line in annotation:
        col = line.split("\t")
        if int(col[2]) == split_coords(args.coords)[1] - 1:
            downstream_intron = (int(col[1]))

            return downstream_intron

def up_intron(args):
    annotation = open(args.annotation, "r").read().splitlines()
    for line in annotation:
        col = line.split("\t")
        if int(col[1]) == split_coords(args.coords)[2]:
            upstream_intron = (int(col[2]))

            return upstream_intron





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


def intron_binding(args, outlist):
    yes_count = 0
    for element in outlist:
        col = element.split("\t")
        if int(col[9]) < split_coords(args.coords)[1] \
                or int(col[9]) > split_coords(args.coords)[2]:
            yes_count += 1

    return yes_count


def write_out(args, outlist, yes_count):
    outfile = open(args.outfile, "w")
    for element in outlist:
        outfile.write(element + "\n")
    outfile.write("Number of intronic Binding sites" + "\t" + str(yes_count) + "\n")
    outfile.write("Number of circRNA Binding sites" + "\t" + str(len(outlist)-yes_count) + "\n")



def main():

    args = argparser()

    downstream_intron = down_intron(args)

    upstream_intron = up_intron(args)

    outlist = ENCORI(args, downstream_intron, upstream_intron)

    yes_count = intron_binding(args, outlist)

    write_out(args, outlist, yes_count)

if __name__ == '__main__':
    main()




