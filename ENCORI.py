import argparse

def argparser():
    parser = argparse.ArgumentParser(description='Ask for path')
    parser.add_argument("-c", "--coords", help = "circRNA coordinates")
    parser.add_argument("-i","--infile", help = "Path to input file")
    parser.add_argument("-o", "--outfile", help = "output file")
    args = parser.parse_args()
    return args

def split_coords(coords):
    chr = coords.split(":")[0]
    start = int(coords.split(":")[1].split("-")[0])
    end = int(coords.split(":")[1].split("-")[1])
    return chr, start, end

def ENCORI(args):
    outlist = []
    infile = open(args.infile, "r").read().splitlines()
    for line in infile:
        if "#" not in line:
            col = line.split("\t")
            if col[5] == split_coords(args.coords)[0] \
                    and int(col[6]) >= split_coords(args.coords)[1] \
                    and int(col[7]) <= split_coords(args.coords)[2]:
                outlist.append(line)

    outfile = open(args.outfile, "w")
    for element in outlist:
        outfile.write(element + "\n")


def main():

    ##args define
    args = argparser()

    ##Run ENCORI
    ENCORI(args)


if __name__ == '__main__':
    main()