import argparse


def argparser():
    parser = argparse.ArgumentParser(description='Ask for path')
    parser.add_argument("-c", "--coords", help = "circRNA coordinates")
    parser.add_argument("-i", "--infile", help = "Path to input file")
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

    return outlist


def writing_out(outlist, args):
    outfile = open(args.outfile, "w")
    for element in outlist:
        outfile.write(element + "\n")


def main():

    ##args define
    args = argparser()

    ##Run ENCORI
    outlist = ENCORI(args)

    print("Counting occurance of each miRNA")
    total_dict = {}
    for i in outlist:
        col = i.split("\t")
        coord = col[6]
        if coord in total_dict:
            total_dict[coord] += 1
        else:
            total_dict[coord] = 1
    for i in (sorted(total_dict.items(), key=lambda x: x[1])):
        print(i)

    # writing_out(outlist, args)


if __name__ == '__main__':
    main()