import argparse


def argparser():
    parser = argparse.ArgumentParser(description='Ask for path')
    parser.add_argument("-i", "--infile", help="Path of input DCC file")
    parser.add_argument("-o", "--outdir", help="Path of output directory")
    args = parser.parse_args()

    return args


def read_file(args):
    infile = open(args.infile, "r").read().splitlines()

    return infile


def count_columns(args):
    with open(args.infile) as infile:
        line = infile.readline()
    total_columns = int(len(line.split()))

    return total_columns


def split_file(args, total_columns, infile):
    out_list = []
    for k in range(3, total_columns):
        for line in infile:
            col = line.split('\t')
            out_list.append(col[0]+":"+col[1]+"-"+col[2]+'\t'+col[k])

        sample_id = (out_list[0].split('\t')[1] + "_" + str(k-2))

        outfile = open(args.outdir + sample_id + '.txt', 'w')
        for element in out_list:
            outfile.write(str(element) + '\n')
        outfile.close()


def main():

    args = argparser()

    infile = read_file(args)

    total_columns = count_columns(args)

    split_file(args, total_columns, infile)


if __name__ == '__main__':
    main()