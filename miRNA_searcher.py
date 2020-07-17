import argparse
from split_coordinates import split_coords


def argparser():
    parser = argparse.ArgumentParser(description='Ask for files')
    parser.add_argument("-A", help="Path to miRNA Annotation file")
    parser.add_argument("-i", help="Path to input circRNA file")
    parser.add_argument("-o", help="Path to output directory")
    args = parser.parse_args()

    return args


def read_gff3(args):
    gff3 = open(args.A).read().splitlines()

    return gff3


def read_circRNA_file(args):
    circRNA_file = open(args.i).read().splitlines()

    return circRNA_file


def miRNA_data(gff3):
    microRNA_list = []
    for line in gff3[1:-1]:
        col = line.split('\t')
        microRNA_chr = col[0]
        microRNA_start = int(col[3])
        microRNA_end = int(col[4])
        microRNA_type = col[2]
        microRNA_ID = col[8].split("Name=")[1].split(";")[0]


        list_n = [microRNA_chr, microRNA_start, microRNA_end, microRNA_type, microRNA_ID]
        microRNA_list.append(list_n)

    return microRNA_list


def circRNA_data(circRNA_file):
    circRNA_list = []
    for line in circRNA_file[1:-1]:
        col = line.split('\t')
        circRNA_chr = split_coords(col[0])[0]
        circRNA_start = split_coords(col[0])[1]
        circRNA_end = split_coords(col[0])[2]
        gene = col[19]

        # Add to list
        list_m = [circRNA_chr, circRNA_start, circRNA_end, gene]

        circRNA_list.append(list_m)

    return circRNA_list


def match(microRNA_list, circRNA_list):
    out_list = []
    for i in microRNA_list:
        for j in circRNA_list:
            if i[0] == j[0] and i[1] >= j[1] and i[2] <= j[2] and i[3] == "miRNA_primary_transcript":
                circRNA = j[0]+":"+str(j[1])+"-"+str(j[2])
                miRNA = i[0]+":"+str(i[1])+"-"+str(i[2])
                list_z = [circRNA, j[3], miRNA, i[4]]
                out_list.append(list_z)

    return out_list


def writing_out(args, out_list):
    outfile = open(args.o, 'w')
    outfile.write('circRNA ID' + "\t" + 'Gene' + '\t' + 'miRNA ID' + '\t' + 'miRNA name' + '\n')
    for i in out_list:
        outfile.write(i[0] + "\t" + i[1] + '\t' + i[2] + '\t' + i[3] + '\n')


def main():

    args = argparser()

    gff3 = read_gff3(args)

    circRNA_file = read_circRNA_file(args)

    microRNA_list = miRNA_data(gff3)

    circRNA_list = circRNA_data(circRNA_file)

    out_list = match(microRNA_list, circRNA_list)

    writing_out(args, out_list)


if __name__ == '__main__':
    main()


