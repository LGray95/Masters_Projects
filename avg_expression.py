# This program was developed by Liam Cheney in 2019-2020
import argparse


def argparser():
    parser = argparse.ArgumentParser(description='Ask for files')
    parser.add_argument("-i", "--infile", help = "input filename")
    parser.add_argument("-o", "--outfile", help = "output filename")
    args = parser.parse_args()

    return args


def read_infile(args):
    infile = open(args.infile, 'r').read().splitlines()

    return infile


def gathering_info(infile):
    print('Gathering Information.')
    total_dict = {}
    count = 0
    for line in infile[1:]:
        if 'Avg_Expression' not in line:
            count = count + 1
            col = line.split('\t')
            circ_name = col[0]
            exp = int(col[29])
            exons = col[16]
            strand = col[5]
            gene = col[14]
            enst = col[15]

            if circ_name not in total_dict.keys():
                total_dict[circ_name] = {'exp':[]}

            total_dict[circ_name]['exp'].append(exp)
            total_dict[circ_name]['exon'] = exons
            total_dict[circ_name]['strand'] = strand
            total_dict[circ_name]['gene'] = gene
            total_dict[circ_name]['enst'] = enst

    return total_dict


def calc_averages(total_dict):
    print('Calculating Averages.')
    average_dict = {}
    for key, value in total_dict.items():
        exp_avg = sum(total_dict[key]['exp']) / len(total_dict[key]['exp'])
        average_dict[key] = {'exp':exp_avg}
        average_dict[key]['strand'] = value['strand']
        average_dict[key]['exon'] = value['exon']
        average_dict[key]['gene'] = value['gene']
        average_dict[key]['enst'] = value['enst']

    return average_dict


def writing_out(args, average_dict):
    print('Writing out.')
    with open(args.outfile, 'w') as out:
        out.write('circRNA_ID' + '\t' + 'Strand' + '\t' + 'Exon Composition' + '\t' + 'Gene' + '\t' +
                  'Ensembl ID' + '\t' + 'Avg exp' + '\n')
        for key in average_dict.keys():
            exp = str(average_dict[key]['exp'])
            strand = average_dict[key]['strand']
            gene = average_dict[key]['gene']
            exon = average_dict[key]['exon']
            enst = average_dict[key]['enst']

            out.write(key + '\t' + strand + '\t' + exon + '\t' + gene + '\t' + enst + '\t' + exp + '\n')


def main():

    # set arguments
    args = argparser()

    # Read infile
    infile = read_infile(args)

    # Gather information
    total_dict = gathering_info(infile)

    # Calculate average CPM
    average_dict = calc_averages(total_dict)

    # Write results to file
    writing_out(args, average_dict)


if __name__ == '__main__':
    main()
