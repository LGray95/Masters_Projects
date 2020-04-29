# This program was developed by Liam Cheney in 2019-2020
import argparse

def argparse():
    parser = argparse.ArgumentParser(description='Ask for files')
    parser.add_argument("-i", "--infile", help = "input filename")
    parser.add_argument("-o", "--outfile", help = "output filename")
    args = parser.parse_args()

    return args

def read_infile(args):
    infile = open(args.infile,'r').read().splitlines()

    return infile


def gathering_info(infile):
    print('Gathering Information.')
    total_dict = {}
    count = 0
    for line in infile[1:]:
        if 'circRNA_ID' not in line:
            count = count + 1
            col = line.split('\t')
            circ_name = col[5]
            cpm = float(col[18].strip('\"'))

            if circ_name not in total_dict.keys():
                total_dict[circ_name] = {'cpm':[]}

            total_dict[circ_name]['cpm'].append(cpm)

    return total_dict


def calc_averages(total_dict):
    print('Calculating Averages.')
    average_dict = {}
    for key in total_dict.keys():
        cpm_avg = sum(total_dict[key]['cpm']) / len(total_dict[key]['cpm'])
        average_dict[key] = {'cpm':cpm_avg}

    return average_dict


def writing_out(args, average_dict):
    print('Writing out.')
    with open(args.o,'w') as out:
        out.write('circRNA_ID' + '\t' + 'Avg tpm' + '\n')
        for key in average_dict.keys():
            cpm = str(average_dict[key]['cpm'])

            out.write(key + '\t' + cpm + '\n')


def main():

    # set arguments
    args = argparse()

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
