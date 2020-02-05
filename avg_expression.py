# This program was developed by Liam Cheney in 2019-2020
import argparse


parser = argparse.ArgumentParser(description='Ask for files')
parser.add_argument("-i", help = "input filename")
parser.add_argument("-o", help = "output filename")
args = parser.parse_args()

infile = open(args.i,'r').read().splitlines()

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


print('Calculating Averages.')
average_dict = {}
for key in total_dict.keys():

    cpm_avg = sum(total_dict[key]['cpm']) / len(total_dict[key]['cpm'])

    average_dict[key] = {'cpm':cpm_avg}

print('Writing out.')
with open(args.o,'w') as out:
    out.write('circRNA_ID' + '\t' + 'Avg tpm' + '\n')
    for key in average_dict.keys():
        cpm = str(average_dict[key]['cpm'])

        out.write(key + '\t' + cpm + '\n')