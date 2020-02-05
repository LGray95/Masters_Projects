import argparse
from split_coordinates import split_coords

parser = argparse.ArgumentParser(description='Ask for files')
parser.add_argument("-A", help="Path to miRNA Annotation file")
parser.add_argument("-i", help="Path to input circRNA list")
# parser.add_argument("-o", help="Path to output file")
args = parser.parse_args()


#reading in both .gff3 and filtered circRNA files
gff3 = open(args.A).read().splitlines()
circRNA_file = open(args.i).read().splitlines()


#creating empty variables for the lists
microRNA_list = []
circRNA_list = []
#looping through .gff3 file
for line in gff3[1:-1]:
    col = line.split('\t')
    microRNA_chr = col[0]
    microRNA_start = int(col[3])
    microRNA_end = int(col[4])
    microRNA_type = col[2]
    microRNA_ID = col[8].split("Name=")[1].split(";")[0]


    #creating a list with chr, start and end. Adding each list to the list variable
    list_n = [microRNA_chr, microRNA_start, microRNA_end, microRNA_type, microRNA_ID]
    microRNA_list.append(list_n)

for line in circRNA_file:
    circRNA_chr = split_coords(line)[0]
    circRNA_start = split_coords(line)[1]
    circRNA_end = split_coords(line)[2]


    #Add in circBase_tissue and Difference
    list_m = [circRNA_chr, circRNA_start, circRNA_end]

    circRNA_list.append(list_m)

yes_count = 0
out_list = []
for i in microRNA_list:
    for j in circRNA_list:
        if i[0] == j[0] and i[1] >= j[1] and i[2] <= j[2] and i[3] == "miRNA_primary_transcript":
            yes_count += 1
            circRNA = j[0]+":"+str(j[1])+"-"+str(j[2])
            miRNA = i[0]+":"+str(i[1])+"-"+str(i[2])
            list_z = [circRNA, miRNA, i[4]]
            out_list.append(list_z)

print(yes_count)
print(out_list)

# # Set column headings for outfile
# outfile = open(args.o, 'w')
# outfile.write('microRNA_coordinates' + ',' + 'microRNA_type' + ',' +
#               'microRNA_ID' + ',' + 'circRNA_ID' + ',' + 'strand' + ',' + 'circRNA_gene' + ',' + 'CPM' + ',' + 'Exon_Count' + '\n')
# for j in out_list:
#     for k in j:
#         outfile.write(str(k) + ',')
#     outfile.write('\n')
# outfile.close()



