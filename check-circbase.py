from split_coordinates import split_coords

circfile =  open("/Users/lachlan/Downloads/infile.txt", "r").read().splitlines()
circBase =  open('/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/circBase/hsa_hg19_circRNA.txt', "r")

circ_list = []
circBase_list = []
for entry in circfile:
    circRNA_coords = entry
    circ_chr = split_coords(circRNA_coords)[0]
    circ_start = int(split_coords(circRNA_coords)[1])
    circ_end = int(split_coords(circRNA_coords)[2])

    list_a = circ_chr, circ_start, circ_end
    circ_list.append(list_a)

for line in circBase:
    # print(line)
    col = line.split("\t")
    if "#" not in line:
        chr = col[0]
        start = int(col[1])
        end = int(col[2])
        list_b = chr, start, end
        circBase_list.append(list_b)

for i in circ_list:
    for j in circBase_list:
        if i[0] == j[0] and (i[1] - j[1]) != 0:
            print(line)
