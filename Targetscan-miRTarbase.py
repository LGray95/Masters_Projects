import argparse

parser = argparse.ArgumentParser(description="Ask for path")
parser.add_argument("-i", help="infile")
args = parser.parse_args()


tscan = open("/Users/lachlan/OneDrive - UNSW/zeeshan/miRNA_binding/miRNA_interactions/Predicted_Targets_Info.default_predictions.txt", "r").read().splitlines()

miRTarbase = open("/Users/lachlan/OneDrive - UNSW/zeeshan/miRNA_binding/miRNA_interactions/hsa_MTI.txt", "r", encoding='latin1').read().splitlines()

infile = open("/Users/lachlan/OneDrive - UNSW/zeeshan/miRNA_binding/miRNA_interactions/"+args.i, "r").read().splitlines()

tscan_list = []
for line in tscan:
    col = line.split("\t")
    miRNA = col[0]
    gene = col[2]

    tscan_list.append(miRNA + "\t" + gene)

mir_list = []
for line in miRTarbase:
    col = line.split("\t")
    miRNA = col[1].replace("hsa-", "")
    gene = col[3]

    mir_list.append(miRNA + "\t" + gene)

outlist = []
for miRNA in infile:
    for element in tscan_list:
        if miRNA.replace("hsa-", "") in element:
            outlist.append(element)

for miRNA in infile:
    for element in mir_list:
        if miRNA.replace("hsa-", "") in element:
            outlist.append(element)

for i in set(outlist):
    print(i)

