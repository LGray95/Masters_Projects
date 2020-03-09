import argparse
import glob

parser = argparse.ArgumentParser(description='Ask for path')
parser.add_argument("-d", help = "Path to File Directory")
parser.add_argument("-i", help = "Path to input gene list")
parser.add_argument("-o", help = "output directory")
args = parser.parse_args()

# Reading in gene list file
targets_file = open(args.i).read().splitlines()
gene_list = []
for line in targets_file:
    col = line.split(',')
    gene = col[4].strip('"')
    gene_list.append(gene)


# Iterate over each input file
for filename in glob.iglob(args.d+"/*/*.gtf"):
    print(filename)
    # Set empty results list
    result_list = []
    # Open the file and split lines
    infile = open(filename,'r').read().splitlines()
    for line in infile:
        # Make sure we are reading in correct lines
        if "#" not in line[0:2]:
            col = line.split("\t")
            # Select only transcripts
            if col[2] == 'transcript':
                # Matching genes
                transcript_gene = col[8].split(";")[0].split()[1].replace('"', "")
                for gene in gene_list:
                    if transcript_gene == gene:
                    # Selecting transcript ID and coverage values
                        transcript_id = col[8].split(";")[1].split()[1].replace('"', "")
                        coverage = col[8].split(";")[3].split()[1].replace('"', "")
                        # Adding values to a list
                        result_list.append(transcript_id + "\t" + coverage)

    #Write out filename
    outfilename = filename.split("/")[-2] + "/" + filename.split("/")[-1].replace(".gtf", ".txt")
    with open(args.o+outfilename, "w") as out:
        for element in result_list:
            out.write(element + "\n")



