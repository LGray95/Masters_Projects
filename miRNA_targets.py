import argparse
import glob
import os


def argparser():
    parser = argparse.ArgumentParser(description='Ask for path')
    parser.add_argument("-d", "--directory", help = "Path to File Directory")
    parser.add_argument("-i", "--infile", help = "Path to input gene list")
    parser.add_argument("-o", "--outdirectory", help = "output directory")
    args = parser.parse_args()
    return args


# Reading in gene list file
def read_genelist(args):
    gene_list = []
    targets_file = open(args.infile).read().splitlines()
    for line in targets_file:
        col = line.split(',')
        gene = col[4].strip('"')
        gene_list.append(gene)
    return gene_list


def miRNA_targets(args, gene_list):
    # Set empty results list
    result_list = []
    # Iterate over each input file
    for filename in glob.iglob(args.directory+"/*/*.gtf"):
        print(filename)
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

        # Write out filename
        print(args.outdirectory + filename.split("/")[-2])
        if os.path.isdir(args.outdirectory + filename.split("/")[-2]) == False:
            os.mkdir(args.outdirectory + filename.split("/")[-2])
        outfilename = filename.split("/")[-2] + "/" + filename.split("/")[-1].replace(".gtf", ".txt")
        with open(args.outdirectory+outfilename, "w") as out:
            for element in result_list:
                out.write(element + "\n")


def main():

    ##args define
    args = argparser()

    ##make genelist
    genelist = read_genelist(args)

    ##Run miRNA_targets
    miRNA_targets(args, genelist)


if __name__ == '__main__':
    main()



