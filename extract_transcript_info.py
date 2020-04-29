import argparse
import glob

parser = argparse.ArgumentParser(description='Ask for path')
parser.add_argument("-d", help = "Path to File Directory")
parser.add_argument("-o", help = "output directory")
args = parser.parse_args()

for filename in glob.iglob(args.d+"/*/*.gtf"):
    print(filename)
    # Set empty results list
    result_list = []
    # Open the file and split lines
    infile = open(filename,'r').read().splitlines()
    for line in infile:
        if "#" not in line[0:2]:
            col = line.split("\t")
            # Select only transcripts
            if col[2] == 'transcript':
                transcript_id = col[8].split(";")[1].split()[1].replace('"', "")
                coverage = col[8].split(";")[5].split()[1].replace('"', "")

                result_list.append(transcript_id + "\t" + coverage)

    outfilename = filename.split("/")[-2] + "/" + filename.split("/")[-1].replace(".gtf", ".txt")
    with open(args.o+outfilename, "w") as out:
        for element in result_list:
            out.write(element + "\n")

