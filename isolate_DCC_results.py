import sys

# Where:
# Infile is the DCC output file
# Outdir is the directory where you want to store the results



# Counts the total columns of the file to parse to range()

with open(sys.argv[1]) as file:
    line = file.readline()
total_columns = int(len(line.split()))

# Open file
in_DCC = open(sys.argv[1]).read().splitlines()
# Read one columns at a time
for k in range(3, total_columns):
    out_list = []
    for i in in_DCC:
        col = i.split('\t')
        chr = col[0]
        start = col[1]
        end = col[2]
        sample = col[k]
        out_list.append(chr+":"+start+"-"+end+'\t'+sample)


    filtered_out_list = []
    sample_ID =(out_list[0].split('\t')[1])
    print(sample_ID)
    # Filter each circRNA for at least two BSJ reads
    for i in out_list[1:-1]:
        if int(i.split('\t')[1]) >= 2:
            filtered_out_list.append(i)

    # Write output to  file
    outfile = open(sys.argv[2]+sample_ID+'.txt', 'w')
    for i in filtered_out_list:
        outfile.write(str(i) + '\n')
    outfile.close()