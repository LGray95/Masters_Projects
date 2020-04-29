with open("/Users/lachlan/OneDrive - UNSW/Masters_Research_2019/Data/linear_analysis/NCAM1/New/MTLE_Hippocampus.txt") as infile:

    count = 0
    for line in infile:
        col = line.split("\t")
        if ".bam" in line:
            print(line.split("-o")[1].split("-A")[0])
            count += 1
        if "NCAM1" in line and col[2] == "transcript" and "NCAM1-AS1" not in line:
            # coords = col[0] + ":" + col[3] + "-" + col[4]
            ID = col[8].split("transcript_id")[1].split(";")[0].replace("\"", "")
            TPM = col[8].split("TPM")[1].replace("\"", "").replace(";", "")
            print(ID + "\t" + TPM)
    print(count)