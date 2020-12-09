import argparse
import pandas as pd


def argparser():
    parser = argparse.ArgumentParser(description='Ask for path')
    parser.add_argument("-c", "--coords", help = "circRNA coordinates")
    parser.add_argument("-i", "--infile", help = "Path to input file")
    parser.add_argument("-o", "--outfile", help = "output file")
    args = parser.parse_args()
    return args


def split_coords(coords):
    chr = coords.split(":")[0]
    start = int(coords.split(":")[1].split("-")[0])
    end = int(coords.split(":")[1].split("-")[1])
    return chr, start, end

def ENCORI(args):
    df = pd.read_csv(args.infile, sep="\t", skiprows=3, header=0)
    df["start"] = pd.to_numeric(df["start"])
    df["end"] = pd.to_numeric(df["end"])
    filtered_df = df[(df["chromosome"] == split_coords(args.coords)[0]) & (df["start"] >= split_coords(args.coords)[1]) & (df["end"] <= split_coords(args.coords)[2])]

    return filtered_df


def writing_out(filtered_df, args):
    filtered_df.to_csv(args.outfile, sep="\t", index=False)


def main():

    # args define
    args = argparser()

    # Run ENCORI
    filtered_df = ENCORI(args)

    # Write out file
    writing_out(filtered_df, args)


if __name__ == '__main__':
    main()