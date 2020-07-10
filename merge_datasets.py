import pandas as pd
import argparse

def argparser():
    parser = argparse.ArgumentParser(description="Ask for path")
    parser.add_argument("-DCC", "--DCC", help="Path to DCC input file")
    parser.add_argument("-CE2", "--CE2", help="Path to CIRCexplorer2 input file")
    parser.add_argument("-CIRI2", "--CIRI2", help="Path to CIRI2 input file")
    parser.add_argument("-o", "--outdirectory", help="output directory")
    args = parser.parse_args()
    return args

def merge_datasets1(args):
    x = pd.read_csv(args.DCC, sep='\t', header=None)
    y = pd.read_csv(args.CE2, sep='\t', header=None)

    # Adding 1 to the circRNA start column of CIRCexplorer2 to compensate for different coordinate systems
    y[1] += 1

    # Making the circRNA_ID column
    y[0] = y.apply(lambda x: '%s:%s' % (x[0], x[1]), axis=1)
    y[0] = y.apply(lambda x: '%s-%s' % (x[0], x[2]), axis=1)

    # Merging The dataframe
    merged_df = pd.merge(x, y, on=0)

    merged_df['4_x'] = pd.to_numeric(merged_df['4_x'])

    merged_df[12] = pd.to_numeric(merged_df[12])

    merged_df['Avg_Expression'] = merged_df[[12, '4_y']].mean(axis=1)

    print(merged_df.shape)

    merged_df.to_csv(args.outdirectory + ("merged_" + args.DCC.split('/')[-1]), sep='\t', index=False)

def merge_datasets2(args):
        x = pd.read_csv(args.CIRI2, sep='\t', header=None)
        y = pd.read_csv(args.CE2, sep='\t', header=None)

        # Making the circRNA_ID column
        y[2] = pd.to_numeric(y[2][1:-1]) - 1
        y[2] = y[2].dropna().astype(int).astype(str)

        # Making the circRNA_ID column
        y[0] = y.apply(lambda x: '%s:%s' % (x[1], x[2]), axis=1)
        y[0] = y.apply(lambda x: '%s-%s' % (x[0], x[3]), axis=1)
        x[0] = x.apply(lambda x: '%s:%s' % (x[0], x[1]), axis=1)
        x[0] = x.apply(lambda x: '%s-%s' % (x[0], x[2]), axis=1)

        # Merging The dataframe
        merged_df = pd.merge(x, y, on=0)

        merged_df['Avg_Expression'] = merged_df[[12, '4_y']].mean(axis=1)

        print(merged_df.shape)

        merged_df.to_csv(args.outdirectory + ("merged_" + args.CIRI2.split('/')[-1].split('_')[-1]), sep='\t', index=False)

def main():

    ## args define
    args = argparser()

    ## Reading in lists
    if args.DCC is not None:
        merge_datasets1(args)

    if args.CIRI2 is not None:
        merge_datasets2(args)


if __name__ == '__main__':
    main()



