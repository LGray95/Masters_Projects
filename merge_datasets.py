import pandas
import argparse

def argparser():
    parser = argparse.ArgumentParser(description="Ask for path")
    parser.add_argument("-DCC_directory", "-DCC_directory", help="Path to DCC directory")
    parser.add_argument("-DCC_list", "--DCC_list", help="Path to DCC file")
    parser.add_argument("-CE2_list", "--CE2_list", help="Path to CIRCexplorer2 list")
    parser.add_argument("-CIRI2", "--CIRI2_list", help="Path to CIRI2 input file")
    parser.add_argument("-o", "--outdirectory", help="output file")
    args = parser.parse_args()
    return args


def DCClist(args):
    if args.DCC_list is not None:
        DCC_list = open(args.DCC_list, "r").read().splitlines()
        return DCC_list


def CE2list(args):
    if args.CE2_list is not None:
        CE2_list = open(args.CE2_list, "r").read().splitlines()
        return CE2_list

def CIRI2list(args):
    if args.CIRI2_list is not None:
        CIRI2_list = open(args.CIRI2_list, "r").read().splitlines()
        return CIRI2_list


def merge_datasets(args, DCC_list, CE2_list, CIRI2_list):
    if args.DCC_list is not None:
        for i, j in zip(DCC_list, CE2_list):
            print(i, j)
            x = pandas.read_csv(args.DCC_directory + i, sep='\t', header=None)
            y = pandas.read_csv(j, sep='\t', header=None)

            # Adding 1 to the circRNA start column of CIRCexplorer2 to compensate for different coordinate systems
            y[1] += 1

            # Making the circRNA_ID column
            y[0] = y.apply(lambda x: '%s:%s' % (x[0], x[1]), axis=1)
            y[0] = y.apply(lambda x: '%s-%s' % (x[0], x[2]), axis=1)

            # Merging The dataframe
            merged_df = pandas.merge(x, y, on=0)

            merged_df['1_x'] = pandas.to_numeric(merged_df['1_x'])

            merged_df[12] = pandas.to_numeric(merged_df[12])

            merged_df[18] = ((merged_df['1_x'] + merged_df[12]) / 2)

            print(merged_df[merged_df.duplicated(subset=0)])

            merged_df.to_csv(args.outdirectory + i, sep='\t', index=False)
    if args.CIRI2_list is not None:
        for i, j in zip(CIRI2_list, CE2_list):
            print(i, j)
            x = pandas.read_csv(args.DCC_directory + i, sep='\t', header=None)
            y = pandas.read_csv(j, sep='\t', header=None)

            # Adding 1 to the circRNA start column of CIRCexplorer2 to compensate for different coordinate systems
            y[1] += 1

            # Making the circRNA_ID column
            y[0] = y[0].replace("|", "-")

            # Merging The dataframe
            merged_df = pandas.merge(x, y, on=0)

            merged_df['1_x'] = pandas.to_numeric(merged_df['1_x'])

            merged_df[12] = pandas.to_numeric(merged_df[12])

            merged_df[18] = ((merged_df['1_x'] + merged_df[12]) / 2)

            print(merged_df[merged_df.duplicated(subset=0)])

            merged_df.to_csv(args.outdirectory + i, sep='\t', index=False)

def main():

    ## args define
    args = argparser()

    ## Reading in lists
    DCC_list = DCClist(args)
    CE2_list = CE2list(args)
    CIRI2_list = CIRI2list(args)

    ## Merge files
    merge_datasets(args, DCC_list, CE2_list, CIRI2_list)


if __name__ == '__main__':
    main()
