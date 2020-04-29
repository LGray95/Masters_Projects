from split_coordinates import split_coords
import argparse


def argparser():
    parser = argparse.ArgumentParser(description='Ask for path')
    parser.add_argument("-c", "--coords", help = "circRNA coordinates")
    parser.add_argument("-i", "--infile", help = "Path to circBase file")
    args = parser.parse_args()
    return args


def read_files(args):
    circBase = open(args.infile, "r").read().splitlines()
    return circBase


def check_circBase(args, circBase):
    match = []
    for line in circBase:
        col = line.split("\t")
        if "#" not in line:
            if split_coords(args.coords)[0] == col[0] and split_coords(args.coords)[1] == int(col[1]) and split_coords(args.coords)[2] == int(col[2]):
                match.append(line)
    return match


def main():

    # Set command line arguments
    args = argparser()

    # Read in circBase file
    circBase = read_files(args)

    # Perform match algorithm
    match = check_circBase(args, circBase)

    # Print all matches in file
    for element in match:
        print(element)


if __name__ == '__main__':
    main()
