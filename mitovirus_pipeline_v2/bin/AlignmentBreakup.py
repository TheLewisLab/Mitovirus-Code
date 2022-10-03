import pandas as pd
import sys
import getopt
import numpy as np

def fileIO(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile="])
    except getopt.GetoptError:
        print('AlignmentBreakup.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('AlignmentBreakup.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
    print('Input file is "', inputfile)
    return inputfile


def main(argv):
    inputFile = fileIO(argv)

    aligns = dict()

    df = pd.read_csv(inputFile, delimiter='\t', header=None, names=['Node', 'Ref', 'NA1', 'NA2', 'NA3', 'NA4', 'NA5', 'NA6', 'NA7', 'NA8', 'NA9', 'NA10'])
    group = df['Ref'].groupby(df['Node'])
    for g in group:
        print(g[0])
        temp = list()
        for s in g[1]:
            spl = s.split('.')
            temp.append(spl[0])
        for s in np.unique(temp):
            if s in aligns.keys():
                aligns[s].append(g[0])
            else:
                aligns[s] = list()
                aligns[s].append(g[0])

    for fam in aligns:
        with open(fam + '.txt', 'w') as f:
            for item in aligns[fam]:
                f.write("%s\n" % item)


if __name__ == "__main__":
    main(sys.argv[1:])