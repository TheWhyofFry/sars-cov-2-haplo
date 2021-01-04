import os
import argparse 

import pandas as pd


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", dest="filelist", nargs="+", help="Input list of csvs")
    parser.add_argument("-o", dest="outfile", help="Output filename")


    args = parser.parse_args()


    csv_list = []

    for filename in args.filelist:
        try:
            csv_list.append(pd.read_csv(filename))
        except pd.errors.EmptyDataError:
            print("File %s is empty"%filename)
            pass
    #print("CSVlist len:",len(csv_list))

    try:
        p = pd.concat(csv_list,axis=0)
        p.to_csv(args.outfile)
    else:
        with open(args.outfile,"w") as empty_file:
            empty_file.write("")



