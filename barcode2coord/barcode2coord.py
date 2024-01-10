#!/usr/bin/env python

import pandas as pd
import gzip
import argparse

def main(args):
    assert args.outfname[-3:] == ".gz", "output must end in .gz"
    assert args.infname[-3:] == ".gz", "input must end in .gz"
    spatial_df = pd.read_csv(args.spatial_barcodes, sep='\t', header=None, names=['barcode', 'x', 'y'])
    barcode_mapping = dict(zip(spatial_df['barcode'], spatial_df.apply(lambda row: str(row['x']) + "x" + str(row['y']), axis=1)))

    with gzip.open(args.infname, 'rt') as inFP, gzip.open(args.outfname, 'wt') as outFP:
        for line in inFP:
            coord = barcode_mapping.get(line.strip(),"NaN")
            outFP.write(coord+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace barcodes in a gzipped TSV file using a mapping from another TSV file.")
    parser.add_argument("infname", help="Input gzipped TSV file containing barcodes")
    parser.add_argument("spatial_barcodes", help="Spatial barcodes TSV file with columns: barcode, column2, column3")
    parser.add_argument("outfname", help="Output file to save the result")
    
    args = parser.parse_args()
    main(args)

