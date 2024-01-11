import h5py
import numpy as np
import sys
import time
from constants import *
from utils import create_datapoint
import pandas as pd

start_time = time.time()

assert sys.argv[1] in ['train', 'test', 'all']
assert sys.argv[2] in ['0', '1', 'all']

h5f = h5py.File(data_dir + 'datafile'
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '_with_gex.h5', 'r')


CHROM = h5f['CHROM'][:]
BARCODE = h5f['BARCODE'][:]
CELLTYPE = h5f['CELLTYPE'][:]
SEQ = h5f['SEQ'][:]
SKIPPED_COUNT = h5f['SKIPPED_COUNT'][:]
INCLUDED_COUNT = h5f['INCLUDED_COUNT'][:]

h5f.close()

# Load gene expression also.
# ccle_gex_file = pd.read_csv(ccle_gex_file)
# Make the column "StrippedName" the index and drop the column.
# ccle_gex_file = ccle_gex_file.set_index('StrippedName')
# print(ccle_gex_file.head())

h5f2 = h5py.File(data_dir + 'dataset'
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '_with_gex.h5', 'w')

CHUNK_SIZE = 1024
for i in range(SEQ.shape[0]//CHUNK_SIZE):
    print(i)
    # Each dataset has CHUNK_SIZE sequences. 
    if (i+1) == SEQ.shape[0]//CHUNK_SIZE:
        NEW_CHUNK_SIZE = CHUNK_SIZE + SEQ.shape[0]%CHUNK_SIZE
    else:
        NEW_CHUNK_SIZE = CHUNK_SIZE

    X_batch = []
    Y_batch = [[] for t in range(1)]
    Z_batch = [] # This is going to be celltype.
    
    for j in range(NEW_CHUNK_SIZE):
        idx = i*CHUNK_SIZE + j
        # celltype = CELLTYPE[idx].decode("utf-8")
        # Strip the double quotes.
        # celltype = celltype.strip('"')
        celltype = CELLTYPE[idx].decode("utf-8")
        # print(celltype)
        # Get gene expression from the df above
        # gex = np.array(ccle_gex_file.loc[celltype].values)
        # print(gex.shape)
        # print(celltype)
        # If the counts are negative, throw error.
        if INCLUDED_COUNT[idx] < 0 or SKIPPED_COUNT[idx] < 0:
            print('Negative count error')
            print(INCLUDED_COUNT[idx], SKIPPED_COUNT[idx])
            print(SEQ[idx], CELLTYPE[idx])
            print(idx)
            sys.exit(1)

        X, Y = create_datapoint(SEQ[idx], CELLTYPE[idx],
                                int(SKIPPED_COUNT[idx]), int(INCLUDED_COUNT[idx]))
        # print(Y[0], Y[1])
        X_batch.extend(X)
        # print("##############")
        # print(Y[0])
        # print("Skipped count:", SKIPPED_COUNT[idx])
        # print("Included count:", INCLUDED_COUNT[idx])
        for t in range(1):
            Y_batch[t].extend([Y])
        Z_batch.extend(celltype)
        # print(np.array(Y_batch).shape)

    X_batch = np.array(X_batch)
    Y_batch = np.array(Y_batch)

    h5f2.create_dataset('X' + str(i), data=X_batch)
    h5f2.create_dataset('Y' + str(i), data=Y_batch)
    h5f2.create_dataset('Z' + str(i), data=np.asarray(Z_batch).astype("|S"))

h5f2.close()

print("--- %s seconds ---" % (time.time() - start_time))