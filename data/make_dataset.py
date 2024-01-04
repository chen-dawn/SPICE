import h5py
import numpy as np
import sys
import time
from constants import *
from utils import create_datapoint

start_time = time.time()

assert sys.argv[1] in ['train', 'test', 'all']
assert sys.argv[2] in ['0', '1', 'all']

h5f = h5py.File(data_dir + 'datafile'
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '.h5', 'r')


CHROM = h5f['CHROM'][:]
BARCODE = h5f['BARCODE'][:]
CELLTYPE = h5f['CELLTYPE'][:]
SEQ = h5f['SEQ'][:]
SKIPPED_COUNT = h5f['SKIPPED_COUNT'][:]
INCLUDED_COUNT = h5f['INCLUDED_COUNT'][:]

h5f.close()

h5f2 = h5py.File(data_dir + 'dataset'
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '.h5', 'w')

CHUNK_SIZE = 400
for i in range(SEQ.shape[0]//CHUNK_SIZE):
    print(i)
    # Each dataset has CHUNK_SIZE sequences. 
    if (i+1) == SEQ.shape[0]//CHUNK_SIZE:
        NEW_CHUNK_SIZE = CHUNK_SIZE + SEQ.shape[0]%CHUNK_SIZE
    else:
        NEW_CHUNK_SIZE = CHUNK_SIZE

    X_batch = []
    Y_batch = [[] for t in range(1)]
    
    for j in range(NEW_CHUNK_SIZE):
        idx = i*CHUNK_SIZE + j
        # celltype = CELLTYPE[idx].decode("utf-8")
        # # Strip the double quotes.
        # celltype = celltype.strip('"')
        # # print 
        # if celltype != "786O":
        #     continue

        X, Y = create_datapoint(SEQ[idx], CELLTYPE[idx],
                                SKIPPED_COUNT[idx], INCLUDED_COUNT[idx])
        # print(Y[0], Y[1])
        X_batch.extend(X)
        # print(Y[0][1])
        for t in range(1):
            Y_batch[t].extend([Y])
        # print(np.array(Y_batch).shape)

    X_batch = np.array(X_batch)
    Y_batch = np.array(Y_batch)

    X_batch = np.asarray(X_batch).astype('int8')
    Y_batch = np.asarray(Y_batch).astype('int8')

    h5f2.create_dataset('X' + str(i), data=X_batch)
    h5f2.create_dataset('Y' + str(i), data=Y_batch)

h5f2.close()

print("--- %s seconds ---" % (time.time() - start_time))