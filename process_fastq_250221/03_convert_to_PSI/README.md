# Convert counts to PSI

We have the counts of the included and exon-skipped reads for each sample. We want to convert these to PSI.

## Merging the counts from the same modes.
Because now I'm outputting the `insert_size` for every index and event mode, I am noticing that sometimes the insert size has some error tolerance. E.g. some has insert size 100, while some has insert size 101. The mapping is probably correct because I'm only using 2 bp edit distance. So first, we would merge the counts from the same index and event modes. 

## Get a list of all events.

I am realizing that every index has many "INCLUDED" and "SKIPPED" events - not all of them are mapped to "0:0:0" which means the coordinates are exactly the same as reference. However, this is not the same for every cell type. So we want to first generate a reference list of all exon skipping events that are observed across all cell types, then calculate the PSI for every event for every cell type. This will be able to generate a comprehensive matrix of PSI values. 