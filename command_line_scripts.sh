# To pretty clean the files.
for file in *umi_dedup_fine_grained_idx.csv; do
    echo $file
    awk -F',' 'NR==1 {print "index\tmode\toffset\tcount"} 
               NR>1 {
                   split($1, id_parts, "_")
                   name=""
                   for (i=1; i<=length(id_parts)-2; i++) {
                       name = (i == 1) ? id_parts[i] : name "_" id_parts[i]
                   }
                   mode = id_parts[length(id_parts)-1]
                   offset = id_parts[length(id_parts)]
                   print name"\t"mode"\t"offset"\t"$2
               }' $file >${file%.csv}_formatted.tsv
done
