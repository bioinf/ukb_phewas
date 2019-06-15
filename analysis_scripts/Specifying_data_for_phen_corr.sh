awk -F\\t  '{print $2"\t"$3}' all_sorted.tsv | sort -k 2,2 | grep -v '6152_6' > data_for_heatmap.tsv
