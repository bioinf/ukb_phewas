dir_name=$1
plink_bin=$2
bim_prefix=$3

#exrtact desired files from directory

for name in $(cd $dir_name; ls -Al | grep 'clump'| awk '{print $9}')
do
	echo $name
	#create output title
	title="$(echo $name | sed 's/\.vcf//' | sed 's/\Common_SNP_jaccard_union_n2_//' | sed 's/\_depict.tsv//' | sed 's/\_to_clump//')"
	# run plink -- beware absolute pathawys
	mkdir -p clumped
	$plink_bin --bfile $bim_prefix  \
	--clump ${dir_name}/$name \
	--clump-field P \
	--clump-p1 5e-08 \
	--clump-r2 0.1 \
	--clump-snp-field SNP \
	--clump-kb 500 \
	--out clumped/$title.clumped
done



