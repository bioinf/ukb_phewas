counter=0
for link in $(cat links.txt)
do
counter=$((counter+1)) 
name=`awk -v var="$counter" 'BEGIN {RS = "" ; FS = "\n" }{print $var}' names.txt`
del_name=`awk -v var="$counter" 'BEGIN {RS = "" ; FS = "\n" }{print $var}' int_names.txt`
wget $link -O - | gunzip > $del_name
awk '$9<5e-08' $del_name > $name
rm -rf $del_name
done

