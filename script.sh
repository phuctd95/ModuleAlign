awk 'NR==FNR{a[$2]=$1;next}($1 in a){print a[$1]"\t"$2"\t"$3}' Node/$1.node $1-$2.blast > temp-1 
awk 'NR==FNR{a[$2]=$1;next}($2 in a){print $1"\t"a[$2]"\t"$3}' Node/$2.node temp-1 > blastNum
rm temp-1

./makeScore $1 $2 blastNum simFile
./align $1 $2 simFile -a $3