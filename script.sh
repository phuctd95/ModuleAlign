mkdir -p Node;
mkdir -p Tree;
mkdir -p Pairs;

#---------------------------------------------------
#Step 1(a) for first network
awk '{print $1"\n"$2}' $1.net > temp-1
awk '!x[$0]++' temp-1 > temp-2
awk '{print NR"\t"$1}' temp-2 > ./Node/$1.node
awk 'NR==FNR{a[$2]=$1;next}($1 in a && $2 in a){print a[$1]"\t"a[$2]}' ./Node/$1.node $1.net > ./Pairs/$1.pairs

#Steo 1(a) for second network
awk '{print $1"\n"$2}' $2.net > temp-1
awk '!x[$0]++' temp-1 > temp-2
awk '{print NR"\t"$1}' temp-2 > ./Node/$2.node

#Step 1(b)
awk 'NR==FNR{a[$2]=$1;next}($1 in a && $2 in a){print a[$1]"\t"a[$2]}' ./Node/$2.node $2.net > ./Pairs/$2.pairs
echo ".pair files of networks is made!\n"

rm temp-1;rm temp-2;

#Step 2
#----------------------------------------------------
./make_static_tree.exe ./Pairs/$1.pairs -n ./Tree/$1
./make_static_tree.exe ./Pairs/$2.pairs -n ./Tree/$2
#----------------------------------------------------

#Step 3
awk 'NR==FNR{a[$2]=$1;next}($1 in a){print a[$1]"\t"$2"\t"$3}' Node/$1.node $1-$2.blast > temp-1 
awk 'NR==FNR{a[$2]=$1;next}($2 in a){print $1"\t"a[$2]"\t"$3}' Node/$2.node temp-1 > blastNum
rm temp-1

#Step 4
g++ makeScore.cpp 
./a.out $1 $2 blastNum simFile

#run the main code
g++ Network.cpp Alignment.cpp ModuleAlign.cpp
./a.out $1  $2  simFile  -a $3
