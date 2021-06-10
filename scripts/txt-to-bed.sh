#!/bin/bash

# USAGE: bash txt-to-bed.sh <inFILE.txt> <outFILE.bed>

awk -F '\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$3,$5,$5+1,$6,$8,$8+1,$4,$7}' $1 > temp.txt && \
awk 'BEGIN{FS=OFS="\t"} {$8=($8=="0") ? "+" : "-" ; $7=($7=="0") ? "+" : "-"}1' temp.txt > temp2.txt
sed '1d' temp2.txt > $2 && \
rm temp*.txt


#head -n -1 $1 > temp.txt && \
#printf '%s\n\n' "$(sed '1d' $1)" > temp.txt && \
#echo -e "$(sed '1d' $1)\n" > temp.txt && \
#awk '{print $3,$5,$5+1,$6,$8,$8+1,$4,$7}' $1 | \
#awk '{$8=($8=="0") ? "+" : "-" ; $7=($7=="0") ? "+" : "-"}1'
#awk -F '\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$3,$5,$5+1,$6,$8,$8+1,$4,$7}'

#awk -F '\t' '{$4=($4=="0") ? "+" : "-" ; $7=($7=="0") ? "+" : "-"}1' rearrangement.txt > temp.txt && \
#awk -F '\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$3,$5,$5+1,$6,$8,$8+1,{$4=($4=="0") ? "+" : "-"},{$7=($7=="0") ? "+" : "-"}}' temp.txt > temp2.txt && \
#sed '1d' temp2.txt > $2

