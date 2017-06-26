#!/bin/bash

MIN_OVERLAP=0
MAX_DISTANCE=100
OUT_FILE_PREFIX=tandem
min_overlap_set=0
max_distance_set=0


bedtools=$(which bedtools)
if [[ -z "$bedtools" ]]; then
    echo "Cannot find bedtools."
    exit
fi

# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash

for i in "$@"
do
case $i in
    --min-overlap=*)
    MIN_OVERLAP="${i#*=}"
    min_overlap_set=1
    shift # past argument=value
    ;;
    --max-distance=*)
    MAX_DISTANCE="${i#*=}"
    max_distance_set=1
    shift # past argument=value
    ;;
    -o=*)
    OUT_FILE_PREFIX="${i#*=}"
    shift
    ;;
esac
done

echo "MIN_OVERLAP=$MIN_OVERLAP"
echo "MAX_DISTANCE=$MAX_DISTANCE"

#flags=$(($min_overlap_set+$max_distance_set))
#echo "FLAGS=$flags"
#if [[ "$flags" -gt 1 ]]; then
#    echo "Set either --min-overlap or --max-distance, not both!"
#    exit
#fi

if [[ -z "$1" ]]; then
    echo "Too few arguments."
    exit
elif [[ -z "$2" ]]; then
    echo "Too few arguments."
    exit
fi


#echo $1 
#echo $2


# https://stackoverflow.com/questions/18709507/maximum-and-minimum-using-awk

grep -F -f $1 $2 > $OUT_FILE_PREFIX.query.gff3
#exit
bedtools closest -d -k 2 -a <(grep $'\t'gene$'\t' $OUT_FILE_PREFIX.query.gff3 | sort -k1,1 -k4,4n) -b <(grep $'\t'gene$'\t' $2 | sort -k1,1 -k4,4n) | awk -v OFS='\t' '{ if ($19 == 0) { if ($13 >= $4) { $19 = -((($5<$14)?$5:$14) - $13 + 1); } else { $19 = -((($5<$14)?5:$14) - $4 + 1); } }; print $0; }' | awk -v OFS='\t' -v min_ovl=$MIN_OVERLAP -v max_dst=$MAX_DISTANCE '{ if (((min_ovl <= -($19))||(min_ovl == 0 && $19 <= max_dst)) && ($4 != $13 || $5 != $14)) { split($9, info1, ";"); split($18, info2, ";"); print $1,$4,$5,$7,$13,$14,$16,substr(info1[1], 4),substr(info2[1], 4),$19; } }' > $OUT_FILE_PREFIX.tmp.tsv


#Triticum_aestivum_CS42_TGACv1_scaffold_061102_1DL	82104	90815	-	80705	81710	+	TRIAE_CS42_1DL_TGACv1_061102_AA0185320	TRIAE_CS42_1DL_TGACv1_061102_AA0185310	394

awk -v OFS='\t' '{ if (($10 >= 0 || (($2 >= $5 && $3 >= $6) || ($5 >= $2 && $6 >= $3))) && ($4 == "+" && $7 == "+")) print $0 }' $OUT_FILE_PREFIX.tmp.tsv > $OUT_FILE_PREFIX.FF.tsv
awk -v OFS='\t' '{ if (($10 >= 0 || (($2 >= $5 && $3 >= $6) || ($5 >= $2 && $6 >= $3))) && ($4 == "-" && $7 == "-")) print $0 }' $OUT_FILE_PREFIX.tmp.tsv > $OUT_FILE_PREFIX.RR.tsv
awk -v OFS='\t' '{ if (($10 >= 0 || (($2 >= $5 && $3 >= $6) || ($5 >= $2 && $6 >= $3))) && ($4 == "+" && $7 == "-")) print $0 }' $OUT_FILE_PREFIX.tmp.tsv > $OUT_FILE_PREFIX.FR.tsv
awk -v OFS='\t' '{ if (($10 >= 0 || (($2 >= $5 && $3 >= $6) || ($5 >= $2 && $6 >= $3))) && ($4 == "-" && $7 == "+")) print $0 }' $OUT_FILE_PREFIX.tmp.tsv > $OUT_FILE_PREFIX.RF.tsv

awk -v OFS='\t' '{ if ($10 < 0 && (($2 >= $5 && $3 <= $6) || ($5 >= $2 && $6 <= $3))) print $0 }' $OUT_FILE_PREFIX.tmp.tsv > $OUT_FILE_PREFIX.full.tsv

rm $OUT_FILE_PREFIX.query.gff3
rm $OUT_FILE_PREFIX.tmp.tsv

