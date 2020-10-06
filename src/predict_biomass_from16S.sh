#!/bin/bash

gn=$1

path=$(readlink -f "$0")
dir=$(dirname "$path")


tmpdir=$(mktemp -d)
cp $gn $tmpdir
cd $tmpdir


gnseqbase=`basename $gn`

# remove empty lines as this caused errors in bedtools
cat $gn | sed '/^[[:space:]]*$/d' | tr -d '\15\32' > $gnseqbase.tmp

# renaming contigs to ensure compatibility with barrnap
awk '/^>/{print ">" ++i; next}{print}' $gnseqbase.tmp > $gnseqbase.tmp2 && mv $gnseqbase.tmp2 $gnseqbase.tmp

# run barrnap to find 16S genes
barrnap --quiet $gnseqbase.tmp | grep 16S > $gnseqbase.gff

# extract 16S gene sequence from genome
bedtools getfasta -s -name -fi $gnseqbase.tmp -bed $gnseqbase.gff -fo $gnseqbase.16S.fasta


fasta="$gnseqbase.16S.fasta"
name="${fasta%.*}"

blastn -db $dir/../dat/seq/16S_db/16S_BMTest.fna -query $fasta -outfmt 6 -qcov_hsp_perc 95 -perc_identity 80 | sort -k12,12nr -k3,3nr | head -11 | awk '$1{print $2}' | grep -Po '=(.*)' > $name.grpred.tmp

posc=`grep -c Gram_pos $name.grpred.tmp`
negc=`grep -c Gram_neg $name.grpred.tmp`
arcc=`grep -c Archaea $name.grpred.tmp`

res=""
if [ "$posc" -gt "$negc" -a "$posc" -gt "$arcc" ]; then
  res="Gram_pos"
elif [ "$negc" -gt "$posc" -a "$negc" -gt "$arcc" ]; then
  res="Gram_neg"
elif [ "$arcc" -gt "$posc" -a "$arcc" -gt "$negc" ]; then
  res="Archaea"
elif [ "$posc" -eq "$negc" ]; then
  res="ambiguous"
fi

if [[ $res == "ambiguous" ]]; then
  if [[ `grep -c "" $name.grpred.tmp` > 0 ]]; then
    restmp=`head -n 1  $name.grpred.tmp`
    res="${restmp#"="}"
  fi
fi

# cleaning up
rm -f -r $tmpdir

echo $res