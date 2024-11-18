#!/bin/bash
zenodoID=10047603

# determine threads
OS=$(uname -s)
if [ "$OS" = "Darwin" -o "$OS" = "FreeBSD" ]; then
    n_threads=$(sysctl hw.ncpu|cut -f2 -d' ')
else
    n_threads=`grep -c ^processor /proc/cpuinfo`
fi

# taxonomy
if [[ -z "$1" ]]; then
    taxonomy="Bacteria"
else taxonomy=$1
fi

# paths
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
export LC_NUMERIC="en_US.UTF-8"
seqpath=$dir/../dat/seq/$taxonomy
seqpath_user=$dir/../dat/seq/$taxonomy/user
mkdir -p $seqpath/rev $seqpath/unrev $seqpath_user $seqpath/rxn

echo Checking updates for $taxonomy $seqpath

url_zenodorecord=https://zenodo.org/records/$zenodoID
url_zenodocurrent=$(curl -Ls -o /dev/null -w %{url_effective} $url_zenodorecord)
url_zenodoseqs=$url_zenodocurrent/files/md5sums.txt
status_zenodo=$(curl -s --head -w %{http_code} $url_zenodoseqs -o /dev/null)

if [[ $status_zenodo -ge 500 ]]; then
    echo "No sequence archive found, manual download needed."
    exit 1
fi

dir_rev=$seqpath/rev
dir_unrev=$seqpath/unrev
dir_rxn=$seqpath/rxn

# get md5 checksums from the online data version
zen_sums=$(mktemp)
wget -q -O $zen_sums $url_zenodoseqs
cat $zen_sums | grep -w $taxonomy > $zen_sums.$taxonomy
rm $zen_sums

# check for mismatches of local md5sums and online version
update_req=0
while read -r md5sum filepath; do
  if [ ! -f "$dir/../dat/seq/$filepath" ]; then
    update_req=1
    break
  fi

  # Calculate the MD5 sum of local file
  calculated_md5=$(md5sum "$dir/../dat/seq/$filepath" | cut -d ' ' -f 1)
  
  # Compare the calculated MD5 sum with the one from the file
  if [ "$md5sum" != "$calculated_md5" ]; then
    update_req=1
    break
  fi
done < $zen_sums.$taxonomy

if [[ $update_req -eq 0 ]]; then
  echo "Reference sequences are up-to-date."
  exit 0
fi

# Download and extract new sequences 
if [[ $update_req -eq 1 ]]; then
  echo "Updating reference sequences to latest version..."
  
  # download+extract taxonomy-specific archive
  allseqs=$(mktemp)
  wget -q -O "$dir/../dat/seq/$taxonomy.tar.gz" $url_zenodocurrent/files/$taxonomy.tar.gz
  tar xzf $dir/../dat/seq/$taxonomy.tar.gz -C $dir/../dat/seq/
  rm $dir/../dat/seq/$taxonomy.tar.gz
  
  # extract rev/unrev/rxn archives
  tar xzf $seqpath/rev/sequences.tar.gz -C $seqpath/rev/
  tar xzf $seqpath/unrev/sequences.tar.gz -C $seqpath/unrev/
  tar xzf $seqpath/rxn/sequences.tar.gz -C $seqpath/rxn/
  
  echo "Reference sequences updated."

  echo "Creating diamond database"
  # concat the fasta files for the query
  touch $seqpath/${taxonomy}_all.fasta
  # add the user files
  for fl in $(ls $seqpath/user/*.fasta); do
      name=$(basename $fl .fasta | tr " " "_")
      sed "s/>/>user_${name}|/g" $fl >> $seqpath/${taxonomy}_all.fasta
  done
  # add the reviewed files
  for fl in $(ls $seqpath/rev/*.fasta); do
      name=$(basename $fl .fasta| tr " " "_")
      sed "s/>/>rev_${name}|/g" $fl >> $seqpath/${taxonomy}_all.fasta
  done
  # add the unreviewed files
  for fl in $(ls $seqpath/unrev/*.fasta| tr " " "_"); do
      name=$(basename $fl .fasta)
      sed "s/>/>unrev_${name}|/g" $fl >> $seqpath/${taxonomy}_all.fasta
  done
  for fl in $(ls $seqpath/rxn/*.fasta| tr " " "_"); do
      name=$(basename $fl .fasta)
      sed "s/>/>rxn_${name}|/g" $fl >> $seqpath/${taxonomy}_all.fasta
  done
  diamond makedb --in $seqpath/${taxonomy}_all.fasta -d $seqpath/${taxonomy}_all --threads $n_threads --quiet
  rm $seqpath/${taxonomy}_all.fasta
  
  # create the subunit table
  Rscript $dir/make_subunitTable
fi

rm $zen_sums.$taxonomy

