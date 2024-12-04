#!/bin/bash

# Porthmeus
# 23.10.24

# keep track of time
start_time=`date +%s`


pathways=""
database="seed"
pwyDatabase="metacyc,custom"
verbose=1
taxonomy="Bacteria"
taxRange="all" # taxonomic range for pathways
bitcutoff=200 # cutoff blast: min bit score
identcutoff=0   # cutoff blast: min identity
identcutoff_exception=70  # min identity for enzymes marked as false friends (hight seq similarity but different function)
covcutoff=75 # cutoff blast: min coverage
subunit_cutoff=50 # more than this % of subunits must be found 
strictCandidates=false
completenessCutoff=66 # consider pathway to be present if other hints (e.g. key enzyme present) are avaiable and pathway completeness is at least as high as completenessCutoff (requires strictCandidates=false)
completenessCutoffNoHints=80 # consider pathway to be present if no hints are avaiable (requires stricCandidates=false)
blast_back=false
noSuperpathways=true
vagueCutoff=0.3 # cutoff for vague reactions. If the amount of vague reactions in a pathways is more then this their influence will not be recognized even with strictCandidates=false
onlyList=false
skipBlast=false
includeSeq=false
use_parallel=true
exhaustive=false
seqSrc=2
anno_genome_cov=false
use_gene_seq=true
stop_on_files_exist=false
update_manually=false
user_temp=false
force_offline=false
input_mode="auto"
output_dir=.
OS=$(uname -s)
if [ "$OS" = "Darwin" -o "$OS" = "FreeBSD" ]; then
    n_threads=$(sysctl hw.ncpu|cut -f2 -d' ')
else
    n_threads=`grep -c ^processor /proc/cpuinfo`
fi

# paths and eariables
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
script_name=$(basename -- "$0")
uniprotIdentity=0.9 # clustered uniprot database (0.5 or 0.9)
metaPwy=$dir/../dat/meta_pwy.tbl
keggPwy=$dir/../dat/kegg_pwy.tbl
seedPwy=$dir/../dat/seed_pwy.tbl
customPwy=$dir/../dat/custom_pwy.tbl
metaRea=$dir/../dat/meta_rea.tbl
reaDB1=$dir/../dat/vmh_reactions.tsv
reaDB2=$dir/../dat/bigg_reactions.tbl
reaDB3=$dir/../dat/seed_reactions_corrected.tsv
reaDB4=$dir/../dat/mnxref_seed.tsv
reaDB5=$dir/../dat/mnxref_seed-other.tsv
reaDB6=$dir/../dat/mnxref_bigg-other.tsv
brenda=$dir/../dat/brenda_ec_edited.csv
seedEC=$dir/../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv
seedEnzymesNames=$dir/../dat/seed_Enzyme_Name_Reactions_Aliases.tsv
altecdb=$dir/../dat/altec.csv
metaGenes=$dir/../dat/meta_genes.csv

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# newly introduced with the diamond search
query_fasta=""
diamond_out=""


usage() {
    echo "Usage"
    echo "$0 -p <keyword> / -e <EC> [-d <database>] [-t <taxonomy>] file.fasta"
    echo "  -p keywords such as pathways or subsystems (for example amino,nucl,cofactor,carbo,polyamine)"
    echo "  -e Search by ec numbers (comma separated)"
    echo "  -r Search by enzyme name (colon separated)"
    echo "  -d Database: vmh or seed (default: $database)"
    echo "  -t Taxonomic range for reference sequences to be used. (Bacteria, Archaea, auto; default: $taxonomy). See Details."
    echo "  -b Bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i Identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c Coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -s Strict candidate reaction handling (do _not_ use pathway completeness, key kenzymes and operon structure to infere if imcomplete pathway could be still present (default: $strictCandidates)"
    echo "  -u Suffix used for output files (default: pathway keyword)"
    echo "  -a blast hits back against uniprot enzyme database"
    echo "  -n Consider superpathways of metacyc database"
    echo "  -l Select the pathway database (MetaCyc, KEGG, SEED, all; default: $pwyDatabase)"
    echo "  -o Only list pathways found for keyword (default: $onlyList)"
    echo "  -x Do not blast only list pathways, reactions and check for available sequences (default: $skipBlast)"
    echo "  -q Include sequences of hits in log files (default: $includeSeq)"
    echo "  -v Verbose level, 0 for nothing, 1 for pathway infos, 2 for full (default: $verbose)"
    echo "  -k Do not use parallel (Deprecated: use '-K 1' instead to disable multi-threading.)"
    echo "  -g Exhaustive search, continue blast even when cutoff is reached (default: $exhaustive)"
    echo "  -z Quality of sequences for homology search: 1:only reviewed (swissprot), 2:unreviewed only if reviewed not available, 3:reviewed+unreviewed, 4:only unreviewed (default: $seqSrc)"
    echo "  -m Limit pathways to taxonomic range (default: $taxRange)"
    echo "  -w Use additional sequences derived from gene names (default: $use_gene_seq)"
    echo "  -y Print annotation genome coverage (default: $anno_genome_cov)"
    echo "  -j Quit if output files already exist (default: $stop_on_files_exist)"
    echo "  -f Path to directory, where output files will be saved (default: current directory)"
    echo "  -U Do not use gapseq sequence archive and update sequences from uniprot manually (very slow) (default: $update_manually)"
    echo "  -T Set user-defined temporary folder (default: $user_temp)"
    echo "  -O Force offline mode (default: $force_offline)"
    echo "  -M Input genome mode. Either 'nucl', 'prot', or 'auto' (default '$input_mode')"
    echo "  -K Number of threads for sequence alignments. If option is not provided, number of available CPUs will be automatically determined."
    echo ""
    echo "Details:"
    echo "\"-t\": if 'auto', gapseq tries to predict if the organism is Bacteria or Archaea based on the provided genome sequence. The prediction is based on the 16S rRNA gene sequence using a classifier that was trained on 16S rRNA genes from organisms with known Gram-staining phenotype. In case no 16S rRNA gene was found, a k-mer based classifier is used instead."

exit 1
}



while getopts "h?p:e:r:d:i:b:c:v:st:nou:al:oxqkgz:m:ywjf:UT:OM:K:" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    p)  
        pathways=$OPTARG
        ;;
    e)  
        ecnumber=$OPTARG
        ;;
    r)  
        reaname="$OPTARG"
        ;;
    d)  
        database=$OPTARG
        ;;
    v)  
        verbose=$OPTARG
        ;;
    b)
        bitcutoff=$OPTARG
        ;;
    i)
        identcutoff=$OPTARG
        ;;
    c)
        covcutoff=$OPTARG
        ;;
    t)  
        taxonomy=$OPTARG
        ;;
    s)
        strictCandidates=true
        ;;
    u)
        output_suffix=$OPTARG
        ;;
    a)
        blast_back=true
        includeSeq=true
        ;;
    n)
        noSuperpathways=false
        ;;
    l)
        pwyDatabase=$OPTARG
        ;;
    o)
        onlyList=true
        ;;
    x)
        skipBlast=true
        ;;
    q)
        includeSeq=true
        ;;
    k)
        use_parallel=false
        n_threads=1
        echo "DEPRECATION NOTICE: Option '-k' is deprecated. To disable multi-threading use '-K 1' instead."
        ;;
    g)
        exhaustive=true
        ;;
    z)
        seqSrc=$OPTARG
        ;;
    m)
        taxRange=$OPTARG
        ;;
    y)
        anno_genome_cov=true
        ;;
    w)
        use_gene_seq=true
        ;;
    j)
        stop_on_files_exist=true
        ;;
    f)
        output_dir=$OPTARG
        ;;
    U)
        update_manually=true
        ;;
    T)
        user_temp=true
        user_temp_folder=$OPTARG
        ;;
    O)
        force_offline=true
        ;;
    M)
        input_mode=$OPTARG
        ;;
    K)
        n_threads=$OPTARG
        if [ $n_threads -eq 1 ]; then
            use_parallel=false
        fi
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

# after parsing arguments, only fasta file should be there - print usage if not
[ "$#" -ne 1 ] && { usage; }
one_dollar=$1


set_output_dir(){
    # set output directory
    case $output_dir in
        /*)
            # absolute path
            output_dir=$output_dir
            ;;
        ~*)
            # relative to $HOME directory
            output_dir="${output_dir/#\~/$HOME}"
            ;;
        *)
            # relative path to current directory
            output_dir=$curdir/$output_dir
            ;;
    esac
    # create path if it does not yet exist
    mkdir -p $output_dir || { echo "$output_dir not writable. Aborting..."; exit 1; }
}


tempdir_orga(){
    # tmp working directory
    fasta=$(readlink -f "$1") # save input file before changing to temporary directory
    tmp_fasta=$(basename "${fasta}" .gz | tr ' ' '_')
    if [[ "$user_temp" = true ]]; then
        mkdir -p $user_temp_folder
        tmpdir=$(mktemp -d $user_temp_folder/"$tmp_fasta"_XXXXXX)
    else
        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT
    fi
    cd $tmpdir
}

fasta_orga(){
    # get fasta file
    if [[ "$fasta" == *.gz ]]; then # in case fasta is in a archive
        gunzip -c "$fasta" > "$tmp_fasta"
        fasta="$tmp_fasta"
    fi
    [[ ! -s "$fasta" ]] && { echo Invalid file: $1; exit 0; }
    tmpvar=$(basename "$fasta")
    fastaID="${tmpvar%.*}"
}

determine_seqtype(){
    # Determine if fasta is nucl or prot
    echo $input_mode
    if [ $input_mode == "auto" ]; then
        n_char=`cat $fasta | grep -v "^>" | awk '{for(i=1;i<=NF;i++)if(!a[$i]++)print $i}' FS="" | wc -l`
        if [ $n_char -ge 15 ]; then
            echo "Protein fasta detected."
            input_mode="prot"
        else
            echo "Nucleotide fasta detected."
            input_mode="nucl"
        fi
    fi
}


translate_nucl(){
    if [ $input_mode == "nucl" ]; then
        echo "Will translate nucleotide to protein seqeuence"
        prodigal -i $fasta -o /dev/null -a ${tmp_fasta}_prodigal.fasta
        fasta=${tmp_fasta}_prodigal.fasta
        input_mode="prot"
    fi
    echo $fasta
    cp $fasta $output_dir/$fasta
}

determine_taxonomy(){
    # determine taxonomy
    if [ $input_mode == "prot" ] && [ $taxonomy == "auto" ]; then
        cp $dir/../dat/seq/hmm/domain.hmm.gz .
        gunzip domain.hmm.gz
        hmmsearch --tblout $fastaID.tblout --cpu $n_threads domain.hmm $fasta > /dev/null
        taxonomy=`Rscript $dir/predict_domain.R "$dir" "$fastaID.tblout"`
        rm domain.hmm
        rm $fastaID.tblout
        
        echo Predicted taxonomy: $taxonomy
    fi
    if [ $input_mode == "nucl" ] && [ "$taxonomy" == "auto" ]; then
        pred_biom=$($dir/predict_biomass_from16S.sh "$fasta")
        if [ "$pred_biom" == "Gram_neg" ] || [ "$pred_biom" == "Gram_pos" ]; then
            taxonomy=Bacteria
        elif [ "$pred_biom" == "Archaea" ]; then
            taxonomy=Archaea
        else
            echo "Taxonomy could be predicted automatically. Assuming default case: Bacteria (Use '-t' parameter to modify it)."
            taxonomy=Bacteria
        fi
        echo Predicted taxonomy: $taxonomy
    fi
    [[ "$taxonomy" == "bacteria" ]] && taxonomy=Bacteria
    [[ "$taxonomy" == "archaea" ]] &&  taxonomy=Archaea

    # Follow taxonomy prediction for pathway tax range if set to "auto"
    if [ "$taxRange" == "auto" ]; then
        taxRange=$taxonomy
    fi
}

define_sequence_directory(){
    # squence directory
    export LC_NUMERIC="en_US.UTF-8"
    seqpath=$dir/../dat/seq/$taxonomy
    seqpath_user=$dir/../dat/seq/$taxonomy/user
    mkdir -p $seqpath/rev $seqpath/unrev $seqpath_user
}

check_for_updates(){
    #check for updates if internet connection is available
    if [[ "$force_offline" = false ]]; then
        $dir/update_sequences.sh $taxonomy
        if [[ ! -f $seqpath/rev/sequences.tar.gz  ]] || [[ ! -f $seqpath/unrev/sequences.tar.gz ]] || [[ ! -f $seqpath/rxn/sequences.tar.gz ]]; then
            echo ERROR: gapseq sequence archives are missing! Sequences needed to be downloaded from uniprot directly!
            exit 1
        fi
    fi
}

create_query_seq(){
    # concat the fasta files for the query
    touch query.fasta
    # add the user files
    for fl in $(ls $seqpath_user/*.fasta); do
        name=$(basename $fl .fasta | tr " " "_")
        sed "s/>/>user_${name}|/g" $fl >> query.fasta
    done
    # add the reviewed files
    for fl in $(ls $seqpath/rev/*.fasta); do
        name=$(basename $fl .fasta| tr " " "_")
        sed "s/>/>rev_${name}|/g" $fl >> query.fasta
    done
    # add the unreviewed files
    for fl in $(ls $seqpath/unrev/*.fasta| tr " " "_"); do
        name=$(basename $fl .fasta)
        sed "s/>/>unrev_${name}|/g" $fl >> query.fasta
    done
    query_fasta=query.fasta
    cp $query_fasta $output_dir/query.fasta
}


create_blastresults(){
    # $fasta should be fasta to the reference genome
    # $query_fasta should be the query fasta
    db_file=$seqpath/${taxonomy}_all
    db_bin=$(basename $fasta)
    diamond_out=${db_bin}_blastresult.tsv
    diamond blastp -q $fasta -d $db_file --outfmt 6 qseqid pident evalue bitscore scovhsp sseqid sstart send qstart qend --threads $n_threads -o $diamond_out 
    cp $diamond_out $output_dir/$diamond_out
}


#echo $n_threads
#parse_opts $@
#echo 1
set_output_dir
#echo 2
tempdir_orga $one_dollar
fasta_orga $one_dollar
#echo 3
determine_seqtype
#echo 4
determine_taxonomy
translate_nucl
define_sequence_directory
#echo 5
check_for_updates
#echo 6
#time_query=$(time create_query_seq)
#echo $time_query
#echo $query_fasta
#echo 7
create_blastresults
#echo 8
echo DEBUG
echo $(pwd $diamond_out)
echo "$(pwd $diamond_out)/$diamond_out"
echo DEBUG
Rscript $dir/gapseq_find_all.R $diamond_out $taxonomy $taxRange $seqSrc $dir $noSuperpathways $identcutoff $identcutoff_exception $database $subunit_cutoff $bitcutoff $covcutoff $completenessCutoff $completenessCutoffNoHints $vagueCutoff

echo $time_blast
ps -q $$ -o %cpu,%mem,args
end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.
