#!/bin/bash
#SBATCH --mem=32g
#SBATCH --nodes=1
#SBATCH --job-name="long"
#SBATCH --ntasks-per-node=1
#SBATCH --time=11:59:00

source ~/.bashrc
conda activate py3

wkdir=/home/gc1a20/hgig_me/lrseq
cd $wkdir
mkdir -p ${1}_gene
cd ${1}_gene

mkdir -p $2
cp /home/gc1a20/bin/make_scores_mat.py .
cat ../$1| while read g;do
    grep -w "$g" ../meta_CADD15.txt |cut -f 1-26 > ${g}_c1
    grep -w "$g" ../meta_CADD15.txt |cut -f 27- | ~/bin/transpose.sh | while read i;
        do
            echo $i | sed 's/ /\n/g' |cut -f 2 -d':' | grep -w -v '\-9' |sort -u >ps
                if [ "$(wc -l ps | awk '{print $1}')" == "1" ]
                then
                echo $i | sed 's/ /\n/g' |cut -f 1 -d'|' | awk '{print"0/"$1}'  >h1
                echo $i | sed 's/ /\n/g' |cut -f 1 -d':' |cut -f 2 -d'|' | awk '{print"0/"$1}' >h2
                paste -sd'\t' h1 >>${g}_1
                paste -sd'\t' h2 >>${g}_2
                else
                echo $i | sed 's/ /\n/g'  |awk '{print"0/0"}'>h0
                paste -sd'\t' h0 >>${g}_1
                paste h0 -sd'\t' >>${g}_2
                fi
                done
    ~/bin/transpose.sh ${g}_1 | sed 's/ /\t/g' > ${g}_1c2
    ~/bin/transpose.sh ${g}_2 | sed 's/ /\t/g' > ${g}_2c2
    cp ../header.meta ${g}_1.meta
    cp ../header.meta ${g}_2.meta
    paste ${g}_c1 ${g}_1c2 | awk -F"\t" '{OFS=FS}{for (i=7;i<=16;i++) if(length($i)<1 || $i==0) $i="3.98e-6"}1' >> ${g}_1.meta
    ~/hgig_me/miniconda3/envs/py3/bin/python make_scores_mat.py --gene ${g}_1 --cadd $2
    paste ${g}_c1 ${g}_2c2 |awk -F"\t" '{OFS=FS}{for (i=7;i<=16;i++) if(length($i)<1 || $i==0) $i="3.98e-6"}1' >> ${g}_2.meta
    ~/hgig_me/miniconda3/envs/py3/bin/python make_scores_mat.py --gene ${g}_2 --cadd $2
    mv ${g}_1_$2_matrix $2/
    mv ${g}_2_$2_matrix $2/
    ls -1 ${g}_* | grep -v 'meta' |xargs rm
done

