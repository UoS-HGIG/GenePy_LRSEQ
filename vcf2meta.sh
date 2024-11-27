#!/bin/bash
#SBATCH --mem=16g
#SBATCH --nodes=1
#SBATCH --job-name="genepy_pre_local"
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00



wkdir=/home/gc1a20/hgig_me/lrseq
cd ${wkdir}


cp header.meta meta_CADDALL.txt
cp header.meta meta_CADD15.txt
cp header.meta meta_CADD20.txt
##variant info
grep -v '#' f.vcf > f1
cut -f 1-8 f1 >p1
cut -f 1-2,4-5 p1 >c1
cut -f 1-2,4-5 p1 |sed 's/\t/\_/g' >c1a

##align the order of alt allele as appears in c1

cut -f 4 c1 >alt
cut -f 8 f1 >c
awk -F";" '{for (i=1;i<=NF;i++) if ($i ~/CSQ\=/) print$i}' p1 |sed 's/CSQ\=//g' >csq
#awk -F";" '{for (i=1;i<=NF;i++) if ($i ~/CSQ\=/) print$i}' p1 | cut -d'|' -f1
paste p1 csq | awk '$5 !~/,/' | cut -f 1-7,9 >p1_s
awk '$5 ~/,/' p1 >p1_m

paste p1 csq | awk '$5 ~/,/' |while read i;
    do
                echo $i | cut -f 5 -d' ' |sed 's/\,/\n/g'>j
                        echo '*' >k
                                echo $i | cut -f 9 -d' ' | sed 's/\,/\n/g' >>k
                                        cat j |while read l
                                                do
                                                                grep -w "$l" k
                                                                        done > x1
                                                                                c1=$(wc -l x1 |cut -f 1 -d' ')
                                                                                        c2=$(wc -l j |cut -f 1 -d' ')
                                                                                                if [ $c1 -eq $c2 ]
                                                                                                            then
                                                                                                                            paste -sd',' x1
                                                                                                                                    else
                                                                                                                                                    paste -sd',' j
                                                                                                                                                            fi
                                                                                                                                                                done >order

paste p1_m order |awk '$9 !~/\|/' |cut -f 5 > alt_re
paste p1_m order |awk '$9 !~/\|/' |cut -f 1-8 > p1_re ###???
paste p1_m order |awk '$9 !~/\|/' | \
    cut -f 8 | \
    awk -F";" '{for (i=1;i<=NF;i++) if ($i ~/CSQ\=/) print$i}' |sed 's/CSQ\=//g' >csq_re
paste p1_m order |awk '$9 ~/\|/' |cut -f 1-7,9 > p1_1 #$9 to $8

paste alt_re csq_re |while read i
    do
        echo $i | cut -f 1 -d' ' |sed 's/\,/\n/g' | awk '{if (length($1)==1 && $1!~/\*/) print"--";else if ($1 ~/\*/) print"-*"; else print$i}'>j
        echo '*' >k
        echo $i | cut -f 2 -d' ' | sed 's/\,/\n/g' >>k
        cat j |while read l
        do
            m=${l:1}
            grep -w "$m" k
        done >x1
        c1=$(wc -l x1 |cut -f 1 -d' ')
        c2=$(wc -l j |cut -f 1 -d' ')
        if [ $c1 -eq $c2 ]
        then
            paste -sd',' x1
        else
            paste -sd',' j
        fi
    done > order_re


paste p1_re order_re |awk '$9 ~/\|/' |cut -f 1-7,9 > p1_2 #$9 to $8

cat p1_s p1_1 p1_2 |\
    sort -k1,1 -k2,2n |\
    awk -F"\t" '{print$1"_"$2"_"$4"_"$5,$6,$7,$8}'  >p1_order

awk 'NR==FNR{a[$1]=$0; next} {print a[$1]}' p1_order c1a >p1_u

cut -f 4 -d' ' p1_u|awk -F"," '{OFS=FS}{for (i=1;i<=NF;i++) if ($i ~/*/) $i="*|*|*|||||||||||||||||||||||"}1' > c_u


##allele funtional consequence
cut -f 2 -d'|' c_u  >c2 #Not used for GenePy

##gene with ensemblID; Note: there are 806 x-genes crossing chunks

cut -f 1-8 f.vcf > f11.vcf

module load biobuilds
bedtools intersect \
    -wao \
    -a f11.vcf \
    -b /mainfs/hgig/public/GENE_DATABASES/gencode.v45.annotation_3.bed |\
    cut -f 1-5,12 >f11.bed

/mainfs/hgig/private/software/datamash-1.8/datamash -g 1,2,3,4,5 collapse 6 <f11.bed |\
    cut -f 6 >c3

perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 |sort -u >gene.lst


##AF
cut -f 3 -d';' c_u |awk -F"|" '{OFS="\t"}{if ($5>0) print$6,$31,$56,$81,$106,$131,$156,$181,$206,$231; else print$17,$42,$67,$92,$117,$142,$167,$192,$217,$242}' >c4


##raw_score_all
cut -f 3 -d';' c_u |awk -F"|" '{OFS="\t"}{print$25,$50,$75,$100,$125,$150,$175,$200,$225,$250}' >c5

##phred_score >=15, which set smaller scores as 0
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<1.387112){$i="";}}1' c5 >c5a

##phred_score >=20
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<2.097252){$i="";}}1' c5 >c5b

#genotype
awk -F"," '{if ($0 !~/*/) print"0"; else for (i=1;i<=NF;i++) if ($i ~/\*/) print i}' alt >x1
cut -f 10- f1 >x2
paste x1 x2 |while read i; do pos=$(echo $i |cut -f 1 -d' '); echo $i | cut -f 2- -d' '| sed "s/${pos}/0/g" | sed 's/ /\t/g';done >c6


##merge;
paste c1 c2 c3 c4 c5 c6 >> meta_CADDALL.txt
paste c1 c2 c3 c4 c5a c6 >> meta_CADD15.txt
paste c1 c2 c3 c4 c5b c6 >> meta_CADD20.txt


#rm c* p1* alt_re order*
rm k j
