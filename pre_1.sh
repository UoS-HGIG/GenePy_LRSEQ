#!/bin/bash
#SBATCH --mem=24g
#SBATCH --nodes=1
#SBATCH --job-name="genesis"
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00


cd /mainfs/hgig/private/gc1a20/lrseq

#bcftools view $1 -R /mainfs/hgig/public/HUMAN_REFS/HG38/target_capture/CCDS_hg38_pad25_sorted.bed -Ov -o f1.vcf
#awk '$1 ~/#/ || $9 ~/PS/' f1.vcf > f2.vcf
#bcftools view -G f2.vcf -Ov -o p1.vcf
#
#awk -F"\t" '$1 ~/#/ || length($4)>1||length($5)>1' p1.vcf | sed 's/^chr//g' > p11.vcf
#
###CADD annoation of the variants
#module unload biobuilds
##Using locally installed conda
#source ~/.bashrc
#conda activate cadd
#
#/mainfs/hgig/private/software/CADD-scripts/CADD.sh \
#    -c 8 \
#    -o lrseq_test_patch.tsv.gz p11.vcf
#
#conda deactivate
#
#
#tabix -p vcf lrseq_test_patch.tsv.gz
#mv lrseq_test_patch* /mainfs/hgig/private/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/
#
#
##Run VEP for each chromosome
#module load ensembl-vep/111.0
#vep -i p1.vcf \
#   		--offline \
#    		--assembly GRCh38 \
#    		--vcf \
#    		--fork 10 \
#    		--cache --force_overwrite \
#    		--pick_allele \
#    		--plugin CADD,/mainfs/hgig/private/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/lrseq_test_patch.tsv.gz,/mainfs/hgig/private/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/wes_202402_patch.tsv.gz,/mainfs/hgig/private/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/whole_genome_SNVs.tsv.gz,/mainfs/hgig/private/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/gnomad.genomes.r3.0.indel.tsv.gz \
#  		--custom /mainfs/hgig/private/software/gnomAD/GRCh38/gnomad.genomes.v3.1.1.RF_flag.vcf.gz,gnomadRF,vcf,exact,,RF_flag \
#   		--af_gnomade --af_gnomadg \
#    	--fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADg_NFE_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_OTH_AF,gnomADg_SAS_AF,gnomADe_AF,gnomADe_NFE_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,CADD_RAW,gnomadRF_RF_flag" \
#    		-o p1.vep.vcf
#module unload ensembl-vep


awk '$1 !~/##|\_/' f2.vcf |cut -f 9-> p2
grep -v '##' p1.vep.vcf >p1
grep '##' p1.vep.vcf >p12.vcf
paste p1 p2 >> p12.vcf

#
awk -F"\t" '$7 ~/PASS/ || $1 ~/#/' p12.vcf >f3.vcf
grep '#' p12.vcf >f.vcf
bcftools query -f '%CHROM %POS %ID %REF %ALT %QUAL %FILTER %INFO [\t%GT:%PS]\n' f3.vcf |sed 's/ /\t/g' |awk -F"\t" '{OFS=FS}{for (i=9;i<=NF;i++) if ($i ~/1\/1/) $i="1\|1:-9";else if ($i ~/\//) $i="0|0:-9"}1' >> f.vcf

