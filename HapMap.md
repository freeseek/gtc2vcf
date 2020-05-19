HapMap
======

How to A tutorial for how to convert HapMap data from Illumina and Affymetrix arrays to a GRCh38 VCF using gtc2vcf

<!--ts-->
   * [Download manifest files](#download-manifest-files)
   * [Download and unpack IDAT and CEL files](#download-and-unpack-idat-and-cel-files)
   * [Create sample maps](#create-sample-maps)
   * [Convert IDATs to GTCs](#convert-idats-to-gtcs)
   * [Convert GTCs to VCF](#convert-gtcs-to-vcf)
   * [Convert CELs to CHPs](#convert-cels-to-chps)
   * [Convert CHPs to VCF](#convert-chps-to-vcf)
<!--te-->

Download manifest files
=======================

Download HumanCNV370v1 manifest and cluster files from <a href="http://support.illumina.com/downloads/humancnv370-duo_v10_product_files.html">Illumina</a> and <a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6986">GEO</a>
```
wget ftp://webdata:webdata@ussd-ftp.illumina.com/downloads/ProductFiles/HumanCNV370/HumanCNV370-Duo/humancnv370v1_c.bpm
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanCNV370/HumanCNV370-Duo/HumanCNV370v1_C.egt
wget -O GSE21091_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE21091&format=file"
tar xvf GSE21091_RAW.tar
gunzip GPL6986_HumanCNV370v1_C.csv.gz
/bin/mv GPL6986_HumanCNV370v1_C.csv HumanCNV370v1_C.csv
/bin/rm GSE21091_RAW.tar
```

Download HumanOmni2.5-4v1 manifest and cluster files from <a href="http://support.illumina.com/downloads/humanomni2-5-quad_product_files.html">Illumina</a>
```
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/94afb35e-7c11-45cc-8a65-d868af527c54/HumanOmni2.5-4v1_H.bpm
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/f003e017-1761-4348-958f-03997a30cf67/HumanOmni2.5-4v1_H.egt
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/d5578cf6-bb3b-4b4b-98d3-21edc5bcbd45/HumanOmni2.5-4v1_H.csv
```

Download HumanOmni25M-8v1-1 manifest and cluster files from <a href="ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/humanomni25">Illumina</a> and <a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL20641">GEO</a>
```
wget https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL20nnn/GPL20641/suppl/GPL20641_HumanOmni2.5M-8v1-1_B.bpm.gz
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/humanomni25/humanomni2-5m-8v1-1_b.egt
wget https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL20nnn/GPL20641/suppl/GPL20641_HumanOmni25M-8v1-1_B.csv.gz
gunzip GPL20641_HumanOmni2.5M-8v1-1_B.bpm.gz
gunzip GPL20641_HumanOmni25M-8v1-1_B.csv.gz
mv GPL20641_HumanOmni2.5M-8v1-1_B.bpm HumanOmni2.5M-8v1-1_B.bpm
mv GPL20641_HumanOmni25M-8v1-1_B.csv HumanOmni25M-8v1-1_B.csv
```

Download GenomeWideEx_6 and GenomeWideSNP_6 library and annotation files from <a href="http://www.affymetrix.com/support/technical/byproduct.affx?product=genomewidesnp_6">Affymetrix</a>
```
wget http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/SNP6_supplemental_axiom_analysis_files.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
unzip -oj genomewidesnp6_libraryfile.zip CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.{cdf,chr{X,Y}probes,specialSNPs}
unzip -o SNP6_supplemental_axiom_analysis_files.zip GenomeWideSNP_6.{generic_prior.txt,apt-probeset-genotype.AxiomGT1.xml,AxiomGT1.sketch}
unzip -o GenomeWideSNP_6.na35.annot.csv.zip GenomeWideSNP_6.na35.annot.csv
/bin/rm genomewidesnp6_libraryfile.zip SNP6_supplemental_axiom_analysis_files.zip GenomeWideSNP_6.na35.annot.csv.zip
```

Re-align flank sequences to GRCh38
```
for chip in HumanCNV370v1_C humanomni25m-8v1-1_b HumanOmni2.5-4v1_H; do
  bcftools +gtc2vcf --csv $chip.csv --fasta-flank | \
    bwa mem -M $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna - | \
    samtools view -bS -o $chip.bam
done
bcftools +affy2vcf --csv GenomeWideSNP_6.na35.annot.csv --fasta-flank | \
  bwa mem -M $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna - | \
  samtools view -bS -o $chip.bam
```

Download and unpack IDAT and CEL files
======================================

```
wget http://bioconductor.org/packages/release/data/annotation/src/contrib/hapmap370k_1.0.1.tar.gz
wget -nH --cut-dirs 2 -r ftp://ftp.ncbi.nlm.nih.gov/hapmap/raw_data/hapmap3_affy6.0/
wget -nH --cut-dirs 5 -r ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_051208.txt
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped

mkdir -p idats
tar xzvf hapmap370k_1.0.1.tar.gz -C idats hapmap370k/inst/idatFiles
tar xzvf hd_genotype_chip/broad_intensities/Omni25_idats_gtcs_2141_samples.tgz -C idats
tar xzvf hd_genotype_chip/sanger_intensities/ALL.wgs.sanger_omni_2_5_8.20130805.snps.genotypes.idats.tar.gz -C idats

mkdir -p cels
for tgz in hapmap3_affy6.0/*.tgz; do tar xzvf $tgz -C cels; done
tar xzvf hd_genotype_chip/coriell_affy6_intensities/Affy60_Coriell_CEL_files.tar.gz -C cels

# one sample is mapped to HG03171 but should be mapped to HG01171, most likely a typo here
/bin/mv "cels/affy6/1000 Genomes phase 1 and 2 cel files/NA18489 .CEL" "cels/affy6/1000 Genomes phase 1 and 2 cel files/NA18489.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG03616.CEL" "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG03616-1.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG03660.CEL" "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG03660-1.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG04149.CEL" "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG04149-1.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG01171.CEL" "cels/affy6/1000 Genomes phase 1 and 2 cel files/HG01171-1.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 3 cel files/HG03616.CEL" "cels/affy6/1000 Genomes phase 3 cel files/HG03616-C1.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 3 cel files/HG03660.CEL" "cels/affy6/1000 Genomes phase 3 cel files/HG03660-C1.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 3 cel files/HG04149.CEL" "cels/affy6/1000 Genomes phase 3 cel files/HG04149-C1.CEL"
/bin/mv "cels/affy6/1000 Genomes phase 3 cel files/HG03171.CEL" "cels/affy6/1000 Genomes phase 3 cel files/HG01171-C1.CEL"
```

Create sample maps
==================

```
awk -F, 'NR>1 {print $5"\t"$1".HumanCNV370v1"}' idats/hapmap370k/inst/idatFiles/samples370k.csv > HapMap.HumanCNV370v1.tsv

awk -F, 'NR>15 {print $2"_"$3"\t"$6".HumanOmni2.5-4v1"}' idats/SampleSheet.csv > HapMap.HumanOmni2.5-4v1.tsv

awk 'NR==FNR {x[$2]=$1} NR>FNR {print $2"\t"x[substr($1,12)]".HumanOmni25M-8v1-1"}' \
  hd_genotype_chip/sanger_intensities/sanger_omni_chip.20130805.internal_to_coriell_id.map \
  idats/omni2.5-8_otgeno_20130805.idats/log.txt > HapMap.HumanOmni25M-8v1-1.tsv

# one sample is mapped to NA19787 but should be mapped to NA19730, most likely a sample swap
# samples mapped to NA21742 and NA21743 are the same individual, most likely a collection issue
cat hapmap3_affy6.0/{passing,excluded}_cels_sample_map.txt | sed 's/.CEL$//' | \
  sed 's/NA19787\tCHEAP_p_HapMapP3Redo2_GenomeWideSNP_6_B09_235604.CEL/NA19730\tCHEAP_p_HapMapP3Redo2_GenomeWideSNP_6_B09_235604.CEL/' | \
  awk '{sm=$1; if (sm in x) sm=sm"-"x[sm]; print $2"\t"sm".GenomeWideEx_6"; x[$1]++}' > HapMap.GenomeWideEx_6.tsv

ls cels/affy6/1000\ Genomes\ phase\ {1\ and\ 2,3}\ cel\ files/*.CEL | sed 's/.CEL$//' | \
  sed 's/.CEL$//' | awk -F/ '{print $4"\t"$4".GenomeWideSNP_6"}' > HapMap.GenomeWideSNP_6.tsv
```

Convert IDATs to GTCs
=====================

declare -A bpm=( ["HumanCNV370v1"]="humancnv370v1_c.bpm"
                 ["HumanOmni2.5-4v1"]="HumanOmni2.5-4v1_H.bpm"
                 ["HumanOmni25M-8v1-1"]="HumanOmni25M-8v1-1_B.bpm" )
declare -A egt=( ["HumanCNV370v1"]="HumanCNV370v1_C.egt"
                 ["HumanOmni2.5-4v1"]="HumanOmni2.5-4v1_H.egt"
                 ["HumanOmni25M-8v1-1"]="humanomni2-5m-8v1-1_b.egt" )
bcftools +gtc2vcf -i $(find idats -iname *.idat) -o gtc2vcf.idat.tsv
mkdir -p HumanCNV370v1 HumanOmni25M-8v1-1 HumanOmni2.5-4v1
for idat in $(cut -f1 gtc2vcf.idat.tsv | grep _Grn.idat$); do
  chip=$(grep ^$idat gtc2vcf.idat.tsv | cut -f16)
  mono $HOME/bin/autoconvert/AutoConvert.exe $(find idats -iname $idat) $chip ${bpm[$chip]} ${egt[$chip]}
done
bcftools +gtc2vcf {HumanCNV370v1,HumanOmni25M-8v1-1,HumanOmni2.5-4v1}/*.gtc -o gtc2vcf.gtc.tsv

Convert	GTCs to VCF
===================

declare -A bpm=( ["HumanCNV370v1"]="humancnv370v1_c.bpm"
                 ["HumanOmni2.5-4v1"]="HumanOmni2.5-4v1_H.bpm"
                 ["HumanOmni25M-8v1-1"]="HumanOmni25M-8v1-1_B.bpm" )
declare -A egt=( ["HumanCNV370v1"]="HumanCNV370v1_C.egt"
                 ["HumanOmni2.5-4v1"]="HumanOmni2.5-4v1_H.egt"
                 ["HumanOmni25M-8v1-1"]="humanomni2-5m-8v1-1_b.egt" )
declare -A csv=( ["HumanCNV370v1"]="HumanCNV370v1_C.csv"
                 ["HumanOmni2.5-4v1"]="HumanOmni2.5-4v1_H.csv"
                 ["HumanOmni25M-8v1-1"]="humanomni25m-8v1-1_b.csv" )
declare -A sam=( ["HumanCNV370v1"]="HumanCNV370v1_C.GRCh38.bam"
                 ["HumanOmni2.5-4v1"]="HumanOmni2.5-4v1_H.GRCh38.bam"
                 ["HumanOmni25M-8v1-1"]="humanomni25m-8v1-1_b.GRCh38.bam" )
for chip in HumanCNV370v1 HumanOmni25M-8v1-1 HumanOmni2.5-4v1; do
  bcftools +gtc2vcf \
    --no-version -Ou \
    --fasta-ref $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --sex hapmap.$chip.sex \
    --do-not-check-bpm \
    --bpm ${bpm[$chip]} \
    --egt ${egt[$chip]} \
    --csv ${csv[$chip]} \
    --sam ${sam[$chip]} \
    --adjust-clusters \
    --gtcs $chip | \
    bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
    bcftools norm --no-version -Ob -o HapMap.$chip.bcf -c x -f $ref && \
    bcftools index -f HapMap.$chip.bcf"
done

Convert CELs to CHPs
====================

(echo cel_files; find -iname *.CEL | grep _) > cels.GenomeWideEx_6.lst
(echo cel_files; find cels/affy6 -iname *.CEL) > cels.GenomeWideSNP_6.lst

for chip in GenomeWideEx_6 GenomeWideSNP_6; do
  mkdir -p $chip
  qsub -l h_rt=72:00:00,h_vmem=32g -N "hapmap.apt-probeset-genotype.$chip" -o uger/ ~/bin/UGERsubmit.sh \
    "apt-probeset-genotype \
    --out-dir $chip \
    --special-snps GenomeWideSNP_6.specialSNPs \
    --read-models-brlmmp GenomeWideSNP_6.generic_prior.txt \
    --chip-type $chip \
    --xml-file GenomeWideSNP_6.apt-probeset-genotype.AxiomGT1.xml \
    --cel-files cels.$chip.lst \
    --table-output false \
    --cc-chp-output \
    --cc-chp-out-dir $chip \
    --write-models"
done

Convert CHPs to VCF
===================

for chip in GenomeWideEx_6 GenomeWideSNP_6; do
  bcftools +affy2vcf \
    --no-version -Ou \
    --fasta-ref HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --csv GenomeWideSNP_6.na35.annot.csv \
    --sam GenomeWideSNP_6.na35.annot.bam \
    --models $chip/AxiomGT1.snp-posteriors.txt \
    --chps $chip | \
    bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
    bcftools norm --no-version -Ob -o HapMap.$chip.bcf -c x -f $ref && \
    bcftools index -f HapMap.$chip.bcf"
done
