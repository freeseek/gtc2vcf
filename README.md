gtc2vcf
=======

A set of tools to convert Illumina and Affymetrix DNA microarray intensity data files into VCF files <b>without</b> using Microsoft Windows. You can use the final output to run the pipeline to detect [mosaic chromosomal alterations](https://github.com/freeseek/mocha). If you use this tool in your publication, please cite this website. For any feedback or questions, contact the [author](mailto:giulio.genovese@gmail.com)

![](gtc2vcf.png)

<!--ts-->
   * [Usage](#usage)
   * [Installation](#installation)
   * [Software Installation](#software-installation)
   * [Identifying chip type for IDAT and CEL files](#identifying-chip-type-for-idat-and-cel-files)
   * [Convert Illumina IDAT files to GTC files](#convert-illumina-idat-files-to-gtc-files)
   * [Convert Illumina GTC files to VCF](#convert-illumina-gtc-files-to-vcf)
   * [Convert Affymetrix CEL files to CHP files](#convert-affymetrix-cel-files-to-chp-files)
   * [Convert Affymetrix CHP files to VCF](#convert-affymetrix-chp-files-to-vcf)
   * [Using an alternative genome reference](#using-an-alternative-genome-reference)
   * [Plot variants](#plot-variants)
   * [Illumina GenCall](#illumina-gencall)
      * [Illumina AutoConvert](#illumina-autoconvert)
      * [Illumina AutoConvert 2.0](#illumina-autoconvert-2-0)
      * [Illumina Array Analysis Platform Genotyping Command Line Interface](#illumina-array-analysis-platform-genotyping-command-line-interface)
      * [Illumina Microarray Analytics Array Analysis Command Line Interface](#illumina-microarray-analytics-array-analysis-command-line-interface)
   * [Acknowledgements](#acknowledgements)
<!--te-->

Usage
=====

Illumina data tool:
```
Usage: bcftools +gtc2vcf [options] [<A.gtc> ...]

Plugin options:
    -l, --list-tags                   list available FORMAT tags with description for VCF output
    -t, --tags LIST                   list of output FORMAT tags [GT,GQ,IGC,BAF,LRR,NORMX,NORMY,R,THETA,X,Y]
    -b, --bpm <file>                  BPM manifest file
    -c, --csv <file>                  CSV manifest file (can be gzip compressed)
    -e, --egt <file>                  EGT cluster file
    -f, --fasta-ref <file>            reference sequence in fasta format
        --set-cache-size <int>        select fasta cache size in bytes
        --gc-window-size <int>        window size in bp used to compute the GC content (-1 for no estimate) [200]
    -g, --gtcs <dir|file>             GTC genotype files from directory or list from file
    -i, --idat                        input IDAT files rather than GTC files
        --capacity <int>              number of variants to read from intensity files per I/O operation [32768]
        --adjust-clusters             adjust cluster centers in (Theta, R) space (requires --bpm and --egt)
        --use-gtc-sample-names        use sample name in GTC files rather than GTC file name
        --do-not-check-bpm            do not check whether BPM and GTC files match manifest file name
        --do-not-check-eof            do not check whether the BPM and EGT readers reach the end of the file
        --genome-studio <file>        input a GenomeStudio final report file (in matrix format)
        --no-version                  do not append version and command line to the header
    -o, --output <file>               write output to a file [standard output]
    -O, --output-type u|b|v|z|t[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF
                                      t: GenomeStudio tab-delimited text output, 0-9: compression level [v]
        --threads <int>               number of extra output compression threads [0]
    -x, --extra <file>                write GTC metadata to a file
    -v, --verbose                     print verbose information
    -W, --write-index[=FMT]           Automatically index the output files [off]

Manifest options:
        --beadset-order               output BeadSetID normalization order (requires --bpm and --csv)
        --fasta-flank                 output flank sequence in FASTA format (requires --csv)
    -s, --sam-flank <file>            input flank sequence alignment in SAM/BAM format (requires --csv)
        --genome-build <assembly>     genome build ID used to update the manifest file [GRCh38]

Examples:
    bcftools +gtc2vcf -i 5434246082_R03C01_Grn.idat
    bcftools +gtc2vcf 5434246082_R03C01.gtc
    bcftools +gtc2vcf -b HumanOmni2.5-4v1_H.bpm -c HumanOmni2.5-4v1_H.csv
    bcftools +gtc2vcf -e HumanOmni2.5-4v1_H.egt
    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv -e GSA-24v3-0_A1_ClusterFile.egt -f human_g1k_v37.fasta -o GSA-24v3-0_A1.vcf
    bcftools +gtc2vcf -c HumanOmni2.5-4v1_H.csv -f human_g1k_v37.fasta 5434246082_R03C01.gtc -o 5434246082_R03C01.vcf
    bcftools +gtc2vcf -f human_g1k_v37.fasta --genome-studio GenotypeReport.txt -o GenotypeReport.vcf

Examples of manifest file options:
    bcftools +gtc2vcf -b GSA-24v3-0_A1.bpm -c GSA-24v3-0_A1.csv --beadset-order
    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv --fasta-flank -o GSA-24v3-0_A1.fasta
    bwa mem -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GSA-24v3-0_A1.fasta -o GSA-24v3-0_A1.sam
    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv --sam-flank GSA-24v3-0_A1.sam -o GSA-24v3-0_A1.GRCh38.csv
```

Affymetrix data tool:
```
Usage: bcftools +affy2vcf [options] --csv <file> --fasta-ref <file> [<A.chp> ...]

Plugin options:
    -l, --list-tags                 list available FORMAT tags with description for  VCF output
    -t, --tags LIST                 list of output FORMAT tags [GT,CONF,BAF,LRR,NORMX,NORMY,DELTA,SIZE]
    -c, --csv <file>                CSV manifest file (can be gzip compressed)
    -f, --fasta-ref <file>          reference sequence in fasta format
        --set-cache-size <int>      select fasta cache size in bytes
        --gc-window-size <int>      window size in bp used to compute the GC content (-1 for no estimate) [200]
        --probeset-ids              tab delimited file with column 'probeset_id' specifying probesets to convert
        --calls <file>              apt-probeset-genotype calls output (can be gzip compressed)
        --confidences <file>        apt-probeset-genotype confidences output (can be gzip compressed)
        --summary <file>            apt-probeset-genotype summary output (can be gzip compressed)
        --snp <file>                apt-probeset-genotype SNP posteriors output (can be gzip compressed)
        --chps <dir|file>           input CHP files rather than tab delimited files
        --cel <file>                input CEL files rather CHP files
        --adjust-clusters           adjust cluster centers in (Contrast, Size) space (requires --snp)
        --no-version                do not append version and command line to the header
    -o, --output <file>             write output to a file [standard output]
    -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
        --threads <int>             number of extra output compression threads [0]
    -x, --extra <file>              write CHP metadata to a file (requires CHP files)
    -v, --verbose                   print verbose information
    -W, --write-index[=FMT]         Automatically index the output files [off]

Manifest options:
        --fasta-flank               output flank sequence in FASTA format (requires --csv)
    -s, --sam-flank <file>          input flank sequence alignment in SAM/BAM format (requires --csv)

Examples:
    bcftools +affy2vcf \
        --csv GenomeWideSNP_6.na35.annot.csv \
        --fasta-ref human_g1k_v37.fasta \
        --chps cc-chp/ \
        --snp AxiomGT1.snp-posteriors.txt \
        --output AxiomGT1.vcf \
        --extra report.tsv
    bcftools +affy2vcf \
        --csv GenomeWideSNP_6.na35.annot.csv \
        --fasta-ref human_g1k_v37.fasta \
        --calls AxiomGT1.calls.txt \
        --confidences AxiomGT1.confidences.txt \
        --summary AxiomGT1.summary.txt \
        --snp AxiomGT1.snp-posteriors.txt \
        --output AxiomGT1.vcf

Examples of manifest file options:
    bcftools +affy2vcf -c GenomeWideSNP_6.na35.annot.csv --fasta-flank -o  GenomeWideSNP_6.fasta
    bwa mem -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GenomeWideSNP_6.fasta -o GenomeWideSNP_6.sam
    bcftools +affy2vcf -c GenomeWideSNP_6.na35.annot.csv -s GenomeWideSNP_6.sam -o GenomeWideSNP_6.na35.annot.GRCh38.csv
```

Installation
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges)
```
sudo apt install wget unzip git g++ zlib1g-dev bwa unzip samtools msitools cabextract mono-devel libgdiplus icu-devtools bcftools
```

Optionally, you can install these libraries to activate further HTSlib features
```
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/GRCh3{7,8} && cd /tmp
```

We recommend compiling the source code but, wherever this is not possible, Linux x86_64 pre-compiled binaries are available for download [here](http://software.broadinstitute.org/software/gtc2vcf). However, notice that you will require BCFtools version 1.20 or newer. You can also download a previous version of the plugin through [bioconda](https://anaconda.org/bioconda/bcftools-gtc2vcf-plugin)

Download latest version of [HTSlib](https://github.com/samtools/htslib) and [BCFtools](https://github.com/samtools/bcftools) (if not downloaded already)
```
wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
tar xjvf bcftools-1.20.tar.bz2
```

Download and compile plugins code (make sure you are using gcc version 5 or newer)
```
cd bcftools-1.20/
/bin/rm -f plugins/{idat2gtc.c,gtc2vcf.{c,h},affy2vcf.c}
wget -P plugins https://raw.githubusercontent.com/freeseek/gtc2vcf/master/{idat2gtc.c,gtc2vcf.{c,h},affy2vcf.c}
make
/bin/cp bcftools plugins/{idat2gtc,gtc2vcf,affy2vcf}.so $HOME/bin/
```

Make sure the directory with the plugins is available to BCFtools
```
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"
```

Install the GRCh37 human genome reference
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > $HOME/GRCh37/human_g1k_v37.fasta
samtools faidx $HOME/GRCh37/human_g1k_v37.fasta
bwa index $HOME/GRCh37/human_g1k_v37.fasta
```

Install the GRCh38 human genome reference (following the suggestion from [Heng Li](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use))
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bwa index $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Affymetrix provides the [Analysis Power Tools (APT)](https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html) for free which allow to call genotypes from raw intensity data using an algorithm derived from [BRLMM-P](http://tools.thermofisher.com/content/sfs/brochures/brlmmp_whitepaper.pdf)
```
mkdir -p $HOME/bin && cd /tmp
wget https://downloads.thermofisher.com/APT/APT_2.11.8/apt_2.11.8_linux_64_x86_binaries.zip
unzip -ojd $HOME/bin apt_2.11.8_linux_64_x86_binaries.zip apt_2.11.8_linux_64_x86_binaries/bin/apt-probeset-genotype
chmod a+x $HOME/bin/apt-probeset-genotype
```

Identifying chip type for IDAT and CEL files
============================================

To convert a pair of green and red IDAT files with raw Illumina intensities into a GTC file with genotype calls you need to provide both a BPM manifest file with the location of the probes and an EGT cluster file with the expected intensities of each genotype cluster. It is important to provide the correct BPM and EGT files otherwise the calling will fail possibly generating a GTC file with meaningless calls. Unfortunately newer IDAT files do not contain information about which BPM manifest file to use. The gtc2vcf bcftools plugin can be used to guess which files to use
```
path_to_idat_folder="..."
bcftools +gtc2vcf \
  -i -g $path_to_idat_folder
```
This will generate a spreadsheet table with information about each IDAT file including a guess for what manifest and cluster files you should use. If a guess is not provided, contact the [author](mailto:giulio.genovese@gmail.com) for troubleshooting

Similarly, you can use the affy2vcf bcftools plugin to extract chip type information from CEL files
```
path_to_cel_folder="..."
bcftools +affy2vcf \
  --cel --chps $path_to_cel_folder
```

Convert Illumina IDAT files to GTC files
========================================

The idat2gtc bcftools plugin can be used to convert Illumina IDAT files to GTC files
```
bpm_manifest_file="..."
egt_cluster_file="..."
bcftools +idat2gtc \
  --bpm $bpm_manifest_file \
  --egt $egt_cluster_file \
  --idats $path_to_idat_folder \
  --output $path_to_gtc_folder
```
The output is equivalent to the output of the Illumina GenCall algorithm while being significantly faster

If you do not have the manifest and cluster files for the Illumina IDAT files you are trying to convert, make sure to check the links [here][Illumina.md]

If you run the command with the option `--autocall-date ""` then the output should be deterministic and using the `--preset` option you can generate output equivalent to the output you obtain with any of the following:

* [Illumina AutoConvert](#autoconvert)
* [Illumina AutoConvert 2.0](#autoconvert-2-0)
* [Illumina Array Analysis Platform Genotyping Command Line Interface](#iaap-cli)
* [Illumina Microarray Analytics Array Analysis Command Line Interface](#array-analysis-cli)

If you similarly patch those tools to make them generate deterministic output, you should be able to verify that you get the same md5sum

Convert Illumina GTC files to VCF
=================================

Specifications for Illumina BPM, EGT, and GTC files were obtained through Illumina's [BeadArrayFiles](https://github.com/Illumina/BeadArrayFiles) library and [GTCtoVCF](https://github.com/Illumina/GTCtoVCF) script. Specifications for IDAT files were obtained through Henrik Bengtsson's [illuminaio](https://github.com/HenrikBengtsson/illuminaio) package
```
bpm_manifest_file="..."
csv_manifest_file="..."
egt_cluster_file="..."
path_to_gtc_folder="..."
ref="$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/GRCh37/human_g1k_v37.fasta"
out_prefix="..."
bcftools +gtc2vcf \
  --no-version -Ou \
  --bpm $bpm_manifest_file \
  --csv $csv_manifest_file \
  --egt $egt_cluster_file \
  --gtcs $path_to_gtc_folder \
  --fasta-ref $ref \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -o $out_prefix.bcf -Ob -c x -f $ref --write-index
```
Heavy random access to the reference will be needed, so it is important that enough extra memory be available for the operating system to cache the reference or else the task can run excruciatingly slowly. Notice that the gtc2vcf bcftools plugin will drop unlocalized variants. The final VCF might contain duplicates. If this is an issue `bcftools norm -d exact` can be used to remove such variants. At least one of the BPM or the CSV manifest files has to be provided. Normalized intensities cannot be computed without the BPM manifest file. Indel alleles cannot be inferred and will be skipped without the CSV manifest file. Information about genotype cluster centers will be included in the VCF if the EGT cluster file is provided. You can use gtc2vcf to convert one GTC file at a time, but we strongly advise to convert multiple files at once as single sample VCF files will consume a lot of storage space. If you convert hundreds of GTC files at once, you can use the `--adjust-clusters` option which will recenter the genotype clusters rather than using those provided in the EGT cluster file and will compute less noisy LRR values. If you use the `--adjust-clusters` option and you are using the output for calling [mosaic chromosomal alterations](https://github.com/freeseek/mocha), then it is safe to turn the median BAF/LRR adjustments off during that step (i.e. use `--adjust-BAF-LRR -1`)

Optionally, between the conversion and the sorting step you can include a `bcftools reheader --samples <file>` command to assign new names to the samples where `<file>` contains `old_name new_name\n` pairs separated by whitespaces, each on a separate line, with `old_name` being the GTC file name without the `.gtc` extension in this case

When running the conversion, the gtc2vcf plugin will double check that the SNP manifest metadata information in the GTC file matches the descriptor file name in the BPM file to make sure you are using the correct manifest file. Sometimes, due to discrepancies between the BPM file name provided by Illumina and the internal descriptor file name, this safety check fails. To turn off this feature in these cases, you can use option `--do-not-check-bpm`

Convert Affymetrix CEL files to CHP files
=========================================

Affymetrix provides a best practice workflow for genotyping data generated using [SNP6](https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-snp6-on-axiom.html) and [Axiom](https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-Axiom-probeset-genotype.html) arrays. As an example, the following command will run the genotyping for the Affymetrix SNP6 array:
```
path_to_output_folder="..."
cel_list_file="..."
apt-probeset-genotype \
  --analysis-files-path . \
  --xml-file GenomeWideSNP_6.apt-probeset-genotype.AxiomGT1.xml \
  --out-dir $path_to_output_folder \
  --cel-files $cel_list_file \
  --special-snps GenomeWideSNP_6.specialSNPs \
  --chip-type GenomeWideEx_6 \
  --chip-type GenomeWideSNP_6 \
  --table-output false \
  --cc-chp-output \
  --write-models \
  --read-models-brlmmp GenomeWideSNP_6.generic_prior.txt
```
Affymetrix provides Library and NetAffx Annotation files for their arrays ([here](http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays), [here](http://media.affymetrix.com/analysis/downloads/lf/genotyping), and [here](https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-annotation-files.html))

As an example, the following commands will obtain the files necessary to run the genotyping for the Affymetrix SNP6 array:
```
wget http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/SNP6_supplemental_axiom_analysis_files.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
unzip -oj genomewidesnp6_libraryfile.zip CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.{cdf,chrXprobes,chrYprobes,specialSNPs}
unzip -o SNP6_supplemental_axiom_analysis_files.zip GenomeWideSNP_6.{generic_prior.txt,apt-probeset-genotype.AxiomGT1.xml,AxiomGT1.sketch}
unzip -o GenomeWideSNP_6.na35.annot.csv.zip GenomeWideSNP_6.na35.annot.csv
```

Note: If the program exits due to different chip types or probe counts with error message such as `Wrong CEL ChipType: expecting: 'GenomeWideSNP_6' and #######.CEL is: 'GenomeWideEx_6'` then make sure you included the option `--chip-type GenomeWideEx_6 --chip-type GenomeWideSNP_6` or `--force` to the command line to solve the problem

Convert Affymetrix CHP files to VCF
===================================

The affy2vcf bcftools plugin can be used to convert Affymetrix CHP files to VCF
```
csv_manifest_file="..." # for example csv_manifest_file="GenomeWideSNP_6.na35.annot.csv"
ref="$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/GRCh37/human_g1k_v37.fasta"
path_to_chp_folder="cc-chp"
path_to_txt_folder="..."
out_prefix="..."
bcftools +affy2vcf \
  --no-version -Ou \
  --csv $csv_manifest_file \
  --fasta-ref $ref \
  --chps $path_to_chp_folder \
  --snp $path_to_txt_folder/AxiomGT1.snp-posteriors.txt \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -o $out_prefix.bcf -Ob -c x -f $ref --write-index
```
Heavy random access to the reference will be needed, so it is important that enough extra memory be available for the operating system to cache the reference or else the task can run excruciatingly slowly. The final VCF might contain duplicates. If this is an issue `bcftools norm -d exact` can be used to remove such variants. There is often no need to use the `--adjust-clusters` option for Affymetrix data as the cluster posteriors are already adjusted using the data processed by the genotype caller

Optionally, between the conversion and the sorting step you can include a `bcftools reheader --samples <file>` command to assign new names to the samples where `<file>` contains `old_name new_name\n` pairs separated by whitespaces, each on a separate line, with `old_name` being the CHP file name without the `.chp` extension

Using an alternative genome reference
=====================================

Illumina provides [GRCh38/hg38](https://support.illumina.com/bulletins/2017/04/infinium-human-genotyping-manifests-and-support-files--with-anno.html) manifests for many of its genotyping arrays. However, if your genotyping array is not supported for the newer reference by Illumina, you can use the `--fasta-flank` and `--sam-flank` options to realign the flank sequences from the manifest files you have and recompute the marker positions. This approach uses [flank sequence](https://support.illumina.com/bulletins/2016/05/infinium-genotyping-manifest-column-headings.html) and [strand](https://support.illumina.com/bulletins/2017/06/how-to-interpret-dna-strand-and-allele-information-for-infinium-.html) information to identify the marker [coordinates](https://support.illumina.com/bulletins/2016/06/-infinium-genotyping-array-manifest-files-what-does-chr-or-mapinfo---mean.html). It will need a sequence aligner such as `bwa` to realign the sequences and it seems to reproduce the coordinates provided from Illumina more than 99.9% of the times. Mapping information will follow the [implicit dbSNP standard](https://github.com/Illumina/GTCtoVCF#manifests). Occasionally the flank sequence provided by Illumina is incorrect and it is impossible to recover the correct marker coordinate from the flank sequence alone

You first have to generate an alignment file for the flank sequences from a CSV manifest file
```
csv_manifest_file="..."
ref="$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/GRCh37/human_g1k_v37.fasta"
bam_alignment_file="..."
bcftools +gtc2vcf \
  -c $csv_manifest_file \
  --fasta-flank | \
  bwa mem -M $ref - | \
  samtools view -bS \
  -o $bam_alignment_file
```
Notice that you need to use the `-M` option to mark shorter split hits as secondary and you should not sort the output BAM file as gtc2vcf expects it to have the sequences in the same order as in the CSV file . Then you load the alignment file while converting your GTC files to VCF including the `-s $bam_alignment_file` option

Some older manifest files from Illumina have thousands of markers with incorrect RefStrand annotations that will lead to incorrect genotypes. While Illumina has not explained why this is the case, it still distributes incorrect manifests. If you are using one of the following manifests
```
Human1M-Duov3_H
Human610-Quadv1_H
Human660W-Quad_v1_H
HumanCytoSNP-12v2-1_Anova
HumanOmni1-Quad_v1-0-Multi_H
HumanOmni1-Quad_v1-0_H
```
We advise to either contact Illumina to demand a fixed version or to use gtc2vcf to realign the flank sequences

Also, Illumina assigns chromosomal positions to indels by first left aligning the flank sequences in an incoherent way (see [here](https://github.com/Illumina/GTCtoVCF/blob/develop/BPMRecord.py)). Apparently this is incoherent enough that Illumina also cannot get the coordinates of homopolymer indels right. For example, chromosome 13 ClinVar indel [rs80359507](https://www.ncbi.nlm.nih.gov/clinvar/variation/37959) is assigned to position 32913838 in the manifest file for the GSA-24v2-0 array, but it is assigned to position 32913837 in the manifest file for GSA-24v3-0 array (GRCh37 coordinates). If you want to trust genotypes at homopolymer indels, we advise to use gtc2vcf to realign the flank sequences

The same functionality exists for the affy2vcf tool to convert Affymetrix data

Plot variants
=============

Install basic tools (Debian/Ubuntu specific if you have admin privileges):
```
sudo apt install r-cran-optparse r-cran-ggplot2 r-cran-data.table r-cran-gridextra
```

Download R scripts
```
/bin/rm -f $HOME/bin/gtc2vcf_plot.R
wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/gtc2vcf/master/gtc2vcf_plot.R
chmod a+x $HOME/bin/gtc2vcf_plot.R
```

Plot variant (for Illumina data)
```
gtc2vcf_plot.R \
  --illumina \
  --vcf input.vcf \
  --chrom 11 \
  --pos 66328095 \
  --png rs1815739.png
```

![](rs1815739.png)

Plot variant (for Affymetrix data)
```
gtc2vcf_plot.R \
  --affymetrix \
  --vcf input.vcf \
  --chrom 1 \
  --pos 196642233 \
  --png rs800292.png
```

![](rs800292.png)

Illumina GenCall
================

To genotype raw Illumina IDAT intensity files using Illumina GenCall algorithms, Illumina over the course of the year has provided several command line interfaces written in the .NET language:
- [AutoConvert](https://support.illumina.com/array/array_software/beeline/downloads.html) (2011)
- [AutoConvert 2.0](https://support.illumina.com/array/array_software/beeline/downloads.html)) (2017)
- [IAAP CLI](https://support.illumina.com/array/array_software/illumina-array-analysis-platform.html) (2019)
- [Array Analysis CLI](http://support.illumina.com/array/array_software/ima-array-analysis-cli/downloads.html) (2023)
We provide instructions to install and run these interfaces. The `sed -i -e ':a' -e 'N' -e '$!ba'` installation commands are used to prevent the interfaces from timestamping the output GTC files by removing the [System.DateTime](https://learn.microsoft.com/en-us/dotnet/api/system.datetime) calls and accesses to the [CreationTime](https://learn.microsoft.com/en-us/dotnet/api/system.io.filesysteminfo.creationtime) property from the binaries, with the goal of making each execution completely reproducible. AutoConvert 2.0, IAAP-CLI, and Array Analysis CLI binaries will both perform version 1.2.0 of the normalization step and seem to produce the exact same results while AutoConvert will only perform version 1.1.2 of the normalization step yielding somewhat different results. If you want to run these binaries but fail to download them, contact the [author](mailto:giulio.genovese@gmail.com) for troubleshooting

Illumina also provides the [Beeline](https://support.illumina.com/array/array_software/beeline.html) software for free and this includes the AutoConvert.exe command line executable which allows to call genotypes from raw intensity data using Illumina's proprietary GenCall algorithm. AutoConvert is almost entirely written in Mono/.Net language, except for one small mathmatical function (findClosestSitesToPointsAlongAxis) which is included within a Windows PE32+ library (MathRoutines.dll). As this is [unmanaged code](http://www.mono-project.com/docs/advanced/embedding/), to be run on Linux with [Mono](https://www.mono-project.com/) it needs to be embedded in an equivalent Linux ELF64 library (libMathRoutines.dll.so) as shown below. This function is run as part of the [normalization](https://doi.org/10.1093/bioinformatics/btm443) of the raw intensities when sampling [400 candidate homozygotes](https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina_gt_normalization.pdf) before calling genotypes.

Illumina AutoConvert
--------------------

To run Illumina AutoConvert (version 1.6.3.1) you will need to fix the hardcoded Windows [backlashes](https://en.wikipedia.org/wiki/Backslash) into UNIX [slashes](https://en.wikipedia.org/wiki/Slash_(punctuation), as shown below
```
mkdir -p $HOME/bin && cd /tmp
wget http://support.illumina.com/content/dam/illumina-support/documents/downloads/software/beeline/autoconvert-software-v1-6-3-installer.zip
wget http://raw.githubusercontent.com/freeseek/gtc2vcf/master/nearest_neighbor.c
unzip -o autoconvert-software-v1-6-3-installer.zip 
msiextract -C Illumina/AutoConvert SetupAutoConvert64_1.6.3.1.msi
msiextract -l SetupAutoConvert64_1.6.3.1.msi | grep DLL$ | while read dll; do mv Illumina/AutoConvert/$dll Illumina/AutoConvert/${dll%DLL}dll; done
gcc -fPIC -shared -O2 -o Illumina/AutoConvert/libMathRoutines.dll.so nearest_neighbor.c
sed -i 's/\x00\x03\\\x00/\x00\x03\/\x00/' Illumina/AutoConvert/AutoCallLib.dll
sed -i 's/G\x00R\x00N\x00.\x00i\x00d\x00a\x00t\x00/G\x00r\x00n\x00.\x00i\x00d\x00a\x00t\x00/' Illumina/AutoConvert/AutoCallLib.dll
sed -i 's/R\x00E\x00D\x00.\x00i\x00d\x00a\x00t\x00/R\x00e\x00d\x00.\x00i\x00d\x00a\x00t\x00/' Illumina/AutoConvert/AutoCallLib.dll
sed -i 's/\\\x00M\x00o\x00d\x00u\x00l\x00e\x00s\x00\\\x00B\x00S\x00G\x00T\x00\\\x00C\x00l\x00u\x00s\x00t\x00e\x00r\x00A\x00l\x00g\x00o\x00r\x00i\x00t\x00h\x00m\x00s\x00\\\x00/\/\x00M\x00o\x00d\x00u\x00l\x00e\x00s\x00\/\x00B\x00S\x00G\x00T\x00\/\x00C\x00l\x00u\x00s\x00t\x00e\x00r\x00A\x00l\x00g\x00o\x00r\x00i\x00t\x00h\x00m\x00s\x00\/\x00/' Illumina/AutoConvert/AutoCallLib.dll
sed -i 's/\\\x00M\x00o\x00d\x00u\x00l\x00e\x00s\x00\\\x00B\x00S\x00G\x00T\x00/\/\x00M\x00o\x00d\x00u\x00l\x00e\x00s\x00\/\x00B\x00S\x00G\x00T\x00/' Illumina/AutoConvert/Modules/BSGT/ClusterAlgorithms/{GoldenGate/GGCA,InfiniumII/I2CA,GenTrain/ILCA}.dll
sed -i 's/\\\x00d\x00a\x00t\x00.\x00b\x00i\x00n\x00/\/\x00d\x00a\x00t\x00.\x00b\x00i\x00n\x00/' Illumina/AutoConvert/Modules/BSGT/ClusterAlgorithms/{GoldenGate/GGCA,InfiniumII/I2CA,GenTrain/ILCA}.dll
sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\x28\xa6\x00\x00\x0a\x13\x40\x12\x40\x28\xa7\x00\x00\x0a\x72\xad\x12\x00\x70\x28\xa6\x00\x00\x0a\x13\x40\x12\x40\x28\xa8\x00\x00\x0a\x28\x23\x00\x00\x0a/\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x7e\x16\x00\x00\x0a\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00/' Illumina/AutoConvert/AutoCallLib.dll
sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\x11\x0e\x6f\xe5\x00\x00\x0a\x13\x11\x12\x11\x28\xe6\x00\x00\x0a\x72\xad\x12\x00\x70\x11\x0e\x6f\xe5\x00\x00\x0a\x13\x12\x12\x12\x28\xe7\x00\x00\x0a\x28\x23\x00\x00\x0a/\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x7e\x16\x00\x00\x0a\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00/' Illumina/AutoConvert/AutoCallLib.dll
rm autoconvert-software-v1-6-3-installer.zip SetupAutoConvert64_1.6.3.1.msi nearest_neighbor.c
mv Illumina/AutoConvert $HOME/bin/
rmdir Illumina
```

You can run Illumina's proprietary GenCall algorithm on a single IDAT file pair
```
mono $HOME/bin/AutoConvert/AutoConvert.exe \
  $idat_green_file \
  $path_to_output_folder \
  $bpm_manifest_file \
  $egt_cluster_file
```
Make sure that the red IDAT file is in the same folder as the green IDAT file. Alternatively you can run on multiple IDAT file pairs
```
mono $HOME/bin/AutoConvert/AutoConvert.exe \
  $path_to_idat_folder \
  $path_to_output_folder \
  $bpm_manifest_file \
  $egt_cluster_file
```

Illumina AutoConvert 2.0
------------------------

To run Illumina AutoConvert 2.0 (version 2.0.1.179) you will need to separately download an additional Mono/.Net library (Heatmap.dll) from [GenomeStudio](https://support.illumina.com/array/array_software/genomestudio.html) or the [polyploid clustering module](https://support.illumina.com/downloads/genomestudio_polyploid_clustering_module_v1-0_software.html) and include it in your binary directory, most likely due to differences in which Mono and .Net resolve library dependencies, as shown below
```
mkdir -p $HOME/bin && cd /tmp
wget http://support.illumina.com/content/dam/illumina-support/documents/downloads/software/beeline/autoconvert-software-v2-0-1-installer.zip
wget http://support.illumina.com/content/dam/illumina-support/documents/downloads/software/genomestudio/genomestudiopolyploidclusteringv1-0.msi
wget http://raw.githubusercontent.com/freeseek/gtc2vcf/master/nearest_neighbor.c
unzip -o autoconvert-software-v2-0-1-installer.zip
msiextract AutoConvertInstaller.msi
msiextract genomestudiopolyploidclusteringv1-0.msi
mv Heatmap.DLL Illumina/AutoConvert\ 2.0/
gcc -fPIC -shared -O2 -o Illumina/AutoConvert\ 2.0/libMathRoutines.dll.so nearest_neighbor.c
sed -i 's/^\(     <AutosomalCallRateThreshold>\)0.97\(<\/AutosomalCallRateThreshold>\r\)$/\10.0\2/' Illumina/AutoConvert\ 2.0/AutoCallConfig.xml
sed -i 's/\\\x00d\x00a\x00t\x00.\x00b\x00i\x00n\x00/\/\x00d\x00a\x00t\x00.\x00b\x00i\x00n\x00/' Illumina/AutoConvert\ 2.0/{GGCA,I2CA,HDCA,ILCA,ILCA3}.dll
sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\x28\xc7\x00\x00\x0a\x13\x3f\x12\x3f\x28\xc8\x00\x00\x0a\x72\xa8\x15\x00\x70\x28\xc7\x00\x00\x0a\x13\x3f\x12\x3f\x28\xc9\x00\x00\x0a\x28\x1f\x00\x00\x0a/\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x7e\x12\x00\x00\x0a\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00/' Illumina/AutoConvert\ 2.0/AutoCallLib.dll
msiextract -l genomestudiopolyploidclusteringv1-0.msi | grep -v Heatmap.DLL | xargs rm
rmdir Modules/BSPC/clusteralgorithms/*
rmdir -p Modules/BSPC/clusteralgorithms
rm autoconvert-software-v2-0-1-installer.zip AutoConvertInstaller.msi genomestudiopolyploidclusteringv1-0.msi nearest_neighbor.c
mv Illumina/AutoConvert\ 2.0 $HOME/bin/
rmdir Illumina
```
We change the autosomal call rate threshold to 0.0 to more aggressively call gender in lower quality samples

If you need to get the Heatmap.dll library from GenomeStudio indtead, you can use the following code
```
wget ftp://webdata2:webdata2@ftp.illumina.com/downloads/software/genomestudio/genomestudio-software-v2-0-4-5-installer.zip
unzip -oj genomestudio-software-v2-0-4-5-installer.zip
cabextract GenomeStudioInstaller.exe
msiextract a0
mv Illumina/GenomeStudio\ 2.0/Heatmap.dll Illumina/AutoConvert\ 2.0/
rm genomestudio-software-v2-0-4-5-installer.zip GenomeStudioInstaller.exe {,a}0 u{0..5} Illumina/GenomeStudio\ 2.0 -r
```

You can run Illumina's proprietary GenCall algorithm on a single IDAT file pair
```
mono $HOME/bin/AutoConvert\ 2.0/AutoConvert.exe \
  $idat_green_file \
  $path_to_output_folder \
  $bpm_manifest_file \
  $egt_cluster_file
```
Make sure that the red IDAT file is in the same folder as the green IDAT file. Alternatively you can run on multiple IDAT file pairs
```
mono $HOME/bin/AutoConvert\ 2.0/AutoConvert.exe \
  $path_to_idat_folder \
  $path_to_output_folder \
  $bpm_manifest_file \
  $egt_cluster_file
```

Make sure that the IDAT files have the same name prefix as the IDAT folder name. The software might require up to 8GB of RAM to run. Illumina provides manifest (BPM) and cluster (EGT) files for their arrays [here](https://support.illumina.com/array/downloads.html). Notice that if you provide the wrong BPM file, you will get an error such as: `Normalization failed!  Unable to normalize!` and if you provide the wrong EGT file, you will get an error such as `System.Exception: Unrecoverable Error...Exiting! Unable to find manifest entry ######## in the cluster file!`

Illumina Array Analysis Platform Genotyping Command Line Interface
------------------------------------------------------------------

Illumina provides the [Illumina Array Analysis Platform Genotyping Command Line Interface](https://support.illumina.com/array/array_software/illumina-array-analysis-platform.html) software for free for research use and this includes the iaap-cli 1.1.0 which runs natively on Linux
```
mkdir -p $HOME/bin && cd /tmp
wget ftp://webdata2:webdata2@ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz
tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz -C $HOME/bin/ iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1
sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\x28\x17\x01\x00\x0a\x13\x07\x12\x07\x72\xdd\x23\x00\x70\x28\x18\x01\x00\x0a/\x00\x00\x00\x00\x00\x00\x00\x00\x00\x7e\x92\x00\x00\x0a\x00\x00\x00\x00\x00/' $HOME/bin/iaap-cli/ArrayAnalysis.NormToGenCall.Services.dll
rm iaap-cli-linux-x64-1.1.0.tar.gz
```

Once iaap-cli is properly installed in your system, run Illumina's proprietary GenCall algorithm on multiple IDAT file pairs
```
CLR_ICU_VERSION_OVERRIDE="$(uconv -V | sed 's/.* //g')" LANG="en_US.UTF-8" $HOME/bin/iaap-cli/iaap-cli \
  gencall \
  $bpm_manifest_file \
  $egt_cluster_file \
  $path_to_output_folder \
  --idat-folder $path_to_idat_folder \
  --output-gtc \
  --gender-estimate-call-rate-threshold 0.0
```
It is important to set the `LANG` environmental variable to `en_US.UTF-8`, if this is set to other values, due to a bug in `iaap-cli` causing malformed GTC files to be generated as a result. Due to another bug in `iaap-cli`, IDAT filenames cannot include more than two `_` characters and should be formatted as `BARCODE_POSITION_(Red|Grn).idat`. When using `iaap_cli` you cannot process old array manifest files with loci data encoded as version 5 or older, such as `HumanHap650Yv3_A.bpm`, as the corresponding code was not carried over and you will get the error `Error in reading file.  Unknown Manifest version`. The AutoConvert command line tool can read older manifest files. We change the autosomal call rate threshold to 0.0 both to more aggressively call gender in lower quality samples and to deal with an implementation issue that causes loci with null cluster scores to be included in the determination of the autosomal call rate threshold

Illumina Microarray Analytics Array Analysis Command Line Interface
-------------------------------------------------------------------

Illumina provides the [Illumina Microarray Analytics Array Analysis Command Line Interface](http://support.illumina.com/array/array_software/ima-array-analysis-cli/downloads.html) software for free for research use and this includes the array-analysis-cli 2.1.0 which runs natively on Linux
```
mkdir -p $HOME/bin && cd /tmp
wget http://support.illumina.com/softwaredownload.html?assetId=72f8a34f-0933-4256-bad6-73d830436c74&assetDetails=IlluminaMicroarrayAnalyticsArrayAnalysisCLIv2.1LinuxInstaller-2.1-array-analysis-cli-linux-x64-v2.1.0.tar.gz
tar xzvf array-analysis-cli-linux-x64-v2.1.0.tar.gz -C $HOME/bin/ --strip-components=1
sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\x28\x89\x00\x00\x0a\x0A\x12\x00\x72\xa3\x15\x00\x70\x28\x8a\x00\x00\x0a/\x00\x00\x00\x00\x00\x00\x00\x00\x72\xfc\x0d\x00\x70\x00\x00\x00\x00\x00/' $HOME/bin/array-analysis-cli//ArrayAnalysis.Core.dll
rm array-analysis-cli-linux-x64-v2.1.0.tar.gz
```

Once array-analysis-cli is properly installed in your system, run Illumina's proprietary GenCall algorithm on multiple IDAT file pairs
```
$HOME/bin/array-analysis-cli/array-analysis-cli \
  genotype call \
  --bpm-manifest $bpm_manifest_file \
  --cluster-file $egt_cluster_file \
  --idat-folder .
```
We cannot change the autosomal call rate threshold to 0.0 both to more aggressively call gender in lower quality samples as the default 0.97 value is hardcoded

Acknowledgements
================

This work is supported by NIH grant [R01 HG006855](http://grantome.com/grant/NIH/R01-HG006855), NIH grant [R01 MH104964](http://grantome.com/grant/NIH/R01-MH104964), NIH grant [R01MH123451](http://grantome.com/grant/NIH/R01-MH123451), US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244), and the Stanley Center for Psychiatric Research
