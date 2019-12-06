gtc2vcf
=======

A set of tools to convert Illumina and Affymetrix array intensity data files into VCF files <b>without</b> using Microsoft Windows. You can use the final output to run the pipeline to detect <a href="https://github.com/freeseek/mocha">mosaic chromosomal alterations</a>. If you use this tool in your publication, please cite this website. For any feedback, send an email to giulio.genovese@gmail.com

![](gtc2vcf.png)

Usage
=====

Illumina tool:
```
Usage: bcftools +gtc2vcf [options] <A.gtc> [...]

Plugin options:
    -l, --list-tags              list available tags with description for VCF output
    -t, --tags LIST              list of output tags [IGC,BAF,LRR,NORMX,NORMY,R,THETA,X,Y]
    -i  --idat <file>            IDAT intensity data file
    -b  --bpm <file>             BPM manifest file
    -c  --csv <file>             CSV manifest file
    -e  --egt <file>             EGT cluster file
    -f, --fasta-ref <file>       reference sequence in fasta format
        --set-cache-size <int>   select fasta cache size in bytes
    -g, --gtc-list <file>        read list of GTC file names from file
        --adjust-clusters        adjust cluster centers in (Theta, R) space (requires --bpm and --egt)
    -x, --sex <file>             output GenCall gender estimate into file
        --do-not-check-bpm       do not check whether BPM and GTC files match manifest file name
        --genome-studio <file>   input a genome studio final report file (in matrix format)
        --beadset-order          output BeadSetID normalization order (requires --bpm and --csv)
        --no-version             do not append version and command line to the header
    -o, --output <file>          write output to a file [standard output]
    -O, --output-type b|u|z|v|g  b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF, g: GenomeStudio [v]
        --threads <int>          number of extra output compression threads [0]

Examples:
    bcftools +gtc2vcf -i 5434246082_R03C01_Grn.idat
    bcftools +gtc2vcf 5434246082_R03C01.gtc
    bcftools +gtc2vcf -b HumanOmni2.5-4v1_H.bpm -c HumanOmni2.5-4v1_H.csv
    bcftools +gtc2vcf -e HumanOmni2.5-4v1_H.egt
    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv -e GSA-24v3-0_A1_ClusterFile.egt -f human_g1k_v37.fasta -o GSA-24v3-0_A1.vcf
    bcftools +gtc2vcf -c HumanOmni2.5-4v1_H.csv -f human_g1k_v37.fasta 5434246082_R03C01.gtc -o 5434246082_R03C01.vcf
    bcftools +gtc2vcf -f human_g1k_v37.fasta --genome-studio GenotypeReport.txt -o GenotypeReport.vcf
```

Affymetrix tool:
```
Usage: bcftools +affy2vcf [options] --fasta-ref <file> --annot <file> --snp-posteriors <file>
                                    --summary <file> --calls <file> --confidences <file>

Plugin options:
    -f, --fasta-ref <file>       reference sequence in fasta format
        --set-cache-size <int>   select fasta cache size in bytes
        --annot <file>           probeset annotation file
        --snp-posteriors <file>  apt-probeset-genotype snp-posteriors output
        --summary <file>         apt-probeset-genotype summary output
        --report <file>          apt-probeset-genotype report output
        --calls <file>           apt-probeset-genotype calls output
        --confidences <file>     apt-probeset-genotype confidences output
        --adjust-clusters        adjust cluster centers in (Contrast, Size) space
    -x, --sex <file>             output apt-probeset-genotype gender estimate into file (requires --report)
        --no-version             do not append version and command line to the header
    -o, --output <file>          write output to a file [standard output]
    -O, --output-type b|u|z|v    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
        --threads <int>          number of extra output compression threads [0]

Example:
    bcftools +affy2vcf \
        --fasta-ref human_g1k_v37.fasta \
        --annot GenomeWideSNP_6.na35.annot.csv \
        --snp-posteriors AxiomGT1.snp-posteriors.txt \
        --summary AxiomGT1.summary.txt \
        --calls AxiomGT1.calls.txt \
        --confidences AxiomGT1.confidences.txt \
        --output AxiomGT1.vcf
```

Installation
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges)
```
sudo apt install wget autoconf zlib1g-dev gzip unzip samtools msitools cabextract mono-devel libgdiplus
```

Optionally, you can install these libraries to activate further HTSlib features:
```
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res && cd /tmp
```

Download latest version of <a href="https://github.com/samtools/htslib">HTSlib</a> and <a href="https://github.com/samtools/bcftools">BCFtools</a> (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Download plugins code
```
/bin/rm -f bcftools/plugins/{gtc,affy}2vcf.c
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/gtc2vcf/master/{gtc,affy}2vcf.c
```

Compile latest version of HTSlib (optionally disable bz2, gcs, and lzma) and BCFtools (make sure you are using gcc version 5 or newer)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{gtc,affy}2vcf.so} $HOME/bin/
```

Make sure the directory with the plugins is available to bcftools
```
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"
```

Install the GRCh37 human genome reference
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > $HOME/res/human_g1k_v37.fasta
samtools faidx $HOME/res/human_g1k_v37.fasta
```

Install the GRCh38 human genome reference (following the suggestion from <a href="http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use">Heng Li</a>)
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Software Installation
=====================

Illumina provides the <a href="https://support.illumina.com/array/array_software/illumina-array-analysis-platform.html">Illumina Array Analysis Platform</a> software for free and this includes the iaap-cli command line executable which runs natively on Linux
```
mkdir -p $HOME/bin && cd /tmp
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz
tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz -C $HOME/bin/ iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1
```
However, notice that in some older Linux machines this approach does not work

Alternatively, Illumina provides the <a href="https://support.illumina.com/array/array_software/beeline.html">Beeline</a> software for free and this includes the AutoConvert.exe command line executable which allows to call genotypes from raw intensity data using Illumina's proprietary GenCall algorithm. AutoConvert is almost entirely written in Mono/.Net language, with the exception of one small mathmatical function (findClosestSitesToPointsAlongAxis) which is contained instead within a Windows PE32+ library (MathRoutines.dll). As this is <a href="http://www.mono-project.com/docs/advanced/embedding/">unmanaged code</a>, to be run on Linux with <a href="https://www.mono-project.com/">Mono</a> it needs to be embedded in an equivalent Linux ELF64 library (libMathRoutines.dll.so) as shown below. This function is run as part of the <a href="http://doi.org/10.1093/bioinformatics/btm443">normalization</a> of the raw intensities when sampling 400 <a href="http://patft.uspto.gov/netacgi/nph-Parser?patentnumber=7035740">candidate homozygotes</a> before calling genotypes (see also <a href="https://patents.google.com/patent/US7035740">here</a>). For some unclear reasons, you will also need to separately download an additional Mono/.Net library (Heatmap.dll) from <a href="https://support.illumina.com/array/array_software/genomestudio.html">GenomeStudio</a> and include it in your binary directory as shown below, most likely due to differences in which Mono and .Net resolve library dependencies
```
mkdir -p $HOME/bin && cd /tmp
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/beeline/autoconvert-software-v2-0-1-installer.zip
unzip -o autoconvert-software-v2-0-1-installer.zip
msiextract AutoConvertInstaller.msi
cp -R Illumina/AutoConvert\ 2.0 $HOME/bin/autoconvert

wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/genomestudio/genomestudio-software-v2-0-4-5-installer.zip
unzip -o genomestudio-software-v2-0-4-5-installer.zip
cabextract GenomeStudio-software-v2-0-4-5-installer/GenomeStudioInstaller.exe
msiextract a0
cp Illumina/GenomeStudio\ 2.0/Heatmap.dll $HOME/bin/autoconvert/

wget https://raw.githubusercontent.com/freeseek/gtc2vcf/master/nearest_neighbor.c
gcc -fPIC -shared -O2 -o $HOME/bin/autoconvert/libMathRoutines.dll.so nearest_neighbor.c
```
If you fail to download the autoconvert software, contact the <a href="mailto:giulio.genovese@gmail.com">author</a> for troubleshooting. Notice that this approach to run AutoConvert on Linux is <strong>not</strong> supported by Illumina

Affymetrix provides the <a href="https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html">Analysis Power Tools (APT)</a> for free which allow to call genotypes from raw intensity data using an algorithm derived from <a href="http://tools.thermofisher.com/content/sfs/brochures/brlmmp_whitepaper.pdf">BRLMM-P</a>
```
mkdir -p $HOME/bin && cd /tmp
wget https://downloads.thermofisher.com/APT/APT2.10.2.2/APT_2.10.2.2_Linux_64_bitx86_binaries.zip
unzip -ojd $HOME/bin APT_2.10.2.2_Linux_64_bitx86_binaries.zip apt-2.10.2.2-x86_64-intel-linux/bin/apt-probeset-genotype
chmod a+x $HOME/bin/apt-probeset-genotype
```

Convert Illumina IDAT files to GTC files
========================================

Once iaap-cli is properly installed in your system, run Illumina's proprietary GenCall algorithm on multiple IDAT file pairs
```
$HOME/bin/iaap-cli/iaap-cli gencall \
  $bpm_manifest_file \
  $egt_cluster_file \
  $path_to_output_folder \
  -f $path_to_idat_folder -g
```

Alternatively, once Mono and AutoConvert are properly installed on your system, run Illumina's proprietary GenCall algorithm on a single IDAT file pair
```
mono $HOME/bin/autoconvert/AutoConvert.exe \
  $idat_green_file \
  $path_to_output_folder \
  $bpm_manifest_file \
  $egt_cluster_file
```
Make sure that the red IDAT file is in the same folder as the green IDAT file. Alternatively you can run on multiple IDAT file pairs
```
mono $HOME/bin/autoconvert/AutoConvert.exe \
  $path_to_idat_folder \
  $path_to_output_folder \
  $bpm_manifest_file \
  $egt_cluster_file
```
Make sure that the IDAT files have the same name prefix as the IDAT folder name. The software might require up to 8GB of RAM to run. Illumina provides manifest (BPM) and cluster (EGT) files for their arrays <a href="https://support.illumina.com/array/downloads.html">here</a>. Notice that if you provide the wrong BPM file, you will get an error such as: `Normalization failed!  Unable to normalize!` and if you provide the wrong EGT file, you will get an error such as `System.Exception: Unrecoverable Error...Exiting! Unable to find manifest entry ######## in the cluster file!`

Some users have encountered an issue when running Mono going along with the following error:
```
System.TypeInitializationException: The type initializer for 'System.Drawing.KnownColors' threw an exception. ---> System.TypeInitializationException: The type initializer for 'System.Drawing.GDIPlus' threw an exception. ---> System.DllNotFoundException: libgdiplus.so.0
```
The problem is related to the fact that you or your system administrator did not install the GDIPlus library. If this is the case, you can manually download an old binary version of the library together with some of its dependencis that should be compatible with your system using the following hack:
```
mkdir -p lib
wget http://old-releases.ubuntu.com/ubuntu/pool/main/libg/libgdiplus/libgdiplus_2.10-2_amd64.deb
ar x libgdiplus_2.10-2_amd64.deb data.tar.gz
tar xzf data.tar.gz -C lib ./usr/lib/libgdiplus.so.0.0.0 --strip-components=3
ln -s libgdiplus.so.0.0.0 lib/libgdiplus.so.0
wget http://old-releases.ubuntu.com/ubuntu/pool/main/t/tiff/libtiff4_3.9.5-1ubuntu1_amd64.deb
ar x libtiff4_3.9.5-1ubuntu1_amd64.deb data.tar.gz
tar xzf data.tar.gz -C lib ./usr/lib/x86_64-linux-gnu/libtiff.so.4.3.4 --strip-components=4
ln -s libtiff.so.4.3.4 lib/libtiff.so.4
wget http://old-releases.ubuntu.com/ubuntu/pool/main/libe/libexif/libexif12_0.6.20-1_amd64.deb
ar x libexif12_0.6.20-1_amd64.deb data.tar.gz
tar xzf data.tar.gz -C lib ./usr/lib/libexif.so.12.3.2 --strip-components=3
ln -s libexif.so.12.3.2 lib/libexif.so.12
wget http://old-releases.ubuntu.com/ubuntu/pool/main/libj/libjpeg8/libjpeg8_8b-1_amd64.deb
ar x libjpeg8_8b-1_amd64.deb data.tar.gz
tar xzf data.tar.gz -C lib ./usr/lib/libjpeg.so.8.0.2 --strip-components=3
ln -s libjpeg.so.8.0.2 lib/libjpeg.so.8
wget http://old-releases.ubuntu.com/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.46-3ubuntu1_amd64.deb
ar x libpng12-0_1.2.46-3ubuntu1_amd64.deb data.tar.gz
tar xzf data.tar.gz -C lib ./lib/x86_64-linux-gnu/libpng12.so.0.46.0 --strip-components=3
ln -s libpng12.so.0.46.0 lib/libpng12.so.0
/bin/rm lib{gdiplus_2.10-2,tiff4_3.9.5-1ubuntu1,exif12_0.6.20-1,jpeg8_8b-1,png12-0_1.2.46-3ubuntu1}_amd64.deb data.tar.gz
```
After downloading the binaries, replace `mono` with `LD_LIBRARY_PATH="lib" mono` when running AutoConvert through Mono

Convert Illumina GTC files to VCF
=================================

Specifications for Illumina BPM, EGT, and GTC files were obtained through Illumina's <a href="https://github.com/Illumina/BeadArrayFiles">BeadArrayFiles</a> library and <a href="https://github.com/Illumina/GTCtoVCF">GTCtoVCF</a> script. Specifications for IDAT files were obtained through Henrik Bengtsson's <a href="https://github.com/HenrikBengtsson/illuminaio">illuminaio</a> package.
```
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
bpm_manifest_file="..."
csv_manifest_file="..."
egt_cluster_file="..."
gtc_list_file="..."
out_prefix="..."
bcftools +gtc2vcf \
  --no-version -Ou \
  -b $bpm_manifest_file \
  -c $csv_manifest_file \
  -e $egt_cluster_file \
  -g $gtc_list_file \
  -f $ref \
  -x $out_prefix.sex | \
  bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  bcftools norm --no-version -Ob -o $out_prefix.bcf -c x -f $ref && \
  bcftools index -f $out_prefix.bcf
```
Notice that the gtc2vcf bcftools plugin will drop unlocalized variants. The final VCF might contain duplicates. If this is an issue `bcftools norm -d` can be used to remove such variants. At least one of the BPM or the CSV manifest files has to be provided. Without the BPM manifest file normalized intensities cannot be computed. Without the CSV manifest file indel alleles cannot be inferred and will be dropped. If the EGT cluster file is provided, information about genotype cluster centers are included in the VCF.

Convert Affymetrix CEL files to genotype calls
==============================================

Affymetrix provides a best practice workflow for genotyping data generated using <a href="https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-snp6-on-axiom.html">SNP6</a> and <a href="https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-Axiom-probeset-genotype.html">Axiom</a> arrays. As an examples, the following command will run the genotyping for the Affymetrix SNP6 array:
```
path_to_output_folder="..."
cel_list_file="..."
apt-probeset-genotype \
  --out-dir $path_to_output_folder \
  --read-models-brlmmp GenomeWideSNP_6.generic_prior.txt \
  --analysis-files-path . \
  --xml-file GenomeWideSNP_6.apt-probeset-genotype.AxiomGT1.xml \
  --cel-files $cel_list_file \
  --summaries \
  --write-models
```
Affymetrix provides Library and NetAffx Annotation files for their arrays <a href="http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays">here</a>

As an example, the following commands will obtain the files necessary to run the genotyping for the Affymetrix SNP6 array:
```
wget http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/SNP6_supplemental_axiom_analysis_files.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
unzip -oj genomewidesnp6_libraryfile.zip CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.{cdf,chrXprobes,chrYprobes}
unzip -o SNP6_supplemental_axiom_analysis_files.zip GenomeWideSNP_6.{generic_prior.txt,apt-probeset-genotype.AxiomGT1.xml,AxiomGT1.sketch}
unzip -o GenomeWideSNP_6.na35.annot.csv.zip GenomeWideSNP_6.na35.annot.csv
```

Note: If the program exits due to different chip types or probe counts with error message such as `Wrong CEL ChipType: expecting: 'GenomeWideSNP_6' and #######.CEL is: 'GenomeWideEx_6'` then add the option `--chip-type GenomeWideEx_6 --chip-type GenomeWideSNP_6` or `--force` to the command line to solve the problem

Convert Affymetrix genotype calls and intensities to VCF
========================================================

The affy2vcf bcftools plugin can be used to convert Affymetrix genotype calls and intensity files to VCF
```
ref="$HOME/res/human_g1k_v37.fasta" # or ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
annot_file="..." # for example annot_file="GenomeWideSNP_6.na35.annot.csv"
path_to_output_folder="..."
out_prefix="..."
bcftools +affy2vcf \
  --no-version -Ou \
  --fasta-ref $ref \
  --annot $annot_file \
  --sex $out_prefix.sex \
  --snp-posteriors $path_to_output_folder/AxiomGT1.snp-posteriors.txt \
  --summary $path_to_output_folder/AxiomGT1.summary.txt \
  --report $path_to_output_folder/AxiomGT1.report.txt \
  --calls $path_to_output_folder/AxiomGT1.calls.txt \
  --confidences $path_to_output_folder/AxiomGT1.confidences.txt | \
  bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  bcftools norm --no-version -Ob -o $out_prefix.bcf -c x -f $ref && \
  bcftools index -f $out_prefix.bcf
```

The final VCF might contain duplicates. If this is an issue `bcftools norm -d` can be used to remove such variants

Plot variants
=============

Install basic tools (Debian/Ubuntu specific if you have admin privileges):
```
sudo apt install r-cran-ggplot2 r-cran-data.table r-cran-gridextra
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
  --chrom 19 \
  --pos 45411941 \
  --png rs429358.png
```

Acknowledgements
================

This work is supported by NIH grant <a href="https://projectreporter.nih.gov/project_info_description.cfm?aid=8852155">R01 HG006855</a> and the Stanley Center for Psychiatric Research and by US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244)
