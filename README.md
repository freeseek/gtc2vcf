gtc2vcf
=======

A set of tools to convert Illumina and Affymetrix DNA microarray intensity data files into VCF files <b>without</b> using Microsoft Windows. You can use the final output to run the pipeline to detect <a href="https://github.com/freeseek/mocha">mosaic chromosomal alterations</a>. If you use this tool in your publication, please cite this website. For any feedback or questions, contact the <a href="mailto:giulio.genovese@gmail.com">author</a>

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
   * [Acknowledgements](#acknowledgements)
<!--te-->

Usage
=====

Illumina tool:
```
Usage: bcftools +gtc2vcf [options] [<A.gtc> ...]

Plugin options:
    -l, --list-tags                 list available FORMAT tags with description for VCF output
    -t, --tags LIST                 list of output FORMAT tags [GT,GQ,IGC,BAF,LRR,NORMX,NORMY,R,THETA,X,Y]
    -b, --bpm <file>                BPM manifest file
    -c, --csv <file>                CSV manifest file
    -e, --egt <file>                EGT cluster file
    -f, --fasta-ref <file>          reference sequence in fasta format
        --set-cache-size <int>      select fasta cache size in bytes
        --gc-window-size <int>      window size in bp used to compute the GC content (-1 for no estimate) [200]
    -g, --gtcs <dir|file>           GTC genotype files from directory or list from file
    -i, --idat                      input IDAT files rather than GTC files
        --capacity <int>            number of variants to read from intensity files per I/O operation [32768]
        --adjust-clusters           adjust cluster centers in (Theta, R) space (requires --bpm and --egt)
        --use-gtc-sample-names      use sample name in GTC files rather than GTC file name
        --do-not-check-bpm          do not check whether BPM and GTC files match manifest file name
        --genome-studio <file>      input a GenomeStudio final report file (in matrix format)
        --no-version                do not append version and command line to the header
    -o, --output <file>             write output to a file [standard output]
    -O, --output-type <b|u|z|v|t>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF
                                    v: uncompressed VCF, t: GenomeStudio tab-delimited text output [v]
        --threads <int>             number of extra output compression threads [0]
    -x, --extra <file>              write GTC metadata to a file
    -v, --verbose                   print verbose information

Manifest options:
        --beadset-order             output BeadSetID normalization order (requires --bpm and --csv)
        --fasta-flank               output flank sequence in FASTA format (requires --csv)
    -s, --sam-flank <file>          input flank sequence alignment in SAM/BAM format (requires --csv)
        --genome-build <assembly>   genome build ID used to update the manifest file [GRCh38]

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
    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv --fasta-source-seq -o GSA-24v3-0_A1.fasta
    bwa mem -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GSA-24v3-0_A1.fasta -o GSA-24v3-0_A1.sam
    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv --sam-source-seq GSA-24v3-0_A1.sam -o GSA-24v3-0_A1.GRCh38.csv
```

Affymetrix tool:
```
Usage: bcftools +affy2vcf [options] --csv <file> --fasta-ref <file> [<A.chp> ...]

Plugin options:
    -l, --list-tags               list available FORMAT tags with description for VCF output
    -t, --tags LIST               list of output FORMAT tags [GT,CONF,BAF,LRR,NORMX,NORMY,DELTA,SIZE]
    -c, --csv <file>              CSV manifest file
    -f, --fasta-ref <file>        reference sequence in fasta format
        --set-cache-size <int>    select fasta cache size in bytes
        --gc-window-size <int>    window size in bp used to compute the GC content (-1 for no estimate) [200]
        --calls <file>            apt-probeset-genotype calls output
        --confidences <file>      apt-probeset-genotype confidences output
        --summary <file>          apt-probeset-genotype summary output
        --snp <file>              apt-probeset-genotype SNP posteriors output
        --chps <dir|file>         input CHP files rather than tab delimited files
        --cel <file>              input CEL files rather CHP files
        --adjust-clusters         adjust cluster centers in (Contrast, Size) space (requires --snp)
        --no-version              do not append version and command line to the header
    -o, --output <file>           write output to a file [standard output]
    -O, --output-type <b|u|z|v>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
        --threads <int>           number of extra output compression threads [0]
    -x, --extra <file>            write CHP metadata to a file (requires CHP files)
    -v, --verbose                 print verbose information

Manifest options:
        --fasta-flank             output flank sequence in FASTA format (requires --csv)
    -s, --sam-flank <file>        input source sequence alignment in SAM/BAM format (requires --csv)

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
    bcftools +affy2vcf -c GenomeWideSNP_6.na35.annot.csv --fasta-flank -o GenomeWideSNP_6.fasta
    bwa mem -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GenomeWideSNP_6.fasta -o GenomeWideSNP_6.sam
    bcftools +affy2vcf -c GenomeWideSNP_6.na35.annot.csv -s GenomeWideSNP_6.sam -o GenomeWideSNP_6.na35.annot.GRCh38.csv
```

Installation
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges)
```
sudo apt install wget git g++ zlib1g-dev bwa unzip samtools msitools cabextract mono-devel libgdiplus libicu66 bcftools
```

Optionally, you can install these libraries to activate further HTSlib features:
```
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res && cd /tmp
```

We recommend compiling the source code but, wherever this is not possible, Linux x86_64 pre-compiled binaries are available for download <a href="http://software.broadinstitute.org/software/gtc2vcf">here</a>. However, notice that you will require a copy of BCFtools 1.10 or newer (available with Ubuntu 20.04)

Download latest version of <a href="https://github.com/samtools/htslib">HTSlib</a> and <a href="https://github.com/samtools/bcftools">BCFtools</a> (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Download plugins code
```
/bin/rm -f bcftools/plugins/{gtc2vcf.{c,h},affy2vcf.c}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/gtc2vcf/master/{gtc2vcf.{c,h},affy2vcf.c}
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
bwa index $HOME/res/human_g1k_v37.fasta
```

Install the GRCh38 human genome reference (following the suggestion from <a href="http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use">Heng Li</a>)
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bwa index $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Software Installation
=====================

Illumina provides the <a href="https://support.illumina.com/array/array_software/illumina-array-analysis-platform.html">Illumina Array Analysis Platform</a> software for free and this includes the iaap-cli command line executable which runs natively on Linux
```
mkdir -p $HOME/bin && cd /tmp
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz
tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz -C $HOME/bin/ iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1
```
However, notice that in some older Linux machines this approach does not work and at the time of this writing iaap-cli is unable to read old BPM manifest files yielding error `Unknown Manifest version`, while the AutoConvert command line tool does not have this limitation

Illumina also provides the <a href="https://support.illumina.com/array/array_software/beeline.html">Beeline</a> software for free and this includes the AutoConvert.exe command line executable which allows to call genotypes from raw intensity data using Illumina's proprietary GenCall algorithm. AutoConvert is almost entirely written in Mono/.Net language, with the exception of one small mathmatical function (findClosestSitesToPointsAlongAxis) which is contained instead within a Windows PE32+ library (MathRoutines.dll). As this is <a href="http://www.mono-project.com/docs/advanced/embedding/">unmanaged code</a>, to be run on Linux with <a href="https://www.mono-project.com/">Mono</a> it needs to be embedded in an equivalent Linux ELF64 library (libMathRoutines.dll.so) as shown below. This function is run as part of the <a href="http://doi.org/10.1093/bioinformatics/btm443">normalization</a> of the raw intensities when sampling <a href="https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina_gt_normalization.pdf">400 candidate homozygotes</a> before calling genotypes. For some unclear reasons, you will also need to separately download an additional Mono/.Net library (Heatmap.dll) from <a href="https://support.illumina.com/array/array_software/genomestudio.html">GenomeStudio</a> and include it in your binary directory as shown below, most likely due to differences in which Mono and .Net resolve library dependencies
```
mkdir -p $HOME/bin && cd /tmp
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/beeline/autoconvert-software-v2-0-1-installer.zip
unzip -o autoconvert-software-v2-0-1-installer.zip
msiextract AutoConvertInstaller.msi
cp -R Illumina/AutoConvert\ 2.0 $HOME/bin/autoconvert

wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/genomestudio/genomestudio-software-v2-0-4-5-installer.zip
unzip -oj genomestudio-software-v2-0-4-5-installer.zip
cabextract GenomeStudioInstaller.exe
msiextract a0
cp Illumina/GenomeStudio\ 2.0/Heatmap.dll $HOME/bin/autoconvert/

wget https://raw.githubusercontent.com/freeseek/gtc2vcf/master/nearest_neighbor.c
gcc -fPIC -shared -O2 -o $HOME/bin/autoconvert/libMathRoutines.dll.so nearest_neighbor.c
```
If you fail to download the autoconvert software, contact the <a href="mailto:giulio.genovese@gmail.com">author</a> for troubleshooting. Notice that this approach to run AutoConvert on Linux is <strong>not</strong> supported by Illumina

Affymetrix provides the <a href="https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html">Analysis Power Tools (APT)</a> for free which allow to call genotypes from raw intensity data using an algorithm derived from <a href="http://tools.thermofisher.com/content/sfs/brochures/brlmmp_whitepaper.pdf">BRLMM-P</a>
```
mkdir -p $HOME/bin && cd /tmp
wget https://downloads.thermofisher.com/APT/2.11.3/apt_2.11.3_linux_64_bit_x86_binaries.zip
unzip -ojd $HOME/bin apt_2.11.3_linux_64_bit_x86_binaries.zip apt_2.11.3_linux_64_bit_x86_binaries/bin/apt-probeset-genotype
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
This will generate a spreadsheet table with information about each IDAT file including a guess for what manifest and cluster files you should use. If a guess is not provided, contact the <a href="mailto:giulio.genovese@gmail.com">author</a> for troubleshooting

Similarly, you can use the affy2vcf bcftools plugin to extract chip type information from CEL files
```
path_to_cel_folder="..."
bcftools +affy2vcf \
  --cel --chps $path_to_cel_folder
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

Specifications for Illumina BPM, EGT, and GTC files were obtained through Illumina's <a href="https://github.com/Illumina/BeadArrayFiles">BeadArrayFiles</a> library and <a href="https://github.com/Illumina/GTCtoVCF">GTCtoVCF</a> script. Specifications for IDAT files were obtained through Henrik Bengtsson's <a href="https://github.com/HenrikBengtsson/illuminaio">illuminaio</a> package
```
bpm_manifest_file="..."
csv_manifest_file="..."
egt_cluster_file="..."
path_to_gtc_folder="..."
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
out_prefix="..."
bcftools +gtc2vcf \
  --no-version -Ou \
  --bpm $bpm_manifest_file \
  --csv $csv_manifest_file \
  --egt $egt_cluster_file \
  --gtcs $path_to_gtc_folder \
  --fasta-ref $ref \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  bcftools norm --no-version -Ob -o $out_prefix.bcf -c x -f $ref && \
  bcftools index -f $out_prefix.bcf
```
Heavy random access to the reference will be needed, so it is important that enough extra memory be available for the operating system to cache the reference or else the task can run excruciatingly slowly. Notice that the gtc2vcf bcftools plugin will drop unlocalized variants. The final VCF might contain duplicates. If this is an issue `bcftools norm -d` can be used to remove such variants. At least one of the BPM or the CSV manifest files has to be provided. Normalized intensities cannot be computed without the BPM manifest file. Indel alleles cannot be inferred and will be skipped without the CSV manifest file. Information about genotype cluster centers will be included in the VCF if the EGT cluster file is provided. You can use gtc2vcf to convert one GTC file at a time, but we strongly advise to convert multiple files at once as single sample VCF files will consume a lot of storage space. If you convert hundreds of GTC files at once, you can use the `--adjust-clusters` option which will recenter the genotype clusters rather than using those provided in the EGT cluster file and will compute less noisy LRR values. If you use the `--adjust-clusters` option and you are using the output for calling <a href="https://github.com/freeseek/mocha">mosaic chromosomal alterations</a>, then it is safe to turn the median BAF/LRR adjustments off during that step (i.e. use `--adjust-BAF-LRR -1`)

When running the conversion, the gtc2vcf plugin will double check that the SNP manifest metadata information in the GTC file matches the descriptor file name in the BPM file to make sure you are using the correct manifest file. Sometimes, due to discrepancies between the BPM file name provided by Illumina and the internal descriptor file name, this safety check fails. To turn off this feature in these cases, you can use `--do-not-check-bpm`
```

``` 
exhibit discrepancies with their internal descriptor file name

Convert Affymetrix CEL files to CHP files
=========================================

Affymetrix provides a best practice workflow for genotyping data generated using <a href="https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-snp6-on-axiom.html">SNP6</a> and <a href="https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-Axiom-probeset-genotype.html">Axiom</a> arrays. As an example, the following command will run the genotyping for the Affymetrix SNP6 array:
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
Affymetrix provides Library and NetAffx Annotation files for their arrays <a href="http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays">here</a> and <a href="http://media.affymetrix.com/analysis/downloads/lf/genotyping">here</a>

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
ref="$HOME/res/human_g1k_v37.fasta" # or ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
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
  bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  bcftools norm --no-version -Ob -o $out_prefix.bcf -c x -f $ref && \
  bcftools index -f $out_prefix.bcf
```
Heavy random access to the reference will be needed, so it is important that enough extra memory be available for the operating system to cache the reference or else the task can run excruciatingly slowly. The final VCF might contain duplicates. If this is an issue `bcftools norm -d` can be used to remove such variants. There is often no need to use the `--adjust-clusters` option for Affymetrix data as the cluster posteriors are already adjusted using the data processed by the genotype caller

Using an alternative genome reference
=====================================

Illumina provides <a href="https://support.illumina.com/bulletins/2017/04/infinium-human-genotyping-manifests-and-support-files--with-anno.html">GRCh38/hg38</a> manifests for many of its genotyping arrays. However, if your genotyping array is not supported for the newer reference by Illumina, you can use the `--fasta-flank` and `--sam-flank` options to realign the source sequences from the manifest files you have and recompute the marker positions. This approach uses <a href="https://support.illumina.com/bulletins/2016/05/infinium-genotyping-manifest-column-headings.html">source sequence</a> and <a href="https://support.illumina.com/bulletins/2017/06/how-to-interpret-dna-strand-and-allele-information-for-infinium-.html">strand</a> information to identify the marker <a href="https://support.illumina.com/bulletins/2016/06/-infinium-genotyping-array-manifest-files-what-does-chr-or-mapinfo---mean.html">coordinates</a>. It will need a sequence aligner such as `bwa` to realign the sequences and it seems to reproduce the coordinates provided from Illumina more than 99.9% of the times. Mapping information will follow the <a href="https://github.com/Illumina/GTCtoVCF#manifests">implicit dbSNP standard</a>. Occasionally the source sequence provided by Illumina is incorrect and it is impossible to recover the correct marker coordinate from the source sequence alone

You first have to generate an alignment file for the source sequences from a CSV manifest file
```
csv_manifest_file="..."
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
bam_alignment_file="..."
bcftools +gtc2vcf \
  -c $csv_manifest_file \
  --fasta-flank | \
  bwa mem -M $ref - | \
  samtools view -bS \
  -o $bam_alignment_file
```
Notice that you need to use the `-M` option to mark shorter split hits as secondary. Then you load the alignment file while converting your GTC files to VCF including the `-s $bam_alignment_file` option

Some older manifest files from Illumina have thousands of markers with incorrect RefStrand annotations that will lead to incorrect genotypes. While Illumina has not explained why this is the case, it still distributes incorrect manifests. If you are using one of the following manifests
```
Human1M-Duov3_H
Human610-Quadv1_H
Human660W-Quad_v1_H
HumanCytoSNP-12v2-1_Anova
HumanOmni1-Quad_v1-0-Multi_H
HumanOmni1-Quad_v1-0_H
```
We advise to either contact Illumina to demand a fixed version or to use gtc2vcf to realign the source sequences

Also, Illumina assigns chromosomal positions to indels by first left aligning the source sequences in an incoherent way (see <a href="https://github.com/Illumina/GTCtoVCF/blob/develop/BPMRecord.py">here</a>). Apparently this is incoherent enough that Illumina also cannot get the coordinates of homopolymer indels right. For example, chromosome 13 ClinVar indel <a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/37959">rs80359507</a> is assigned to position 32913838 in the manifest file for the GSA-24v2-0 array, but it is assigned to position 32913837 in the manifest file for GSA-24v3-0 array (GRCh37 coordinates). If you want to trust genotypes at homopolymer indels, we advise to use gtc2vcf to realign the source sequences

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
  --png	rs800292.png
```

![](rs800292.png)

Acknowledgements
================

This work is supported by NIH grant <a href="http://grantome.com/grant/NIH/R01-HG006855">R01 HG006855</a>, NIH grant <a href="http://grantome.com/grant/NIH/R01-MH104964">R01 MH104964</a>, US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244), and the Stanley Center for Psychiatric Research
