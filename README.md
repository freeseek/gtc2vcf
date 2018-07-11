gtc2vcf
=======

A set of tools to convert Illumina files containing intensity data into VCF files <b>without</b> using Microsoft Windows. You can use the final output to run the pipeline to detect <a href="https://github.com/freeseek/mocha">mosaic chromosomal alterations</a>. If you use this tool in your publication, please cite this website. For any feedback, send an email to giulio.genovese@gmail.com

![](gtc2vcf.png)

Usage
=====

```
Usage: bcftools +gtc2vcf [options] <A.gtc> [...]

Plugin options:
    -l, --list-tags                    list available tags with description for VCF output
    -t, --tags LIST                    list of output tags [IGC,BAF,LRR]
    -i  --idat <file>                  IDAT file
    -b  --bpm <file>                   BPM manifest file
    -e  --egt <file>                   EGT cluster file
    -f, --fasta-ref <file>             reference sequence in fasta format
    -g, --gtc-list <file>              read GTC file names from file
        --do-not-check-bpm             do not check whether BPM and GTC files match manifest file name
        --genome-studio                input a genome studio file
        --no-version                   do not append version and command line to the header
    -o, --output <file>                write output to a file [standard output]
    -O, --output-type b|u|z|v|g        b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF, g GenomeStudio [v]
        --threads <int>                number of extra output compression threads [0]
```

Installation
============

Install basic tools (Debian/Ubuntu specific)
```
sudo apt-get install wget liblzma-dev libbz2-dev libgsl0-dev gzip samtools unzip wine64 mono-devel
```

Preparation steps
```
mkdir -p $HOME/bin && cd /tmp
```

Download latest version of `htslib` and `bcftools` (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Add patches and code for plugin
```
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/gtc2vcf/master/{gtc2vcf.c,fixref.patch}
cd bcftools/plugins && patch < fixref.patch && cd ../..
```

Compile latest version of `htslib` (optionally disable `bz2` and `lzma`) and `bcftools`
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{fixref,gtc2vcf}.so} $HOME/bin/
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

Illumina provides the <a href="https://support.illumina.com/array/array_software/beeline.html">Beeline</a> software for free and this includes the AutoConvert executable which allows to call genotypes from raw intensity data using Illumina's proprietary GenCall algorithm. AutoConvert is almost entirely written in Mono/.Net language, with the exception of one small mathmatical function (findClosestSitesToPointsAlongAxis) which is contained instead within a Windows PE32+ library (MathRoutines.dll). As this is <a href="http://www.mono-project.com/docs/advanced/embedding/">unmanaged code</a>, to be run on Linux with Mono it needs to be embedded in an equivalent Linux ELF64 library (libMathRoutines.dll.so) as shown below. This function is run as part of the <a href="http://doi.org/10.1093/bioinformatics/btm443">normalization</a> of the raw intensities when sampling 400 <a href="http://patft.uspto.gov/netacgi/nph-Parser?patentnumber=7035740">candidate homozygotes</a> before calling genotypes. For some unclear reasons, you will also need to separately download an additional Mono/.Net library (Heatmap.dll) from <a href="https://support.illumina.com/array/array_software/genomestudio.html">GenomeStudio</a> and include it in your binary directory as shown below, most likely due to differences in which Mono and .Net resolve library dependencies.
```
mkdir -p $HOME/bin && cd /tmp
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/beeline/autoconvert-software-v2-0-1-installer.zip
unzip -o autoconvert-software-v2-0-1-installer.zip
wine64 msiexec /i AutoConvertInstaller.msi WIXUI_DONTVALIDATEPATH="1"
wget http://download.microsoft.com/download/E/2/1/E21644B5-2DF2-47C2-91BD-63C560427900/NDP452-KB2901907-x86-x64-AllOS-ENU.exe
wine64 NDP452-KB2901907-x86-x64-AllOS-ENU.exe
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/genomestudio/genomestudio-software-v2-0-4-5-installer.zip
unzip -o genomestudio-software-v2-0-4-5-installer.zip
wine64 GenomeStudio-software-v2-0-4-5-installer/GenomeStudioInstaller.exe
wget https://raw.githubusercontent.com/freeseek/gtc2vcf/master/nearest_neighbor.c
cp -R $HOME/.wine/drive_c/Program\ Files/Illumina/AutoConvert\ 2.0 $HOME/bin/autoconvert
cp $HOME/.wine/drive_c/Program\ Files/Illumina/GenomeStudio\ 2.0/Heatmap.dll $HOME/bin/autoconvert/
gcc -fPIC -shared -O2 -o $HOME/bin/autoconvert/libMathRoutines.dll.so nearest_neighbor.c
```
Notice that this approach to run AutoConvert on Linux is <strong>not</strong> supported by Illumina

Convert Illumina IDAT files to GTC files
========================================

Once Mono and AutoConvert are properly installed on your system, run Illumina's proprietary GenCall algorithm
```
mono $HOME/bin/autoconvert/AutoConvert.exe path_to_idat_folder path_to_output_folder manifest_file egt_file
```
Make sure that the IDAT files have the same name prefix as the IDAT folder name. The software might require up to 8GB of RAM to run. Illumina provides manifest (BPM) and cluster (EGT) files for their arrays <a href="https://support.illumina.com/array/downloads.html">here</a>

Convert Illumina GTC files to VCF
=================================

Specifications for Illumina BPM, EGT, and GTC files were obtained through Illumina's <a href="https://github.com/Illumina/BeadArrayFiles">BeadArrayFiles</a> library and <a href="https://github.com/Illumina/GTCtoVCF">GTCtoVCF</a> script. Specifications for IDAT files were obtained through Henrik Bengtsson's <a href="https://github.com/HenrikBengtsson/illuminaio">illuminaio</a> package. Reference strand determination is performed using Illumina's <a href="https://www.illumina.com/documents/products/technotes/technote_topbot.pdf">TOP/BOT</a> strand assignment in the manifest file. The resulting bcftools plugin is hundreds of times faster than Illumina's script and can be used to convert GTC files to VCF
```
$HOME/bin/bcftools +$HOME/bin/gtc2vcf.so --no-version -Ou -b manifest_file -e egt_file -g $gtcs -f $ref | \
  $HOME/bin/bcftools sort -Ou -T . | \
  $HOME/bin/bcftools +$HOME/bin/fixref.so --no-version -Ou -- -f $ref -m top -b | \
  $HOME/bin/bcftools norm --no-version -Ob -o $out.bcf -c x -f $ref && \
  $HOME/bin/bcftools index -f $out.bcf
```
Notice that this will drop unlocalized variants and indels

Convert Illumina GenomeStudio final report to VCF
=================================================

Alternatively, if a GenomeStudio final report in <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2-0/genomestudio-genotyping-module-v2-user-guide-11319113-01.pdf#page=67">matrix format</a> is provided instead, this can be converted to VCF
```
$HOME/bin/bcftools  +$HOME/bin/gtc2vcf.so --no-version -Ou --genome-studio genome_studio_file -f $ref | \
  $HOME/bin/bcftools sort -Ou -T . | \
  $HOME/bin/bcftools +$HOME/bin/fixref.so --no-version -Ou -e 'REF="N" || ALT="N"' -- -f $ref -m top -b | \
  $HOME/bin/bcftools norm --no-version -Ob -o $out.bcf -c x -f $ref && \
  $HOME/bin/bcftools index -f $out.bcf
```
Notice that this will drop unlocalized variants, monomorphic variants, and indels
