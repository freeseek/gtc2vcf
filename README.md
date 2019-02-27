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
    -l, --list-tags                    list available tags with description for VCF output
    -t, --tags LIST                    list of output tags [IGC,BAF,LRR]
    -i  --idat <file>                  IDAT file
    -b  --bpm <file>                   BPM manifest file
    -c  --csv <file>                   CSV manifest file
    -e  --egt <file>                   EGT cluster file
    -f, --fasta-ref <file>             reference sequence in fasta format
    -g, --gtc-list <file>              read GTC file names from file
        --adjust-clusters              adjust cluster centers in (Theta, R) space
    -x, --sex <file>                   output GenCall gender estimate into file
        --do-not-check-bpm             do not check whether BPM and GTC files match manifest file name
        --genome-studio                input a genome studio final report file (in matrix format)
        --no-version                   do not append version and command line to the header
    -o, --output <file>                write output to a file [standard output]
    -O, --output-type b|u|z|v|g        b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF, g GenomeStudio [v]
        --threads <int>                number of extra output compression threads [0]
```

Affymetrix tool:
```
Usage: bcftools +affy2vcf [options] --fasta-ref <fasta> --annot <file> --summary <file>
                            --snp-posteriors <file> --confidences <file> --calls <file>

Plugin options:
    -f, --fasta-ref <file>                     reference sequence in fasta format
        --annot <file>                         probeset annotation file
        --summary <summary.txt>                apt-probeset-genotype summary output
        --snp-posteriors <snp-posteriors.txt>  apt-probeset-genotype snp-posteriors output
        --report <report.txt>                  apt-probeset-genotype report output
        --confidences <confidences.txt>        apt-probeset-genotype confidences output
        --calls <calls.txt>                    apt-probeset-genotype calls output
        --adjust-clusters                      adjust cluster centers in (Contrast, Size) space
    -x, --sex <file>                           output apt-probeset-genotype gender estimate into file
        --no-version                           do not append version and command line to the header
    -o, --output <file>                        write output to a file [standard output]
    -O, --output-type b|u|z|v                  b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF
        --threads <int>                        number of extra output compression threads [0]
```

Installation
============

Install basic tools (Debian/Ubuntu specific)
```
sudo apt install wget gzip unzip samtools wine64 mono-devel libgdiplus
```

Optionally, you can install these libraries to activate further bcftools features:
```
sudo apt install liblzma-dev libbz2-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res && cd /tmp
```

Download latest version of `htslib` and `bcftools` (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Add patch (to allow the fixref plugin to flip BAF values) and code for plugins
```
/bin/rm -f bcftools/plugins/{gtc2vcf.c,affy2vcf.c,fixref.patch}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/gtc2vcf/master/{gtc2vcf.c,affy2vcf.c,fixref.patch}
cd bcftools/plugins && patch < fixref.patch && cd ../..
```
If for any reason the patch fails with an error message, contact the <a href="mailto:giulio.genovese@gmail.com">author</a> for a fix.

Compile latest version of `htslib` (optionally disable `bz2` and `lzma`) and `bcftools` (make sure you are using gcc version 5 or newer)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{fixref,gtc2vcf,affy2vcf}.so} $HOME/bin/
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

Illumina provides the <a href="https://support.illumina.com/array/array_software/beeline.html">Beeline</a> software for free and this includes the AutoConvert executable which allows to call genotypes from raw intensity data using Illumina's proprietary GenCall algorithm. AutoConvert is almost entirely written in Mono/.Net language, with the exception of one small mathmatical function (findClosestSitesToPointsAlongAxis) which is contained instead within a Windows PE32+ library (MathRoutines.dll). As this is <a href="http://www.mono-project.com/docs/advanced/embedding/">unmanaged code</a>, to be run on Linux with <a href="https://www.mono-project.com/">Mono</a> it needs to be embedded in an equivalent Linux ELF64 library (libMathRoutines.dll.so) as shown below. This function is run as part of the <a href="http://doi.org/10.1093/bioinformatics/btm443">normalization</a> of the raw intensities when sampling 400 <a href="http://patft.uspto.gov/netacgi/nph-Parser?patentnumber=7035740">candidate homozygotes</a> before calling genotypes. For some unclear reasons, you will also need to separately download an additional Mono/.Net library (Heatmap.dll) from <a href="https://support.illumina.com/array/array_software/genomestudio.html">GenomeStudio</a> and include it in your binary directory as shown below, most likely due to differences in which Mono and .Net resolve library dependencies.
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
cp -R $HOME/.wine/drive_c/Program\ Files/Illumina/AutoConvert\ 2.0 $HOME/bin/autoconvert
cp $HOME/.wine/drive_c/Program\ Files/Illumina/GenomeStudio\ 2.0/Heatmap.dll $HOME/bin/autoconvert/
wget https://raw.githubusercontent.com/freeseek/gtc2vcf/master/nearest_neighbor.c
gcc -fPIC -shared -O2 -o $HOME/bin/autoconvert/libMathRoutines.dll.so nearest_neighbor.c
```
If you do not succeed at installing AutoConvert and GenomeStudio using wine following the instructions above, you can always resort to install them using a Windows machine, copy the files on Linux, and then add the required Linux ELF64 library as described in the last steps above. Notice that this approach to run AutoConvert on Linux is <strong>not</strong> supported by Illumina

Affymetrix provides the <a href="http://www.affymetrix.com/support/developer/powertools/changelog/index.html">Analysis Power Tools (APT)</a> for free which allow to call genotypes from raw intensity data using an algorithm derived from <a href="http://tools.thermofisher.com/content/sfs/brochures/brlmmp_whitepaper.pdf">BRLMM-P</a>.
```
mkdir -p $HOME/bin && cd /tmp
wget https://downloads.thermofisher.com/Affymetrix_Softwares/APT_2.10.0/apt-2.10.0-x86_64-intel-linux.zip
unzip -ojd $HOME/bin apt-2.10.0-x86_64-intel-linux.zip apt-2.10.0-x86_64-intel-linux/bin/apt-probeset-genotype
chmod a+x $HOME/bin/apt-probeset-genotype
```

Convert Illumina IDAT files to GTC files
========================================

Once Mono and AutoConvert are properly installed on your system, run Illumina's proprietary GenCall algorithm on a single IDAT file pair
```
mono $HOME/bin/autoconvert/AutoConvert.exe $idat_green_file $path_to_output_folder $manifest_file $egt_file
```
Make sure that the red IDAT file is in the same folder as the green IDAT file. Alternatively you can run on multiple IDAT file pairs
```
mono $HOME/bin/autoconvert/AutoConvert.exe $path_to_idat_folder $path_to_output_folder $manifest_file $egt_file
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
After downloading the binaries, replace `mono` with `LD_LIBRARY_PATH="lib" mono` when running AutoConvert through Mono.

Convert Illumina GTC files to VCF
=================================

Specifications for Illumina BPM, EGT, and GTC files were obtained through Illumina's <a href="https://github.com/Illumina/BeadArrayFiles">BeadArrayFiles</a> library and <a href="https://github.com/Illumina/GTCtoVCF">GTCtoVCF</a> script. Specifications for IDAT files were obtained through Henrik Bengtsson's <a href="https://github.com/HenrikBengtsson/illuminaio">illuminaio</a> package. Reference strand determination is performed using Illumina's <a href="https://www.illumina.com/documents/products/technotes/technote_topbot.pdf">TOP/BOT</a> strand assignment in the manifest file. The gtc2vcf bcftools plugin is hundreds of times faster than Illumina's script and can be used to convert GTC files to VCF
```
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
manifest_file="..."
egt_file="..."
gtc_list="..."
out="..."
$HOME/bin/bcftools +$HOME/bin/gtc2vcf.so --no-version -Ou -b $manifest_file -e $egt_file -g $gtc_list -f $ref -x $out.sex | \
  $HOME/bin/bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  $HOME/bin/bcftools +$HOME/bin/fixref.so --no-version -Ou -- -f $ref -m top --flip-baf | \
  $HOME/bin/bcftools norm --no-version -Ob -o $out.bcf -c x -f $ref && \
  $HOME/bin/bcftools index -f $out.bcf
```
Notice that this will drop unlocalized variants and indels

Convert Illumina GenomeStudio final report to VCF
=================================================

Alternatively, if a GenomeStudio final report in <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2-0/genomestudio-genotyping-module-v2-user-guide-11319113-01.pdf#page=67">matrix format</a> is provided instead, and it follows the following convention (Illumina does not share specifications for this file format):
```
Chromosome	Position	IlmnStrand	SNP	Name	SM.GType	SM.Score	SM.Theta	SM.R	SM.B Allele Freq	SM.Log R Ratio
9	139906359	BOT	[T/C]	200003	AA	0.9299	0.029	1.300	0.0027	0.2150
9	139926402	TOP	[A/G]	200006	AB	0.7877	0.435	2.675	0.4742	0.3024
2	220084902	BOT	[T/C]	200047	AA	0.8612	0.083	0.476	0.0532	-0.1173
2	220089685	TOP	[C/G]	200050	BB	0.8331	0.995	1.499	1.0000	0.3068
```

It can be converted into VCF format with the following command
```
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
genome_studio_file="..."
out="..."
$HOME/bin/bcftools +$HOME/bin/gtc2vcf.so --no-version -Ou --genome-studio $genome_studio_file -f $ref | \
  $HOME/bin/bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  $HOME/bin/bcftools +$HOME/bin/fixref.so --no-version -Ou -e 'REF="N" || ALT="N"' -- -f $ref -m top --flip-baf | \
  $HOME/bin/bcftools norm --no-version -Ob -o $out.bcf -c x -f $ref && \
  $HOME/bin/bcftools index -f $out.bcf
```
Notice that this will drop unlocalized variants, monomorphic variants, and indels

Convert Affymetrix CEL files to genotype calls
==============================================

Affymetrix provides a best practice workflow for genotyping data generated using <a href="https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-snp6-on-axiom.html">SNP6</a> and <a href="https://www.affymetrix.com/support/developer/powertools/changelog/VIGNETTE-Axiom-probeset-genotype.html">Axiom</a> arrays. As an examples, the following command will run the genotyping for the Affymetrix SNP6 array:
```
dir="..."
cel_list="..."
$HOME/bin/apt-probeset-genotype \
  --out-dir $dir \
  --read-models-brlmmp GenomeWideSNP_6.generic_prior.txt \
  --analysis-files-path . \
  --xml-file GenomeWideSNP_6.apt-probeset-genotype.AxiomGT1.xml \
  --cel-files $cel_list \
  --summaries \
  --write-models
```
Affymetrix provides Library and NetAffx Annotation files for their arrays <a href="http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays">here</a>.

As an example, the following commands will obtain the files necessary to run the genotyping for the Affymetrix SNP6 array:
```
wget http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/SNP6_supplemental_axiom_analysis_files.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
unzip -oj genomewidesnp6_libraryfile.zip CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.{cdf,chrXprobes,chrYprobes}
unzip -o SNP6_supplemental_axiom_analysis_files.zip GenomeWideSNP_6.{generic_prior.txt,apt-probeset-genotype.AxiomGT1.xml,AxiomGT1.sketch}
unzip -o GenomeWideSNP_6.na35.annot.csv.zip GenomeWideSNP_6.na35.annot.csv
```

Note: If the program exits due to different chip types or probe counts with error message such as `Wrong CEL ChipType: expecting: 'GenomeWideSNP_6' and #######.CEL is: 'GenomeWideEx_6'` then add the option `--chip-type GenomeWideEx_6 --chip-type GenomeWideSNP_6` or `--force` to the command line to solve the problem.

Convert Affymetrix genotype calls and intensities to VCF
========================================================

The affy2vcf bcftools plugin can be used to convert Affymetrix genotype calls and intensity files to VCF
```
ref="$HOME/res/human_g1k_v37.fasta" # or ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
dir="..."
out="..."
$HOME/bin/bcftools +$HOME/bin/affy2vcf.so --no-version -Ou --fasta-ref $ref --annot $annot_file --sex $out.sex \
  --snp-posteriors $dir/AxiomGT1.snp-posteriors.txt \
  --summary $dir/AxiomGT1.summary.txt \
  --report $dir/AxiomGT1.report.txt \
  --calls $dir/AxiomGT1.calls.txt \
  --confidences $dir/AxiomGT1.confidences.txt | \
  $HOME/bin/bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  $HOME/bin/bcftools +$HOME/bin/fixref.so --no-version -Ou -e 'REF="N" || ALT="N"' -- -f $ref -m swap --flip-baf | \
  $HOME/bin/bcftools norm --no-version -Ob -o $out.bcf -c x -f $ref && \
  $HOME/bin/bcftools index -f $out.bcf
```

Acknowledgements
================

This work is supported by NIH grant <a href="https://projectreporter.nih.gov/project_info_description.cfm?aid=8852155">R01 HG006855</a> and the Stanley Center for Psychiatric Research and by US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244)
