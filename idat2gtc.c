/* The MIT License

   Copyright (c) 2024 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

// the code in this file reimplements functionalities and ideas present in:
// - AutoConvert (v1.6.3.1)
// - GTCtoVCF
// - BeadArrayFiles
// these resources were provided by Illumina without license restrictions

// the code in this file can be used as a replacement of the Illumina AutoCall software to convert IDAT intensity files
// into GTC genotype files for Infinium arrays which was implemented over time in different proprietary software:
// - AutoConvert (v1.6.3.1) - http://support.illumina.com/downloads/beeline_software_v10.html
// - AutoConvert 2.0 (v2.0.1.179) - http://support.illumina.com/downloads/beeline-software-2-0.html
// - IAAP CLI (v1.1) - http://support.illumina.com/downloads/iaap-genotyping-cli.html
// - Array Analysis CLI (v2.1) -
// http://support.illumina.com/downloads/illumina-microarray-analytics-array-analysis-cli-v2-installers.html

// the Illumina AutoCall software performs three main steps:
// - Normalization
// - Genotyping
// - Gender Estimation
// if AutoConvert and AutoConvert 2.0 are run without an input cluster file, only the normalization will be performed

// the normalization, clustering, and genotype calling functionalities of Illumina AutoCall were covered by the
// following patents:
// - http://patents.google.com/patent/US7035740 - covers normalization algorithm (2024-05-05)
// - http://patents.google.com/patent/US7467117 - divisional, covers clustering and genotyping (2024-03-24)
// - http://patents.google.com/patent/US20050216207 - same as US7035740
// - http://patents.google.com/patent/US20060224529 - same as US7467117

// GenCall GenTrain 2.0 uses the following algorithms:
// - Normalization algorithm (version 1.1.2)
// - Clustering algorithm (version 6.3.1)
// - Genotyping algorithm (version 6.3.0)
// GenCall GenTrain 3.0 uses the following algorithms:
// - Normalization algorithm version 1.2.0
// - Clustering algorithm version 7.0.0
// - Genotyping algorithm version 7.0.0

// the Illumina GenCall Source Code (http://support.illumina.com/downloads/gencall_software.html) includes:
// - NormalizationGoldenGate.cs - normalization routines (version 1.1.0)
// - NormalizationInfinium.cs - normalization routines (version 1.1.2)
// - GenTrain60.cs - clustering (version 6.3.1) and genotyping (6.3.0) routines
// - Utils.cs - closest points to axis, MATLAB robust fit, and other MATLAB routines

// the InfiniumIDATParser Java implementation of the normalization algorithm (version 1.1.2) by Jay Carey includes:
// - InfiniumIDATParser.java - IDAT parsing routines (2010-02-25)
// - InfiniumNormalization.java - normalization routines (version 1.1.2) (2010-01-07)
// - InfiniumUtils.java - closest points to axis, MATLAB robust fit, and other MATLAB routines (2010-01-08)
// this software was used in the 1000 Genomes project (Supplementary chapter 5.3 of http://doi.org/10.1038/nature15394)
// as part of the intensity rank sum test (IRS test) in the Genome STRiP software

// the differences between the normalization algorithm version 1.1.2 and version 1.2.0 are:
// - the original implementation of the madsigma function for robust line fitting is updated as it was updated in MATLAB
// - HandleScale will not use loci with missing data anymore for sub-bead pool bins with less than 192 loci
// - NormalizeSingleBinSingleChannel handles Infinium I (A/T and C/G) probes for sub-bead pool bins with less than 192
// loci
//   for which version 1.1.2 would previously not attempt to compute a background intensity offset

// each AutoCall software determines gender in a slightly different way:
// - AutoConvert (v1.6.3.1) - only uses X chromosome heterozygosity and checks whether it is higher than 0.1
// - AutoConvert 2.0 (v2.0.1.179) - checks whether Y chromosome intensity R values are higher than 0.3 if autosomal call
// rate is higher than 0.97
// - IAAP CLI (v1.1) - same as above but there is a bug in the determination of the autosomal call rate that includes
// loci with null cluster scores as missing
// - Array Analysis CLI (v2.1) - same as above but with the bug removed
// we follow the approach of AutoConvert 2.0 and Array Analysis CLI as default and allow the user to use the approach of
// AutoConvert if requested for inexplicable reasons, AutoConvert 2.0, IAAP CLI, and Array Analysis CLI downsample to
// 10000 random autosomal loci to estimate the autosomal call rate this behavior can be suppressed by setting the
// autosomal call rate threshold from 0.97 to 0.0. However, this cannot be done with Array Analysis CLI

// to replicate the functionality for interoperability purposes, the following bugs were reimplemented:
// matlab_robustfit0 deviates from the original MATLAB implementation (statrobustfit) to match Illumina implementation
// (robustLineFit) when input option addconst/calcoffset is false by erroneously summing the vector into a scalar and
// causing the adjfactor variable to be always equal to 100.0 normalization IDs are allowed to overflow beyond 255,
// which happens with some probes in the Omni5 arrays, which can cause some Infinium I (G/C) probes to be normalized
// together with some Infinium II probes probe pairs with missing values are still used in the normalization step as
// probes with zero values the additional code included in GenTrain 3.0 in the Illumina implementation
// (NormalizeSingleBinSingleChannel) calls MATLAB function trimmean on an array where some values are artificially set
// to zero for no good reasons while other values are left out when determining scale_x with GenTrain 2.0 for
// normalization bins with less than 192 loci we include failed loci as AA loci

/****************************************
 * LITERATURE MENTIONING NORMALIZATION  *
 ****************************************/

// http://doi.org/10.1101/sqb.2003.68.69
// Fan,J.B. et al. (2003) Highly parallel SNP genotyping. Cold Spring Harb Symp Quant Biol, 68, 69–78
// first document that mentions GenCall and GenTrain

// http://patents.google.com/patent/US7035740
// Kermani 2005, Artificial intelligence and global normalization methods for genotyping
// explains how normalization works

// http://patents.google.com/patent/US7467117
// Kermani 2006, Artificial intelligence and global normalization methods for genotyping
// also explains how normaliation works(???)

// http://www.illumina.com/Documents/products/technotes/technote_gencall_data_analysis_software.pdf
// Illumina 2005, Illumina GenCall Data Analysis Software
// it does not describe the normalization but it refers to it

// http://doi.org/10.1016/j.mrfmmm.2004.07.022
// Shen 2005, High-throughput SNP genotyping on universal bead arrays
// introduces the GenTrain algorithm. It explains the GenScores are computed using fuzzy logic

// http://doi.org/10.1038/sj.ejhg.5201528;
// Moorhead et al. 2006, Optimal genotype determination in highly multiplexed SNP data
// in the supplement a normalization procedure very similar to Illumina's is proposed

// http://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina_gt_normalization.pdf
// http://dnatech.genomecenter.ucdavis.edu/documents/illumina_gt_normalization.pdf
// Illumina 2006, Illumina’s Genotyping Data Normalization Methods
// has color versions of the patent figures with details that are missing from the patent including the use of 400
// homozygotes

// http://doi.org/10.1101/gr.5402306
// Peiffer et al. 2006, High-resolution genomic profiling of chromosomal aberrations using Infinium whole-genome
// genotyping explains Illumina normalization with minimum details

// http://www.illumina.com/documents/products/technotes/technote_cnv_algorithms.pdf
// Illumina 2007, DNA Copy Number and Loss of Heterozygosity Analysis Algorithms
// explains how LRR and BAF behave over CNVs

// http://doi.org/10.1093/bioinformatics/btm443
// Teo et al. 2007, A genotype calling algorithm for the Illumina BeadArray platform
// explains Illumina normalization with details that are missing from the patent including the use of 400 homozygotes
// (paper about Illuminus caller)

// http://doi.org/10.1101/gr.5686107
// Oosting et al. 2007, High-resolution copy number analysis of paraffin-embedded archival tissue using SNP BeadArrays
// explains an alternative normalization strategy

// http://doi.org/10.1101/gr.6861907
// Wang et al. 2007, PennCNV: An integrated hidden Markov model designed for high-resolution copy number variation
// detection in whole-genome SNP genotyping data explains Illumina normalization with minimum details

// http://doi.org/10.1093/bioinformatics/btn386
// Giannoulatou et al. 2008 GenoSNP: a variational Bayes within-sample SNP genotyping algorithm that does not require a
// reference population explains an alternative normalization strategy still based on beadpools (paper about GenoSNP
// caller)

// http://doi.org/10.1186/1471-2105-9-409
// Staaf et al. 2008 Normalization of Illumina Infinium whole-genome SNP data improves copy number estimates and allelic
// intensity ratios explains Illumina normalization with minimum details

// http://www.illumina.com/documents/products/technotes/technote_gentrain2.pdf
// Illumina 2009, Improved Cluster Generation with Gentrain2
// explains Gentrain 2.0

// http://doi.org/10.1093/bioinformatics/btp470
// Ritchie et al. 2009 R/Bioconductor software for Illumina’s Inﬁnium whole-genome genotyping BeadChips
// explains an alternative normalization strategy

// http://doi.org/10.1093/nar/gkp552
// LaFramboise et al. 2009 Single nucleotide polymorphism arrays: a decade of biological, computational and
// technological advances explains Illumina normalization with minimum details but defines it as "The computational
// workhorse in the Illumina protocol"

// http://support.illumina.com/documents/products/technotes/technote_array_analysis_workflows.pdf
// Illumina 2011, Microarray Data Analysis Workflows
// explains how IDAT are converted to GTC with AutoCall

// http://doi.org/10.1186/1471-2105-12-68
// Ritchie et al. 2011 Comparing genotyping algorithms for Illumina’s Infinium whole-genome SNP BeadChips
// explains Illumina normalization with minimum details (paper comparing GenCall GenTrain 1.0, Infinium, GenoSNP, CRLMM)

// http://doi.org/10.1007/978-1-61779-555-8_29
// Teo 2011 Genotype Calling for the Illumina Platform
// explains Illumina normalization with details that are missing from the patent including the use of 400 homozygotes

// http://doi.org/10.1093/bioinformatics/bts47
// Goldstein et al. 2012 zCall: a rare variant caller for array-based genotyping
// uses Illumina normalization but no details provided

// http://doi.org/10.1093/bioinformatics/btr673
// Li et al. 2012, M3 : an improved SNP calling algorithm for Illumina BeadArray data
// explains Illumina normalization with minimum details (paper about M3 caller)

// http://doi.org/10.1093/bioinformatics/bts180
// Shah et al. 2012, optiCall: a robust genotype-calling algorithm for rare, low-frequency and common variants
// explains Illumina normalization with minimum details (paper about optiCall caller which uses Illumina normalization)

// http://doi.org/10.1093/bioinformatics/btu107
// Zhou et al. 2014, iCall: a genotype-calling algorithm for rare, low-frequency and common variants on the Illumina
// exome array paper about iCall which uses Illumina normalization

// http://web.stat.tamu.edu/sheather/PDF/WZhou_MSProject.pdf
// Zhou 2014, Segmentation-Based Detection of Mosaic Chromosomal Abnormality in Bladder Cancer Cells Using Whole Genome
// SNP Array includes explanation of the normalization following Illumina's technical note

// http://doi.org/10.1111/pbi.12183
// Wang,S. et al. (2014) Characterization of polyploid wheat genomic diversity using a high-density 90,000 single
// nucleotide polymorphism array. Plant Biotechnol J, 12, 787–796. introduces the polyploid clustering algorithm
// released by Illumina on 2013-10-07

// http://emea.illumina.com/content/dam/illumina-marketing/documents/products/technotes/gentrain3-technical-note-370-2016-015.pdf
// Illumina 2016, Improved Genotype Clustering with GenTrain 3.0
// explains that with less than 192 loci in a single normalization bin it will perform an affine normalization with two
// degrees of freedom rather than six

// http://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/gentrain-tech-note-m-gl-01258/gentrain-tech-note-m-gl-01258.pdf
// Illumina 2023, Genotype clustering with GenTrain 3.0
// explains that with less than 192 loci in a single normalization bin it will perform an affine normalization with two
// degrees of freedom rather than six

#include <ctype.h>
#include <getopt.h>
#include <errno.h>
#include <time.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <htslib/hts.h>
#include <htslib/hfile.h>
#include <htslib/khash.h>
#include <htslib/ksort.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#define IDAT2GTC_VERSION "2024-09-27"

#define AUTOCALL_DATE_FORMAT_DFLT "%m/%d/%y %#I:%M %p" // equivalent to "MM/dd/yyyy h:mm tt"
#define AUTOCALL_VERSION_DFLT "3.0.0"

KSORT_INIT_GENERIC(float)
KSORT_INIT_GENERIC(int)

// void error(const char *format, ...)
//{
//     va_list ap;
//     va_start(ap, format);
//     vfprintf(stderr, format, ap);
//     va_end(ap);
//     exit(-1);
// }
//
// static inline int iupac2bitmask(char iupac)
//{
//     const int A = 1;
//     const int C = 2;
//     const int G = 4;
//     const int T = 8;
//     if ( iupac >= 97 ) iupac -= 32;
//     if ( iupac == 'A' ) return A;
//     if ( iupac == 'C' ) return C;
//     if ( iupac == 'G' ) return G;
//     if ( iupac == 'T' ) return T;
//     if ( iupac == 'M' ) return A|C;
//     if ( iupac == 'R' ) return A|G;
//     if ( iupac == 'W' ) return A|T;
//     if ( iupac == 'S' ) return C|G;
//     if ( iupac == 'Y' ) return C|T;
//     if ( iupac == 'K' ) return G|T;
//     if ( iupac == 'V' ) return A|C|G;
//     if ( iupac == 'H' ) return A|C|T;
//     if ( iupac == 'D' ) return A|G|T;
//     if ( iupac == 'B' ) return C|G|T;
//     if ( iupac == 'N' ) return A|C|G|T;
//     return -1;
// }
//
///**
// *  mkdir_p() - create new directory for a file $fname
// *  @fname:   the file name to create the directory for, the part after last "/" is ignored
// */
// void mkdir_p(const char *fmt, ...)
//{
//    va_list ap;
//    va_start(ap, fmt);
//    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
//    va_end(ap);
//
//    char *path = (char*)malloc(n);
//    va_start(ap, fmt);
//    vsnprintf(path, n, fmt, ap);
//    va_end(ap);
//
//    char *tmp = strdup(path), *p = tmp+1;
//    while (*p)
//    {
//        while (*p && *p!='/') p++;
//        if ( !*p ) break;
//        char ctmp = *p;
//        *p = 0;
//        int ret = mkdir(tmp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//        if ( ret!=0 && errno!=EEXIST ) error("Error creating directory %s: %s\n", path,strerror(errno));
//        *p = ctmp;
//        while ( *p && *p=='/' ) p++;
//    }
//    free(tmp);
//    free(path);
//}

/****************************************
 * hFILE READING FUNCTIONS              *
 ****************************************/

static inline ssize_t HTS_RESULT_USED md5_hread(hFILE *fp, void *buffer, size_t nbytes, hts_md5_context *md5) {
    ssize_t ret = hread(fp, buffer, nbytes);
    if (md5 && ret > 0) hts_md5_update(md5, buffer, ret);
    return ret;
}

static inline int md5_hgetc(hFILE *fp, hts_md5_context *md5) {
    int c = hgetc(fp);
    if (md5 && c != EOF) hts_md5_update(md5, &c, 1);
    return c;
}

// read or skip a fixed number of bytes
static void read_bytes(hFILE *hfile, void *buffer, size_t nbytes, hts_md5_context *md5) {
    if (buffer) {
        if (md5_hread(hfile, buffer, nbytes, md5) < nbytes) {
            error("Failed to read %ld bytes from stream\n", nbytes);
        }
    } else {
        int i, c = 0;
        for (i = 0; i < nbytes; i++) c = md5_hgetc(hfile, md5);
        if (c == EOF) error("Failed to reposition stream forward %ld bytes\n", nbytes);
    }
}

// tests the end-of-file indicator for an hFILE
static int heof(hFILE *hfile) {
    if (hgetc(hfile) == EOF) return 1;
    hfile->begin--;
    return 0;
}

// read or skip a fixed length array
static void read_array(hFILE *hfile, void **arr, size_t *m_arr, size_t nmemb, size_t size, size_t term,
                       hts_md5_context *md5) {
    if (arr) {
        if (!m_arr) {
            *arr = malloc((nmemb + term) * size);
            if (!*arr) error("Failed to allocate memory for array\n");
        } else if (*m_arr < nmemb + term) {
            void *tmp = realloc(*arr, (nmemb + term) * size);
            if (!tmp) error("Failed to allocate memory for array\n");
            *arr = tmp;
            *m_arr = nmemb + term;
        }
        if (md5_hread(hfile, *arr, nmemb * size, md5) < nmemb * size) {
            error("Failed to read %ld bytes from stream\n", nmemb * size);
        }
    } else {
        int i, c = 0;
        for (i = 0; i < nmemb * size; i++) c = md5_hgetc(hfile, md5);
        if (c == EOF) error("Failed to reposition stream forward %ld bytes\n", nmemb * size);
    }
}

// read or skip a length-prefixed string
// http://en.wikipedia.org/wiki/LEB128#Decode_unsigned_integer
static void read_pfx_string(hFILE *hfile, char **str, size_t *m_str, hts_md5_context *md5) {
    uint8_t byte;
    size_t n = 0, shift = 0;
    while (1) {
        if (md5_hread(hfile, (void *)&byte, 1, md5) < 1) {
            error("Failed to read 1 byte from stream\n");
        }
        n |= (size_t)(byte & 0x7F) << shift;
        if (!(byte & 0x80)) break;
        shift += 7;
    }
    if (n || m_str) {
        read_array(hfile, (void **)str, m_str, n, 1, 1, md5);
        if (str) (*str)[n] = '\0';
    }
}

// check whether file is compressed with gzip
static int is_gzip(hFILE *hfile) {
    uint8_t buffer[2];
    if (hpeek(hfile, (void *)buffer, 2) < 2) error("Failed to read 2 bytes from stream\n");
    return (buffer[0] == 0x1f && buffer[1] == 0x8b);
}

static inline int hwrite_uint16(hFILE *hfile, uint16_t num) { return hwrite(hfile, &num, sizeof(uint16_t)); }

static inline int hwrite_int32(hFILE *hfile, int32_t num) { return hwrite(hfile, &num, sizeof(int32_t)); }

// http://en.wikipedia.org/wiki/LEB128#Encode_unsigned_integer
static int hwrite_pfx_string(hFILE *hfile, const char *str) {
    if (!str) {
        hputc(0, hfile);
        return 0;
    }
    size_t n = strlen(str);
    size_t value = n;
    int ret = n;
    do {
        uint8_t byte = value & 0x7f;
        value >>= 7;
        if (value) byte ^= 0x80;
        if (hputc(byte, hfile) == EOF) return -1;
        ret++;
    } while (value);
    if (hwrite(hfile, str, n) < 0) return -1;
    return ret;
}

/****************************************
 * IDAT FILE IMPLEMENTATION             *
 ****************************************/

// http://github.com/snewhouse/glu-genetics/blob/master/glu/lib/illumina.py
// http://github.com/HenrikBengtsson/illuminaio/blob/master/R/readIDAT.R
// /humgen/cnp04/sandbox/bobh/idat_parser/src/edu/mit/broad/gapcore/apps/infinium_idat_parser/InfiniumIDATParser.java

#define NUM_SNPS_READ 1000 // ID_N_CORES
// #define ... 100 // ID_BACKGROUNDS - not used
// #define ... 101 // ID_BACKGROUND_DEVS - not used
#define ILLUMINA_ID 102 // ID_BEAD_TYPES
#define SD 103          // ID_DEVS
#define MEAN 104        // ID_MEANS
// #define ... 105 // ID_MEDIANS - not used
// #define ... 106 // ID_N_BEADS - not used
#define NBEADS 107 // ID_N_GOOD_BEADS
// #define ... 108 // ID_TRIMMED_MEANS - not used
#define MID_BLOCK 200         // ID_ILLUMICODES
#define RUN_INFO 300          // ID_PROCESS_HISTORY
#define RED_GREEN 400         // ID_TENTH_PERCENTILE
#define IDAT_SNP_MANIFEST 401 // ID_SAMPLE_BEADSET
#define SENTRIX_BARCODE 402   // ID_BARCODE
#define CHIP_TYPE 403         // ID_SENTRIX_FORMAT
#define SENTRIX_POSITION 404  // ID_SECTION_LABEL
#define BEADSET 405           // ID_BEADSET
#define IDAT_SAMPLE_NAME 406  // ID_DNA
#define DESCRIPTION 407       // ID_OPA
#define IDAT_SAMPLE_PLATE 408 // ID_DNA_PLATE
#define IDAT_SAMPLE_WELL 409  // ID_WELL
#define IDAT_SAMPLE_COUNT 410 // ID_SAMPLE_COUNT
// #define ... 411 // ID_DX - not used
#define IDAT_VLN 510 // ID_VLN

typedef struct {
    const char *chip_type;
    int num_snps;
    int num_mid_blocks;
    const char *chip_type_guess;
} chip_type_t;

static chip_type_t chip_types[] = {
    {"1-95um_multi-swath_for_4x5M", 4568350, 4568350, "HumanOmni5-4-v1-0"},
    {"1-95um_multi-swath_for_4x5M", 4640213, 4640213, "HumanOmni5-4v1-1"},
    {"1-95um_multi-swath_for_4x5M", 4685673, 4685673, "InfiniumOmni5-4v1-2"},
    {"1-95um_multi-swath_for_4x5M", 4696316, 4696316, "HumanOmni5-4-v1-0"},
    {"1-95um_multi-swath_for_8x2-5M", 2266191, 2266191, "Multi-EthnicGlobal"},
    {"1-95um_multi-swath_for_8x2-5M", 2266367, 2266367, "Multi-EthnicGlobal"},
    {"1-95um_multi-swath_for_8x2-5M", 2266404, 2266404, "Multi-EthnicGlobal"},
    {"1-95um_multi-swath_for_8x2-5M", 2266406, 2266406, "Multi-EthnicGlobal"},
    {"1-95um_multi-swath_for_8x2-5M", 2268676, 2268676, "MEGAEx_BioVU_15075710"},
    {"1-95um_multi-swath_for_8x2-5M", 2315574, 2315574, "Multi-EthnicGlobal"},
    {"1-95um_multi-swath_for_8x2-5M", 2389000, 2389000, "CCPMBiobankMEGA2_20002558X345183"},
    {"1-95um_multi-swath_for_8x2-5M", 2508689, 2508689, "GDA-8v1-0"},
    {"1-95um_multi-swath_for_8x2-5M", 2550870, 2550870, "HumanOmni2.5-8v1"},
    {"1-95um_multi-swath_for_8x2-5M", 2563064, 2563064, "HumanOmni25M-8v1-1"},
    {"1-95um_multi-swath_for_8x2-5M", 2575219, 2575219, "HumanOmni2.5-8v1"},
    {"1-95um_multi-swath_for_8x2-5M", 2605775, 2605775, "HumanOmni25M-8v1-1"},
    {"BeadChip 12x1", 55300, 55300, "humanmethylation27_270596_v1-2 ???"},
    {"BeadChip 12x1Q", 191668, 191668, "CanineHD"},
    {"BeadChip 12x1Q", 299260, 299260, "HumanCytoSNP-12v2-1"},
    {"BeadChip 12x8", 301084, 301084, "HumanCore-12v1-0"},
    {"BeadChip 12x8", 304138, 304138, "HumanExome-12v1-1"},
    {"BeadChip 12x8", 567727, 567727, "HumanCoreExome-12-v1-0"},
    {"BeadChip 12x8", 569060, 569060, "HumanCoreExome-12-v1-0"},
    {"BeadChip 12x8", 573012, 573012, "HumanCoreExome-12-v1-1"},
    {"BeadChip 12x8", 576769, 576769, "HumanCoreExome-12-v1-1"},
    {"BeadChip 12x8", 622399, 622399, "humanmethylation450_15017482_v-1-2 ???"},
    {"BeadChip 12x8", 722405, 722405, "HumanOmniExpress-12-v1-1"},
    {"BeadChip 12x8", 734889, 734889, "HumanOmniExpress-12-v1-0"},
    {"BeadChip 12x8", 736136, 736136, "HumanOmniExpress-12-v1-0"},
    {"BeadChip 1x12", 577085, 8627, "HumanHap550v3"},
    {"BeadChip 1x12", 661182, 49163, "HumanHap650Yv3"},
    {"BeadChip 1x40", 1129736, 57373, "Human1Mv1"},
    {"BeadChip 1x40 66", 1078890, 52497, "Human1Mv1"},
    {"BeadChip 24x1x4", 306776, 306776, "InfiniumCore-24v1-2"},
    {"BeadChip 24x1x4", 527136, 527136, "OncoArray-500K"},
    {"BeadChip 24x1x4", 577781, 577781, "HumanCoreExome-24v1-0"},
    {"BeadChip 24x1x4", 581261, 581261, "HumanCoreExome-24v1-2"},
    {"BeadChip 24x1x4", 582684, 582684, "HumanCoreExome-24v1-1"},
    {"BeadChip 24x1x4", 611866, 611866, "HumanCoreExome-24v1-4"},
    {"BeadChip 24x1x4", 623302, 623302, "PsychChip_15048346"},
    {"BeadChip 24x1x4", 623513, 623513, "InfiniumPsychArray-24v1-1"},
    {"BeadChip 24x1x4", 638714, 638714, "PsychChip_v1-1_15073391"},
    {"BeadChip 24x1x4", 647864, 647864, "InfiniumPsychArray-24v1-3"},
    {"BeadChip 24x1x4", 663209, 663209, "GSA-24v1-0"},
    {"BeadChip 24x1x4", 704215, 704215, "GSA-24v3-0"},
    {"BeadChip 24x1x4", 708013, 708013, "DeCodeGenetics_V1_20012591"},
    {"BeadChip 24x1x4", 710576, 710576, "GSAMD-24v1-0_20011747"},
    {"BeadChip 24x1x4", 710606, 710606, "GSAMD-24v1-0_20011747"},
    {"BeadChip 24x1x4", 710608, 710608, "GSAMD-24v1-0_20011747"},
    {"BeadChip 24x1x4", 715653, 715653, "HumanOmniExpress-24v1-1"},
    {"BeadChip 24x1x4", 716279, 716279, "InfiniumOmniExpress-24v1-2"},
    {"BeadChip 24x1x4", 718963, 718963, "HumanOmniExpress-24-v1-0"},
    {"BeadChip 24x1x4", 719234, 719234, "HumanOmniExpress-24-v1-0"},
    {"BeadChip 24x1x4", 729110, 729110, "ASA-24v1-0"},
    {"BeadChip 24x1x4", 733354, 733354, "GSA-24v2-0"},
    {"BeadChip 24x1x4", 749019, 749019, "DeCodeGenetics_V3_20032937X331991"},
    {"BeadChip 24x1x4", 751614, 751614, "GSAMD-24v3-0-EA_20034606"},
    {"BeadChip 24x1x4", 766804, 766804, "JSA-24v1-0"},
    {"BeadChip 24x1x4", 776509, 776509, "ASA-24v1-0"},
    {"BeadChip 24x1x4", 780343, 780343, "GSAMD-24v2-0_20024620"},
    {"BeadChip 24x1x4", 780509, 780509, "GSAMD-24v2-0_20024620"},
    {"BeadChip 24x1x4", 818205, 818205, "GSA-24v2-0"},
    {"BeadChip 2x10", 321354, 37161, "HumanHap300v2"},
    {"BeadChip 2x12", 381079, 29275, "HumanCNV370v1"},
    {"BeadChip 2x20", 561686, 54936, "HumanHap550v3"},
    {"BeadChip 2x6Q", 1224000, 180026, "Human1M-Duov3"},
    {"BeadChip 2x6Q", 1224629, 180026, "Human1M-Duov3"},
    {"BeadChip 48x4", 730546, 730546, "GSA-MD-48v4-0_20098041"},
    {"BeadChip 4x10", 2623923, 1300482, "HumanOmni2.5-4v1"},
    {"BeadChip 4x10", 2623923, 1323441, "HumanOmni2.5-4v1"},
    {"BeadChip 4x10", 2624666, 1300941, "HumanOmni2.5-4v1"},
    {"BeadChip 4x10", 2624666, 1323725, "HumanOmni2.5-4v1"},
    {"BeadChip 4x10", 2624671, 1323726, "HumanOmni2.5-4v1"},
    {"BeadChip 4x10", 2655594, 1354653, "HumanOmni2.5-4v1"},
    {"BeadChip 4X1X14", 1186430, 1186430, "HumanOmni1-Quad_v1-0"},
    {"BeadChip 4x2Q", 376216, 186490, "HumanCNV370-Quadv3"},
    {"BeadChip 4x3Q", 626122, 208778, "Human610-Quadv1"},
    {"BeadChip 4x3Q", 667447, 208778, "Human660W-Quad_v1"},
    {"BeadChip 8x5", 1052641, 1052641, "infinium-methylationepic-v-1-0 ???"},
    {"BeadChip 8x5", 867478, 867478, "CytoSNP-850K"},
    {"BeadChip 8x5", 988240, 988240, "HumanOmniExpressExome-8-v1-1"},
    {"BeadChip 8x5", 989536, 989536, "HumanOmniExpressExome-8-v1-1"},
    {"BeadChip 8x5", 992824, 992824, "HumanOmniExpressExome-8-v1-4"},
    {"BeadChip 8x5", 996003, 996003, "HumanOmniExpressExome-8-v1-2"},
    {"BeadChip 8x5", 996055, 996055, "HumanOmniExpressExome-8-v1-2"},
    {"SLIDE.15028542.24x1x3", 307984, 307984, "HumanCore-24v1-0"},
    {"SLIDE.15028542.24x1x3", 311460, 311460, "HumanCore-24v1-0"},
    {NULL, 0, 0, NULL}};

typedef struct {
    char *run_time;
    char *block_type;
    char *block_pars;
    char *block_code;
    char *code_version;
} RunInfo;

typedef struct {
    char *fn;
    hFILE *hfile;
    int64_t version;
    int32_t number_toc_entries;
    uint16_t *id;
    int64_t *toc;
    int32_t num_snps;
    int32_t num_mid_blocks;
    int32_t *ilmn_id;
    uint16_t *sd;
    uint16_t *mean;
    uint8_t *nbeads;
    const uint16_t *trimmed_mean; // only used for historical purposes
    uint8_t *mid_block;
    uint8_t red_green[4];
    char *snp_manifest;
    char *sentrix_barcode;
    char *chip_type;
    char *sentrix_position;
    char *beadset;
    char *sample_name;
    char *description;
    char *sample_plate;
    char *sample_well;
    int32_t sample_count;
    char *vln;
    RunInfo *run_infos;
    int32_t m_run_infos;
    const char *chip_type_guess;
    const char *imaging_date;
    const char *scanner_data;
    void *ilmn_id2index;
} idat_t;

KHASH_MAP_INIT_INT(32, int32_t)

static int idat_read(idat_t *idat, uint16_t id) {
    int i;
    for (i = 0; i < idat->number_toc_entries && id != idat->id[i]; i++);
    if (i == idat->number_toc_entries) return -1;
    if (hseek(idat->hfile, idat->toc[i], SEEK_SET) < 0)
        error("Fail to seek to position %ld in IDAT %s file\n", idat->toc[i], idat->fn);

    switch (id) {
    case NUM_SNPS_READ:
        read_bytes(idat->hfile, (void *)&idat->num_snps, sizeof(int32_t), NULL);
        break;
    case ILLUMINA_ID:
        idat->ilmn_id = (int32_t *)malloc(idat->num_snps * sizeof(int32_t));
        read_bytes(idat->hfile, (void *)idat->ilmn_id, idat->num_snps * sizeof(int32_t), NULL);
        int ret;
        idat->ilmn_id2index = kh_init(32);
        khash_t(32) *hash = (khash_t(32) *)idat->ilmn_id2index;
        for (i = 0; i < idat->num_snps; i++) {
            khiter_t k = kh_put(32, hash, idat->ilmn_id[i], &ret);
            if (ret < 0) error("Unable to insert Illumina ID %d in hash table\n", idat->ilmn_id[i]);
            if (ret > 0)
                kh_val(hash, k) = kh_size(hash) - 1;
            else
                error("Duplicate Illumina ID %d in hash table\n", idat->ilmn_id[i]);
        }
        break;
    case SD:
        idat->sd = (uint16_t *)malloc(idat->num_snps * sizeof(uint16_t));
        read_bytes(idat->hfile, (void *)idat->sd, idat->num_snps * sizeof(uint16_t), NULL);
        break;
    case MEAN:
        idat->mean = (uint16_t *)malloc(idat->num_snps * sizeof(uint16_t));
        read_bytes(idat->hfile, (void *)idat->mean, idat->num_snps * sizeof(uint16_t), NULL);
        idat->trimmed_mean = idat->mean;
        break;
    case NBEADS:
        idat->nbeads = (uint8_t *)malloc(idat->num_snps * sizeof(uint8_t));
        read_bytes(idat->hfile, (void *)idat->nbeads, idat->num_snps * sizeof(uint8_t), NULL);
        break;
    case MID_BLOCK:
        read_bytes(idat->hfile, (void *)&idat->num_mid_blocks, sizeof(int32_t), NULL);
        idat->mid_block = (uint8_t *)malloc(idat->num_mid_blocks * sizeof(uint8_t));
        read_bytes(idat->hfile, (void *)idat->mid_block, idat->num_mid_blocks * sizeof(uint8_t), NULL);
        break;
    case RED_GREEN:
        read_bytes(idat->hfile, (void *)&idat->red_green, 4 * sizeof(uint8_t), NULL);
        break;
    case IDAT_SNP_MANIFEST:
        read_pfx_string(idat->hfile, &idat->snp_manifest, NULL, NULL);
        break;
    case SENTRIX_BARCODE:
        read_pfx_string(idat->hfile, &idat->sentrix_barcode, NULL, NULL);
        break;
    case CHIP_TYPE:
        read_pfx_string(idat->hfile, &idat->chip_type, NULL, NULL);
        break;
    case SENTRIX_POSITION:
        read_pfx_string(idat->hfile, &idat->sentrix_position, NULL, NULL);
        break;
    case BEADSET:
        read_pfx_string(idat->hfile, &idat->beadset, NULL, NULL);
        break;
    case IDAT_SAMPLE_NAME:
        read_pfx_string(idat->hfile, &idat->sample_name, NULL, NULL);
        break;
    case DESCRIPTION:
        read_pfx_string(idat->hfile, &idat->description, NULL, NULL);
        break;
    case IDAT_SAMPLE_PLATE:
        read_pfx_string(idat->hfile, &idat->sample_plate, NULL, NULL);
        break;
    case IDAT_SAMPLE_WELL:
        read_pfx_string(idat->hfile, &idat->sample_well, NULL, NULL);
        break;
    case IDAT_SAMPLE_COUNT:
        read_bytes(idat->hfile, (void *)&idat->sample_count, sizeof(int32_t), NULL);
        break;
    case IDAT_VLN:
        read_pfx_string(idat->hfile, &idat->vln, NULL, NULL);
        break;
    case RUN_INFO:
        read_bytes(idat->hfile, (void *)&idat->m_run_infos, sizeof(int32_t), NULL);
        idat->run_infos = (RunInfo *)malloc(idat->m_run_infos * sizeof(RunInfo));
        for (i = 0; i < idat->m_run_infos; i++) {
            read_pfx_string(idat->hfile, &idat->run_infos[i].run_time, NULL, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].block_type, NULL, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].block_pars, NULL, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].block_code, NULL, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].code_version, NULL, NULL);
        }
        break;
    default:
        error("IDAT file format does not support TOC entry %d\n", id);
        break;
    }
    return 0;
}

static idat_t *idat_init(const char *fn, int load_arrays) {
    idat_t *idat = (idat_t *)calloc(1, sizeof(idat_t));
    idat->fn = strdup(fn);
    idat->hfile = hopen(idat->fn, "rb");
    if (idat->hfile == NULL) error("Could not open %s: %s\n", idat->fn, strerror(errno));
    if (is_gzip(idat->hfile)) error("File %s is gzip compressed and currently cannot be sought\n", idat->fn);

    int i;
    uint8_t buffer[4];
    if (hread(idat->hfile, (void *)buffer, 4) < 4) error("Failed to read magic number from %s file\n", idat->fn);
    if (memcmp(buffer, "IDAT", 4) != 0) error("IDAT file %s format identifier is bad\n", idat->fn);

    read_bytes(idat->hfile, (void *)&idat->version, sizeof(int64_t), NULL);
    if (idat->version < 3)
        error("Cannot read IDAT file %s. Unsupported IDAT file format version: %ld\n", idat->fn, idat->version);

    read_bytes(idat->hfile, (void *)&idat->number_toc_entries, sizeof(int32_t), NULL);
    idat->id = (uint16_t *)malloc(idat->number_toc_entries * sizeof(uint16_t));
    idat->toc = (int64_t *)malloc(idat->number_toc_entries * sizeof(int64_t));
    for (i = 0; i < idat->number_toc_entries; i++) {
        read_bytes(idat->hfile, (void *)&idat->id[i], sizeof(uint16_t), NULL);
        read_bytes(idat->hfile, (void *)&idat->toc[i], sizeof(int64_t), NULL);
    }

    for (i = 0; i < idat->number_toc_entries; i++) {
        if (!load_arrays && idat->id[i] <= MID_BLOCK) {
            if (idat->id[i] == MID_BLOCK) {
                if (hseek(idat->hfile, idat->toc[i], SEEK_SET) < 0)
                    error("Fail to seek to position %ld in IDAT %s file\n", idat->toc[i], idat->fn);
                read_bytes(idat->hfile, (void *)&idat->num_mid_blocks, sizeof(int32_t), NULL);
            }
            continue;
        }
        idat_read(idat, idat->id[i]);
    }

    if (idat->chip_type) {
        const chip_type_t *ptr;
        for (ptr = chip_types; ptr->chip_type; ptr++) {
            if (strcmp(idat->chip_type, ptr->chip_type) == 0 && ptr->num_snps == idat->num_snps
                && ptr->num_mid_blocks == idat->num_mid_blocks)
                idat->chip_type_guess = ptr->chip_type_guess;
        }
    }

    for (i = 0; i < idat->m_run_infos; i++) {
        if (strcmp(idat->run_infos[i].block_type, "Scan") != 0) continue;
        idat->imaging_date = idat->run_infos[i].run_time;
        idat->scanner_data = idat->run_infos[i].block_pars;
    }

    return idat;
}

static void idat_destroy(idat_t *idat) {
    if (!idat) return;
    if (hclose(idat->hfile) < 0) error("Error closing IDAT file %s\n", idat->fn);
    free(idat->fn);
    free(idat->id);
    free(idat->toc);
    free(idat->snp_manifest);
    free(idat->sentrix_barcode);
    free(idat->chip_type);
    free(idat->sentrix_position);
    free(idat->beadset);
    free(idat->sample_name);
    free(idat->description);
    free(idat->sample_plate);
    free(idat->sample_well);
    free(idat->vln);
    int i;
    for (i = 0; i < idat->m_run_infos; i++) {
        free(idat->run_infos[i].run_time);
        free(idat->run_infos[i].block_type);
        free(idat->run_infos[i].block_pars);
        free(idat->run_infos[i].block_code);
        free(idat->run_infos[i].code_version);
    }
    free(idat->run_infos);
    free(idat->ilmn_id);
    free(idat->sd);
    free(idat->mean);
    free(idat->nbeads);
    free(idat->mid_block);
    if (idat->ilmn_id2index) kh_destroy(32, idat->ilmn_id2index);
    free(idat);
}

static void idat_to_csv(const idat_t *idat, FILE *stream, int verbose) {
    int i;
    fprintf(stream, "Illumina, Inc.\n");
    fprintf(stream, "[Heading]\n");
    fprintf(stream, "Descriptor File Name,%s\n", strrchr(idat->fn, '/') ? strrchr(idat->fn, '/') + 1 : idat->fn);
    fprintf(stream, "IDAT file version,%ld\n", idat->version);
    fprintf(stream, "Number of TOC entries,%d\n", idat->number_toc_entries);
    fprintf(stream, "Probes Count,%d\n", idat->num_snps);
    fprintf(stream, "Mid Blocks Count,%d\n", idat->num_mid_blocks);
    fprintf(stream, "Red Green,%02x %02x %02x %02x\n", idat->red_green[0], idat->red_green[1], idat->red_green[2],
            idat->red_green[3]);
    fprintf(stream, "SNP Manifest,%s\n", idat->snp_manifest ? idat->snp_manifest : "");
    fprintf(stream, "Sentrix Barcode,%s\n", idat->sentrix_barcode);
    fprintf(stream, "Chip Type,%s\n", idat->chip_type);
    fprintf(stream, "Sentrix Position,%s\n", idat->sentrix_position);
    fprintf(stream, "BeadSet,%s\n", idat->beadset ? idat->beadset : "");
    fprintf(stream, "Sample Name,%s\n", idat->sample_name ? idat->sample_name : "");
    fprintf(stream, "Description,%s\n", idat->description ? idat->description : "");
    fprintf(stream, "Sample Plate,%s\n", idat->sample_plate ? idat->sample_plate : "");
    fprintf(stream, "Sample Well,%s\n", idat->sample_well ? idat->sample_well : "");
    fprintf(stream, "Sample Count,%d\n", idat->sample_count);
    fprintf(stream, "Vln,%s\n", idat->vln ? idat->vln : "");
    fprintf(stream, "Chip Prefix (Guess),%s\n", idat->chip_type_guess ? idat->chip_type_guess : "Unknown");
    fprintf(stream, "[Assay]\n");
    fprintf(stream, "IlmnID,Sd,Mean,Nbeads\n");
    if (verbose) {
        for (i = 0; i < idat->num_snps; i++)
            fprintf(stream, "%d,%d,%d,%d\n", idat->ilmn_id[i], idat->sd[i], idat->mean[i], idat->nbeads[i]);
        fprintf(stream, "[Mid Blocks]\n");
        for (i = 0; i < idat->num_mid_blocks; i++) fprintf(stream, "%d\n", idat->mid_block[i]);
    } else {
        fprintf(stream, "... use --verbose to visualize Assay data ...\n");
        fprintf(stream, "[Mid Blocks]\n");
        fprintf(stream, "... use --verbose to visualize Mid Blocks data ...\n");
    }
    fprintf(stream, "[Run Infos]\n");
    for (i = 0; i < idat->m_run_infos; i++) {
        fprintf(stream, "%s\t%s\t%s\t%s\t%s\n", idat->run_infos[i].run_time, idat->run_infos[i].block_type,
                idat->run_infos[i].block_pars, idat->run_infos[i].block_code, idat->run_infos[i].code_version);
    }
}

static void idats_to_tsv(idat_t **idats, int n, FILE *stream) {
    fprintf(stream,
            "idat\tnumber_probes\tnumber_mid_blocks\tred_green\tmanifest_file\tsentrix_"
            "barcode\tchip_type\t"
            "sentrix_position\tbeadset\tsample_name\tdescription\tsample_plate\tsample_"
            "well\tsample_count\tvln\t"
            "chip_type_guess\tscan_date\tscanner_data\n");
    int i;
    for (i = 0; i < n; i++) {
        idat_t *idat = idats[i];
        fprintf(stream,
                "%s\t%d\t%d\t%02x %02x %02x "
                "%02x\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n",
                strrchr(idat->fn, '/') ? strrchr(idat->fn, '/') + 1 : idat->fn, idat->num_snps, idat->num_mid_blocks,
                idat->red_green[0], idat->red_green[1], idat->red_green[2], idat->red_green[3],
                idat->snp_manifest ? idat->snp_manifest : "", idat->sentrix_barcode, idat->chip_type,
                idat->sentrix_position, idat->beadset ? idat->beadset : "", idat->sample_name ? idat->sample_name : "",
                idat->description ? idat->description : "", idat->sample_plate ? idat->sample_plate : "",
                idat->sample_well ? idat->sample_well : "", idat->sample_count, idat->vln ? idat->vln : "",
                idat->chip_type_guess ? idat->chip_type_guess : "Unknown", idat->imaging_date ? idat->imaging_date : "",
                idat->scanner_data ? idat->scanner_data : "");
    }
}

/****************************************
 * GTC FILE IMPLEMENTATION              *
 ****************************************/

// http://github.com/broadinstitute/picard/blob/master/src/main/java/picard/arrays/illumina/InfiniumGTCFile.java
// http://github.com/Illumina/BeadArrayFiles/blob/develop/docs/GTC_File_Format_v5.pdf
// http://github.com/Illumina/BeadArrayFiles/blob/develop/module/GenotypeCalls.py

#define NUM_SNPS 1
#define PLOIDY 2      // AutoConvert 2.0
#define PLOIDY_TYPE 3 // AutoConvert 2.0
#define GTC_SAMPLE_NAME 10
#define GTC_SAMPLE_PLATE 11
#define GTC_SAMPLE_WELL 12
#define CLUSTER_FILE 100
#define GTC_SNP_MANIFEST 101
#define IMAGING_DATE 200
#define AUTOCALL_DATE 201
#define AUTOCALL_VERSION 300
#define NORMALIZATION_TRANSFORMS 400
#define CONTROLS_X 500
#define CONTROLS_Y 501
#define RAW_X 1000
#define RAW_Y 1001
#define GENOTYPES 1002
#define BASE_CALLS 1003
#define GENOTYPE_SCORES 1004
#define SCANNER_DATA 1005
#define CALL_RATE 1006
#define GENDER 1007
#define LOGR_DEV 1008
#define GC10 1009
#define DX 1010
#define SAMPLE_DATA 1011
#define B_ALLELE_FREQS 1012   // AutoConvert 2.0
#define LOGR_RATIOS 1013      // AutoConvert 2.0
#define PERCENTILES_X 1014    // AutoConvert 2.0
#define PERCENTILES_Y 1015    // AutoConvert 2.0
#define SLIDE_IDENTIFIER 1016 // AutoConvert 2.0

// static const char *code2genotype[] = {
//     "NC",       "AA",       "AB",       "BB",       "NULL",     "A",        "B",        "AAA",
//     "AAB",      "ABB",      "BBB",      "AAAA",     "AAAB",     "AABB",     "ABBB",     "BBBB",
//     "AAAAA",    "AAAAB",    "AAABB",    "AABBB",    "ABBBB",    "BBBBB",    "AAAAAA",   "AAAAAB",
//     "AAAABB",   "AAABBB",   "AABBBB",   "ABBBBB",   "BBBBBB",   "AAAAAAA",  "AAAAAAB",  "AAAAABB",
//     "AAAABBB",  "AAABBBB",  "AABBBBB",  "ABBBBBB",  "BBBBBBB",  "AAAAAAAA", "AAAAAAAB", "AAAAAABB",
//     "AAAAABBB", "AAAABBBB", "AAABBBBB", "AABBBBBB", "ABBBBBBB", "BBBBBBBB"};

typedef struct {
    int32_t version;
    float offset_x;
    float offset_y;
    float scale_x;
    float scale_y;
    float shear;
    float theta;
    float cvx;
    float cvy;
    float nn12;
    float rr12;
    float taa;
    float tbb;
} XForm;

typedef char BaseCall[2];

typedef struct {
    char *scanner_name;
    int32_t pmt_green;
    int32_t pmt_red;
    char *scanner_version;
    char *imaging_user;
} ScannerData;

typedef struct {
    float p50gc;
    int32_t num_calls;
    int32_t num_no_calls;
    int32_t num_intensity_only;
} SampleData;

typedef uint16_t Percentiles[3];

typedef struct {
    char *fn;
    hFILE *hfile;
    int32_t version;
    int32_t number_toc_entries;
    uint16_t *id;
    int32_t *toc;
    int32_t num_snps;
    int32_t ploidy;
    int32_t ploidy_type;
    char *sample_name;
    char *sample_plate;
    char *sample_well;
    char *cluster_file;
    char *snp_manifest;
    char *imaging_date;
    char *autocall_date;
    char *autocall_version;
    XForm *normalization_transforms;
    size_t m_normalization_transforms;
    uint16_t *controls_x;
    size_t m_controls_x;
    uint16_t *controls_y;
    size_t m_controls_y;
    ScannerData scanner_data;
    float call_rate;
    char gender;
    float logr_dev;
    float p10gc;
    int32_t dx;
    SampleData sample_data;
    Percentiles percentiles_x;
    Percentiles percentiles_y;
    char *sentrix_id;

    char *display_name;
    float *sin_theta; // precomputed sine transforms
    float *cos_theta; // precomputed cosine transforms

    uint16_t *raw_x;
    size_t m_raw_x;
    uint16_t *raw_y;
    size_t m_raw_y;
    uint8_t *genotypes;
    size_t m_genotypes;
    BaseCall *base_calls;
    size_t m_base_calls;
    float *genotype_scores;
    size_t m_genotype_scores;
    float *b_allele_freqs;
    size_t m_b_allele_freqs;
    float *logr_ratios;
    size_t m_logr_ratios;
} gtc_t;

// returns the length of a string including the variable-length prefix encoding the number of characters
static int leb128_strlen(const char *s) {
    if (!s) return 1;
    size_t n = strlen(s);
    size_t value = n++;
    while (value >>= 7) n++;
    return n;
}

static int gtc_write(const gtc_t *gtc, const char *fn, int gtc_file_version) {
    hFILE *hfile = hopen(fn, "wb");
    if (hfile == NULL) error("Could not open %s: %s\n", fn, strerror(errno));
    const uint8_t header[4] = {'g', 't', 'c', gtc_file_version};
    if (hwrite(hfile, header, 4) < 0) return -1;
    int32_t number_toc_entries = gtc_file_version == 3 ? 24 : 31;
    if (hwrite_int32(hfile, number_toc_entries) < 0) return -1;
    int offset = 4 + sizeof(int32_t) + number_toc_entries * (sizeof(uint16_t) + sizeof(int32_t));
    if (hwrite_uint16(hfile, NUM_SNPS) < 0) return -1;
    if (hwrite_int32(hfile, gtc->num_snps) < 0) return -1;
    if (gtc_file_version != 3) {
        if (hwrite_uint16(hfile, PLOIDY) < 0) return -1;
        if (hwrite_int32(hfile, gtc->ploidy) < 0) return -1;
        if (hwrite_uint16(hfile, PLOIDY_TYPE) < 0) return -1;
        if (hwrite_int32(hfile, gtc->ploidy_type) < 0) return -1;
    }
    if (hwrite_uint16(hfile, GTC_SAMPLE_NAME) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->sample_name);
    if (hwrite_uint16(hfile, GTC_SAMPLE_PLATE) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->sample_plate);
    if (hwrite_uint16(hfile, GTC_SAMPLE_WELL) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->sample_well);
    if (hwrite_uint16(hfile, CLUSTER_FILE) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->cluster_file);
    if (hwrite_uint16(hfile, GTC_SNP_MANIFEST) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->snp_manifest);
    if (hwrite_uint16(hfile, IMAGING_DATE) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->imaging_date);
    if (hwrite_uint16(hfile, AUTOCALL_DATE) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->autocall_date);
    if (hwrite_uint16(hfile, AUTOCALL_VERSION) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->autocall_version);
    if (hwrite_uint16(hfile, NORMALIZATION_TRANSFORMS) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->m_normalization_transforms * sizeof(XForm) + sizeof(int32_t);
    if (hwrite_uint16(hfile, CONTROLS_X) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->m_controls_x * sizeof(uint16_t) + sizeof(int32_t);
    if (hwrite_uint16(hfile, CONTROLS_Y) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->m_controls_y * sizeof(uint16_t) + sizeof(int32_t);
    if (hwrite_uint16(hfile, RAW_X) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->num_snps * sizeof(uint16_t) + sizeof(int32_t);
    if (hwrite_uint16(hfile, RAW_Y) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->num_snps * sizeof(uint16_t) + sizeof(int32_t);
    if (hwrite_uint16(hfile, GENOTYPES) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->num_snps * sizeof(uint8_t) + sizeof(int32_t);
    if (hwrite_uint16(hfile, BASE_CALLS) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->num_snps * sizeof(BaseCall) + sizeof(int32_t);
    if (hwrite_uint16(hfile, GENOTYPE_SCORES) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += gtc->num_snps * sizeof(float) + sizeof(int32_t);
    if (hwrite_uint16(hfile, SCANNER_DATA) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += leb128_strlen(gtc->scanner_data.scanner_name) + sizeof(float) + sizeof(float)
              + leb128_strlen(gtc->scanner_data.scanner_version) + leb128_strlen(gtc->scanner_data.imaging_user);
    if (hwrite_uint16(hfile, CALL_RATE) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += sizeof(float);
    if (hwrite_uint16(hfile, GENDER) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += sizeof(char);
    if (hwrite_uint16(hfile, LOGR_DEV) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += sizeof(float);
    if (hwrite_uint16(hfile, GC10) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += sizeof(float);
    if (hwrite_uint16(hfile, DX) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += sizeof(int32_t);
    if (hwrite_uint16(hfile, SAMPLE_DATA) < 0) return -1;
    if (hwrite_int32(hfile, offset) < 0) return -1;
    offset += sizeof(SampleData);
    if (gtc_file_version != 3) {
        if (hwrite_uint16(hfile, B_ALLELE_FREQS) < 0) return -1;
        if (hwrite_int32(hfile, offset) < 0) return -1;
        offset += gtc->num_snps * sizeof(float) + sizeof(int32_t);
        if (hwrite_uint16(hfile, LOGR_RATIOS) < 0) return -1;
        if (hwrite_int32(hfile, offset) < 0) return -1;
        offset += gtc->num_snps * sizeof(float) + sizeof(int32_t);
        if (hwrite_uint16(hfile, PERCENTILES_X) < 0) return -1;
        if (hwrite_int32(hfile, offset) < 0) return -1;
        offset += sizeof(Percentiles);
        if (hwrite_uint16(hfile, PERCENTILES_Y) < 0) return -1;
        if (hwrite_int32(hfile, offset) < 0) return -1;
        offset += sizeof(Percentiles);
        if (hwrite_uint16(hfile, SLIDE_IDENTIFIER) < 0) return -1;
        if (hwrite_int32(hfile, offset) < 0) return -1;
        offset += leb128_strlen(gtc->sentrix_id);
    }

    if (hwrite_pfx_string(hfile, gtc->sample_name) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->sample_plate) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->sample_well) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->cluster_file) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->snp_manifest) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->imaging_date) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->autocall_date) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->autocall_version) < 0) return -1;
    if (hwrite(hfile, (const void *)&gtc->m_normalization_transforms, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->normalization_transforms, gtc->m_normalization_transforms * sizeof(XForm)) < 0)
        return -1;
    if (hwrite(hfile, (const void *)&gtc->m_controls_x, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->controls_x, gtc->m_controls_x * sizeof(uint16_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)&gtc->m_controls_y, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->controls_y, gtc->m_controls_y * sizeof(uint16_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)&gtc->num_snps, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->raw_x, gtc->num_snps * sizeof(uint16_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)&gtc->num_snps, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->raw_y, gtc->num_snps * sizeof(uint16_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)&gtc->num_snps, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->genotypes, gtc->num_snps * sizeof(uint8_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)&gtc->num_snps, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->base_calls, gtc->num_snps * sizeof(BaseCall)) < 0) return -1;
    if (hwrite(hfile, (const void *)&gtc->num_snps, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, (const void *)gtc->genotype_scores, gtc->num_snps * sizeof(float)) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->scanner_data.scanner_name) < 0) return -1;
    if (hwrite(hfile, &gtc->scanner_data.pmt_green, sizeof(float)) < 0) return -1;
    if (hwrite(hfile, &gtc->scanner_data.pmt_red, sizeof(float)) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->scanner_data.scanner_version) < 0) return -1;
    if (hwrite_pfx_string(hfile, gtc->scanner_data.imaging_user) < 0) return -1;
    if (hwrite(hfile, &gtc->call_rate, sizeof(float)) < 0) return -1;
    if (hwrite(hfile, &gtc->gender, sizeof(char)) < 0) return -1;
    if (hwrite(hfile, &gtc->logr_dev, sizeof(float)) < 0) return -1;
    if (hwrite(hfile, &gtc->p10gc, sizeof(float)) < 0) return -1;
    if (hwrite(hfile, &gtc->dx, sizeof(int32_t)) < 0) return -1;
    if (hwrite(hfile, &gtc->sample_data, sizeof(SampleData)) < 0) return -1;
    if (gtc_file_version != 3) {
        if (hwrite(hfile, (const void *)&gtc->num_snps, sizeof(int32_t)) < 0) return -1;
        if (hwrite(hfile, (const void *)gtc->b_allele_freqs, gtc->num_snps * sizeof(float)) < 0) return -1;
        if (hwrite(hfile, (const void *)&gtc->num_snps, sizeof(int32_t)) < 0) return -1;
        if (hwrite(hfile, (const void *)gtc->logr_ratios, gtc->num_snps * sizeof(float)) < 0) return -1;
        if (hwrite(hfile, (const void *)gtc->percentiles_x, sizeof(Percentiles)) < 0) return -1;
        if (hwrite(hfile, (const void *)gtc->percentiles_y, sizeof(Percentiles)) < 0) return -1;
        if (hwrite_pfx_string(hfile, gtc->sentrix_id) < 0) return -1;
    }
    if (hclose(hfile) < 0) error("Error closing GTC file %s\n", fn);
    return 0;
}

static void gtc_destroy(gtc_t *gtc) {
    if (!gtc) return;
    if (gtc->hfile && hclose(gtc->hfile) < 0) error("Error closing GTC file %s\n", gtc->fn);
    free(gtc->fn);
    free(gtc->id);
    free(gtc->toc);
    free(gtc->sample_name);
    free(gtc->sample_plate);
    free(gtc->sample_well);
    free(gtc->cluster_file);
    free(gtc->snp_manifest);
    free(gtc->imaging_date);
    free(gtc->autocall_date);
    free(gtc->autocall_version);
    free(gtc->normalization_transforms);
    free(gtc->controls_x);
    free(gtc->controls_y);
    free(gtc->scanner_data.scanner_name);
    free(gtc->scanner_data.scanner_version);
    free(gtc->scanner_data.imaging_user);
    free(gtc->sentrix_id);

    free(gtc->display_name);
    free(gtc->sin_theta);
    free(gtc->cos_theta);

    free(gtc->raw_x);
    free(gtc->raw_y);
    free(gtc->genotypes);
    free(gtc->base_calls);
    free(gtc->genotype_scores);
    free(gtc->b_allele_freqs);
    free(gtc->logr_ratios);
    free(gtc);
}

/****************************************
 * BPM FILE IMPLEMENTATION              *
 ****************************************/

// http://github.com/snewhouse/glu-genetics/blob/master/glu/lib/illumina.py
// http://github.com/Illumina/BeadArrayFiles/blob/develop/module/BeadPoolManifest.py

typedef struct {
    int32_t version;
    uint8_t norm_id; // Normalization lookups from manifest. This indexes into list of
                     // normalization transforms read from GTC file
    char *ilmn_id;   // IlmnID (probe identifier) of locus
    char *name;      // Name (variant identifier) of locus
    int32_t index;
    char *ilmn_strand; // TOP BOT PLUS MINUS or Top Bot P M
    char *snp;         // SNP value for locus (e.g., [A/C])
    char *chrom;       // Chromosome for the locus (e.g., XY)
    char *ploidy;
    char *species;
    char *map_info; // Mapping location of locus
    char *customer_strand;
    int32_t address_a;        // AddressA ID of locus
    char *allele_a_probe_seq; // CSV files or BPM files with version 4 data block
    int32_t address_b;        // AddressB ID of locus (0 if none)
    char *allele_b_probe_seq; // CSV files or BPM files with version 4 data block (empty if
                              // none)
    char *genome_build;
    char *source;
    char *source_version;
    char *source_strand;
    char *source_seq;      // CSV files or BPM files with version 4 data block
    char *top_genomic_seq; // CSV files or BPM files with version 4 data block
    int32_t beadset_id;    // CSV files
    uint8_t exp_clusters;
    uint8_t intensity_only;
    uint8_t assay_type; // Identifies type of assay (0 - Infinium II, 1 - Infinium I (A/T),
                        // 2 - Infinium I (G/C)
    uint8_t assay_type_csv;
    float frac_a;
    float frac_c;
    float frac_g;
    float frac_t;
    char *ref_strand; // RefStrand annotation
} LocusEntry;

// retrieve assay type following (allele_a_probe_seq, source_seq) -> assay_type map
// (...W., ...W[./.]W...) -> 1
// (...S., ...S[./.]S...) -> 2
// (...S., ...S[./.]W...) -> 1
// (...S., ...W[./.]S...) -> 1
// (...W., ...S[./.]W...) -> 2
// (...W., ...W[./.]S...) -> 2
static uint8_t get_assay_type(const char *allele_a_probe_seq, const char *allele_b_probe_seq, const char *source_seq) {
    if (!allele_a_probe_seq || !source_seq) return 0xFF;
    if (!allele_b_probe_seq) return 0;
    const char *left = strchr(source_seq, '[');
    const char *right = strchr(source_seq, ']');
    if (!left || !right) error("Source sequence is malformed: %s\n", source_seq);
    char trail_left = toupper(*(left - 1));
    char trail_right = toupper(*(right + 1));
    if ((trail_left == 'A' || trail_left == 'T') && (trail_right == 'A' || trail_right == 'T')) return 1;
    if ((trail_left == 'C' || trail_left == 'G') && (trail_right == 'C' || trail_right == 'G')) return 2;
    int i = 2;
    while (!(iupac2bitmask(allele_a_probe_seq[strlen(allele_a_probe_seq) - i])
             & iupac2bitmask(allele_b_probe_seq[strlen(allele_b_probe_seq) - i])))
        i++;
    char trail_a_probe_seq = toupper(allele_a_probe_seq[strlen(allele_a_probe_seq) - i]);
    if (trail_a_probe_seq == 'C' || trail_a_probe_seq == 'G' || trail_a_probe_seq == 'S') return 1;
    if (trail_a_probe_seq == 'A' || trail_a_probe_seq == 'T' || trail_a_probe_seq == 'W') return 2;
    // these weird rule were deduced from manifests for array GDA_PGx-8v1-0_20042614
    if (trail_a_probe_seq == 'Y' && trail_right == 'G') return 1;
    if (trail_a_probe_seq == 'Y' && trail_right == 'T') return 1;
    if (trail_a_probe_seq == 'Y' && trail_right == 'A') return 2;
    if (trail_a_probe_seq == 'K' && trail_right == 'C') return 1;
    if (trail_a_probe_seq == 'K' && trail_right == 'A') return 2;
    if (trail_a_probe_seq == 'M' && trail_right == 'G') return 1;
    if (trail_a_probe_seq == 'M' && trail_right == 'T') return 2;
    if (trail_a_probe_seq == 'R' && trail_right == 'C') return 1;
    if (trail_a_probe_seq == 'R' && trail_right == 'T') return 2;
    fprintf(stderr, "Warning: Unable to retrieve assay type: %s %s %s\n", allele_a_probe_seq, allele_b_probe_seq,
            source_seq);
    return 0xFF;
}

static void locusentry_read(LocusEntry *locus_entry, hFILE *hfile, hts_md5_context *md5) {
    locus_entry->norm_id = 0xFF;
    read_bytes(hfile, (void *)&locus_entry->version, sizeof(int32_t), md5);
    if (locus_entry->version < 4 || locus_entry->version == 5 || locus_entry->version > 8)
        error("Locus version %d in manifest file not supported\n", locus_entry->version);
    read_pfx_string(hfile, &locus_entry->ilmn_id, NULL, md5);
    read_pfx_string(hfile, &locus_entry->name, NULL, md5);
    read_pfx_string(hfile, NULL, NULL, md5);
    read_pfx_string(hfile, NULL, NULL, md5);
    read_pfx_string(hfile, NULL, NULL, md5);
    read_bytes(hfile, (void *)&locus_entry->index, sizeof(int32_t), md5);
    read_pfx_string(hfile, NULL, NULL, md5);
    read_pfx_string(hfile, &locus_entry->ilmn_strand, NULL, md5);
    read_pfx_string(hfile, &locus_entry->snp, NULL, md5);
    read_pfx_string(hfile, &locus_entry->chrom, NULL, md5);
    read_pfx_string(hfile, &locus_entry->ploidy, NULL, md5);
    read_pfx_string(hfile, &locus_entry->species, NULL, md5);
    read_pfx_string(hfile, &locus_entry->map_info, NULL, md5);
    read_pfx_string(hfile, &locus_entry->top_genomic_seq, NULL, md5); // only version 4
    read_pfx_string(hfile, &locus_entry->customer_strand, NULL, md5);
    read_bytes(hfile, (void *)&locus_entry->address_a, sizeof(int32_t), md5);
    read_bytes(hfile, (void *)&locus_entry->address_b, sizeof(int32_t), md5);
    read_pfx_string(hfile, &locus_entry->allele_a_probe_seq, NULL, md5); // only version 4
    read_pfx_string(hfile, &locus_entry->allele_b_probe_seq, NULL, md5); // only version 4
    read_pfx_string(hfile, &locus_entry->genome_build, NULL, md5);
    read_pfx_string(hfile, &locus_entry->source, NULL, md5);
    read_pfx_string(hfile, &locus_entry->source_version, NULL, md5);
    read_pfx_string(hfile, &locus_entry->source_strand, NULL, md5);
    read_pfx_string(hfile, &locus_entry->source_seq, NULL, md5); // only version 4
    if (locus_entry->source_seq) {
        char *ptr = strchr(locus_entry->source_seq, '-');
        if (ptr && *(ptr - 1) == '/') {
            *ptr = *(ptr - 2);
            *(ptr - 2) = '-';
        }
    }

    if (locus_entry->version >= 6) {
        read_bytes(hfile, NULL, 1, md5);
        read_bytes(hfile, (void *)&locus_entry->exp_clusters, sizeof(int8_t), md5);
        read_bytes(hfile, (void *)&locus_entry->intensity_only, sizeof(int8_t), md5);
        read_bytes(hfile, (void *)&locus_entry->assay_type, sizeof(uint8_t), md5);

        if (locus_entry->assay_type < 0 || locus_entry->assay_type > 2)
            error("Format error in reading assay type from locus entry\n");
        if (locus_entry->address_b == 0 && locus_entry->assay_type != 0)
            error("Manifest format error: Assay type is inconsistent with address B\n");
        if (locus_entry->address_b != 0 && locus_entry->assay_type == 0)
            error("Manifest format error: Assay type is inconsistent with address B\n");
    } else {
        locus_entry->assay_type =
            get_assay_type(locus_entry->allele_a_probe_seq, locus_entry->allele_b_probe_seq, locus_entry->source_seq);
    }

    if (locus_entry->version >= 7) {
        read_bytes(hfile, &locus_entry->frac_a, sizeof(float), md5);
        read_bytes(hfile, &locus_entry->frac_c, sizeof(float), md5);
        read_bytes(hfile, &locus_entry->frac_t, sizeof(float), md5);
        read_bytes(hfile, &locus_entry->frac_g, sizeof(float), md5);
    }
    if (locus_entry->version >= 8) read_pfx_string(hfile, &locus_entry->ref_strand, NULL, md5);
}

typedef struct {
    char *fn;
    hFILE *hfile; // bpm file
    htsFile *fp;  // csv file
    int32_t version;
    char *manifest_name;  // Name of manifest
    char *control_config; // Control description from manifest
    int32_t num_loci;     // Number of loci in manifest
    int32_t *indexes;
    char **names; // Names of loci from manifest
    void *names2index;
    uint8_t *norm_ids;
    LocusEntry *locus_entries;
    uint8_t *norm_lookups;
    char **header;
    size_t m_header;
    char unsigned md5_buf[16];
} bpm_t;

static uint8_t *bpm_norm_lookups(bpm_t *bpm) {
    int i;
    uint8_t sorted_norm_ids[256];
    for (i = 0; i < 256; i++) sorted_norm_ids[i] = 0xFF;
    for (i = 0; i < bpm->num_loci; i++) {
        int norm_id = bpm->locus_entries[i].norm_id;
        sorted_norm_ids[norm_id] = norm_id;
    }
    int j = 0;
    for (i = 0; i < 256; i++)
        if (sorted_norm_ids[i] != 0xFF) sorted_norm_ids[j++] = sorted_norm_ids[i];
    uint8_t *norm_lookups = (uint8_t *)malloc(256 * sizeof(uint8_t *));
    memset((void *)norm_lookups, 0xFF, 256 * sizeof(uint8_t *));
    for (i = 0; i < j; i++) norm_lookups[sorted_norm_ids[i]] = i;
    return norm_lookups;
}

static bpm_t *bpm_init(const char *fn, int eof_check, int make_dict, int checksum) {
    bpm_t *bpm = (bpm_t *)calloc(1, sizeof(bpm_t));
    bpm->fn = strdup(fn);
    bpm->hfile = hopen(bpm->fn, "rb");
    if (bpm->hfile == NULL) error("Could not open %s: %s\n", bpm->fn, strerror(errno));
    if (is_gzip(bpm->hfile)) error("File %s is gzip compressed and currently cannot be sought\n", bpm->fn);

    hts_md5_context *md5 = checksum ? hts_md5_init() : NULL;

    int i;
    uint8_t buffer[4];
    if (md5_hread(bpm->hfile, (void *)buffer, 4, md5) < 4) error("Failed to read magic number from %s file\n", bpm->fn);
    if (memcmp(buffer, "BPM", 3) != 0) error("BPM file %s format identifier is bad\n", bpm->fn);
    if (buffer[3] != 1) error("BPM file %s version is unknown\n", bpm->fn);

    read_bytes(bpm->hfile, (void *)&bpm->version, sizeof(int32_t), md5);
    if (bpm->version & 0x1000) bpm->version ^= 0x1000;
    if (bpm->version > 5 || bpm->version < 3) error("BPM file %s version %d is unsupported\n", bpm->fn, bpm->version);
    read_pfx_string(bpm->hfile, &bpm->manifest_name, NULL, md5);

    if (bpm->version > 1) read_pfx_string(bpm->hfile, &bpm->control_config, NULL, md5);

    read_bytes(bpm->hfile, (void *)&bpm->num_loci, sizeof(int32_t), md5);
    read_array(bpm->hfile, (void **)&bpm->indexes, NULL, bpm->num_loci, sizeof(int32_t), 0, md5);
    bpm->names = (char **)malloc(bpm->num_loci * sizeof(char *));
    for (i = 0; i < bpm->num_loci; i++) read_pfx_string(bpm->hfile, &bpm->names[i], NULL, md5);
    if (make_dict) {
        bpm->names2index = khash_str2int_init();
        for (i = 0; i < bpm->num_loci; i++) {
            if (khash_str2int_has_key(bpm->names2index, bpm->names[i]))
                error("Illumina probe %s present multiple times in file %s\n", bpm->names[i], fn);
            khash_str2int_inc(bpm->names2index, bpm->names[i]);
        }
    }
    read_array(bpm->hfile, (void **)&bpm->norm_ids, NULL, bpm->num_loci, sizeof(uint8_t), 0, md5);

    bpm->locus_entries = (LocusEntry *)malloc(bpm->num_loci * sizeof(LocusEntry));
    LocusEntry locus_entry;
    for (i = 0; i < bpm->num_loci; i++) {
        memset(&locus_entry, 0, sizeof(LocusEntry));
        locusentry_read(&locus_entry, bpm->hfile, md5);
        int idx = locus_entry.index - 1;
        if (idx < 0 || idx >= bpm->num_loci) error("Locus entry index %d is out of boundaries\n", locus_entry.index);
        if (bpm->norm_ids[idx] > 100)
            error("Manifest format error: read invalid normalization ID %d\n", bpm->norm_ids[idx]);
        // To mimic the flawed byte-wrapping behavior from GenomeStudio, AutoCall, and
        // IAAP, this value is allowed to overflow beyond 255, which happens with some
        // probes in the Omni5 arrays
        bpm->norm_ids[idx] += 100 * locus_entry.assay_type;
        locus_entry.norm_id = bpm->norm_ids[idx];
        memcpy(&bpm->locus_entries[idx], &locus_entry, sizeof(LocusEntry));
    }
    bpm->norm_lookups = bpm_norm_lookups(bpm);
    for (i = 0; i < bpm->num_loci; i++) {
        if (i != bpm->locus_entries[i].index - 1)
            error("Manifest format error: read invalid number of assay entries\n");
    }
    if (bpm->locus_entries[0].version < 8)
        fprintf(stderr, "Warning: RefStrand annotation missing from manifest file %s\n", bpm->fn);

    read_bytes(bpm->hfile, (void *)&bpm->m_header, sizeof(int32_t), md5);
    bpm->header = (char **)malloc(bpm->m_header * sizeof(char *));
    for (i = 0; i < bpm->m_header; i++) read_pfx_string(bpm->hfile, &bpm->header[i], NULL, md5);

    if (!heof(bpm->hfile)) {
        if (eof_check)
            error(
                "BPM reader did not reach the end of file %s at position %ld\nUse --do-not-check-eof to suppress this "
                "check\n",
                bpm->fn, htell(bpm->hfile));
        if (checksum)
            while (md5_hgetc(bpm->hfile, md5) != EOF);
    }

    if (md5) {
        hts_md5_final(bpm->md5_buf, md5);
        hts_md5_destroy(md5);
    }

    return bpm;
}

static void bpm_destroy(bpm_t *bpm) {
    if (!bpm) return;
    int i;
    if (bpm->hfile && hclose(bpm->hfile) < 0) error("Error closing BPM file %s\n", bpm->fn);
    free(bpm->fn);
    if (bpm->fp && hts_close(bpm->fp) < 0) error("Error closing CSV file %s\n", bpm->fp->fn);
    free(bpm->manifest_name);
    free(bpm->control_config);
    free(bpm->indexes);
    if (bpm->names) {
        for (i = 0; i < bpm->num_loci; i++) free(bpm->names[i]);
        free(bpm->names);
    }
    khash_str2int_destroy(bpm->names2index);
    free(bpm->norm_ids);
    for (i = 0; i < bpm->num_loci; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        free(locus_entry->ilmn_id);
        free(locus_entry->name);
        free(locus_entry->ilmn_strand);
        free(locus_entry->snp);
        free(locus_entry->chrom);
        free(locus_entry->ploidy);
        free(locus_entry->species);
        free(locus_entry->map_info);
        free(locus_entry->customer_strand);
        free(locus_entry->allele_a_probe_seq);
        free(locus_entry->allele_b_probe_seq);
        free(locus_entry->genome_build);
        free(locus_entry->source);
        free(locus_entry->source_version);
        free(locus_entry->source_strand);
        free(locus_entry->source_seq);
        free(locus_entry->top_genomic_seq);
        free(locus_entry->ref_strand);
    }
    free(bpm->locus_entries);
    free(bpm->norm_lookups);
    for (i = 0; i < bpm->m_header; i++) free(bpm->header[i]);
    free(bpm->header);
    free(bpm);
}

/****************************************
 * EGT FILE IMPLEMENTATION              *
 ****************************************/

// http://github.com/broadinstitute/picard/blob/master/src/main/java/picard/arrays/illumina/InfiniumEGTFile.java
// http://github.com/Illumina/BeadArrayFiles/blob/develop/module/ClusterFile.py

typedef struct {
    int32_t N;        // Number of samples assigned to cluster during training
    float r_dev;      // R (intensity) std deviation value
    float r_mean;     // R (intensity) mean value
    float theta_dev;  // Theta std devation value
    float theta_mean; // Theta mean value
} ClusterStats;

typedef struct {
    float cluster_separation; // A score measure the separation between genotype clusters
    float total_score;        // The GenTrain score
    float original_score;     // The original score before editing this cluster
    uint8_t edited;           // Whether this cluster has been manually manipulated
} ClusterScore;

typedef struct {
    ClusterStats aa_cluster_stats; // Describes AA genotype cluster
    ClusterStats ab_cluster_stats; // Describes AB genotype cluster
    ClusterStats bb_cluster_stats; // Describes BB genotype cluster
    float intensity_threshold;     // Intensity threshold for no-call
    ClusterScore cluster_score;    // Various scores for cluster
    int32_t address;               // Bead type identifier for probe A
    float r_mean;                  // precomputed clusters mean
} ClusterRecord;

typedef struct {
    char *fn;
    hFILE *hfile;
    int32_t version;
    char *gencall_version;       // The GenCall version
    char *cluster_version;       // The clustering algorithm version
    char *call_version;          // The genotyping algorithm version
    char *normalization_version; // The normalization algorithm version
    char *date_created;          // The date the cluster file was created (e.g., 3/9/2017 2:18:30 PM)
    uint8_t is_wgt;
    int32_t data_block_version;
    char *opa;
    char *manifest_name; // The manifest name used to build this cluster file
    int32_t num_records;
    ClusterRecord *cluster_records;
    char **names; // Names of records from manifest
    void *names2index;
    char unsigned md5_buf[16];
} egt_t;

static void clusterscore_read(ClusterScore *clusterscore, hFILE *hfile, hts_md5_context *md5) {
    read_bytes(hfile, (void *)&clusterscore->cluster_separation, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterscore->total_score, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterscore->original_score, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterscore->edited, sizeof(uint8_t), md5);
}

static void clusterrecord_read(ClusterRecord *clusterrecord, hFILE *hfile, int32_t data_block_version,
                               hts_md5_context *md5) {
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.N, sizeof(int32_t), md5);
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.N, sizeof(int32_t), md5);
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.N, sizeof(int32_t), md5);
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.r_dev, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.r_dev, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.r_dev, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.r_mean, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.r_mean, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.r_mean, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.theta_dev, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.theta_dev, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.theta_dev, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.theta_mean, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.theta_mean, sizeof(float), md5);
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.theta_mean, sizeof(float), md5);
    if (data_block_version >= 7) {
        read_bytes(hfile, (void *)&clusterrecord->intensity_threshold, sizeof(float), md5);
        read_bytes(hfile, NULL, 14 * sizeof(float), md5);
    } else {
        clusterrecord->intensity_threshold = NAN;
    }
}

static egt_t *egt_init(const char *fn, int eof_check, int checksum) {
    int i;
    egt_t *egt = (egt_t *)calloc(1, sizeof(egt_t));
    egt->fn = strdup(fn);
    egt->hfile = hopen(egt->fn, "rb");
    if (egt->hfile == NULL) error("Could not open %s: %s\n", egt->fn, strerror(errno));
    if (is_gzip(egt->hfile)) error("File %s is gzip compressed and currently cannot be sought\n", egt->fn);

    hts_md5_context *md5 = checksum ? hts_md5_init() : NULL;

    read_bytes(egt->hfile, (void *)&egt->version, sizeof(int32_t), md5);
    if (egt->version != 3) error("EGT cluster file version %d not supported\n", egt->version);

    read_pfx_string(egt->hfile, &egt->gencall_version, NULL, md5);
    read_pfx_string(egt->hfile, &egt->cluster_version, NULL, md5);
    read_pfx_string(egt->hfile, &egt->call_version, NULL, md5);
    read_pfx_string(egt->hfile, &egt->normalization_version, NULL, md5);
    read_pfx_string(egt->hfile, &egt->date_created, NULL, md5);

    read_bytes(egt->hfile, (void *)&egt->is_wgt, sizeof(uint8_t), md5);
    if (egt->is_wgt != 1) error("Only WGT cluster file version supported\n");

    read_pfx_string(egt->hfile, &egt->manifest_name, NULL, md5);

    read_bytes(egt->hfile, (void *)&egt->data_block_version, sizeof(int32_t), md5);
    if (egt->data_block_version < 5 || egt->data_block_version == 6 || egt->data_block_version > 9)
        error("Data block version %d in cluster file not supported\n", egt->data_block_version);
    read_pfx_string(egt->hfile, &egt->opa, NULL, md5);

    read_bytes(egt->hfile, (void *)&egt->num_records, sizeof(int32_t), md5);
    egt->cluster_records = (ClusterRecord *)malloc(egt->num_records * sizeof(ClusterRecord));
    for (i = 0; i < egt->num_records; i++)
        clusterrecord_read(&egt->cluster_records[i], egt->hfile, egt->data_block_version, md5);
    for (i = 0; i < egt->num_records; i++) clusterscore_read(&egt->cluster_records[i].cluster_score, egt->hfile, md5);

    // toss useless strings such as aa_ab_bb/aa_ab/aa_bb/ab_bb
    for (i = 0; i < egt->num_records; i++) read_pfx_string(egt->hfile, NULL, NULL, md5);

    egt->names = (char **)malloc(egt->num_records * sizeof(char *));
    egt->names2index = khash_str2int_init();
    for (i = 0; i < egt->num_records; i++) {
        read_pfx_string(egt->hfile, &egt->names[i], NULL, md5);
        if (khash_str2int_has_key(egt->names2index, egt->names[i]))
            error("Illumina probe %s present multiple times in file %s\n", egt->names[i], fn);
        khash_str2int_inc(egt->names2index, egt->names[i]);
    }
    for (i = 0; i < egt->num_records; i++)
        read_bytes(egt->hfile, (void *)&egt->cluster_records[i].address, sizeof(int32_t), md5);

    int32_t aa_n, ab_n, bb_n;
    for (i = 0; i < egt->num_records; i++) {
        read_bytes(egt->hfile, (void *)&aa_n, sizeof(int32_t), md5);
        read_bytes(egt->hfile, (void *)&ab_n, sizeof(int32_t), md5);
        read_bytes(egt->hfile, (void *)&bb_n, sizeof(int32_t), md5);
        if (egt->cluster_records[i].aa_cluster_stats.N != aa_n || egt->cluster_records[i].ab_cluster_stats.N != ab_n
            || egt->cluster_records[i].bb_cluster_stats.N != bb_n)
            error("Cluster counts don't match with EGT cluster file %s\n", egt->fn);
    }

    if (egt->data_block_version == 9) read_bytes(egt->hfile, NULL, egt->num_records * sizeof(float), md5);
    if (eof_check && !heof(egt->hfile))
        error(
            "EGT reader did not reach the end of file %s at position %ld\nUse --do-not-check-eof to suppress this "
            "check\n",
            egt->fn, htell(egt->hfile));
    if (!heof(egt->hfile)) {
        if (eof_check)
            error(
                "EGT reader did not reach the end of file %s at position %ld\nUse --do-not-check-eof to suppress this "
                "check\n",
                egt->fn, htell(egt->hfile));
        if (checksum)
            while (md5_hgetc(egt->hfile, md5) != EOF);
    }

    if (md5) {
        hts_md5_final(egt->md5_buf, md5);
        hts_md5_destroy(md5);
    }

    for (i = 0; i < egt->num_records; i++) {
        ClusterStats *aa = &egt->cluster_records[i].aa_cluster_stats;
        ClusterStats *ab = &egt->cluster_records[i].ab_cluster_stats;
        ClusterStats *bb = &egt->cluster_records[i].bb_cluster_stats;
        egt->cluster_records[i].r_mean =
            (aa->N * aa->r_mean + ab->N * ab->r_mean + bb->N * bb->r_mean) / (aa->N + ab->N + bb->N);
    }
    return egt;
}

static void egt_destroy(egt_t *egt) {
    if (!egt) return;
    int i;
    if (hclose(egt->hfile) < 0) error("Error closing EGT file %s\n", egt->fn);
    free(egt->fn);
    free(egt->gencall_version);
    free(egt->cluster_version);
    free(egt->call_version);
    free(egt->normalization_version);
    free(egt->date_created);
    free(egt->opa);
    free(egt->manifest_name);
    free(egt->cluster_records);
    for (i = 0; i < egt->num_records; i++) free(egt->names[i]);
    free(egt->names);
    khash_str2int_destroy(egt->names2index);
    free(egt);
}

// static void egt_to_csv(const egt_t *egt, FILE *stream, int verbose) {
//     fprintf(stream, "Illumina, Inc.\n");
//     fprintf(stream, "[Heading]\n");
//     fprintf(stream, "Descriptor File Name,%s\n", strrchr(egt->fn, '/') ? strrchr(egt->fn, '/') + 1 : egt->fn);
//     fprintf(stream, "GenCall version,%s\n", egt->gencall_version);
//     fprintf(stream, "Clustering algorithm version,%s\n", egt->cluster_version);
//     fprintf(stream, "Genotyping algorithm version,%s\n", egt->call_version);
//     fprintf(stream, "Normalization algorithm version,%s\n", egt->normalization_version);
//     fprintf(stream, "Date Manufactured,%s\n", egt->date_created);
//     fprintf(stream, "Manifest name used to build this cluster file,%s\n", egt->manifest_name);
//     fprintf(stream, "OPA,%s\n", egt->opa ? egt->opa : "");
//     fprintf(stream, "Loci Count,%d\n", egt->num_records);
//     fprintf(stream, "[Assay]\n");
//     fprintf(stream,
//             "Name,AA.N,AA.R_dev,AA.R_mean,AA.Theta_dev,AA.Theta_mean,AB.N,AB.R_dev,AB.R_mean,AB."
//             "Theta_dev,AB.Theta_mean,BB.N,BB.R_dev,BB.R_mean,BB.Theta_dev,BB.Theta_mean,Intensity "
//             "Threshold,Cluster Separation,GenTrain Score,Original Score,Edited,Address\n");
//     if (verbose) {
//         int i;
//         for (i = 0; i < egt->num_records; i++) {
//             ClusterRecord *cluster_record = &egt->cluster_records[i];
//             fprintf(stream, "%s,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\n", egt->names[i],
//                     cluster_record->aa_cluster_stats.N, cluster_record->aa_cluster_stats.r_dev,
//                     cluster_record->aa_cluster_stats.r_mean, cluster_record->aa_cluster_stats.theta_dev,
//                     cluster_record->aa_cluster_stats.theta_mean, cluster_record->ab_cluster_stats.N,
//                     cluster_record->ab_cluster_stats.r_dev, cluster_record->ab_cluster_stats.r_mean,
//                     cluster_record->ab_cluster_stats.theta_dev, cluster_record->ab_cluster_stats.theta_mean,
//                     cluster_record->bb_cluster_stats.N, cluster_record->bb_cluster_stats.r_dev,
//                     cluster_record->bb_cluster_stats.r_mean, cluster_record->bb_cluster_stats.theta_dev,
//                     cluster_record->bb_cluster_stats.theta_mean, cluster_record->intensity_threshold,
//                     cluster_record->cluster_score.cluster_separation, cluster_record->cluster_score.total_score,
//                     cluster_record->cluster_score.original_score, cluster_record->cluster_score.edited,
//                     cluster_record->address);
//         }
//     } else {
//         fprintf(stream, "... use --verbose to visualize Assay data ...\n");
//     }
// }

/****************************************
 * MATLAB ROBUST FIT ROUTINES           *
 ****************************************/

// the code for these routines was derived from the Statistics and Machine Learning Toolbox in MATLAB
// Illumina implemented the whole robustfit() function despite the fact that they only needed
// the one dimensional case of it that could easily do away with matrices
// the implementation here reimplements one dimensional linear regression to solve the linear least squares problem
// and so doing away with the need for matrix routines for computing the QR matrix factorization
// the original implementation of robustfit() from Tom Lane was the one adopted in GenTrain 2.0:
// http://github.com/iarsenal95/computer_vision/blob/master/final_project/MATLAB/boosting/weightedstats/private/statrobustfit.m
// in 2002 Tom Lane realized that the MATLAB implementation of the madsigma() sub-routine was problematic:
// http://groups.google.com/g/comp.soft-sys.matlab/c/Raf-VYUh9yY/m/gIi16wAR4VQJ
// this must have led to the new version of madsigma() being adopted in GenTrain 3.0:
// http://github.com/stephane-on/Spectral_analysis/blob/master/statrobustfit.m

inline static double sqr(double x) { return x * x; }

inline static float sqrf(float x) { return x * x; }

// equivalent to MATLAB linsolve(x,y)
// http://www.mathworks.com/help/matlab/ref/linsolve.html
static int matlab_linsolve0(int n, const float *x, const float *y, double *m) {
    int i;
    double sumx2 = 0.0;
    double sumxy = 0.0;
    for (i = 0; i < n; i++) {
        sumx2 += sqr((double)x[i]);
        sumxy += (double)x[i] * (double)y[i];
    }
    if (sumx2 == 0) return 1;
    *m = sumxy / sumx2;
    return 0;
}

// equivalent to MATLAB linsolve([ones(n,1), x],y)
// http://www.mathworks.com/help/matlab/ref/linsolve.html
static int matlab_linsolve1(int n, const float *x, const float *y, double *b, double *m) {
    int i;
    double sumx2 = 0.0;
    double sumxy = 0.0;
    double sumx = 0.0;
    double sumy = 0.0;
    for (i = 0; i < n; i++) {
        sumx2 += sqr((double)x[i]);
        sumxy += (double)y[i] * (double)x[i];
        sumx += (double)x[i];
        sumy += (double)y[i];
    }
    double denom = (double)n * sumx2 - sumx * sumx;
    if (denom == 0) return 1;
    *m = (n * sumxy - sumx * sumy) / denom;
    *b = (sumy * sumx2 - sumx * sumxy) / denom;
    return 0;
}

// equivalent to MATLAB wfit(y,x,w) which is equivalent to linsolve(diag(sqrt(w))*x,diag(sqrt(w))*y)
// stats/private/statrobustfit.m
static int matlab_wfit0(int n, const float *y, const float *x, const double *w, double *m) {
    int i;
    double wsumx2 = 0.0;
    double wsumxy = 0.0;
    for (i = 0; i < n; i++) {
        wsumx2 += w[i] * sqr((double)x[i]);
        wsumxy += w[i] * (double)x[i] * (double)y[i];
    }
    if (wsumx2 == 0) return 1;
    *m = wsumxy / wsumx2;
    return 0;
}

// equivalent to MATLAB wfit(y,[ones(n,1),x],w) which is equivalent to
// linsolve(diag(sqrt(w))*[ones(n,1),x],diag(sqrt(w))*y) stats/private/statrobustfit.m
static int matlab_wfit1(int n, const float *y, const float *x, const double *w, double *b, double *m) {
    int i;
    double wsumx2 = 0.0;
    double wsumxy = 0.0;
    double wsumx = 0.0;
    double wsumy = 0.0;
    double wsum = 0.0;
    for (i = 0; i < n; i++) {
        wsumx2 += w[i] * sqr((double)x[i]);
        wsumxy += w[i] * (double)x[i] * (double)y[i];
        wsumx += w[i] * (double)x[i];
        wsumy += w[i] * (double)y[i];
        wsum += w[i];
    }
    double denom = wsum * wsumx2 - wsumx * wsumx;
    if (denom == 0) return 1;
    *m = (wsum * wsumxy - wsumx * wsumy) / denom;
    *b = (wsumy * wsumx2 - wsumx * wsumxy) / denom;
    return 0;
}

// http://www.mathworks.com/help/stats/nanmean.html
static float matlab_nanmean(int n, const float *vals) {
    if (n == 0) return NAN;
    int i, j;
    double sum = 0.0;
    for (i = 0, j = 0; i < n; i++) {
        if (!isnan(vals[i])) {
            sum += vals[i];
            j++;
        }
    }
    return (float)(sum / (double)j);
}

// http://www.mathworks.com/help/matlab/ref/mean.html
static float matlab_mean(int n, const float *vals) {
    if (n == 0) return NAN;
    int i;
    double sum = 0.0;
    for (i = 0; i < n; i++) sum += vals[i];
    return (float)(sum / (double)n);
}

// the input array does not need to be sorted
// http://www.mathworks.com/help/matlab/ref/median.html
static float matlab_median(int n, float *vals) {
    if (n == 0) return 0.0f;
    ks_introsort_float((size_t)n, vals);
    if (n % 2 == 1) return vals[n / 2];
    return (vals[n / 2 - 1] + vals[n / 2]) * 0.5f;
}

// stats/private/statrobustfit.m
// function s = madsigma(r,p)
// %MADSIGMA    Compute sigma estimate using MAD of residuals from 0
// rs = sort(abs(r));
// s = median(rs(max(1,p):end)) / 0.6745; % 0.6745 ~ qnorm(0.75)
static double matlab_madsigma_new(int n, const double *r, int p) {
    int i;
    float *rs = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) rs[i] = (float)fabs(r[i]);
    ks_introsort_float((size_t)n, rs);
    double s = (double)matlab_median(n - (p - 1), rs + (p - 1)) / 0.6745;
    if (s == 0.0) s = 0.5 * (double)matlab_mean(n, rs);
    free(rs);
    return s;
}

// a separate implementation from Illumina can be found in function madsigma in file Utils.cs
// the code follows the original implementation from Tom Lane in 2000
// stats/private/statrobustfit.m
// function s = madsigma(r,p);
// %MADSIGMA    Compute sigma estimate using MAD of residuals
// m = median(r);
// rs = sort(abs(r-m));
// if (abs(m) > rs(end))
//     % Unexpectedly all residuals are very small
//     rs = sort(abs(r));
// end
// s = median(rs(p:end)) / 0.6745; % 0.6745 ~ qnorm(0.75)
// if (s==0), s = .5*mean(rs); end
static double matlab_madsigma_old(int n, const double *r, int p) {
    int i;
    float *rs = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) rs[i] = (float)r[i];
    float m = matlab_median((size_t)n, rs);
    for (i = 0; i < n; i++) rs[i] = fabsf(rs[i] - m);
    ks_introsort_float((size_t)n, rs);
    if (fabsf(m) > rs[n - 1]) {
        for (i = 0; i < n; i++) rs[i] = fabsf((float)r[i]);
        ks_introsort_float((size_t)n, rs);
    }
    double s = (double)matlab_median(n - (p - 1), rs + (p - 1)) / 0.6745;
    if (s == 0.0) s = 0.5 * (double)matlab_mean(n, rs);
    free(rs);
    return s;
}

// roughly equivalent to MATLAB robustfit(x,y,'bisquare',4.685,'off')
// http://www.mathworks.com/help/stats/robustfit.html
// stats/private/statrobustfit.m
// stats/private/statrobustwfun.m
static void matlab_robustfit0(int n, const float *x, const float *y, double (*madsigma)(int, const double *, int),
                              float *out_m) {
    int i;
    double *r = (double *)malloc(n * sizeof(double));
    double *w = (double *)malloc(n * sizeof(double));
    double m, m0 = 0.0;
    if (matlab_linsolve0(n, x, y, &m)) error("Error while running linsolve0\n");
    // [Q,R] = qr(x,0);
    // R = [sqrt(sum(x.^2)]
    // E = X/R = [x/sqrt(sum(x.^2)]
    // h = min(.9999, sum(E.*E,2)) = min(.9999, x.^2 / sum(x.^2))
    // adjfactor = 1 ./ sqrt(1-h)
    // as GenCall messed up the implementation, here we use instead
    // h = min(.9999, sum(E.*E)) = min(.9999, sum(x.^2) / sum(x.^2)) = .9999
    // adjfactor = 1 / sqrt(1 - 0.9999) ~ 100;
    double adjfactor = 100.0 + 24832 * DBL_EPSILON;

    int iter = 0;
    do {
        // as Illumina messed up the implementation, here we use adjfactor instead of adjfactor[i]
        for (i = 0; i < n; i++) r[i] = ((double)y[i] - m * (double)x[i]) * adjfactor;
        double s = madsigma(n, r, 1);
        if (s == 0.0) s = 1.0;
        for (i = 0; i < n; i++) {
            r[i] *= 1.0 / (s * 4.685);
            w[i] = fabs(r[i]) < 1 ? sqr(1.0 - sqr(r[i])) : 0.0;
        }
        m0 = m;
        if (matlab_wfit0(n, y, x, w, &m)) error("Error while running wfit0\n");
        iter++;
    } while (iter < 50 && fabs(m - m0) > 1e-6 * (double)fmaxf((float)fabs(m), (float)fabs(m0)));

    free(r);
    free(w);
    *out_m = (float)m;
}

// roughly equivalent to MATLAB robustfit(x,y,'bisquare',4.685,'on')
// http://www.mathworks.com/help/stats/robustfit.html
// stats/private/statrobustfit.m
// stats/private/statrobustwfun.m
static void matlab_robustfit1(int n, const float *x, const float *y, double (*madsigma)(int, const double *, int),
                              float *out_b, float *out_m) {
#ifdef __GNUC__
    if (n <= 0) __builtin_unreachable(); // to prevent an unnecessary "may be used uninitialized" warning
#endif
    int i;
    double *adjfactor = (double *)malloc(n * sizeof(double));
    double *r = (double *)malloc(n * sizeof(double));
    double *w = (double *)malloc(n * sizeof(double));
    double b, m, b0 = 0.0, m0 = 0.0;
    if (matlab_linsolve1(n, x, y, &b, &m))
        error(
            "Error while running linsolve1\nFailed to normalize and gencall\nThis typically happens when the wrong "
            "manifest file is used\n");
    // [Q,R] = qr([ones(n,1),x],0);
    // R = [-sqrt(n), -sum(x)/sqrt(n); 0, sqrt(sum(x.^2)-sum(x)^2/n)]
    // E = X/R = [-ones(n,1)/sqrt(n), (sum(x)/n-x)/sqrt(sum(x.^2)-sum(x)^2/n)]
    // h = min(.9999, sum(E.*E,2)) = min(.9999, (n*x.^2 - 2*sum(x)*x + sum(x.^2))/(n*sum(x.^2) - sum(x)^2))
    double sumx = 0.0;
    double sumx2 = 0.0;
    for (i = 0; i < n; i++) {
        sumx += (double)x[i];
        sumx2 += sqr((double)x[i]);
    }
    double denom = (double)n * sumx2 - sqr(sumx);
    for (i = 0; i < n; i++) {
        double h = fmin(.9999, ((double)n * sqr((double)x[i]) - 2.0 * sumx * (double)x[i] + sumx2) / denom);
        adjfactor[i] = 1.0 / sqrt(1.0 - h);
    }

    int iter = 0;
    do {
        for (i = 0; i < n; i++) r[i] = ((double)y[i] - b - m * (double)x[i]) * adjfactor[i];
        double s = madsigma(n, r, 2);
        if (s == 0.0) s = 1.0;
        for (i = 0; i < n; i++) {
            r[i] *= 1.0 / (s * 4.685);
            w[i] = fabs(r[i]) < 1 ? sqr(1.0 - sqr(r[i])) : 0.0;
        }
        b0 = b;
        m0 = m;
        if (matlab_wfit1(n, y, x, w, &b, &m)) error("Error while running wfit1\n");
        iter++;
    } while (iter < 50
             && (fabs(b - b0) > 1e-6 * (double)fmaxf((float)fabs(b), (float)fabs(b0))
                 || fabs(m - m0) > 1e-6 * (double)fmaxf((float)fabs(m), (float)fabs(m0))));

    free(adjfactor);
    free(r);
    free(w);
    *out_b = (float)b;
    *out_m = (float)m;
}

/****************************************
 * NEAREST NEIGHBOR ROUTINES            *
 ****************************************/

// a separate implementation from Illumina of these functions in GenCall can be found in file Utils.cs
// It seems like Illumina at first used a function with O(n^2) complexity for the same task and then when they switched
// from GoldenGate to larger Infinium arrays this solution did not scale anymore. This led to a reimplementation in C as
// the C# version was not fast enough. For this reason AutoConvert, an almost entirely C# executable, requires this
// specific function as unmanaged C code while IAAP and ACLI have their equivalent version in C#, maybe because by then
// computers had become fast enough

int elementsInBin[12];
int *binData[12];
int elementsInShiftedBin[11];
int *binDataShifted[11];

// a separate implementation from Illumina of this function can be found in function ClosestPointsB
int findClosestSitesToPointsAlongAxis(int n_raw, float *raw_x, float *raw_y, int n_axis, float *axis_x, float *axis_y,
                                      int *ret) {
    int i;
    float *raw_a = NULL;
    float *raw_b = NULL;
    float *axis_a = NULL;
    float axis_max_val;
    float bin_width;
    int bin_idx;
    float quotient;
    float reminder;
    int *curr_bin_data;
    int curr_bin_size;
    float curr_axis_x;
    float curr_axis_y;
    float x_dist;
    float y_dist;
    double best_val;
    int best_idx;
    int j;
    int curr_idx;
    double sq_dist;
    double axis_max_dist;
    int use_y = 1;
    int use_x = 1;

    for (i = 0; i < n_axis; i++) {
        if (axis_x[i] > 0.0001) {
            use_y = 0;
            break;
        }
    }

    for (i = 0; i < n_axis; i++) {
        if (axis_y[i] > 0.0001) {
            use_x = 0;
            break;
        }
    }

    if (use_y) {
        raw_a = raw_y;
        raw_b = raw_x;
        axis_a = axis_y;
    } else if (use_x) {
        raw_a = raw_x;
        raw_b = raw_y;
        axis_a = axis_x;
    } else {
        return -1;
    }

    axis_max_val = axis_a[n_axis - 1];
    bin_width = axis_max_val / 12.0f;
    axis_max_dist = (double)bin_width;

    for (i = 0; i < n_raw; i++) {
        if ((double)raw_b[i] > axis_max_dist) continue;
        bin_idx = (int)(raw_a[i] / bin_width);
        if (bin_idx < 0) bin_idx = 0;
        if (bin_idx > 11) bin_idx = 11;
        elementsInBin[bin_idx]++;
        bin_idx = (int)(raw_a[i] / bin_width - 0.5f);
        if (bin_idx < 0) bin_idx = 0;
        if (bin_idx > 10) bin_idx = 10;
        elementsInShiftedBin[bin_idx]++;
    }

    for (i = 0; i <= 11; i++) {
        binData[i] = (int *)malloc((size_t)elementsInBin[i] * sizeof(int));
        elementsInBin[i] = 0;
        if (i == 11) continue;
        binDataShifted[i] = (int *)malloc((size_t)elementsInShiftedBin[i] * sizeof(int));
        elementsInShiftedBin[i] = 0;
    }

    for (i = 0; i < n_raw; i++) {
        if ((double)raw_b[i] > axis_max_dist) continue;
        bin_idx = (int)(raw_a[i] / bin_width);
        if (bin_idx < 0) bin_idx = 0;
        if (bin_idx > 11) bin_idx = 11;
        binData[bin_idx][elementsInBin[bin_idx]] = i;
        elementsInBin[bin_idx]++;
        bin_idx = (int)(raw_a[i] / bin_width - 0.5f);
        if (bin_idx < 0) bin_idx = 0;
        if (bin_idx > 10) bin_idx = 10;
        binDataShifted[bin_idx][elementsInShiftedBin[bin_idx]] = i;
        elementsInShiftedBin[bin_idx]++;
    }

    for (i = 0; i < n_axis; i++) {
        quotient = axis_a[i] / bin_width;
        bin_idx = (int)quotient;
        reminder = quotient - (float)bin_idx;
        curr_bin_data = NULL;
        curr_bin_size = 0;
        if (bin_idx < 0) bin_idx = 0;
        if (bin_idx > 11) bin_idx = 11;

        if (0.25f <= reminder && reminder <= 0.75f) {
            curr_bin_data = binData[bin_idx];
            curr_bin_size = elementsInBin[bin_idx];
        } else {
            if (reminder < 0.25f) {
                if (bin_idx == 0) {
                    curr_bin_data = binData[bin_idx];
                    curr_bin_size = elementsInBin[bin_idx];
                } else {
                    curr_bin_data = binDataShifted[bin_idx - 1];
                    curr_bin_size = elementsInShiftedBin[bin_idx - 1];
                }
            } else if (bin_idx == 11) {
                curr_bin_data = binData[bin_idx];
                curr_bin_size = elementsInBin[bin_idx];
            } else {
                curr_bin_data = binDataShifted[bin_idx];
                curr_bin_size = elementsInShiftedBin[bin_idx];
            }
        }

        curr_axis_x = axis_x[i];
        curr_axis_y = axis_y[i];
        best_val = 1e20;
        best_idx = -1;

        for (j = 0; j < curr_bin_size; j++) {
            curr_idx = curr_bin_data[j];
            x_dist = raw_x[curr_idx] - curr_axis_x;
            y_dist = raw_y[curr_idx] - curr_axis_y;
            sq_dist = (double)(x_dist * x_dist + y_dist * y_dist);
            if (sq_dist < best_val) {
                best_val = sq_dist;
                best_idx = curr_idx;
            }
        }

        ret[i] = best_idx;
    }

    for (i = 0; i <= 11; i++) {
        free((void *)binData[i]);
        elementsInBin[i] = 0;
        if (i > 10) continue;
        free((void *)binDataShifted[i]);
        elementsInShiftedBin[i] = 0;
    }

    return 0;
}

// a separate implementation from Illumina of this function can be found in function ClosestPointsSlow
// as explained in the patent, this approach is slow as it runs in O(n^2)
static int *closest_points_slow(int nref, const float *xref, const float *yref, int n, float *x, float *y) {
    int i, j, *closest_sites = (int *)malloc(n * sizeof(int));
    for (i = 0; i < n; i++) {
        float xv = x[i];
        float yv = y[i];
        double mindist = (xv - xref[0]) * (xv - xref[0]) + (yv - yref[0]) * (yv - yref[0]);
        int mini = 0;
        for (j = 1; j < nref; j++) {
            double dist = (xv - xref[j]) * (xv - xref[j]) + (yv - yref[j]) * (yv - yref[j]);
            if (dist < mindist) {
                mindist = dist;
                mini = j;
            }
        }
        closest_sites[i] = mini;
    }
    return closest_sites;
}

#define SAMPLE 2000 // mentioned in file Utils.cs

// a separate implementation from Illumina of this function can be found in function ClosestPoints
static int *closest_points(int nref, float *xref, float *yref, int n, float *x, float *y) {
    if (nref < SAMPLE) return closest_points_slow(nref, xref, yref, n, x, y);
    int *closest_sites = (int *)malloc(n * sizeof(int));
    findClosestSitesToPointsAlongAxis(nref, xref, yref, n, x, y, closest_sites);
    return closest_sites;
}

/****************************************
 * MATLAB UTILS ROUTINES                *
 ****************************************/

// a separate implementation from Illumina of these functions in GenCall can be found in file Utils.cs

// the input array does need to be sorted
static float percentile(int n, const float *vals, int percentile) {
    if (n == 0) return NAN;
    int i1 = n * percentile / 100;
    float f = (float)(n * percentile) / 100.0f - (float)i1;
    if (f < 0.5f) {
        i1--;
    }
    if (i1 < 0) {
        return vals[0];
    }
    if (i1 >= n - 1) {
        return vals[n - 1];
    }

    float x1 = 100.0f * ((float)i1 + 0.5f) / (float)n;
    float x2 = 100.0f * ((float)(i1 + 1) + 0.5f) / (float)n;
    float y1 = (float)vals[i1];
    float y2 = (float)vals[i1 + 1];
    float m = (y2 - y1) / (x2 - x1);
    return y1 + m * ((float)percentile - x1);
}

// http://www.mathworks.com/help/matlab/ref/iqr.html
static float matlab_iqr(int n, const float *vals) {
    int i;
    float *vs = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) vs[i] = vals[i];
    ks_introsort_float((size_t)n, vs);
    float iqr = percentile(n, vs, 75) - percentile(n, vs, 25);
    free(vs);
    return iqr;
}

// http://www.mathworks.com/help/stats/trimmean.html
static float matlab_trimmean(int n, float *vals, int percent) {
    ks_introsort_float((size_t)n, vals);
    float high = percentile(n, vals, 100 - percent / 2);
    float low = percentile(n, vals, percent / 2);
    double sum = 0.0;
    int i, count = 0;
    for (i = 0; i < n; i++) {
        if (vals[i] >= low && vals[i] <= high) {
            sum += (double)vals[i];
            count++;
        }
    }
    return (float)sum / (float)count;
}

// http://www.mathworks.com/help/matlab/ref/linspace.html
static float *matlab_linspace(int n, float minv, float maxv) {
    int i;
    float *vals = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) vals[i] = minv + (maxv - minv) * (float)i / (float)(n - 1);
    return vals;
}

// the input array does not need to be sorted
// http://www.mathworks.com/help/matlab/ref/unique.html
static int matlab_unique(int n, int *indices) {
    int i, j;
    ks_introsort_int((size_t)n, indices);
    for (i = 0; indices[i] == -1; i++);
    indices[0] = indices[i++];
    for (j = 1; i < n; i++)
        if (indices[i] != indices[i - 1]) indices[j++] = indices[i];
    return j;
}

// the input arrays do not need to be sorted
// http://www.mathworks.com/help/matlab/ref/union.html
static int *matlab_union(int na, const int *a, int nb, const int *b, int *n) {
    int i, *c = (int *)malloc((na + nb) * sizeof(int));
    for (i = 0; i < na; i++) c[i] = a[i];
    for (i = 0; i < nb; i++) c[i + na] = b[i];
    *n = matlab_unique(na + nb, c);
    return c;
}

// http://www.mathworks.com/help/matlab/ref/min.html
static float matlab_min(int n, const float *vals) {
    int i;
    float minval = FLT_MAX;
    for (i = 0; i < n; i++) {
        if (isnan(vals[i])) continue;
        if (vals[i] < minval) minval = vals[i];
    }
    return minval;
}

// http://www.mathworks.com/help/matlab/ref/max.html
static float matlab_max(int n, const float *vals) {
    int i;
    float maxval = -FLT_MAX;
    for (i = 0; i < n; i++) {
        if (isnan(vals[i])) continue;
        if (vals[i] > maxval) maxval = vals[i];
    }
    return maxval;
}

/****************************************
 * NORMALIZATION ROUTINES               *
 ****************************************/

// a thorough explanation of the normalization steps can be found in the document
// Kermani, B. G. Artificial intelligence and global normalization methods for genotyping. U.S. Patent No. 7,035,740
// (2005-09-29) http://patents.google.com/patent/US7035740 Peiffer, D. A. et al. High-resolution genomic profiling of
// chromosomal aberrations using Infinium whole-genome genotyping. Genome Res., 16, 1136–1148 (2006-08-09)
// http://doi-org.ezp-prod1.hul.harvard.edu/10.1101/gr.5402306
// Illumina, Inc. Illumina’s Genotyping Data Normalization Methods. Pub. No. 970-2006-010 (2006-09-26)
// http://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina_gt_normalization.pdf
// Illumina, Inc. Improved Genotype Clustering with GenTrain 3.0. Pub. No. 370-2016-015-A (2016)
// http://emea.illumina.com/content/dam/illumina-marketing/documents/products/technotes/gentrain3-technical-note-370-2016-015.pdf
// http://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/gentrain-tech-note-m-gl-01258/gentrain-tech-note-m-gl-01258.pdf
// a separate implementation from Illumina of these functions in GenCall can be found in file NormalizationInfinium.cs
// http://support.illumina.com/downloads/gencall_software.html

// Illumina software includes three normalization protocols for Infinium arrays:
// 1.1.0 Normalization10 not used
// 1.1.2 Normalization111+Normalization10 used in AutoConvert and IAAP Genotyping CLI with option --gentrain-id 2
// 1.2.0 NormalizationDragonfish+Normalization111_Dragonfish+Normalization10_Dragonfish used in AutoConvert 2.0, IAAP
// Genotyping CLI, and Array Analysis CLI we implement version 1.1.2 and 1.2.0 for interoperability purposes with
// existing Illumina cluster files

// Peiffer, D. A. et al. High-resolution genomic profiling of chromosomal aberrations using Infinium whole-genome
// genotyping. Genome Res., 16, 1136–1148 (2006-08-09) The data for each BeadChip is self-normalized using infor- mation
// contained within the array. This normalization algo- rithm removes outliers, adjusts for channel-dependent back-
// ground and global intensity differences, and also scales the data.
// The X and Y color channels undergo an affine coordinate trans-
// formation to make the data appear as canonical as possible with
// the homozygotes lying along the transformed x- and y-axes. The
// following five steps are applied: (1) outlier removal; (2) a trans-
// lation correction in which the asymptotes are fitted to candidate
// AA and BB homozygotes; the intersection of these fit lines de-
// fines the translated origin; (3) rotational correction: the angle of
// the AA homozygote asymptote with respect to the translated
// X-axis is used to define the rotational correction; (4) shear cor-
// rection: the angle of the BB homozygote asymptote with respect
// to the translated and rotated y-axis is used to define the shear
// correction; (5) scaling correction: statistical centroids are com-
// puted for the candidate AA homozygotes to define an x-axis scal-
// ing parameter, and for candidate BB homozygotes to define a
// y-axis scaling parameter. The translated, rotated, shear-corrected
// data are normalized to a scale of ∼1 using the scaling parameters

#define SAMPLING 400 // mentioned in Illumina’s Genotyping Data Normalization Methods. Pub. No. 970-2006-010 (2006)
#define ROBUST_THRESHOLD                                                                                               \
    192 // mentioned in Improved Genotype Clustering with GenTrain 3.0. Pub. No. 370-2016-015-A (2016)

// a separate implementation from Illumina can be found in functions RemoveOutliers from classes
// Normalization10_Dragonfish and Normalization10 Illumina, Inc. Illumina’s Genotyping Data Normalization Methods. Pub.
// No. 970-2006-010 (2006-09-26) Outlier SNPs are removed from consideration during normalization parameter estimation.
// These SNPs are only considered outliers during the normalization process and are not excluded from downstream
// analysis. A SNP is considered an outlier if its intensity meets any of the following criteria:
// - Its value of x, y, or x/(x+y) is smaller than either the 5th smallest or the 1st percentile (whichever is smaller)
// of those values across all SNPs.
// - Its value of x, y, or x/(x+y) is larger than either the 5th largest or the 99 th percentile (whichever is larger)
// of those values across all SNPs.
static void remove_outliers(int *n, float *x, float *y) {
    if (*n < SAMPLING) return;
    int i, j;
    float *xs = (float *)malloc(*n * sizeof(float));
    float *ys = (float *)malloc(*n * sizeof(float));
    float *ts = (float *)malloc(*n * sizeof(float));
    for (i = 0; i < *n; i++) {
        xs[i] = x[i];
        ys[i] = y[i];
        ts[i] = y[i] / (FLT_MIN * FLT_EPSILON + x[i] + y[i]);
    }
    ks_introsort_float((size_t)(*n), xs);
    ks_introsort_float((size_t)(*n), ys);
    ks_introsort_float((size_t)(*n), ts);

    int M = 5;
    int Nb = 1;

    float tcut1a = ts[M - 1];
    float tcut2a = ts[*n - M + 1 - 1];

    float xcut1a = xs[M - 1];
    float xcut2a = xs[*n - M + 1 - 1];

    float ycut1a = ys[M - 1];
    float ycut2a = ys[*n - M + 1 - 1];

    float tcut1b = percentile(*n, ts, Nb);
    float tcut2b = percentile(*n, ts, 100 - Nb);

    float xcut1b = percentile(*n, xs, Nb);
    float xcut2b = percentile(*n, xs, 100 - Nb);

    float ycut1b = percentile(*n, ys, Nb);
    float ycut2b = percentile(*n, ys, 100 - Nb);

    float tcut1 = fminf(tcut1a, tcut1b);
    float tcut2 = fmaxf(tcut2a, tcut2b);

    float xcut1 = fminf(xcut1a, xcut1b);
    float xcut2 = fmaxf(xcut2a, xcut2b);

    float ycut1 = fminf(ycut1a, ycut1b);
    float ycut2 = fmaxf(ycut2a, ycut2b);

    for (i = 0, j = 0; i < *n; i++) {
        if (y[i] <= ycut1 || x[i] <= xcut1 || y[i] >= ycut2 || x[i] >= xcut2) {
            continue;
        }
        double t = y[i] / (double)(y[i] + x[i]);
        if (t <= tcut1 || t >= tcut2) {
            continue;
        }
        x[j] = x[i];
        y[j] = y[i];
        j++;
    }
    *n = j;

    free(xs);
    free(ys);
    free(ts);
}

// a separate implementation from Illumina can be found in function RemoveOffset from class Normalization10_Dragonfish
// and Normalization10 Illumina, Inc. Illumina’s Genotyping Data Normalization Methods. Pub. No. 970-2006-010
// (2006-09-26) a. An x-sweep is performed by sampling 400 points along the x-axis, from the smallest x value to the
// largest. The closest SNP to each sampled point along the axis is added to the set of candidate homozygote As. b. The
// same analysis is performed along the y-axis to find the candidate homozygote Bs. c. A straight line is fit into
// candidate homozygote A alleles. d. A straight line is fit into candidate homozygote B alleles. e. The intercept of
// the two lines is computed, and this coordinate corresponds to offset_x and offset_y.
static void remove_offset(int n, float *x, float *y, int *naa, int **iaa, int *nbb, int **ibb,
                          double (*madsigma)(int, const double *, int), float *offset_x, float *offset_y) {
    if (n < ROBUST_THRESHOLD) {
        *offset_x = 0.0f;
        *offset_y = 0.0f;
        return;
    }

    int i;
    float mx = matlab_min(n, x);
    float my = matlab_min(n, y);
    float *xt = (float *)malloc(n * sizeof(float));
    float *yt = (float *)malloc(n * sizeof(float));

    for (i = 0; i < n; i++) {
        xt[i] = x[i] - mx;
        yt[i] = y[i] - my;
    }
    float *xsweep = matlab_linspace(SAMPLING, 0.0f, matlab_max(n, xt));
    float *ysweep = matlab_linspace(SAMPLING, 0.0f, matlab_max(n, yt));
    float *zeros = (float *)calloc(SAMPLING, sizeof(float));

    *iaa = closest_points(n, xt, yt, SAMPLING, xsweep, zeros);
    *ibb = closest_points(n, xt, yt, SAMPLING, zeros, ysweep);
    *naa = matlab_unique(SAMPLING, *iaa);
    *nbb = matlab_unique(SAMPLING, *ibb);

    float *xaa = (float *)malloc(*naa * sizeof(float));
    float *yaa = (float *)malloc(*naa * sizeof(float));
    for (i = 0; i < *naa; i++) {
        xaa[i] = xt[(*iaa)[i]];
        yaa[i] = yt[(*iaa)[i]];
    }

    float *xbb = (float *)malloc(*nbb * sizeof(float));
    float *ybb = (float *)malloc(*nbb * sizeof(float));
    for (i = 0; i < *nbb; i++) {
        xbb[i] = xt[(*ibb)[i]];
        ybb[i] = yt[(*ibb)[i]];
    }

    float baa, maa;
    float bbb, mbb;
    matlab_robustfit1(*naa, xaa, yaa, madsigma, &baa, &maa);
    matlab_robustfit1(*nbb, ybb, xbb, madsigma, &bbb, &mbb);

    float ox = (bbb + mbb * baa) / (1.0f - mbb * maa);
    float oy = baa + maa * ox;

    *offset_x = ox + mx;
    *offset_y = oy + my;

    for (i = 0; i < n; i++) {
        x[i] -= *offset_x;
        y[i] -= *offset_y;
    }

    free(xt);
    free(yt);
    free(xsweep);
    free(ysweep);
    free(zeros);
    free(xaa);
    free(yaa);
    free(xbb);
    free(ybb);
}

// a separate implementation from Illumina can be found in function HandleRotation from class Normalization10_Dragonfish
// and Normalization10 Illumina, Inc. Illumina’s Genotyping Data Normalization Methods. Pub. No. 970-2006-010
// (2006-09-26) a. The points are corrected for translation and another x-sweep is performed to determine a set of
// control points. b. A straight line is fit into the control points. The angle between this line and the x-axis defines
// the amount of rotation in the data. This angle corresponds to the theta parameter.
static void handle_rotation(int n, float *x, float *y, int *naa, int **iaa, int *nbb, int **ibb,
                            double (*madsigma)(int, const double *, int), float *theta) {
    if (n < ROBUST_THRESHOLD) {
        *theta = 0.0f;
        return;
    }

    int i;
    float *xsweep = matlab_linspace(SAMPLING, matlab_min(n, x), matlab_max(n, x));
    float *ysweep = matlab_linspace(SAMPLING, matlab_min(n, y), matlab_max(n, y));
    float *zeros = (float *)calloc(SAMPLING, sizeof(float));
    int *tiaa = closest_points(n, x, y, SAMPLING, xsweep, zeros);
    int *tibb = closest_points(n, x, y, SAMPLING, zeros, ysweep);

    int naa_in = *naa, *iaa_in = *iaa;
    int nbb_in = *nbb, *ibb_in = *ibb;
    *iaa = matlab_union(naa_in, iaa_in, SAMPLING, tiaa, naa);
    *ibb = matlab_union(nbb_in, ibb_in, SAMPLING, tibb, nbb);

    float *tx = (float *)malloc(*naa * sizeof(float));
    float *ty = (float *)malloc(*naa * sizeof(float));
    for (i = 0; i < *naa; i++) {
        tx[i] = x[(*iaa)[i]];
        ty[i] = y[(*iaa)[i]];
    }

    float m;
    matlab_robustfit0(*naa, tx, ty, madsigma, &m);

    double taa = atan((double)m);
    double ct = cos(taa);
    double st = sin(taa);

    for (i = 0; i < n; i++) {
        float tmp = x[i];
        x[i] = (float)(ct * (double)x[i] + st * (double)y[i]);
        y[i] = (float)((0.0 - st) * (double)tmp + ct * (double)y[i]);
    }

    *theta = (float)taa;

    free(xsweep);
    free(ysweep);
    free(zeros);
    free(iaa_in);
    free(ibb_in);
    free(tiaa);
    free(tibb);
    free(tx);
    free(ty);
}

// a separate implementation from Illumina can be found in function HandleShear from class Normalization10_Dragonfish
// and Normalization10 Illumina, Inc. Illumina’s Genotyping Data Normalization Methods. Pub. No. 970-2006-010
// (2006-09-26) a. The points are corrected for rotation and another y-sweep is performed to determine a set of control
// points. b. A straight line is fit to these control points. The angle of this line identifies the shear parameter
static void handle_shear(int n, float *x, float *y, int *nbb, int **ibb, double (*madsigma)(int, const double *, int),
                         float *shear) {
    if (n < ROBUST_THRESHOLD) {
        *shear = 0.0f;
        return;
    }

    int i;
    float *ysweep = matlab_linspace(SAMPLING, 0.0f, matlab_max(n, y));
    float *zeros = (float *)calloc(SAMPLING, sizeof(float));
    int *tibb = closest_points(n, x, y, SAMPLING, zeros, ysweep);
    int nbb_in = *nbb, *ibb_in = *ibb;
    *ibb = matlab_union(nbb_in, ibb_in, SAMPLING, tibb, nbb);

    float *tx = (float *)malloc(*nbb * sizeof(float));
    float *ty = (float *)malloc(*nbb * sizeof(float));
    for (i = 0; i < *nbb; i++) {
        tx[i] = x[(*ibb)[i]];
        ty[i] = y[(*ibb)[i]];
    }

    float m;
    matlab_robustfit0(*nbb, ty, tx, madsigma, &m);

    double tbb = atan((double)m);
    double shy = tan(tbb);

    for (i = 0; i < n; i++) x[i] = (float)((double)x[i] - shy * (double)y[i]);

    *shear = (float)shy;

    free(ibb_in);
    free(ysweep);
    free(zeros);
    free(tibb);
    free(tx);
    free(ty);
}

// a separate implementation from Illumina can be found in function HandleScale from classes Normalization10_Dragonfish
// and Normalization10 0.7413 ~ 1/(2*qnorm(0.75))
static void base_handle_scale(int n, float *x, float *y, int gentrain_version, float *scale_x, float *scale_y) {
    int i, naa, nbb, *iaa, *ibb;
    // this should never happen
    for (i = 0; i < n; i++) {
        if (x[i] < 0.0f) x[i] = 0.0f;
        if (y[i] < 0.0f) y[i] = 0.0f;
    }

    if (n < ROBUST_THRESHOLD) {
        float *t = (float *)malloc(n * sizeof(float));
        // for GenTrain 2.0 we replicate the bug by allowing failed probes as AA points
        for (i = 0; i < n; i++)
            t[i] = x[i] > 0.0f || y[i] > 0.0f || gentrain_version == 2
                       ? (float)(180.0 * M_1_PI * atan2((double)y[i], (double)x[i]))
                       : NAN;

        naa = 0;
        nbb = 0;
        for (i = 0; i < n; i++) {
            if (t[i] < 10.0f) naa++;
            if (t[i] > 80.0f) nbb++;
        }
        iaa = (int *)malloc(naa * sizeof(int));
        ibb = (int *)malloc(nbb * sizeof(int));
        naa = 0;
        nbb = 0;
        for (i = 0; i < n; i++) {
            if (t[i] < 10.0f) iaa[naa++] = i;
            if (t[i] > 80.0f) ibb[nbb++] = i;
        }

        free(t);
    } else {
        float *xsweep = matlab_linspace(SAMPLING, 0.0f, matlab_max(n, x));
        float *ysweep = matlab_linspace(SAMPLING, 0.0f, matlab_max(n, y));
        float *zeros = (float *)calloc(SAMPLING, sizeof(float));

        iaa = closest_points(n, x, y, SAMPLING, xsweep, zeros);
        ibb = closest_points(n, x, y, SAMPLING, zeros, ysweep);
        naa = matlab_unique(SAMPLING, iaa);
        nbb = matlab_unique(SAMPLING, ibb);

        free(xsweep);
        free(ysweep);
        free(zeros);
    }

    float *xaa = (float *)malloc(naa * sizeof(float));
    float *ybb = (float *)malloc(nbb * sizeof(float));
    for (i = 0; i < naa; i++) xaa[i] = x[iaa[i]];
    for (i = 0; i < nbb; i++) ybb[i] = y[ibb[i]];

    if (n < ROBUST_THRESHOLD) {
        *scale_x = matlab_trimmean(naa, xaa, 20);
        *scale_y = matlab_trimmean(nbb, ybb, 20);
    } else {
        *scale_x = 0.5f * matlab_trimmean(naa, xaa, 50) + 0.7413f * matlab_iqr(naa, xaa);
        *scale_y = 0.5f * matlab_trimmean(nbb, ybb, 50) + 0.7413f * matlab_iqr(nbb, ybb);
    }

    for (i = 0; i < n; i++) {
        x[i] /= *scale_x;
        y[i] /= *scale_y;
    }

    free(iaa);
    free(ibb);
    free(xaa);
    free(ybb);
}

// a separate implementation from Illumina can be found in function HandleScale from class Normalization111_Dragonfish
// and Normalization111 Illumina, Inc. Illumina’s Genotyping Data Normalization Methods. Pub. No. 970-2006-010
// (2006-09-26) a. The points are corrected for shear, and another x-sweep is performed to identify a set of virtual
// points. b. A statistical robust measure of the mean of these control points is used to determine scale_x. c. A
// Y-sweep is done, and some virtual points are identified via triangulation. A statistical robust measure of the mean
// of these control points is used to determine scale_y.
static void handle_scale(int n, float *x, float *y, int gentrain_version, float *scale_x, float *scale_y) {
    if (n < ROBUST_THRESHOLD) {
        base_handle_scale(n, x, y, gentrain_version, scale_x, scale_y);
        return;
    }

    int i;
    int naa = 0;
    int nbb = 0;
    float xthrsh = 0.1f * percentile(n, x, 99);
    float ythrsh = 0.1f * percentile(n, y, 99);
    for (i = 0; i < n; i++) {
        if (x[i] > 5.0f * y[i] && x[i] > xthrsh) naa++;
        if (y[i] > 5.0f * x[i] && y[i] > ythrsh) nbb++;
    }
    int *iaa = (int *)malloc(naa * sizeof(int));
    int *ibb = (int *)malloc(nbb * sizeof(int));

    naa = 0;
    nbb = 0;
    for (i = 0; i < n; i++) {
        if (x[i] > 5.0f * y[i] && x[i] > xthrsh) iaa[naa++] = i;
        if (y[i] > 5.0f * x[i] && y[i] > ythrsh) ibb[nbb++] = i;
    }

    float *xaa = (float *)malloc(naa * sizeof(float));
    float *ybb = (float *)malloc(nbb * sizeof(float));
    for (i = 0; i < naa; i++) xaa[i] = x[iaa[i]];
    for (i = 0; i < nbb; i++) ybb[i] = y[ibb[i]];

    float xscale = matlab_trimmean(naa, xaa, 50);
    float yscale = matlab_trimmean(nbb, ybb, 50);

    for (i = 0; i < n; i++) {
        x[i] /= xscale;
        y[i] /= yscale;
    }

    *scale_x = (float)xscale;
    *scale_y = (float)yscale;

    free(iaa);
    free(ibb);
    free(xaa);
    free(ybb);
}

static void get_nn12_rr12(int n, const float *x, const float *y, float *nn12, float *rr12) {
    int i, j;
    float *xs = (float *)malloc(n * sizeof(float));
    float *ys = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) {
        xs[i] = x[i];
        ys[i] = y[i];
    }
    ks_introsort_float((size_t)n, xs);
    ks_introsort_float((size_t)n, ys);
    float xthrsh = 0.2f * percentile(n, xs, 95);
    float ythrsh = 0.2f * percentile(n, ys, 95);
    int count = 0;
    for (i = 0; i < n; i++)
        if (x[i] < xthrsh && y[i] < ythrsh) count++;
    *nn12 = (float)count / (float)n;
    if (count) {
        float *xy = (float *)malloc(n * sizeof(float));
        for (i = 0; i < n; i++) xy[i] = x[i] + y[i];
        float *xy12 = (float *)malloc(count * sizeof(float));
        for (i = 0, j = 0; i < n; i++)
            if (x[i] < xthrsh && y[i] < ythrsh) xy12[j++] = xy[i];
        float mean_xy12 = matlab_nanmean(count, xy12);
        float mean_xy = matlab_nanmean(n, xy);
        *rr12 = mean_xy12 / mean_xy;
        free(xy);
        free(xy12);
    } else {
        *rr12 = 1.0f;
    }
    free(xs);
    free(ys);
}

// a separate implementation from Illumina can be found in function NormalizeSingleBin from class
// Normalization111_Dragonfish
static void normalize_single_bin(int n, float *x, float *y, int gentrain_version, XForm *xform) {
    int naa, *iaa = NULL, nbb, *ibb = NULL;
    double (*madsigma)(int, const double *, int) = gentrain_version == 3 ? matlab_madsigma_new : matlab_madsigma_old;
    xform->version = 1;
    remove_outliers(&n, x, y);
    remove_offset(n, x, y, &naa, &iaa, &nbb, &ibb, madsigma, &xform->offset_x, &xform->offset_y);
    get_nn12_rr12(n, x, y, &xform->nn12, &xform->rr12);
    handle_rotation(n, x, y, &naa, &iaa, &nbb, &ibb, madsigma, &xform->theta);
    free(iaa);
    handle_shear(n, x, y, &nbb, &ibb, madsigma, &xform->shear);
    free(ibb);
    handle_scale(n, x, y, gentrain_version, &xform->scale_x, &xform->scale_y);

    xform->taa = (float)((double)(xform->theta * 180.0f) * M_1_PI);
    xform->tbb = (float)((atan((double)xform->shear) - (double)xform->theta) * 180.0 * M_1_PI);
}

// a separate implementation from Illumina can be found in function MirrorData from class NormalizationDragonfish
static void mirror_data(int n, float *x, float *y) {
    int i;
    for (i = 0; i < n; i++) {
        if (y[i] > x[i]) {
            float tmp = x[i];
            x[i] = y[i];
            y[i] = tmp;
        }
    }
}

// a separate implementation from Illumina can be found in function GetAA_Values from class NormalizationDragonfish
static int *get_aa_values(int n, const float *r, const float *t, int *naa) {
    int i, j;
    int *iaa = (int *)malloc(n * sizeof(int));
    for (i = 0, j = 0; i < n; i++)
        if (t[i] < 0.1f && !isnan(r[i]) && r[i] != FLT_MIN * FLT_EPSILON) iaa[j++] = i;
    *naa = j;
    return iaa;
}

// a separate implementation from Illumina can be found in functions RectToPolar from classes
// Normalization111_Dragonfish and Normalization111
static void rect_to_polar(int n, float *x, float *y) {
    int i;
    float *r = x;
    float *t = y;
    for (i = 0; i < n; i++) {
        if (x[i] == 0.0f && y[i] == 0.0f) {
            r[i] = NAN;
            t[i] = NAN;
            continue;
        }
        float tmp = x[i];
        r[i] = x[i] < 0.0f && y[i] < 0.0f ? FLT_MIN * FLT_EPSILON : fabsf(x[i]) + fabsf(y[i]);
        t[i] = (float)(atan2((double)y[i], (double)tmp) * M_2_PI);
    }
}

// a separate implementation from Illumina can be found in function NormalizeSingleBinSingleChannel from class
// NormalizationDragonfish Illumina, Inc. Improved Genotype Clustering with GenTrain 3.0. Pub. No. 370-2016-015-A (2016)
// In the sample intensity normalization process, specific
// groups of loci are normalized together in “normalization
// bins.” Due to differences in probe design, Infinium I loci
// (two probes per locus) and Infinium II loci (one probe per
// locus) are normalized in separate bins. If the number of loci
// in a normalization bin is small (< 192 loci), the normalization
// process can be negatively impacted. With the low bead
// pool complexity supported on the Infinium XT platform,
// the occurrence of small normalization bins may be more
// prevalent, especially with normalization bins consisting
// of Infinium I loci. With the GenTrain 2.0 algorithm, small
// normalization bin size negatively impacts the normalization
// of intensity data for the given locus (Figure 1A).
// The GenTrain 3.0 algorithm improves the normalization
// of small bins by taking advantage of the special nature
// of Infinium I assay data, where the signal intensity for
// both alleles originates in the same color channel. This
// affords the possibility to fit a normalization model with
// only two free parameters, instead of six. When applied to
// the same data mishandled by GenTrain 2.0, GenTrain 3.0
// improves the performance of the intensity normalization
// and generates tight clusters (Figure 1B). The GenTrain 3.0
// algorithm applies the improved normalization model for any
// normalization bin containing fewer than 192 Infinium I loci.
static void normalize_single_bin_single_channel(int n, float *x, float *y, XForm *xform) {
    int i, j, k;
    mirror_data(n, x, y);
    float *aux = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) aux[i] = y[i];
    ks_introsort_float((size_t)n, aux);
    float ythrsh = percentile(n, aux, 2);
    for (i = 0; i < n; i++) {
        x[i] -= ythrsh;
        y[i] -= ythrsh;
        aux[i] = y[i];
    }
    int naa;
    rect_to_polar(n, x, y);
    float *r = x;
    float *t = y;
    int *iaa = get_aa_values(n, r, t, &naa);
    for (i = 0; i < naa; i++) aux[i] = aux[iaa[i]];
    float ymean = matlab_trimmean(naa, aux, 20) - 50.0f;
    // here we replicate the bug in the normalization protocol
    for (i = 0, j = 0; i < n; i++)
        if (!isnan(r[i])) j++;
    for (i = 0, k = 0; i < j; i++)
        if (!isnan(r[i])) r[k++] = r[i];
    for (i = k; i < j; i++) r[i] = 0.0f;
    float rmean = matlab_trimmean(j, r, 20) - 2.0f * ymean;
    ythrsh += ymean;

    xform->version = 1;
    xform->offset_x = ythrsh;
    xform->offset_y = ythrsh;
    xform->theta = 0.0f;
    xform->shear = 0.0f;
    xform->scale_x = rmean;
    xform->scale_y = rmean;
    xform->rr12 = 1.0f;
    free(iaa);
    free(aux);
}

// a separate implementation from Illumina can be found in function Normalize from class NormalizationDragonfish
static XForm *normalize(int n, const uint16_t *xin, const uint16_t *yin, const uint8_t *norm_ids, int gentrain_version,
                        size_t *n_xforms) {
    int i, j, max_count = 0, counts[256];
    *n_xforms = 0;
    memset(counts, 0, 256 * sizeof(int));
    int *aux = (int32_t *)malloc(n * sizeof(int));

    // count size of sub-bead pool bins and sort coordinates by bin
    for (i = 0; i < n; i++) {
        aux[i] = (norm_ids[i] << 23) + i;
        counts[norm_ids[i]]++;
    }
    ks_introsort_int((size_t)n, aux);

    // compute number of sub-bead pool bins and size of the largest bin
    for (i = 0, j = 0; i < 256; i++) {
        if (counts[i]) {
            if (counts[i] > max_count) max_count = counts[i];
            counts[j++] = counts[i];
        }
    }
    *n_xforms = j;
    XForm *xform = (XForm *)calloc(*n_xforms, sizeof(XForm));
    float *x = (float *)malloc(max_count * sizeof(float));
    float *y = (float *)malloc(max_count * sizeof(float));

    // compute the normalization transform one sub-bead pool bin at a time
    int k = 0;
    for (i = 0; i < *n_xforms; i++) {
        if (counts[i] < 10) error("Error in normalization. Not enough good loci. Found %d\n", counts[i]);
        uint8_t norm_id = aux[k] >> 23;
        for (j = 0; j < counts[i]; j++) {
            int idx = aux[k] & 0x7FFFFF;
            x[j] = (float)xin[idx];
            y[j] = (float)yin[idx];
            k++;
        }
        if (gentrain_version == 3 && norm_id >= 100 && counts[i] < ROBUST_THRESHOLD) { // GenTrain 3.0 and Infinium I
            normalize_single_bin_single_channel(counts[i], x, y, &xform[i]);
        } else { // GenTrain 2.0 or Infinium II
            normalize_single_bin(counts[i], x, y, gentrain_version, &xform[i]);
        }
    }

    free(aux);
    free(x);
    free(y);
    return xform;
}

/****************************************
 * MATH ROUTINES                        *
 ****************************************/

// a separate implementation from Illumina of these functions in GenCall can be found in file Utils.cs

// http://www.mathworks.com/help/fuzzy/zmf.html
static float matlab_zmf(float x, float a, float b) {
    if (a >= b) error("Invalid arguments for zmf (a >= b)");
    if (x <= a) return 1;
    if (a < x && x <= (a + b) / 2.0f) return 1.0f - 2.0f * sqrf((x - a) / (b - a));
    if ((a + b) / 2.0f < x && x <= b) return 2.0f * sqrf((x - b) / (b - a));
    return 0;
}

// http://www.mathworks.com/help/fuzzy/smf.html
static float matlab_smf(float x, float a, float b) {
    if (a >= b) return x >= (a + b) / 2.0f ? 1.0f : 0.0f;
    if (x <= a) return 0;
    if (a < x && x <= (a + b) / 2.0f) return 2.0f * sqrf((x - a) / (b - a));
    if ((a + b) / 2.0f < x && x <= b) return 1.0f - 2.0f * sqrf((x - b) / (b - a));
    return 1;
}

// http://www.mathworks.com/help/stats/normpdf.html
static double matlab_normpdf_vleft(float x, float mu, float sigma) {
    if (sigma <= 0.0f) return NAN;
    if (x < mu) return 0.5 * M_2_SQRTPI * M_SQRT1_2 / (double)sigma;
    return exp(-0.5 * sqr((double)((x - mu) / sigma))) * 0.5 * M_2_SQRTPI * M_SQRT1_2 / (double)sigma;
}

// http://www.mathworks.com/help/stats/normpdf.html
static double matlab_normpdf_vright(float x, float mu, float sigma) {
    if (sigma <= 0.0f) return NAN;
    if (x > mu) return 0.5 * M_2_SQRTPI * M_SQRT1_2 / (double)sigma;
    return exp(-0.5 * sqr((double)((x - mu) / sigma))) * 0.5 * M_2_SQRTPI * M_SQRT1_2 / (double)sigma;
}

// http://www.mathworks.com/help/stats/normpdf.html
static double matlab_normpdf(float x, float mu, float sigma) {
    if (sigma <= 0.0f) return NAN;
    return exp(-0.5 * sqr((double)((x - mu) / sigma))) * 0.5 * M_2_SQRTPI * M_SQRT1_2 / (double)sigma;
}

/****************************************
 * GENOTYPE CALLING ROUTINES            *
 ****************************************/

// compute normalized intensities (http://github.com/Illumina/BeadArrayFiles/blob/develop/docs/GTC_File_Format_v5.pdf)
// a separate implementation from Illumina can be found in function Transform from class NormalizationTransform
static inline void raw_x_y2norm_x_y(uint16_t raw_x, uint16_t raw_y, float offset_x, float offset_y, float cos_theta,
                                    float sin_theta, float shear, float scale_x, float scale_y, float *norm_x,
                                    float *norm_y) {
    float temp_x = (float)raw_x - offset_x;
    float temp_y = (float)raw_y - offset_y;
    float temp_x2 = cos_theta * temp_x + sin_theta * temp_y;
    float temp_y2 = -sin_theta * temp_x + cos_theta * temp_y;
    float temp_x3 = temp_x2 - shear * temp_y2;
    *norm_x = temp_x3 < 0.0f ? 0.0f : temp_x3 / scale_x;
    *norm_y = temp_y2 < 0.0f ? 0.0f : temp_y2 / scale_y;
}

// compute Theta and R from raw intensities
static inline void norm_x_y2ilmn_theta_r(float norm_x, float norm_y, float *ilmn_theta, float *ilmn_r) {
    if (norm_x == 0.0f && norm_y == 0.0f) {
        *ilmn_theta = (float)NAN;
        *ilmn_r = (float)NAN;
        return;
    }
    *ilmn_theta = (float)(atan2((double)norm_y, (double)norm_x) * M_2_PI);
    if (norm_x < 0.0f && norm_y < 0.0f) {
        *ilmn_r = FLT_MIN * FLT_EPSILON;
    } else {
        *ilmn_r = fabsf(norm_x) + fabsf(norm_y);
    }
}

// http://stackoverflow.com/questions/23392321/most-efficient-way-to-find-median-of-three-integers
static inline float median3(float a, float b, float c) { return fmaxf(fminf(a, b), fminf(fmaxf(a, b), c)); }

// a separate implementation from Illumina can be found in function gen_std_flair from class GenTrain62 or in file
// GenTrain60.cs
static ClusterRecord *gen_std_flair(const ClusterRecord *cluster_record) {
    ClusterRecord *out_cluster = (ClusterRecord *)malloc(sizeof(ClusterRecord));
    memcpy(out_cluster, cluster_record, sizeof(ClusterRecord));

    int Mtight = 3;
    int Mloose = 3;

    float z1 = 0.5f * (cluster_record->ab_cluster_stats.theta_mean - cluster_record->aa_cluster_stats.theta_mean)
               / (cluster_record->ab_cluster_stats.theta_dev + cluster_record->aa_cluster_stats.theta_dev);
    float z2 = 0.5f * (cluster_record->bb_cluster_stats.theta_mean - cluster_record->ab_cluster_stats.theta_mean)
               / (cluster_record->bb_cluster_stats.theta_dev + cluster_record->ab_cluster_stats.theta_dev);
    float mz = fminf(z1, z2);
    float alpha = fmaxf((1.0f / (float)Mtight) * mz, 1.0f);
    float beta = fminf((1.0f / (float)Mloose) * mz, 1.0f);
    float eta = alpha * beta;

    out_cluster->aa_cluster_stats.theta_dev *= eta;
    out_cluster->ab_cluster_stats.theta_dev *= eta;
    out_cluster->bb_cluster_stats.theta_dev *= eta;

    float min_dispersion_t = 0.02f;
    if (out_cluster->aa_cluster_stats.theta_dev < min_dispersion_t)
        out_cluster->aa_cluster_stats.theta_dev = min_dispersion_t;
    if (out_cluster->ab_cluster_stats.theta_dev < min_dispersion_t)
        out_cluster->ab_cluster_stats.theta_dev = min_dispersion_t;
    if (out_cluster->bb_cluster_stats.theta_dev < min_dispersion_t)
        out_cluster->bb_cluster_stats.theta_dev = min_dispersion_t;

    float min_dispersion_r_cte = 0.2f;
    int M = 7;
    float min_dispersion_r = min_dispersion_r_cte;
    // compute median of the three values
    float med = median3(cluster_record->aa_cluster_stats.theta_dev, cluster_record->ab_cluster_stats.theta_dev,
                        cluster_record->bb_cluster_stats.theta_dev);
    if (min_dispersion_r < med) min_dispersion_r = med;

    float min_dispersion_r_aa = cluster_record->aa_cluster_stats.r_mean / (float)M;
    float min_dispersion_r_ab = cluster_record->ab_cluster_stats.r_mean / (float)M;
    float min_dispersion_r_bb = cluster_record->bb_cluster_stats.r_mean / (float)M;
    if (min_dispersion_r_aa < min_dispersion_r) min_dispersion_r_aa = min_dispersion_r;
    if (min_dispersion_r_ab < min_dispersion_r) min_dispersion_r_ab = min_dispersion_r;
    if (min_dispersion_r_bb < min_dispersion_r) min_dispersion_r_bb = min_dispersion_r;

    if (out_cluster->aa_cluster_stats.r_dev < min_dispersion_r_aa)
        out_cluster->aa_cluster_stats.r_dev = min_dispersion_r_aa;
    if (out_cluster->ab_cluster_stats.r_dev < min_dispersion_r_ab)
        out_cluster->ab_cluster_stats.r_dev = min_dispersion_r_ab;
    if (out_cluster->bb_cluster_stats.r_dev < min_dispersion_r_bb)
        out_cluster->bb_cluster_stats.r_dev = min_dispersion_r_bb;

    return out_cluster;
}

// a separate implementation from Illumina can be found in function modilik from class GenTrain62 or in file
// GenTrain60.cs this function computes the likelihood for each cluster
static void modilik(ClusterRecord *c, float t, float r, double *Laa, double *Lab, double *Lbb) {
    double alpha = 100.0; // what is the relevance of this?
    *Laa = 0.0;
    *Lab = 0.0;
    *Lbb = 0.0;

    // computes the Mahalanobis distance
    double daa =
        alpha * (double)fabsf(t - c->aa_cluster_stats.theta_mean) + (double)fabsf(r - c->aa_cluster_stats.r_mean);
    double dab =
        alpha * (double)fabsf(t - c->ab_cluster_stats.theta_mean) + (double)fabsf(r - c->ab_cluster_stats.r_mean);
    double dbb =
        alpha * (double)fabsf(t - c->bb_cluster_stats.theta_mean) + (double)fabsf(r - c->bb_cluster_stats.r_mean);

    int bCovered = 0;
    if (daa <= dbb && dab <= dbb && c->aa_cluster_stats.r_mean <= c->ab_cluster_stats.r_mean && !isnan(t)) {
        *Laa = matlab_normpdf_vleft(t, c->aa_cluster_stats.theta_mean, c->aa_cluster_stats.theta_dev)
               * matlab_normpdf_vleft(r, c->aa_cluster_stats.r_mean, c->aa_cluster_stats.r_dev);
        *Lab = matlab_normpdf_vright(t, c->ab_cluster_stats.theta_mean, c->ab_cluster_stats.theta_dev)
               * matlab_normpdf_vright(r, c->ab_cluster_stats.r_mean, c->ab_cluster_stats.r_dev);
        bCovered = 1;
    }
    if (daa <= dab && dbb <= dab && c->aa_cluster_stats.r_mean <= c->bb_cluster_stats.r_mean && !isnan(t)) {
        *Laa = matlab_normpdf_vleft(t, c->aa_cluster_stats.theta_mean, c->aa_cluster_stats.theta_dev)
               * matlab_normpdf_vleft(r, c->aa_cluster_stats.r_mean, c->aa_cluster_stats.r_dev);
        *Lbb = matlab_normpdf_vright(t, c->bb_cluster_stats.theta_mean, c->bb_cluster_stats.theta_dev)
               * matlab_normpdf_vright(r, c->bb_cluster_stats.r_mean, c->bb_cluster_stats.r_dev);
        bCovered = 1;
    }
    if (dab <= daa && dbb <= daa && c->ab_cluster_stats.r_mean <= c->bb_cluster_stats.r_mean && !isnan(t)) {
        *Lab = matlab_normpdf_vleft(t, c->ab_cluster_stats.theta_mean, c->ab_cluster_stats.theta_dev)
               * matlab_normpdf_vleft(r, c->ab_cluster_stats.r_mean, c->ab_cluster_stats.r_dev);
        *Lbb = matlab_normpdf_vright(t, c->bb_cluster_stats.theta_mean, c->bb_cluster_stats.theta_dev)
               * matlab_normpdf_vright(r, c->bb_cluster_stats.r_mean, c->bb_cluster_stats.r_dev);
        bCovered = 1;
    }
    if (daa <= dbb && dab <= dbb && c->aa_cluster_stats.r_mean > c->ab_cluster_stats.r_mean && !isnan(t)) {
        *Laa = matlab_normpdf_vleft(t, c->aa_cluster_stats.theta_mean, c->aa_cluster_stats.theta_dev)
               * matlab_normpdf_vright(r, c->aa_cluster_stats.r_mean, c->aa_cluster_stats.r_dev);
        *Lab = matlab_normpdf_vright(t, c->ab_cluster_stats.theta_mean, c->ab_cluster_stats.theta_dev)
               * matlab_normpdf_vleft(r, c->ab_cluster_stats.r_mean, c->ab_cluster_stats.r_dev);
        bCovered = 1;
    }
    if (dab <= daa && dbb <= daa && c->ab_cluster_stats.r_mean > c->bb_cluster_stats.r_mean && !isnan(t)) {
        *Lab = matlab_normpdf_vleft(t, c->ab_cluster_stats.theta_mean, c->ab_cluster_stats.theta_dev)
               * matlab_normpdf_vright(r, c->ab_cluster_stats.r_mean, c->ab_cluster_stats.r_dev);
        *Lbb = matlab_normpdf_vright(t, c->bb_cluster_stats.theta_mean, c->bb_cluster_stats.theta_dev)
               * matlab_normpdf_vleft(r, c->bb_cluster_stats.r_mean, c->bb_cluster_stats.r_dev);
        bCovered = 1;
    }
    if (daa <= dab && dbb <= dab && c->aa_cluster_stats.r_mean > c->bb_cluster_stats.r_mean && !isnan(t)) {
        *Laa = matlab_normpdf_vleft(t, c->aa_cluster_stats.theta_mean, c->aa_cluster_stats.theta_dev)
               * matlab_normpdf_vright(r, c->aa_cluster_stats.r_mean, c->aa_cluster_stats.r_dev);
        *Lbb = matlab_normpdf_vright(t, c->bb_cluster_stats.theta_mean, c->bb_cluster_stats.theta_dev)
               * matlab_normpdf_vleft(r, c->bb_cluster_stats.r_mean, c->bb_cluster_stats.r_dev);
        bCovered = 1;
    }

    if (!bCovered) {
        *Laa = matlab_normpdf(t, c->aa_cluster_stats.theta_mean, c->aa_cluster_stats.theta_dev)
               * matlab_normpdf(r, c->aa_cluster_stats.r_mean, c->aa_cluster_stats.r_dev);
        *Lab = matlab_normpdf(t, c->ab_cluster_stats.theta_mean, c->ab_cluster_stats.theta_dev)
               * matlab_normpdf(r, c->ab_cluster_stats.r_mean, c->ab_cluster_stats.r_dev);
        *Lbb = matlab_normpdf(t, c->bb_cluster_stats.theta_mean, c->bb_cluster_stats.theta_dev)
               * matlab_normpdf(r, c->bb_cluster_stats.r_mean, c->bb_cluster_stats.r_dev);
    }
}

// a separate implementation from Illumina can be found in function computeScoreCallPrelim from class GenTrain62 or in
// file GenTrain60.cs Illumina, Inc. Illumina GenCall Data Analysis Software. Pub. No. 370-2004-009 (2004) To call
// genotypes for an individual’s DNA, the calling algorithm takes the DNA’s intensity values and the information
// generated by the clustering algorithm; subsequently, it then identifies to which cluster the data for any specific
// locus (of the DNA of interest) corre- sponds. The DNA data is first normalized (using the same procedure as for the
// clustering algorithm). The calling operation (classification) is performed using a Bayesian model. The score for each
// call (GenCall Score) is the product of the GenTrain Score and a data-to-model fit score. After scoring all the loci
// in the DNA of interest, the application computes a composite score for that DNA (DNA Score). Subsequently, the
// GenCall score of each locus for this DNA is further penalized by the DNA Score. Shen,R. et al. High-throughput SNP
// genotyping on universal bead arrays. Mutat Res, 573, 70–82 (2005-06-03) A quality score, the GenCall score, is
// calculated for each SNP call, reflecting the degree of separation be- tween homozygote and heterozygote clusters for
// that SNP and the placement of the individual call within a cluster. To make a genotype call, the software looks at
// many factors but one of the first is the distribution of beads of the same type and in this way outliers are rejected
// to ensure genotyping accuracy. The GenCall score is composed of various sub-scores, of which the most important one
// is the clustering score. This score is a locus-specific score, and is computed by a fuzzy logic inference system. It
// varies from 0.0 to 1.0, and correlates with accuracy of the genotype call. GenCall scores have been shown to
// correlate with the accuracy of the genotyping call.
static float compute_score_call_prelim(float r, float t, const ClusterRecord *cluster_record, uint8_t *iAPmax) {
    if (r < cluster_record->intensity_threshold) {
        *iAPmax = (uint8_t)0;
        return (float)NAN;
    }
    double omega = 1.0;
    double Den = 1.0 + 3.0 * omega;

    ClusterRecord *c = gen_std_flair(cluster_record);

    // likelihoods
    double Laa;
    double Lab;
    double Lbb;
    if (isnan(t)) {
        Laa = Lab = Lbb = NAN;
    } else {
        modilik(c, t, r, &Laa, &Lab, &Lbb);
        if (isnan(Laa)) Laa = 0.0;
        if (isnan(Lab)) Lab = 0.0;
        if (isnan(Lbb)) Lbb = 0.0;
    }

    int N = c->aa_cluster_stats.N + c->ab_cluster_stats.N + c->bb_cluster_stats.N;

    // priors
    double Paa = ((double)((float)c->aa_cluster_stats.N / (float)N) + omega) / Den;
    double Pab = ((double)((float)c->ab_cluster_stats.N / (float)N) + omega) / Den;
    double Pbb = ((double)((float)c->bb_cluster_stats.N / (float)N) + omega) / Den;

    double Evidence = Laa * Paa + Lab * Pab + Lbb * Pbb;

    // posteriors
    double APaa = Laa * Paa / Evidence;
    double APab = Lab * Pab / Evidence;
    double APbb = Lbb * Pbb / Evidence;

    //
    if (APaa >= APab && APaa >= APbb)
        *iAPmax = (uint8_t)1; // AA
    else if (APab >= APaa && APab >= APbb)
        *iAPmax = (uint8_t)2; // AB
    else if (APbb >= APaa && APbb >= APab)
        *iAPmax = (uint8_t)3; // BB
    else
        *iAPmax = (uint8_t)0; // NC

    double mx = 0.0;
    double scndmx = 0.0;
    if (APaa > APab) {
        mx = APaa;
        scndmx = APab;
    } else {
        mx = APab;
        scndmx = APaa;
    }
    if (APbb > mx) {
        scndmx = mx;
        mx = APbb;
    } else if (APbb > scndmx) {
        scndmx = APbb;
    }

    double ap_ratio = mx / (DBL_MIN * DBL_EPSILON + scndmx);
    double ap_lod = log10(ap_ratio);

    double score_ap = matlab_smf((float)ap_lod, 0.0f, 2.0f);

    float score_r1 = matlab_smf(r, 0.0f, 0.1f);
    float score_r2 = 0.0f;
    float score_r3 = 0.0f;
    float score_r4 = 0.0f;

    int numClusters = 0;
    if (c->aa_cluster_stats.N > 0) numClusters++;
    if (c->ab_cluster_stats.N > 0) numClusters++;
    if (c->bb_cluster_stats.N > 0) numClusters++;
    float score_misclust;
    if (numClusters == 1 && c->ab_cluster_stats.N == 0)
        score_misclust = 0.7f;
    else if (numClusters != 3)
        score_misclust = 0.95f;
    else
        score_misclust = 1.0f;

    float score_t = 1.0f;

    float RdropBegin = 6.0f;
    float RdropEnd = 12.0f;
    switch (*iAPmax) {
    case 1: // AA
        score_t = matlab_zmf(t, c->aa_cluster_stats.theta_mean + 2.0f * c->aa_cluster_stats.theta_dev,
                             c->aa_cluster_stats.theta_mean + 6.0f * c->aa_cluster_stats.theta_dev);
        score_r2 = matlab_smf(r, 0.0f, c->aa_cluster_stats.r_mean / 10.0f);
        score_r3 = matlab_smf(r, c->aa_cluster_stats.r_mean - 6.0f * c->aa_cluster_stats.r_dev,
                              c->aa_cluster_stats.r_mean - 2.0f * c->aa_cluster_stats.r_dev);
        score_r4 = matlab_zmf(r, c->aa_cluster_stats.r_mean + RdropBegin * c->aa_cluster_stats.r_dev,
                              c->aa_cluster_stats.r_mean + RdropEnd * c->aa_cluster_stats.r_dev);
        break;
    case 2: // AB
        score_t = matlab_zmf(fabsf((t - c->ab_cluster_stats.theta_mean) / c->ab_cluster_stats.theta_dev), 2.0f, 6.0f);
        score_r2 = matlab_smf(r, 0.0f, c->ab_cluster_stats.r_mean / 10.0f);
        score_r3 = matlab_smf(r, c->ab_cluster_stats.r_mean - 6.0f * c->ab_cluster_stats.r_dev,
                              c->ab_cluster_stats.r_mean - 2.0f * c->ab_cluster_stats.r_dev);
        score_r4 = matlab_zmf(r, c->ab_cluster_stats.r_mean + RdropBegin * c->ab_cluster_stats.r_dev,
                              c->ab_cluster_stats.r_mean + RdropEnd * c->ab_cluster_stats.r_dev);
        break;
    case 3: // BB
        score_t = matlab_smf(t, c->bb_cluster_stats.theta_mean - 6.0f * c->bb_cluster_stats.theta_dev,
                             c->bb_cluster_stats.theta_mean - 2.0f * c->bb_cluster_stats.theta_dev);
        score_r2 = matlab_smf(r, 0.0f, c->bb_cluster_stats.r_mean / 10.0f);
        score_r3 = matlab_smf(r, c->bb_cluster_stats.r_mean - 6.0f * c->bb_cluster_stats.r_dev,
                              c->bb_cluster_stats.r_mean - 2.0f * c->bb_cluster_stats.r_dev);
        score_r4 = matlab_zmf(r, c->bb_cluster_stats.r_mean + RdropBegin * c->bb_cluster_stats.r_dev,
                              c->bb_cluster_stats.r_mean + RdropEnd * c->bb_cluster_stats.r_dev);
        break;
    }
    float score_r = score_r1 * score_r2 * score_r3 * score_r4;

    float score_n = 1.0f;

    if (isnan(cluster_record->cluster_score.total_score)) score_t = score_r = (float)NAN;
    if (isnan(t)) score_t = score_r = score_n = (float)NAN;

    float score_call_prelim =
        (float)score_ap * cluster_record->cluster_score.total_score * score_t * score_r * score_n * score_misclust;

    free(c);
    return score_call_prelim;
}

// http://www.mathworks.com/help/fuzzy/gbellmf.html
static double matlab_gbellmf(double x, double a, double b, double c) {
    double tmp = sqr((x - c) / a);
    if (tmp == 0.0 && b == 0.0) return 0.5;
    return 1.0 / (1.0 + pow(tmp, b));
}

// a separate implementation from Illumina can be found in function gencall_score_map from class GenTrain62 or in file
// GenTrain60.cs 0.35 = 0.5 * 0.7 0.504 = 0.8 × 0.7 × 0.9 1.71 = 0.9 * 1.9 1.08 = 0.9 * 1.2
static double gencall_score_map(double x) { return pow(x, 0.35) * matlab_gbellmf(x, 0.504, 1.71, 1.08); }

static inline char rev_allele(char allele) {
    static const char allele_complement[128] = {
        0, 0,   0, 0,   0,   0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0,   0, 0,   0,   0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 'T', 0, 'G', 'D', 0, 0, 'C', 0, 'I', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    if (allele > 95) return 0;
    return allele_complement[(int)allele];
}

// a separate implementation from Illumina can be found in function GetBaseCall from class AutoCallPollerThread
static void get_base_call(const char *snp, const char *ilmn_strand, uint8_t genotype, BaseCall *base_call) {
    char a = toupper(ilmn_strand[0]) == 'T' ? snp[1] : rev_allele(snp[1]);
    char b = toupper(ilmn_strand[0]) == 'T' ? snp[3] : rev_allele(snp[3]);
    switch (genotype) {
    case 1:
        (*base_call)[0] = a;
        (*base_call)[1] = a;
        return;
    case 2:
        (*base_call)[0] = a;
        (*base_call)[1] = b;
        return;
    case 3:
        (*base_call)[0] = b;
        (*base_call)[1] = b;
        return;
    }
    (*base_call)[0] = '-';
    (*base_call)[1] = '-';
}

// a separate implementation from Illumina can be found in function MakeCalls from class AutoCallPollerThread
static void make_calls(gtc_t *gtc, const bpm_t *bpm, const egt_t *egt, float gencall_cutoff) {
    int i, n = bpm->num_loci;
    gtc->sample_data.num_calls = 0;
    gtc->sample_data.num_intensity_only = 0;
    for (i = 0; i < n; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        uint16_t raw_x = gtc->raw_x[i];
        uint16_t raw_y = gtc->raw_y[i];
        float norm_x = -NAN;
        float norm_y = -NAN;
        float ilmn_theta = -NAN;
        float ilmn_r = -NAN;
        if (raw_x || raw_y) {
            int norm_id = bpm->norm_lookups[bpm->norm_ids[i]];
            XForm *xform = &gtc->normalization_transforms[norm_id];
            raw_x_y2norm_x_y(gtc->raw_x[i], gtc->raw_y[i], xform->offset_x, xform->offset_y, gtc->cos_theta[norm_id],
                             gtc->sin_theta[norm_id], xform->shear, xform->scale_x, xform->scale_y, &norm_x, &norm_y);
            norm_x_y2ilmn_theta_r(norm_x, norm_y, &ilmn_theta, &ilmn_r);

            int idx;
            int ret = khash_str2int_get(egt->names2index, locus_entry->name, &idx);
            if (ret < 0) error("Illumina probe %s not found in cluster file\n", locus_entry->name);
            ClusterRecord *cluster_record = &egt->cluster_records[idx];
            float min_dispersion_r = 0.1f;
            if (cluster_record->aa_cluster_stats.r_dev < min_dispersion_r)
                cluster_record->aa_cluster_stats.r_dev = min_dispersion_r;
            if (cluster_record->ab_cluster_stats.r_dev < min_dispersion_r)
                cluster_record->ab_cluster_stats.r_dev = min_dispersion_r;
            if (cluster_record->bb_cluster_stats.r_dev < min_dispersion_r)
                cluster_record->bb_cluster_stats.r_dev = min_dispersion_r;
            uint8_t genotype;
            float score_call_prelim = compute_score_call_prelim(ilmn_r, ilmn_theta, cluster_record, &genotype);
            float score_call = (float)gencall_score_map(score_call_prelim);
            gtc->genotype_scores[i] = isnan(score_call) ? 0.0f : score_call;
            gtc->genotypes[i] = genotype;
        } else {
            gtc->genotype_scores[i] = 0.0f;
            gtc->genotypes[i] = 0;
        }

        if (gtc->genotype_scores[i] < gencall_cutoff) gtc->genotypes[i] = 0;
        if (locus_entry->intensity_only) {
            gtc->genotypes[i] = 0;
            gtc->sample_data.num_intensity_only++;
        }
        if (gtc->genotypes[i]) gtc->sample_data.num_calls++;
        get_base_call(locus_entry->snp, locus_entry->ilmn_strand, gtc->genotypes[i], &gtc->base_calls[i]);
    }

    gtc->sample_data.num_no_calls = gtc->num_snps - gtc->sample_data.num_intensity_only - gtc->sample_data.num_calls;
    gtc->call_rate = (float)gtc->sample_data.num_calls
                     / ((float)gtc->num_snps - (float)gtc->sample_data.num_intensity_only + FLT_MIN * FLT_EPSILON);
}

typedef struct {
    int version;
    int min_loci;
    int max_loci;
    int min_x_loci;
    int min_y_loci;
    float call_rate_threshold;
    float y_threshold;
    float x_threshold;
    float x_het_rate_threshold;
} gender_t;

// a separate implementation from Illumina can be found in function EstimateGender from class AutoCallPollerThread
// TODO what happened here to gender->max_loci?
static void estimate_gender(gtc_t *gtc, const bpm_t *bpm, const egt_t *egt, const gender_t *gender) {
    int i, n = bpm->num_loci;
    int x_count = 0;
    int x_hets_count = 0;
    int x_non_missing_count = 0;
    int y_count = 0;
    int auto_count = 0;
    int auto_non_missing_count = 0;
    float *r_x = (float *)malloc(n * sizeof(float));
    float *r_y = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        const char *chrom =
            strncasecmp(locus_entry->chrom, "CHR", 3) == 0 ? locus_entry->chrom + 3 : locus_entry->chrom;
        int idx;
        int ret = khash_str2int_get(egt->names2index, locus_entry->name, &idx);
        if (ret < 0) error("Illumina probe %s not found in cluster file\n", locus_entry->name);
        ClusterRecord *cluster_record = &egt->cluster_records[idx];
        int norm_id = bpm->norm_lookups[bpm->norm_ids[i]];
        XForm *xform = &gtc->normalization_transforms[norm_id];
        float norm_x, norm_y, t;
        if (cluster_record->cluster_score.total_score != 0.0f) {
            if (strcmp(chrom, "X") == 0) {
                raw_x_y2norm_x_y(gtc->raw_x[i], gtc->raw_y[i], xform->offset_x, xform->offset_y,
                                 gtc->cos_theta[norm_id], gtc->sin_theta[norm_id], xform->shear, xform->scale_x,
                                 xform->scale_y, &norm_x, &norm_y);
                norm_x_y2ilmn_theta_r(norm_x, norm_y, &t, &r_x[x_count]);
                if (gtc->genotypes[i] == 2) x_hets_count++;
                if (gtc->genotypes[i]) x_non_missing_count++;
                x_count++;
            } else if (strcmp(chrom, "Y") == 0) {
                raw_x_y2norm_x_y(gtc->raw_x[i], gtc->raw_y[i], xform->offset_x, xform->offset_y,
                                 gtc->cos_theta[norm_id], gtc->sin_theta[norm_id], xform->shear, xform->scale_x,
                                 xform->scale_y, &norm_x, &norm_y);
                norm_x_y2ilmn_theta_r(norm_x, norm_y, &t, &r_y[y_count]);
                y_count++;
            } else if (strcmp(chrom, "XY") != 0 && strcmp(chrom, "MT") != 0) {
                auto_count++;
                if (gtc->genotypes[i]) auto_non_missing_count++;
            }
        }
    }

    gtc->gender = 'U';
    if (gender->version == 1 || y_count < gender->min_y_loci || auto_count < gender->min_loci) {
        if (x_non_missing_count > gender->min_x_loci) {
            gtc->gender = !((double)((float)x_hets_count / (float)x_non_missing_count) > gender->x_het_rate_threshold)
                              ? 'M'
                              : 'F';
        }
    } else if (auto_count > 0 && (double)auto_non_missing_count / (double)auto_count > gender->call_rate_threshold) {
        for (i = 0; i < y_count; i++)
            if (isnan(r_y[i]) || isinf(r_y[i])) r_y[i] = 0.0f;
        float y_med = matlab_median(y_count, r_y);
        if ((double)y_med > gender->y_threshold) {
            gtc->gender = 'M';
        } else if (x_count < gender->min_x_loci) {
            gtc->gender = 'F';
        } else {
            for (i = 0; i < x_count; i++)
                if (isnan(r_x[i]) || isinf(r_x[i])) r_x[i] = 0.0f;
            float x_med = matlab_median(x_count, r_x);
            gtc->gender = (double)x_med > gender->x_threshold ? 'F' : 'U';
        }
    }
    free(r_x);
    free(r_y);
}

// compute BAF and LRR from Theta and R as explained in Peiffer, D. A. et al. High-resolution genomic profiling of
// chromosomal aberrations using Infinium whole-genome genotyping. Genome Res. 16, 1136–1148 (2006)
// Peiffer, D. A. et al. High-resolution genomic profiling of chromosomal aberrations using Infinium whole-genome
// genotyping. Genome Res., 16, 1136–1148 (2006-08-09)
static inline void get_baf_lrr(float ilmn_theta, float ilmn_r, float aa_theta, float ab_theta, float bb_theta,
                               float aa_r, float ab_r, float bb_r, float r_mean, float *baf, float *lrr) {
    float r_ref;
    if (ilmn_theta == ab_theta) {
        r_ref = ab_r;
        *baf = 0.5f;
    } else if (ilmn_theta < ab_theta) {
        r_ref = aa_r + (ilmn_theta - aa_theta) * (aa_r - ab_r) / (aa_theta - ab_theta);
        *baf = (ilmn_theta - aa_theta) / (ab_theta - aa_theta) * 0.5f;
    } else if (ilmn_theta > ab_theta) {
        r_ref = ab_r + (ilmn_theta - ab_theta) * (bb_r - ab_r) / (bb_theta - ab_theta);
        *baf = 0.5f + (ilmn_theta - ab_theta) / (bb_theta - ab_theta) * 0.5f;
    } else {
        *lrr = -NAN;
        *baf = -NAN;
        return;
    }
    *lrr = ilmn_r != 0.0f ? (float)log2(ilmn_r / (isnan(r_mean) ? r_ref : r_mean)) : -FLT_MAX;
}

// a separate implementation from Illumina can be found in functions CalculateLogRDev and CalculateBAlleleFreq from
// class AutoCallPollerThread
static void calculate_baf_lrr(gtc_t *gtc, const bpm_t *bpm, const egt_t *egt) {
    int i, count = 0, n = bpm->num_loci;
    double sum = 0.0;
    double sum2 = 0.0;
    for (i = 0; i < n; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        int idx;
        int ret = khash_str2int_get(egt->names2index, locus_entry->name, &idx);
        if (ret < 0) error("Illumina probe %s not found in cluster file\n", locus_entry->name);
        ClusterRecord *c = &egt->cluster_records[idx];
        float baf = -NAN;
        float lrr = -NAN;
        if ((gtc->raw_x[i] || gtc->raw_y[i]) && c) {
            int norm_id = bpm->norm_lookups[bpm->norm_ids[i]];
            XForm *xform = &gtc->normalization_transforms[norm_id];
            float norm_x, norm_y, t, r;

            raw_x_y2norm_x_y(gtc->raw_x[i], gtc->raw_y[i], xform->offset_x, xform->offset_y, gtc->cos_theta[norm_id],
                             gtc->sin_theta[norm_id], xform->shear, xform->scale_x, xform->scale_y, &norm_x, &norm_y);
            norm_x_y2ilmn_theta_r(norm_x, norm_y, &t, &r);
            get_baf_lrr(t, r, c->aa_cluster_stats.theta_mean, c->ab_cluster_stats.theta_mean,
                        c->bb_cluster_stats.theta_mean, c->aa_cluster_stats.r_mean, c->ab_cluster_stats.r_mean,
                        c->bb_cluster_stats.r_mean, locus_entry->intensity_only ? c->r_mean : NAN, &baf, &lrr);
        }
        gtc->b_allele_freqs[i] = baf < 0.0 ? 0.0f : baf > 1.0 ? 1.0f : (float)baf;
        gtc->logr_ratios[i] = (float)lrr;

        char start_chrom = strncasecmp(locus_entry->chrom, "CHR", 3) == 0 ? toupper(locus_entry->chrom[3])
                                                                          : toupper(locus_entry->chrom[0]);
        if (!locus_entry->intensity_only && (start_chrom != 'X' && start_chrom != 'Y' && start_chrom != 'M')
            && !isinf(lrr) && !isnan(lrr)) {
            sum += (double)lrr;
            sum2 += sqr((double)lrr);
            count++;
        }
    }
    gtc->logr_dev = (float)sqrt(sum2 / (double)count - sqr(sum / (double)count));
}

// a separate implementation from Illumina can be found in function CalculateIntensityPercentiles from class
// AutoCallPollerThread
static void calculate_intensity_percentiles(gtc_t *gtc) {
    int i, n = gtc->num_snps;
    float *xs = (float *)malloc(n * sizeof(float));
    float *ys = (float *)malloc(n * sizeof(float));
    for (i = 0; i < n; i++) {
        xs[i] = (float)gtc->raw_x[i];
        ys[i] = (float)gtc->raw_y[i];
    }
    ks_introsort_float((size_t)n, xs);
    ks_introsort_float((size_t)n, ys);
    gtc->percentiles_x[0] = (uint16_t)percentile(n, xs, 5);
    gtc->percentiles_x[1] = (uint16_t)percentile(n, xs, 50);
    gtc->percentiles_x[2] = (uint16_t)percentile(n, xs, 95);
    gtc->percentiles_y[0] = (uint16_t)percentile(n, ys, 5);
    gtc->percentiles_y[1] = (uint16_t)percentile(n, ys, 50);
    gtc->percentiles_y[2] = (uint16_t)percentile(n, ys, 95);
    free(xs);
    free(ys);
}

// a separate implementation from Illumina can be found in function ComputeSampleStats from class AutoCallPollerThread
// Illumina, Inc. Illumina GenCall Data Analysis Software. Pub. No. 370-2004-009 (2004)
// GenCall Scores may be averaged among DNAs and
// among loci for purposes of evaluating the quality of the
// genotyping within a particular DNA or locus. For example,
// we often evaluate “GC10” and “GC50” scores that are calcu-
// lated by taking the 10th percentile and the 50th percentile
// (median) of the GenCall Scores for a certain locus, respec-
// tively. Using GC10 and GC50 Scores, a user may choose
// to fail particularly poor performing loci, for instance,
// by discarding loci with GC10 of 0.1 or lower. Also, a series
// of aggregate statistics (i.e., average) of the GC10 or GC50
// scores for each DNA can be used to identify low-quality
// DNAs (for instance, a user may discard DNA samples with
// average GC10 scores of 0.2 or lower). The GenCall Score
// can also be used in situations where users have a mini-
// mum required call rate. This rate translates to making
// calls on a certain percentile of the data. Users can sort
// all their genotypes based on the GenCall Score, and then
// choose the top (Nth) percentile of interest for their study.
static void compute_sample_stats(gtc_t *gtc, const bpm_t *bpm, float gencall_cutoff) {
    int i, j, n = gtc->num_snps;
    float *gs = (float *)malloc(n * sizeof(float));
    for (i = 0, j = 0; i < n; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        if (gtc->genotype_scores[i] > gencall_cutoff && !locus_entry->intensity_only)
            gs[j++] = (float)gtc->genotype_scores[i];
    }
    ks_introsort_float((size_t)j, gs);
    gtc->p10gc = percentile(j, gs, 10);
    gtc->sample_data.p50gc = percentile(j, gs, 50);
    free(gs);
}

/****************************************
 * CREATE NEW GTC STRUCTURE             *
 ****************************************/

// a separate implementation from Illumina can be found in class MD5ChecksumFile of the Array Analysis CLI
static char *basename(const char *fn, const char unsigned *md5_buf) {
    const char str[] = "(MD5Checksum=";
    char *ptr = strrchr(fn, '/');
    if (ptr)
        ptr++;
    else
        ptr = (char *)fn;
    char *ret;
    if (md5_buf) {
        int len = strlen(ptr);
        ret = (char *)malloc((len + 47) * sizeof(char));
        memcpy((void *)ret, (void *)ptr, (size_t)len);
        ptr = ret + len;
        memcpy((void *)ptr, &str, sizeof(str) - 1);
        ptr += sizeof(str) - 1;
        hts_md5_hex(ptr, md5_buf);
        ptr += 32;
        *(ptr++) = ')';
        *ptr = 0;
    } else {
        ret = strdup(ptr);
    }
    return ret;
}

// TODO this should be done once only for the BPM structure
static int32_t *get_control_addresses(const char *str, int *n_addresses) {
    int i, j;
    int moff = 0, *off = NULL;
    int moff2 = 0, *off2 = NULL;
    int32_t *addresses = NULL;
    int m_addresses = 0;
    *n_addresses = 0;

    char *s = strdup(str);
    int noff = ksplit_core(s, '\n', &moff, &off);
    for (i = 0; i < noff; i++) {
        char *ptr = strchr(&s[off[i]], ',');
        *ptr = '\0';
        int noff2 = ksplit_core(&s[off[i]], ':', &moff2, &off2);
        hts_expand(int32_t, *n_addresses + noff2, m_addresses, addresses);
        for (j = 0; j < noff2; j++) {
            char *endptr;
            addresses[*n_addresses + j] = (int32_t)strtol(&s[off[i] + off2[j]], &endptr, 10);
        }
        *n_addresses += noff2;
    }
    free(s);
    free(off);
    free(off2);
    return addresses;
}

static char *get_string_parameter(const char *str, const char *id) {
    const char *ptr = strstr(str, id);
    if (!ptr) return NULL;
    ptr += strlen(id);
    if (*ptr != '=') return NULL;
    ptr++;
    const char *ptr2 = strchr(ptr, '|');
    return strndup(ptr, ptr2 ? ptr2 - ptr : strlen(ptr));
}

static int32_t get_int32_parameter(const char *str, const char *id) {
    const char *ptr = strstr(str, id);
    if (!ptr) return 0;
    ptr += strlen(id);
    if (*ptr != '=') return 0;
    ptr++;
    char *endptr;
    return (int32_t)strtol(ptr, &endptr, 10);
}

// a separate implementation from Illumina can be found in function LoadSampleSection from class SampleData
// AutoConvert used the creation time of the IDAT file for the imaging date field of the GTC file:
// imagingDate = fileInfo.CreationTime.ToLongDateString() + " " + fileInfo.CreationTime.ToLongTimeString();
// this was later updated to instead use the imaging date field of the last Scan entry in the IDAT file
static void load_sample_section(gtc_t *gtc, const idat_t *idat, int imaging_date) {
    int i;
    RunInfo *run_info = NULL;
    for (i = 0; i < idat->m_run_infos; i++)
        if (strcmp(idat->run_infos[i].block_type, "Scan") == 0) run_info = &idat->run_infos[i];
    if (run_info) {
        gtc->imaging_date = imaging_date ? strdup(run_info->run_time) : NULL;
        gtc->scanner_data.scanner_name = get_string_parameter(run_info->block_pars, "sherlockID");
        gtc->scanner_data.pmt_green = get_int32_parameter(run_info->block_pars, "PMTGainCY3");
        gtc->scanner_data.pmt_red = get_int32_parameter(run_info->block_pars, "PMTGainCY5");
        gtc->scanner_data.scanner_version = strdup(run_info->code_version);
        gtc->scanner_data.imaging_user = get_string_parameter(run_info->block_pars, "Username");
    }
}

static int32_t get32_index(void *dict, int32_t key) {
    khash_t(32) *hash = (khash_t(32) *)dict;
    khiter_t k = kh_get(32, hash, key);
    if (k == kh_end(hash)) return -1;
    return kh_val(hash, k);
}

// a separate implementation from Illumina can be found in function fillArray from class SampleData
static void fill_array(const idat_t *grn_idat, const idat_t *red_idat, const bpm_t *bpm, gtc_t *gtc) {
    int i;
    int32_t idx1, idx2;
    for (i = 0; i < bpm->num_loci; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        if (locus_entry->assay_type == 0) { // 0 - Infinium II probes
            idx1 = get32_index(red_idat->ilmn_id2index, locus_entry->address_a);
            idx2 = get32_index(grn_idat->ilmn_id2index, locus_entry->address_a);
            if (idx1 == -1 || idx2 == -1) continue; // warning?
            if (red_idat->nbeads[idx1] >= 2 && grn_idat->nbeads[idx2] >= 2) {
                gtc->raw_x[i] = red_idat->trimmed_mean[idx1];
                gtc->raw_y[i] = grn_idat->trimmed_mean[idx2];
            }
        } else if (locus_entry->assay_type == 1) { // 1 - Infinium I (A/T) probes
            idx1 = get32_index(red_idat->ilmn_id2index, locus_entry->address_a);
            idx2 = get32_index(red_idat->ilmn_id2index, locus_entry->address_b);
            if (idx1 == -1 || idx2 == -1) continue; // warning?
            if (red_idat->nbeads[idx1] >= 2 && red_idat->nbeads[idx2] >= 2) {
                gtc->raw_x[i] = red_idat->trimmed_mean[idx1];
                gtc->raw_y[i] = red_idat->trimmed_mean[idx2];
            }
        } else if (locus_entry->assay_type == 2) { // 2 - Infinium I (G/C) probes
            idx1 = get32_index(grn_idat->ilmn_id2index, locus_entry->address_a);
            idx2 = get32_index(grn_idat->ilmn_id2index, locus_entry->address_b);
            if (idx1 == -1 || idx2 == -1) continue; // warning?
            if (grn_idat->nbeads[idx1] >= 2 && grn_idat->nbeads[idx2] >= 2) {
                gtc->raw_x[i] = grn_idat->trimmed_mean[idx1];
                gtc->raw_y[i] = grn_idat->trimmed_mean[idx2];
            }
        } else {
            error("Assay type %d for probe %s not valid\n", locus_entry->assay_type, locus_entry->ilmn_id);
        }
    }
}

// a separate implementation from Illumina can be found in function fillControlsArray from class SampleData
static void fill_controls_array(const idat_t *grn_idat, const idat_t *red_idat, const bpm_t *bpm, gtc_t *gtc) {
    int i, n_controls;
    int32_t idx1, idx2;
    int *control_addresses = get_control_addresses(bpm->control_config, &n_controls);
    gtc->m_controls_x = n_controls;
    gtc->m_controls_y = n_controls;
    gtc->controls_x = (uint16_t *)calloc(n_controls, sizeof(uint16_t));
    gtc->controls_y = (uint16_t *)calloc(n_controls, sizeof(uint16_t));
    for (i = 0; i < n_controls; i++) {
        idx1 = get32_index(red_idat->ilmn_id2index, control_addresses[i]);
        idx2 = get32_index(grn_idat->ilmn_id2index, control_addresses[i]);
        if (idx1 == -1 || idx2 == -1) continue; // warning?
        gtc->controls_x[i] = red_idat->trimmed_mean[idx1];
        gtc->controls_y[i] = grn_idat->trimmed_mean[idx2];
    }
    free(control_addresses);
}

// a separate implementation from Illumina can be found in function Process from class AutoCallPollerThread
static gtc_t *gtc_init(const idat_t *grn_idat, const idat_t *red_idat, const bpm_t *bpm, const egt_t *egt,
                       int gentrain_version, int gtc_file_version, float gencall_cutoff, int sample_name, int checksums,
                       int imaging_date, const char *autocall_date_format, const char *autocall_version,
                       const gender_t *gender) {
    if (!grn_idat || !red_idat || !bpm) return NULL;

    gtc_t *gtc = (gtc_t *)calloc(1, sizeof(gtc_t));
    gtc->version = gtc_file_version;
    gtc->ploidy = 2;
    gtc->ploidy_type = 1;
    gtc->sample_name = red_idat->sample_name ? strdup(red_idat->sample_name) : NULL;
    char *ptr, *ptr2;
    if (!gtc->sample_name && sample_name) {
        ptr = strrchr(grn_idat->fn, '/');
        if (ptr)
            ptr++;
        else
            ptr = grn_idat->fn;
        ptr2 = strstr(ptr, "_Grn.idat");
        gtc->sample_name = strndup(ptr, ptr2 - ptr);
    }
    gtc->sample_plate = red_idat->sample_plate ? strdup(red_idat->sample_plate) : NULL;
    gtc->sample_well = red_idat->sample_well ? strdup(red_idat->sample_well) : NULL;
    gtc->sentrix_id = red_idat->sentrix_barcode ? strdup(red_idat->sentrix_barcode) : NULL;
    if (egt) gtc->cluster_file = basename(egt->fn, checksums && egt ? egt->md5_buf : NULL);
    gtc->snp_manifest = basename(bpm->fn, checksums ? bpm->md5_buf : NULL);

    time_t timer;
    char buffer[26];
    struct tm *tm_info;
    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime(buffer, 26, autocall_date_format, tm_info);
    gtc->autocall_date = strdup(buffer);
    gtc->autocall_version = strdup(autocall_version);

    gtc->num_snps = bpm->num_loci;
    gtc->raw_x = (uint16_t *)calloc(gtc->num_snps, sizeof(uint16_t));
    gtc->raw_y = (uint16_t *)calloc(gtc->num_snps, sizeof(uint16_t));
    gtc->genotypes = (uint8_t *)calloc(gtc->num_snps, sizeof(uint8_t));
    gtc->base_calls = (BaseCall *)malloc(gtc->num_snps * sizeof(BaseCall));
    memset(gtc->base_calls, '-', gtc->num_snps * sizeof(BaseCall));
    gtc->genotype_scores = (float *)calloc(gtc->num_snps, sizeof(float));
    gtc->b_allele_freqs = (float *)calloc(gtc->num_snps, sizeof(float));
    gtc->logr_ratios = (float *)calloc(gtc->num_snps, sizeof(float));

    fill_array(grn_idat, red_idat, bpm, gtc);

    fprintf(stderr, "Normalizing...\n");
    gtc->normalization_transforms = normalize(gtc->num_snps, gtc->raw_x, gtc->raw_y, bpm->norm_ids, gentrain_version,
                                              &gtc->m_normalization_transforms);

    gtc->sin_theta = (float *)malloc(gtc->m_normalization_transforms * sizeof(float));
    gtc->cos_theta = (float *)malloc(gtc->m_normalization_transforms * sizeof(float));
    int i;
    for (i = 0; i < gtc->m_normalization_transforms; i++) {
        gtc->sin_theta[i] = (float)sin((double)gtc->normalization_transforms[i].theta);
        gtc->cos_theta[i] = (float)cos((double)gtc->normalization_transforms[i].theta);
    }

    if (egt) {
        fprintf(stderr, "Calling...\n");
        make_calls(gtc, bpm, egt, gencall_cutoff);
        fprintf(stderr, "Call rate: %.7f\n", gtc->call_rate);
        estimate_gender(gtc, bpm, egt, gender);
        fprintf(stderr, "Gender: %s\n", gtc->gender == 'M' ? "Male" : gtc->gender == 'F' ? "Female" : "Unknown");
        calculate_baf_lrr(gtc, bpm, egt);
        compute_sample_stats(gtc, bpm, gencall_cutoff);
    }

    calculate_intensity_percentiles(gtc);

    fill_controls_array(grn_idat, red_idat, bpm, gtc);

    load_sample_section(gtc, red_idat, imaging_date);

    return gtc;
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Convert Illumina IDAT files for Infinium arrays to GTC files.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: convert Illumina IDAT files for Infinium arrays to GTC files.\n"
           "(version " IDAT2GTC_VERSION
           " http://github.com/freeseek/idat2vcf)\n"
           "[ Kermani, B. G. Artificial intelligence and global normalization methods for\n"
           "  genotyping. U.S. Patents No. 7,035,740 (2005-09-29) and 7,467,117 (2006-10-05) ]\n"
           "[ Peiffer, D. A. et al. High-resolution genomic profiling of chromosomal aberrations\n"
           "  using Infinium whole-genome genotyping. Genome Res., 16, 1136–1148 (2006-08-09) ]\n"
           "[ Illumina, Inc. Illumina GenCall Data Analysis Software. Pub. No. 370-2004-009 (2004) ]\n"
           "[ Illumina, Inc. Illumina’s Genotyping Data Normalization Methods. Pub. No. 970-2006-010 (2006-09-26) ]\n"
           "[ Illumina, Inc. Improved Cluster Generation with Gentrain2. Pub. No. 037-2009-015 (2009-01-26)]\n"
           "[ Illumina, Inc. Improved Genotype Clustering with GenTrain 3.0. Pub. No. 370-2016-015-A (2016) ]\n"
           "Usage: bcftools +idat2gtc --bpm <file> [options]\n"
           "\n"
           "Plugin options:\n"
           "    -b, --bpm <file>                    BPM manifest file\n"
           "    -e, --egt <file>                    EGT cluster file\n"
           "    -i, --idats <dir>                   IDAT files from directory\n"
           "    -g, --grn-idats <file>              file with list of green IDATs\n"
           "    -r, --red-idats <file>              file with list of red IDATs\n"
           "    -o, --output <dir>                  write output to a directory\n"
           "    -v, --gentrain-version <int>        whether to use GenTrain 2.0 (2) or GenTrain 3.0 (3) for "
           "normalization [3]\n"
           "    -c, --gencall-cutoff <int>          cutoff score for GenCall algorithm [0.15]\n"
           "        --snp-map <file>                create SNP map file\n"
           "        --do-not-check-eof              do not check whether the BPM and EGT readers reach the end of the "
           "file\n"
           "        --preset <int>                  Illumina AutoCall software to emulate [4]\n"
           "                                        AutoConvert (1), AutoConvert 2.0 (2), IAAP CLI (3), Array Analysis "
           "CLI (4)\n"
           "GTC output files options:\n"
           "        --gtc-version <int>             whether use the old (3) or the new (5) GTC file format [5]\n"
           "        --no-sample-name                leave sample name empty if missing from IDAT files\n"
           "        --no-checksums                  do not include cluster and manifest files checksums\n"
           "        --no-imaging-date               do not include imaging date\n"
           "        --autocall-date <string>        AutoCall date format to use [" AUTOCALL_DATE_FORMAT_DFLT
           "]\n"
           "        --autocall-version <string>     AutoCall version label to use [" AUTOCALL_VERSION_DFLT
           "]\n"
           "\n"
           "Gender estimation options:\n"
           "        --gender-version <int>          whether to only use heterozygosity (1) or also intensities (2) "
           "[2]\n"
           "        --min-loci <int>                minimum number of autosomal loci for gender estimation [100]\n"
           "        --max-loci <int>                maximum number of autosomal loci for gender estimation [10000]\n"
           "        --min-x-loci <int>              minimum number of X loci for gender estimation [20]\n"
           "        --min-y-loci <int>              minimum number of Y loci for gender estimation [20]\n"
           "        --call-rate-threshold <float>   threshold for autosomal call rate for gender estimation [0.0]\n"
           "        --y-threshold <float>           threshold for Y intensity for gender estimation [0.3]\n"
           "        --x-threshold <float>           threshold for X intensity for gender estimation [0.9]\n"
           "        --x-het-rate-threshold <float>  threshold for X Het Rate for gender estimation [0.1]\n"
           "\n"
           "Examples:\n"
           "    bcftools +idat2gtc --bpm GSA-24v3-0_A1.bpm --egt GSA-24v3-0_A1_ClusterFile.egt \\\n"
           "      5434246082_R03C01_Grn.idat 5434246082_R03C01_Red.idat\n"
           "    bcftools +idat2gtc --bpm GSA-24v3-0_A1.bpm --egt GSA-24v3-0_A1_ClusterFile.egt --snp-map "
           "GSA-24v3-0_A1.bpm.csv\n"
           "    bcftools +idat2gtc --bpm GSA-24v3-0_A1.bpm --egt GSA-24v3-0_A1_ClusterFile.egt --gentrain-version 2 "
           "--gtc-version 3 \\\n"
           "      --no-sample-name --no-checksums --no-imaging-date --autocall-date \"\" --autocall-version 1.6.3.1 "
           "--gender-version 1\n"
           "    bcftools +idat2gtc --bpm GSA-24v3-0_A1.bpm --egt GSA-24v3-0_A1_ClusterFile.egt --no-sample-name \\\n"
           "      --no-checksums --autocall-date \"\" --autocall-version 2.0.1.179 --min-loci 10 --max-loci 100\n"
           "    bcftools +idat2gtc --bpm GSA-24v3-0_A1.bpm --egt GSA-24v3-0_A1_ClusterFile.egt --no-sample-name \\\n"
           "      --no-checksums --autocall-date \"\"\n"
           "\n";
}

static inline FILE *get_file_handle(const char *str) {
    if (!str) return NULL;
    FILE *ret;
    if (strcmp(str, "-") == 0) {
        ret = stdout;
    } else {
        ret = fopen(str, "w");
        if (!ret) error("Failed to open %s: %s\n", str, strerror(errno));
    }
    return ret;
}

// to recapitulate the .NET behavior of ToString() I need to add a small value
// http://stackoverflow.com/questions/2085449
// http://stackoverflow.com/questions/11085052
// http://stackoverflow.com/questions/14325214
static double round_adjust(double x) {
    double y = 5e-8, z;
    if (x > 1.0) {
        z = 1.0;
        while (x > z) {
            y *= 10.0;
            z *= 10.0;
        }
    } else {
        z = 0.1;
        while (x < z) {
            y *= 0.1;
            z *= 0.1;
        }
    }
    return x + y;
}

// this is the same file that can be generated by AutoConvert, AutoConvert 2.0 or by Picard
// BpmToNormalizationManifestCsv most likely this file was generated by AutoConvert to allow other software such as
// Illuminus, GenoSNP, Birdseed, optiCall, zCall, and iCall, to normalize intensities across sub-bead pools
// http://gatk.broadinstitute.org/hc/en-us/articles/360057440631-BpmToNormalizationManifestCsv-Picard
// a separate implementation from Illumina can be found in function CreateSNPMapFile from class AutoCallPollerThread
static void snp_map_write(const bpm_t *bpm, const egt_t *egt, const char *fn) {
    int i;
    FILE *out_txt = get_file_handle(fn);
    fprintf(out_txt, "Index,Name,Chromosome,Position,GenTrain Score,SNP,ILMN Strand,Customer Strand,NormID\n");
    for (i = 0; i < bpm->num_loci; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        const char *chrom =
            strncasecmp(locus_entry->chrom, "CHR", 3) == 0 ? locus_entry->chrom + 3 : locus_entry->chrom;
        double gentrain_score = NAN;
        if (egt) {
            int idx;
            int ret = khash_str2int_get(egt->names2index, locus_entry->name, &idx);
            if (ret < 0) error("Illumina probe %s not found in cluster file\n", locus_entry->name);
            ClusterRecord *cluster_record = &egt->cluster_records[idx];
            gentrain_score = round_adjust(cluster_record->cluster_score.total_score);
        }
        fprintf(out_txt, "%d,%s,%s,%s,%.4f,%s,%s,%s,%d\n", locus_entry->index, locus_entry->name, chrom,
                locus_entry->map_info, gentrain_score, locus_entry->snp, locus_entry->ilmn_strand,
                locus_entry->source_strand, locus_entry->norm_id);
    }
    if (out_txt != stdout && out_txt != stderr) fclose(out_txt);
}

void mkdir_p(const char *fmt, ...);

int run(int argc, char *argv[]) {
    const char *bpm_fname = NULL;
    const char *egt_fname = NULL;
    const char *snp_map_fname = NULL;
    const char *idat_pathname = NULL;
    const char *grn_idat_fname = NULL;
    const char *red_idat_fname = NULL;
    const char *output_pathname = ".";
    const char *autocall_date_format = AUTOCALL_DATE_FORMAT_DFLT;
    const char *autocall_version = AUTOCALL_VERSION_DFLT;
    char *tmp;
    int gentrain_version = 3;
    float gencall_cutoff = 0.15;
    int eof_check = 1;
    int preset = 4;
    int gtc_file_version = 5;
    int sample_name = 1;
    int checksums = 1;
    int imaging_date = 1;
    gender_t gender;
    gender.version = 2;                // 1 in AutoConvert
    gender.min_loci = 100;             // version 2
    gender.max_loci = 10000;           // version 2 for downsampling
    gender.min_x_loci = 20;            // shared between version 1 and 2
    gender.min_y_loci = 20;            // version 2
    gender.call_rate_threshold = 0.0;  // changed from 0.97
    gender.y_threshold = 0.3;          // version 2
    gender.x_threshold = 0.9;          // version 2
    gender.x_het_rate_threshold = 0.1; // version 1

    static struct option loptions[] = {{"bpm", required_argument, NULL, 'b'},
                                       {"egt", required_argument, NULL, 'e'},
                                       {"idats", required_argument, NULL, 'i'},
                                       {"grn-idats", required_argument, NULL, 'g'},
                                       {"red-idats", required_argument, NULL, 'r'},
                                       {"output", required_argument, NULL, 'o'},
                                       {"gentrain-version", required_argument, NULL, 'v'},
                                       {"gencall-cutoff", required_argument, NULL, 'c'},
                                       {"snp-map", required_argument, NULL, 1},
                                       {"do-not-check-eof", no_argument, NULL, 2},
                                       {"preset", required_argument, NULL, 3},
                                       {"gtc-version", required_argument, NULL, 4},
                                       {"no-sample-name", no_argument, NULL, 5},
                                       {"no-cheksums", no_argument, NULL, 6},
                                       {"no-imaging-date", no_argument, NULL, 7},
                                       {"autocall-date", required_argument, NULL, 8},
                                       {"autocall-version", required_argument, NULL, 9},
                                       {"gender-version", no_argument, NULL, 10},
                                       {"min-loci", no_argument, NULL, 11},
                                       {"max-loci", no_argument, NULL, 12},
                                       {"min-x-loci", no_argument, NULL, 13},
                                       {"min-y-loci", no_argument, NULL, 14},
                                       {"call-rate-threshold", no_argument, NULL, 15},
                                       {"y-threshold", no_argument, NULL, 16},
                                       {"x-threshold", no_argument, NULL, 17},
                                       {"x-het-rate-threshold", no_argument, NULL, 18},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?b:e:i:g:r:o:v:c:", loptions, NULL)) >= 0) {
        switch (c) {
        case 'b':
            bpm_fname = optarg;
            break;
        case 'e':
            egt_fname = optarg;
            break;
        case 'i':
            idat_pathname = optarg;
            break;
        case 'g':
            grn_idat_fname = optarg;
            break;
        case 'r':
            red_idat_fname = optarg;
            break;
        case 'o':
            output_pathname = optarg;
            break;
        case 'v':
            gentrain_version = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --gentrain-version %s\n", optarg);
            if (gentrain_version != 2 && gentrain_version != 3)
                error("The --gentrain-version option only allows values 2, and 3\n%s", usage_text());
            break;
        case 'c':
            gencall_cutoff = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --gencall-cutoff %s\n", optarg);
            break;
        case 1:
            snp_map_fname = optarg;
            break;
        case 2:
            eof_check = 0;
            break;
        case 3:
            preset = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --preset %s\n", optarg);
            if (preset < 1 || preset > 4)
                error("The --preset option only allows values 1, 2, 3, and 4\n%s", usage_text());
            switch (preset) {
            case 1:
                gentrain_version = 2;
                gender.version = 1;
                gtc_file_version = 3;
                sample_name = 0;
                checksums = 0;
                imaging_date = 0;
                autocall_version = "1.6.3.1";
                break;
            case 2:
                gentrain_version = 3;
                gender.version = 2;
                gender.call_rate_threshold = 0.97;
                gtc_file_version = 5;
                sample_name = 0;
                checksums = 0;
                imaging_date = 1;
                autocall_version = "2.0.1.179";
                break;
            case 3:
                gentrain_version = 3;
                gender.version = 2;
                // we did not reimplement the bug of estimating the autosomal call rate including loci with 0 cluster
                // scores as missing
                gender.call_rate_threshold = 0.97;
                gtc_file_version = 5;
                sample_name = 0;
                checksums = 0;
                imaging_date = 1;
                autocall_version = "3.0.0";
                break;
            case 4:
                gentrain_version = 3;
                gender.version = 2;
                gender.call_rate_threshold = 0.97;
                gtc_file_version = 5;
                sample_name = 1;
                checksums = 1;
                imaging_date = 1;
                autocall_version = "3.0.0";
                break;
            }
            break;
        case 4:
            gtc_file_version = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --gtc-version %s\n", optarg);
            if (gtc_file_version != 3 && gtc_file_version != 5)
                error("The --gtc-version option only allows values 3, and 5\n%s", usage_text());
            break;
        case 5:
            sample_name = 0;
            break;
        case 6:
            checksums = 0;
            break;
        case 7:
            imaging_date = 0;
            break;
        case 8:
            autocall_date_format = optarg;
            break;
        case 9:
            autocall_version = optarg;
            break;
        case 10:
            gender.version = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --gender-version %s\n", optarg);
            if (gender.version != 1 && gender.version != 2)
                error("The --gender-version option only allows values 1 and 2\n%s", usage_text());
            break;
        case 11:
            gender.min_loci = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --min-loci %s\n", optarg);
            break;
        case 12:
            gender.max_loci = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --max-loci %s\n", optarg);
            break;
        case 13:
            gender.min_x_loci = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --min-x-loci %s\n", optarg);
            break;
        case 14:
            gender.min_y_loci = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --min-y-loci %s\n", optarg);
            break;
        case 15:
            gender.call_rate_threshold = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --call-rate-threshold %s\n", optarg);
            break;
        case 16:
            gender.y_threshold = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --y-threshold %s\n", optarg);
            break;
        case 17:
            gender.x_threshold = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --x-threshold %s\n", optarg);
            break;
        case 18:
            gender.x_het_rate_threshold = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --x-het-rate-threshold %s\n", optarg);
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
        }
    }
    if (bpm_fname == NULL) error("The --bpm option is required\n%s", usage_text());
    if (idat_pathname != NULL && (grn_idat_fname != NULL || red_idat_fname != NULL))
        error("Cannot use option --idats with either option --grn-idats or --red-idats\n%s", usage_text());
    if (grn_idat_fname != NULL && red_idat_fname == NULL)
        error("Option --grn-idats requires option --red-idats\n%s", usage_text());
    if (grn_idat_fname == NULL && red_idat_fname != NULL)
        error("Option --red-idats requires option --grn-idats\n%s", usage_text());
    if (idat_pathname == NULL && grn_idat_fname == NULL && red_idat_fname == NULL) {
        if (snp_map_fname == NULL && argc - optind == 0) error("No IDAT files provided as input\n%s", usage_text());
        if (argc - optind % 2 == 1)
            error(
                "If options --idats/--grn-idats/--red-idats are not used, input an alternating list of green and red "
                "IDATs\n%s",
                usage_text());
    }

    fprintf(stderr, "idat2gtc " IDAT2GTC_VERSION " http://github.com/freeseek/gtc2vcf\n");
    fprintf(stderr, "Using normalization algorithm version %s\n", gentrain_version == 2 ? "1.1.2" : "1.2.0");

    if (strcmp(output_pathname, ".") != 0) mkdir_p("%s/", output_pathname);

    // read SNP manifest file
    fprintf(stderr, "Reading BPM file %s\n", bpm_fname);
    bpm_t *bpm = bpm_init(bpm_fname, eof_check, 0, checksums);

    // read cluster file
    egt_t *egt = NULL;
    if (egt_fname) {
        fprintf(stderr, "Reading EGT file %s\n", egt_fname);
        egt = egt_init(egt_fname, eof_check, checksums);
        if (!strcmp(egt->normalization_version, "1.2.0")) {
            if (gentrain_version != 3)
                fprintf(stderr, "Normalization algorithm version %s for cluster file %s corresponds to GenTrain 3.0\n",
                        egt->normalization_version, egt->fn);
        } else if (!strcmp(egt->normalization_version, "1.1.2")) {
            if (gentrain_version != 2)
                fprintf(stderr, "Normalization algorithm version %s for cluster file %s corresponds to GenTrain 2.0\n",
                        egt->normalization_version, egt->fn);
        } else if (!strcmp(egt->normalization_version, "1.1.0")) {
            if (gentrain_version != 1)
                fprintf(stderr, "Normalization algorithm version %s for cluster file %s corresponds to GenTrain 1.0\n",
                        egt->normalization_version, egt->fn);
        } else {
            fprintf(stderr, "Normalization algorithm version %s for cluster file %s is not recognized\n",
                    egt->normalization_version, egt->fn);
        }
    } else {
        fprintf(stderr, "No cluster file specified or forcing no cluster use\n");
        if (!gentrain_version) gentrain_version = 3;
    }

    // write SNP map file if requested
    if (snp_map_fname) snp_map_write(bpm, egt, snp_map_fname);

    // generate lists of green and red IDATs to process
    int i, n = 0;
    char **grn_idats = NULL;
    char **red_idats = NULL;
    if (idat_pathname != NULL) {
        // this code for now does not recursively looks for IDAT files
        DIR *d = opendir(idat_pathname);
        if (!d) error("Failed to open directory %s\n", idat_pathname);
        struct dirent *dir;
        int m_grn = 0;
        int m_red = 0;
        int p = strlen(idat_pathname);
        grn_idats = NULL;
        red_idats = NULL;
        while ((dir = readdir(d))) {
            char *ptr = strstr(dir->d_name, "_Grn.idat");
            if (!ptr) continue;
            hts_expand0(char *, n + 1, m_grn, grn_idats);
            hts_expand0(char *, n + 1, m_red, red_idats);
            int q = strlen(dir->d_name);
            grn_idats[n] = (char *)malloc((p + q + 2) * sizeof(char));
            memcpy(grn_idats[n], idat_pathname, p);
            grn_idats[n][p] = '/';
            memcpy(&grn_idats[n][p + 1], dir->d_name, q + 1);
            dir->d_name[q - 8] = 'R';
            dir->d_name[q - 7] = 'e';
            dir->d_name[q - 6] = 'd';
            red_idats[n] = (char *)malloc((p + q + 2) * sizeof(char));
            memcpy(red_idats[n], idat_pathname, p);
            red_idats[n][p] = '/';
            memcpy(&red_idats[n][p + 1], dir->d_name, q + 1);
            n++;
        }
        closedir(d);

    } else if (grn_idat_fname != NULL && red_idat_fname != NULL) {
        grn_idats = hts_readlines(grn_idat_fname, &n);
        int n_check;
        red_idats = hts_readlines(red_idat_fname, &n_check);
        if (n != n_check)
            error("File %s contains %d filenames while file %s contains %d filenames\n", grn_idat_fname, n,
                  red_idat_fname, n_check);
    } else if (argc > optind) {
        n = (argc - optind) / 2;
        grn_idats = (char **)malloc(n * sizeof(char *));
        red_idats = (char **)malloc(n * sizeof(char *));
        for (i = 0; i < n; i++) {
            grn_idats[i] = argv[optind++];
            red_idats[i] = argv[optind++];
        }
    }

    if (n > 0) {
        if (egt) {
            fprintf(stderr, "Using genotyping algorithm version %s\n", gentrain_version == 2 ? "6.3.0" : "7.0.0");
            fprintf(stderr, "Gender estimation parameters\n");
            fprintf(stderr, "\tVersion: %d\n", gender.version);
            fprintf(stderr, "\tMinX_Loci: %d\n", gender.min_x_loci);
            fprintf(stderr, "\tX_HetRateThreshold: %f\n", gender.x_het_rate_threshold);
            fprintf(stderr, "\tMinAutosomalLoci: %d\n", gender.min_loci);
            fprintf(stderr, "\tMaxAutosomalLoci: %d\n", gender.max_loci);
            fprintf(stderr, "\tMinY_Loci: %d\n", gender.min_y_loci);
            fprintf(stderr, "\tAutosomalCallRateThreshold: %f\n", gender.call_rate_threshold);
            fprintf(stderr, "\tX_IntensityThreshold: %f\n", gender.x_threshold);
            fprintf(stderr, "\tY_IntensityThreshold: %f\n", gender.y_threshold);
        }
        DIR *d = opendir(output_pathname);
        if (!d) error("Failed to open directory %s\n", output_pathname);
        kstring_t gtc_fname = {0, 0, NULL};
        for (i = 0; i < n; i++) {
            fprintf(stderr, "Reading GRN IDAT file %s\n", grn_idats[i]);
            idat_t *grn_idat = idat_init(grn_idats[i], 1);
            fprintf(stderr, "Reading RED IDAT file %s\n", red_idats[i]);
            idat_t *red_idat = idat_init(red_idats[i], 1);
            gtc_t *gtc =
                gtc_init(grn_idat, red_idat, bpm, egt, gentrain_version, gtc_file_version, gencall_cutoff, sample_name,
                         checksums, imaging_date, autocall_date_format, autocall_version, &gender);
            const char *ptr = strstr(grn_idats[i], "_Grn.idat");
            if (!ptr) ptr = strstr(grn_idats[i], ".idat");
            const char *ptr2 = strrchr(grn_idats[i], '/');
            if (ptr2)
                ptr2++;
            else
                ptr2 = grn_idats[i];
            ksprintf(&gtc_fname, "%s/%.*s.gtc", output_pathname, (int)(ptr ? ptr - ptr2 : strlen(ptr2)), ptr2);
            idat_destroy(grn_idat);
            idat_destroy(red_idat);
            fprintf(stderr, "Writing GTC file %s\n", gtc_fname.s);
            if (gtc_write(gtc, gtc_fname.s, gtc_file_version) < 0) error("Failed to write GTC file: %s\n", gtc_fname.s);
            gtc_destroy(gtc);
            gtc_fname.l = 0;
        }
        free(gtc_fname.s);
        closedir(d);
    }

    if (idat_pathname != NULL || grn_idat_fname != NULL || red_idat_fname != NULL) {
        for (i = 0; i < n; i++) {
            free(grn_idats[i]);
            free(red_idats[i]);
        }
    }
    free(grn_idats);
    free(red_idats);

    bpm_destroy(bpm);
    egt_destroy(egt);

    return 0;
}
