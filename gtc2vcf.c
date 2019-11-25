/* The MIT License

   Copyright (c) 2018-2019 Giulio Genovese

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

#include <getopt.h>
#include <errno.h>
#include <ctype.h>
#include <sys/resource.h>
#include <htslib/hfile.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "tsv2vcf.h"

#define GTC2VCF_VERSION "2019-11-24"

#define FT_GS (1<<4)

static inline char revnt(char nt)
{
    if ( nt=='A' ) return 'T';
    if ( nt=='C' ) return 'G';
    if ( nt=='G' ) return 'C';
    if ( nt=='T' ) return 'A';
    return -1;
}

// tests the end-of-file indicator for an hFILE
static inline int heof(hFILE *fp)
{
    if ( hgetc(fp) == EOF ) return 1;
    fp->begin--;
    return 0;
}

// read or skip a fixed number of bytes
static inline void read_bytes(hFILE *fp,
                              void *buffer,
                              size_t nbytes)
{
    if (buffer)
    {
        if ( hread(fp, buffer, nbytes) < nbytes )
        {
            error("Failed to read %ld bytes from stream\n", nbytes);
        }
    }
    else
    {
        int c = 0;
        for (int i=0; i<nbytes; i++) c = hgetc(fp);
        if ( c == EOF ) error("Failed to reposition stream forward %ld bytes\n", nbytes);
    }
}

// read or skip a fixed length array
static inline void read_array(hFILE *fp,
                              void **arr,
                              size_t *m_arr,
                              size_t nmemb,
                              size_t size,
                              size_t term)
{
    if (arr)
    {
        if (!m_arr)
        {
            *arr = malloc((nmemb + term) * size);
            if (!*arr) error("Failed to allocate memory for array\n");
        }
        else if (*m_arr < nmemb + term)
        {
            void *tmp = realloc(*arr, (nmemb + term) * size);
            if (!tmp) error("Failed to allocate memory for array\n");
            *arr = tmp;
            *m_arr = nmemb + term;
        }
        if ( hread(fp, *arr, nmemb * size) < nmemb * size )
        {
            error("Failed to read %ld bytes from stream\n", nmemb * size);
        }
    }
    else
    {
        int c = 0;
        for (int i=0; i<nmemb * size; i++) c = hgetc(fp);
        if ( c == EOF ) error("Failed to reposition stream forward %ld bytes\n", nmemb * size);
    }
}

// read or skip a length-prefixed array
static inline void read_pfx_array(hFILE *fp,
                                  void **arr,
                                  size_t *m_arr,
                                  size_t item_size)
{
    int32_t n;
    if ( hread(fp, (void *)&n, 4) < 4 )
    {
        error("Failed to read 4 bytes from stream\n");
    }
    read_array(fp, arr, m_arr, n, item_size, 0);
}

// read or skip a length-prefixed string
// https://en.wikipedia.org/wiki/LEB128#Decode_unsigned_integer
static inline void read_pfx_string(hFILE *fp,
                                   char **str,
                                   size_t *m_str)
{
    uint8_t byte;
    size_t n = 0, shift = 0;
    while (1)
    {
        if ( hread(fp, (void *)&byte, 1) < 1 )
        {
            error("Failed to read 1 byte from stream\n");
        }
        n |= (size_t)(byte & 0x7F) << shift;
        if ( !(byte & 0x80) ) break;
        shift += 7;
    }
    read_array(fp, (void **)str, m_str, n, 1, 1);
    if (str) (*str)[n] = '\0';
}

static inline int is_gzip(hFILE *fp)
{
    uint8_t buffer[2];
    if ( hpeek(fp, (void *)buffer, 2) < 2 )
    {
        error("Failed to read 2 bytes from stream\n");
    }
    return (buffer[0] == 0x1f && buffer[1] == 0x8b);
}

/****************************************
 * BUFFER ARRAY IMPLEMENTATION          *
 ****************************************/

typedef struct
{
    hFILE *fp;
    off_t offset;
    int32_t item_num;
    int32_t item_offset;
    size_t item_capacity;
    size_t item_size;
    char *buffer;
}
buffer_array_t;

static buffer_array_t *buffer_array_init(hFILE *fp,
                                         size_t capacity,
                                         size_t item_size)
{
    buffer_array_t *arr = (buffer_array_t *)malloc(1 * sizeof(buffer_array_t));
    arr->fp = fp;
    read_bytes(fp, (void *)&arr->item_num, sizeof(int32_t));
    arr->offset = htell(arr->fp);
    arr->item_offset = 0;
    arr->item_capacity = (capacity <= 0) ? 32768 : capacity;
    arr->item_size = item_size;
    arr->buffer = (char *)malloc(arr->item_capacity * item_size);
    read_bytes(fp, (void *)arr->buffer, (arr->item_num < arr->item_capacity ? arr->item_num : arr->item_capacity) * item_size);
    return arr;
}

static inline int get_element(buffer_array_t *arr,
                              void *dst,
                              size_t item_idx)
{
    if ( !arr || item_idx >= arr->item_num) return -1;
    else if (item_idx - arr->item_offset < arr->item_capacity)
    {
        memcpy(dst, (void *)(arr->buffer + (item_idx - arr->item_offset) * arr->item_size), arr->item_size);
        return 0;
    }
    arr->item_offset = item_idx;
    if ( hseek(arr->fp, arr->offset + item_idx * arr->item_size, SEEK_SET) < 0 )
    {
        error("Fail to seek to position %ld in file\n", arr->offset + item_idx * arr->item_size);
    }
    read_bytes(arr->fp, (void *)arr->buffer, ((arr->item_num - arr->item_offset) < arr->item_capacity ?
        (arr->item_num - arr->item_offset) : arr->item_capacity) * arr->item_size);
    memcpy(dst, (void *)arr->buffer, arr->item_size);
    return 0;
}

static void buffer_array_destroy(buffer_array_t *arr)
{
    if (!arr) return;
    free(arr->buffer);
    free(arr);
}

/****************************************
 * IDAT FILE IMPLEMENTATION             *
 ****************************************/

// see readIDAT_nonenc.R from https://github.com/HenrikBengtsson/illuminaio

#define NUM_SNPS_READ 1000
#define ILLUMINA_ID 102
#define SD 103
#define MEAN 104
#define NBEADS 107
#define MID_BLOCK 200
#define RED_GREEN 400
#define MOSTLY_NULL 401
#define SENTRIX_BARCODE 402
#define CHIP_TYPE 403
#define SENTRIX_POSITION 404
#define UNKNOWN_1 405
#define UNKNOWN_2 406
#define UNKNOWN_3 407
#define UNKNOWN_4 408
#define UNKNOWN_5 409
#define UNKNOWN_6 410
#define UNKNOWN_7 510
#define RUN_INFO 300

typedef struct
{
    char *run_time;
    char *block_type;
    char *block_pars;
    char *block_code;
    char *code_version;
}
RunInfo;

typedef struct
{
    char *fn;
    hFILE *fp;
    int64_t version;
    int32_t number_toc_entries;
    uint16_t *id;
    int64_t *toc;
    int32_t num_snps;
    buffer_array_t *ilmn_id;
    buffer_array_t *sd;
    buffer_array_t *mean;
    buffer_array_t *nbeads;
    buffer_array_t *mid_block;
    int32_t red_green;
    char *mostly_null;
    char *sentrix_barcode;
    char *chip_type;
    char *sentrix_position;
    char *unknown1;
    char *unknown2;
    char *unknown3;
    char *unknown4;
    char *unknown5;
    char *unknown6;
    char *unknown7;
    RunInfo *run_infos;
    int32_t m_run_infos;
}
idat_t;

static int idat_read(idat_t *idat, uint16_t id)
{
    int i;
    for (i=0; i<idat->number_toc_entries && id != idat->id[i]; i++);
    if (i == idat->number_toc_entries) return -1;
    if ( hseek(idat->fp, idat->toc[i], SEEK_SET) < 0 ) error("Fail to seek to position %ld in IDAT %s file\n", idat->toc[i], idat->fn);

    switch (id)
    {
        case NUM_SNPS_READ: read_bytes(idat->fp, (void *)&idat->num_snps, sizeof(int32_t)); break;
        case ILLUMINA_ID: idat->ilmn_id = buffer_array_init(idat->fp, 0, sizeof(int32_t)); break;
        case SD: idat->sd = buffer_array_init(idat->fp, 0, sizeof(uint16_t)); break;
        case MEAN: idat->mean = buffer_array_init(idat->fp, 0, sizeof(uint16_t)); break;
        case NBEADS: idat->nbeads = buffer_array_init(idat->fp, 0, sizeof(uint8_t)); break;
        case MID_BLOCK: idat->mid_block = buffer_array_init(idat->fp, 0, sizeof(uint8_t)); break;
        case RED_GREEN: read_bytes(idat->fp, (void *)&idat->red_green, sizeof(int32_t)); break;
        case MOSTLY_NULL: read_pfx_string(idat->fp, &idat->mostly_null, NULL); break;
        case SENTRIX_BARCODE: read_pfx_string(idat->fp, &idat->sentrix_barcode, NULL); break;
        case CHIP_TYPE: read_pfx_string(idat->fp, &idat->chip_type, NULL); break;
        case SENTRIX_POSITION: read_pfx_string(idat->fp, &idat->sentrix_position, NULL); break;
        case UNKNOWN_1: read_pfx_string(idat->fp, &idat->unknown1, NULL); break;
        case UNKNOWN_2: read_pfx_string(idat->fp, &idat->unknown2, NULL); break;
        case UNKNOWN_3: read_pfx_string(idat->fp, &idat->unknown3, NULL); break;
        case UNKNOWN_4: read_pfx_string(idat->fp, &idat->unknown4, NULL); break;
        case UNKNOWN_5: read_pfx_string(idat->fp, &idat->unknown5, NULL); break;
        case UNKNOWN_6: read_pfx_string(idat->fp, &idat->unknown6, NULL); break;
        case UNKNOWN_7:read_pfx_string(idat->fp, &idat->unknown7, NULL);  break;
        case RUN_INFO:
            read_bytes(idat->fp, (void *)&idat->m_run_infos, sizeof(int32_t));
            idat->run_infos = (RunInfo *)malloc(idat->m_run_infos * sizeof(RunInfo));
            for (int i=0; i<idat->m_run_infos; i++)
            {
                read_pfx_string(idat->fp, &idat->run_infos[i].run_time, NULL);
                read_pfx_string(idat->fp, &idat->run_infos[i].block_type, NULL);
                read_pfx_string(idat->fp, &idat->run_infos[i].block_pars, NULL);
                read_pfx_string(idat->fp, &idat->run_infos[i].block_code, NULL);
                read_pfx_string(idat->fp, &idat->run_infos[i].code_version, NULL);
            }
            break;
        default:
            error("IDAT file format does not support TOC entry %d\n", id);
            break;
    }
    return 0;
}

static idat_t *idat_init(const char *fn)
{
    idat_t *idat = (idat_t *)calloc(1, sizeof(idat_t));
    idat->fn = strdup(fn);
    idat->fp = hopen(idat->fn, "rb");
    if ( idat->fp == NULL ) error("Could not open %s: %s\n", idat->fn, strerror(errno));
    if ( is_gzip(idat->fp) ) error("File %s is gzip compressed and currently cannot be sought\n", idat->fn);

    uint8_t buffer[4];
    if ( hread(idat->fp, (void *)buffer, 4) < 4 ) error("Failed to read magic number from %s file\n", idat->fn);
    if ( memcmp(buffer, "IDAT", 4) != 0 ) error("IDAT file %s format identifier is bad\n", idat->fn);

    read_bytes(idat->fp, (void *)&idat->version, sizeof(int64_t));
    if ( idat->version<3 ) error("Cannot read IDAT file %s. Unsupported IDAT file format version: %ld\n", idat->fn, idat->version);

    read_bytes(idat->fp, (void *)&idat->number_toc_entries, sizeof(int32_t));
    idat->id = (uint16_t *)malloc(idat->number_toc_entries * sizeof(uint16_t));
    idat->toc = (int64_t *)malloc(idat->number_toc_entries * sizeof(int64_t));
    for (int i=0; i<idat->number_toc_entries; i++)
    {
        read_bytes(idat->fp, (void *)&idat->id[i], sizeof(uint16_t));
        read_bytes(idat->fp, (void *)&idat->toc[i], sizeof(int64_t));
    }

    for (int i=0; i<idat->number_toc_entries; i++) idat_read(idat, idat->id[i]);

    return idat;
}

static void idat_destroy(idat_t *idat)
{
    if (!idat) return;
    free(idat->fn);
    if ( hclose(idat->fp) < 0 ) error("Error closing IDAT file\n");
    free(idat->id);
    free(idat->toc);
    free(idat->mostly_null);
    free(idat->sentrix_barcode);
    free(idat->chip_type);
    free(idat->sentrix_position);
    free(idat->unknown1);
    free(idat->unknown2);
    free(idat->unknown3);
    free(idat->unknown4);
    free(idat->unknown5);
    free(idat->unknown6);
    free(idat->unknown7);
    for (int i=0; i<idat->m_run_infos; i++)
    {
        free(idat->run_infos[i].run_time);
        free(idat->run_infos[i].block_type);
        free(idat->run_infos[i].block_pars);
        free(idat->run_infos[i].block_code);
        free(idat->run_infos[i].code_version);
    }
    free(idat->run_infos);

    buffer_array_destroy(idat->ilmn_id);
    buffer_array_destroy(idat->sd);
    buffer_array_destroy(idat->mean);
    buffer_array_destroy(idat->nbeads);
    buffer_array_destroy(idat->mid_block);
    free(idat);
}

static void idat_summary(const idat_t *idat, FILE *stream)
{
    fprintf(stream, "IDAT file version = %ld\n", idat->version);
    fprintf(stream, "Number of TOC entries = %d\n", idat->number_toc_entries);
    fprintf(stream, "Number of SNPs = %d\n", idat->num_snps);
    fprintf(stream, "Number of mid blocks = %d\n", idat->mid_block->item_num);
    fprintf(stream, "Red Green = %d\n", idat->red_green);
    fprintf(stream, "Mostly Null = %s\n", idat->mostly_null);
    fprintf(stream, "Sentrix Barcode = %s\n", idat->sentrix_barcode);
    fprintf(stream, "Chip Type = %s\n", idat->chip_type);
    fprintf(stream, "Sentrix Position = %s\n", idat->sentrix_position);
}

static void idat_to_csv(const idat_t *idat, FILE *stream)
{

    fprintf(stream, "Illumina, Inc.\n");
    fprintf(stream, "[Heading]\n");
    fprintf(stream, "Descriptor File Name,%s\n", strrchr(idat->fn, '/') ? strrchr(idat->fn, '/') + 1 : idat->fn);
    fprintf(stream, "Sentrix Barcode,%s\n", idat->sentrix_barcode);
    fprintf(stream, "Chip Type,%s\n", idat->chip_type);
    fprintf(stream, "Sentrix Position,%s\n",idat->sentrix_position);
    fprintf(stream, "Probes Count ,%d\n", idat->num_snps);
    fprintf(stream, "Mid Blocks Count ,%d\n", idat->mid_block->item_num);
    fprintf(stream, "[Assay]\n");
    fprintf(stream, "IlmnID,Sd,Mean,Nbeads\n");
    for (int i=0; i<idat->num_snps; i++)
    {
        int32_t ilmn_id;
        get_element(idat->ilmn_id, (void *)&ilmn_id, i);
        int16_t sd;
        get_element(idat->sd, (void *)&sd, i);
        int16_t mean;
        get_element(idat->mean, (void *)&mean, i);
        int8_t nbeads;
        get_element(idat->nbeads, (void *)&nbeads, i);
        fprintf(stream, "%d,%d,%d,%d\n", ilmn_id, sd, mean, nbeads);
    }
    fprintf(stream, "[Mid Blocks]\n");
    for (int i=0; i<idat->mid_block->item_num; i++)
    {
        int8_t mid_block;
        get_element(idat->mid_block, (void *)&mid_block, i);
        fprintf(stream, "%d\n", mid_block);
    }
    fprintf(stream, "[Run Infos]\n");
    for (int i=0; i<idat->m_run_infos; i++)
    {
        fprintf(stream, "%s\t%s\t%s\t%s\t%s\n", idat->run_infos[i].run_time, idat->run_infos[i].block_type,
        idat->run_infos[i].block_pars, idat->run_infos[i].block_code, idat->run_infos[i].code_version);
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

/****************************************
 * BPM FILE IMPLEMENTATION              *
 ****************************************/

// see BeadPoolManifest.py from https://github.com/Illumina/BeadArrayFiles
// see also https://support.illumina.com/bulletins/2016/05/infinium-genotyping-manifest-column-headings.html

typedef struct
{
    int32_t version;
    uint8_t norm_id; // Normalization lookups from manifest. This indexes into list of normalization transforms read from GTC file
    char *ilmn_id; // IlmnID (probe identifier) of locus
    char *name; // Name (variant identifier) of locus
    int32_t index;
    char *ilmn_strand;
    char *snp; // SNP value for locus (e.g., [A/C])
    char *chrom; // Chromosome for the locus (e.g., XY)
    char *ploidy;
    char *species;
    char *map_info; // Mapping location of locus
    char *unknown_strand;
    int32_t address_a; // AddressA ID of locus
    char *allele_a_probe_seq; // CSV files
    int32_t address_b; // AddressB ID of locus (0 if none)
    char *allele_b_probe_seq; // CSV files (empty if none)
    char *genome_build;
    char *source;
    char *source_version;
    char *source_strand;
    char *source_seq; // CSV files
    char *top_genomic_seq; // CSV files
    int32_t beadset_id; // CSV files
    uint8_t exp_clusters;
    uint8_t intensity_only;
    uint8_t assay_type; // Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)
    float frac_a;
    float frac_c;
    float frac_g;
    float frac_t;
    char *ref_strand; // RefStrand annotation
}
LocusEntry;

static void locusentry_read(LocusEntry *locus_entry, hFILE *fp)
{
    locus_entry->norm_id = 0xFF;
    read_bytes(fp, (void *)&locus_entry->version, sizeof(int32_t));
    read_pfx_string(fp, &locus_entry->ilmn_id, NULL);
    read_pfx_string(fp, &locus_entry->name, NULL);
    read_pfx_string(fp, NULL, NULL);
    read_pfx_string(fp, NULL, NULL);
    read_pfx_string(fp, NULL, NULL);
    read_bytes(fp, (void *)&locus_entry->index, sizeof(int32_t));
    read_pfx_string(fp, NULL, NULL);
    read_pfx_string(fp, &locus_entry->ilmn_strand, NULL);
    read_pfx_string(fp, &locus_entry->snp, NULL);
    read_pfx_string(fp, &locus_entry->chrom, NULL);
    read_pfx_string(fp, &locus_entry->ploidy, NULL);
    read_pfx_string(fp, &locus_entry->species, NULL);
    read_pfx_string(fp, &locus_entry->map_info, NULL);
    read_pfx_string(fp, NULL, NULL);
    read_pfx_string(fp, &locus_entry->unknown_strand, NULL);
    read_bytes(fp, (void *)&locus_entry->address_a, sizeof(int32_t));
    read_bytes(fp, (void *)&locus_entry->address_b, sizeof(int32_t));
    read_pfx_string(fp, NULL, NULL);
    read_pfx_string(fp, NULL, NULL);
    read_pfx_string(fp, &locus_entry->genome_build, NULL);
    read_pfx_string(fp, &locus_entry->source, NULL);
    read_pfx_string(fp, &locus_entry->source_version, NULL);
    read_pfx_string(fp, &locus_entry->source_strand, NULL);
    read_pfx_string(fp, NULL, NULL);
    read_bytes(fp, NULL, 1);
    read_bytes(fp, (void *)&locus_entry->exp_clusters, sizeof(int8_t));
    read_bytes(fp, (void *)&locus_entry->intensity_only, sizeof(int8_t));
    read_bytes(fp, (void *)&locus_entry->assay_type, sizeof(uint8_t));

    // if ( strcmp(locus_entry->unknown_strand, locus_entry->source_strand) ) fprintf(stderr, "%s %s\n", locus_entry->unknown_strand, locus_entry->source_strand);

    if ( locus_entry->assay_type<0 || locus_entry->assay_type>2 ) error("Format error in reading assay type from locus entry\n");
    if ( locus_entry->address_b==0 && locus_entry->assay_type!=0 ) error("Manifest format error: Assay type is inconsistent with address B\n");
    if ( locus_entry->address_b!=0 && locus_entry->assay_type==0 ) error("Manifest format error: Assay type is inconsistent with address B\n");

    if (locus_entry->version>=7)
    {
        read_bytes(fp, &locus_entry->frac_a, sizeof(float));
        read_bytes(fp, &locus_entry->frac_c, sizeof(float));
        read_bytes(fp, &locus_entry->frac_t, sizeof(float));
        read_bytes(fp, &locus_entry->frac_g, sizeof(float));
    }
    if (locus_entry->version==8) read_pfx_string(fp, &locus_entry->ref_strand, NULL);
}

typedef struct
{
    char *fn;
    hFILE *fp;
    int32_t version;
    char *manifest_name; // Name of manifest
    char *control_config; // Control description from manifest
    int32_t num_loci; // Number of loci in manifest
    int32_t *indexes;
    char **names; // Names of loci from manifest
    uint8_t *norm_ids;
    LocusEntry *locus_entries;
    uint8_t *norm_lookups;
    char **header;
    size_t m_header;
}
bpm_t;

static uint8_t *bpm_norm_lookups(bpm_t *bpm)
{
    uint8_t sorted_norm_ids[256];
    for (int i=0; i<256; i++) sorted_norm_ids[i] = 0xFF;
    for (int i=0; i<bpm->num_loci; i++)
    {
        int norm_id = bpm->locus_entries[i].norm_id;
        sorted_norm_ids[norm_id] = norm_id;
    }
    int j=0;
    for (int i=0; i<256; i++) if ( sorted_norm_ids[i] != 0xFF ) sorted_norm_ids[j++] = sorted_norm_ids[i];
    uint8_t *norm_lookups = (uint8_t *)malloc(256 * sizeof(uint8_t *));
    memset((void *)norm_lookups, 0xFF, 256 * sizeof(uint8_t *));
    for (int i=0; i<j; i++) norm_lookups[sorted_norm_ids[i]] = i;
    return norm_lookups;
}

static bpm_t *bpm_init(const char *fn)
{
    bpm_t *bpm = (bpm_t *)calloc(1, sizeof(bpm_t));
    bpm->fn = strdup(fn);
    bpm->fp = hopen(bpm->fn, "rb");
    if ( bpm->fp == NULL ) error("Could not open %s: %s\n", bpm->fn, strerror(errno));
    if ( is_gzip(bpm->fp) ) error("File %s is gzip compressed and currently cannot be sought\n", bpm->fn);

    uint8_t buffer[4];
    if ( hread(bpm->fp, (void *)buffer, 4) < 4 ) error("Failed to read magic number from %s file\n", bpm->fn);
    if ( memcmp(buffer, "BPM", 3) != 0 ) error("BPM file %s format identifier is bad\n", bpm->fn);
    if ( buffer[3]!=1 ) error("BPM file %s version is unknown\n", bpm->fn);

    read_bytes(bpm->fp, (void *)&bpm->version, sizeof(int32_t));
    if ( bpm->version & 0x1000 ) bpm->version ^= 0x1000;
    if ( bpm->version > 5 || bpm->version < 3 ) error("BPM file %s version %d is unsupported\n", bpm->fn, bpm->version);
    read_pfx_string(bpm->fp, &bpm->manifest_name, NULL);

    if ( bpm->version > 1 ) read_pfx_string(bpm->fp, &bpm->control_config, NULL);

    read_bytes(bpm->fp, (void *)&bpm->num_loci, sizeof(int32_t));
    read_array(bpm->fp, (void **)&bpm->indexes, NULL, bpm->num_loci, sizeof(int32_t), 0);
    bpm->names = (char **)malloc(bpm->num_loci * sizeof(char *));
    for (int i=0; i<bpm->num_loci; i++) read_pfx_string(bpm->fp, &bpm->names[i], NULL);
    read_array(bpm->fp, (void **)&bpm->norm_ids, NULL, bpm->num_loci, sizeof(uint8_t), 0);

    bpm->locus_entries = (LocusEntry *)malloc(bpm->num_loci * sizeof(LocusEntry));
    LocusEntry locus_entry;
    for (int i=0; i<bpm->num_loci; i++)
    {
        memset(&locus_entry, 0, sizeof(LocusEntry));
        locusentry_read(&locus_entry, bpm->fp);
        int idx = locus_entry.index - 1;
        if ( idx < 0 || idx >= bpm->num_loci ) error("Locus entry index %d is out of boundaries\n", locus_entry.index);
        if ( bpm->norm_ids[idx] > 100 ) error("Manifest format error: read invalid normalization ID %d\n", bpm->norm_ids[idx]);
        locus_entry.norm_id = bpm->norm_ids[idx] + 100 * locus_entry.assay_type;
        memcpy(&bpm->locus_entries[ idx ], &locus_entry, sizeof(LocusEntry));
    }
    bpm->norm_lookups = bpm_norm_lookups(bpm);
    for (int i=0; i<bpm->num_loci; i++)
    {
        if ( i != bpm->locus_entries[i].index - 1 ) error("Manifest format error: read invalid number of assay entries\n");
    }

    read_bytes(bpm->fp, (void *)&bpm->m_header, sizeof(int32_t));
    bpm->header = (char **)malloc(bpm->m_header * sizeof(char *));
    for (int i=0; i<bpm->m_header; i++) read_pfx_string(bpm->fp, &bpm->header[i], NULL);

    if ( !heof(bpm->fp) ) error("BPM reader did not reach the end of file %s at position %ld\n", bpm->fn, htell(bpm->fp));

    return bpm;
}

static void bpm_destroy(bpm_t *bpm)
{
    if (!bpm) return;
    free(bpm->fn);
    if ( hclose(bpm->fp) < 0 ) error("Error closing BPM file\n");
    free(bpm->manifest_name);
    free(bpm->control_config);
    free(bpm->indexes);
    if ( bpm->names )
    {
        for (int i=0; i<bpm->num_loci; i++) free(bpm->names[i]);
        free(bpm->names);
    }
    free(bpm->norm_ids);
    for (int i=0; i<bpm->num_loci; i++)
    {
        free(bpm->locus_entries[i].ilmn_id);
        free(bpm->locus_entries[i].name);
        free(bpm->locus_entries[i].ilmn_strand);
        free(bpm->locus_entries[i].snp);
        free(bpm->locus_entries[i].chrom);
        free(bpm->locus_entries[i].ploidy);
        free(bpm->locus_entries[i].species);
        free(bpm->locus_entries[i].map_info);
        free(bpm->locus_entries[i].unknown_strand);
        free(bpm->locus_entries[i].allele_a_probe_seq);
        free(bpm->locus_entries[i].allele_b_probe_seq);
        free(bpm->locus_entries[i].genome_build);
        free(bpm->locus_entries[i].source);
        free(bpm->locus_entries[i].source_version);
        free(bpm->locus_entries[i].source_strand);
        free(bpm->locus_entries[i].source_seq);
        free(bpm->locus_entries[i].top_genomic_seq);
        free(bpm->locus_entries[i].ref_strand);
    }
    free(bpm->locus_entries);
    free(bpm->norm_lookups);
    for (int i=0; i<bpm->m_header; i++) free(bpm->header[i]);
    free(bpm->header);
    free(bpm);
}

static void bpm_summary(const bpm_t *bpm, FILE *stream)
{
    fprintf(stream, "BPM manifest file version = %d\n", bpm->version);
    fprintf(stream, "Name of manifest = %s\n", bpm->manifest_name);
    fprintf(stream, "Number of loci = %d\n", bpm->num_loci);
}

static void bpm_to_csv(const bpm_t *bpm, FILE *stream)
{
    for (int i=0; i<bpm->m_header; i++) fprintf(stream, "%s\n", bpm->header[i]);
    fprintf(stream, "Index,NormID,IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID,Exp_Clusters,Intensity_Only,Assay_Type,Frac A,Frac C,Frac G,Frac T,RefStrand\n");
    for (int i=0; i<bpm->num_loci; i++)
    {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        fprintf(stream, "%d,%d,%s,%s,%s,%s,%010d,%-s,%010d,%-s,%s,%s,%s,%s,%s,%s,%s,%s,%-s,%-s,%d,%d,%d,%d,%f,%f,%f,%f,%s\n",
            locus_entry->index,
            locus_entry->norm_id,
            locus_entry->ilmn_id,
            locus_entry->name,
            locus_entry->ilmn_strand,
            locus_entry->snp,
            locus_entry->address_a,
            locus_entry->allele_a_probe_seq,
            locus_entry->address_b,
            locus_entry->allele_b_probe_seq,
            locus_entry->genome_build,
            locus_entry->chrom,
            locus_entry->map_info,
            locus_entry->ploidy,
            locus_entry->species,
            locus_entry->source,
            locus_entry->source_version,
            locus_entry->source_strand,
            locus_entry->source_seq,
            locus_entry->top_genomic_seq,
            locus_entry->beadset_id,
            locus_entry->exp_clusters,
            locus_entry->intensity_only,
            locus_entry->assay_type,
            locus_entry->frac_a,
            locus_entry->frac_c,
            locus_entry->frac_g,
            locus_entry->frac_t,
            locus_entry->ref_strand);
    }
    fprintf(stream, "[Controls]\n");
    fprintf(stream, "%s\n", bpm->control_config);
    if (stream != stdout && stream != stderr) fclose(stream);
}

/****************************************
 * CSV FILE IMPLEMENTATION              *
 ****************************************/

int tsv_read_uint8(tsv_t *tsv,
                   bcf1_t *rec,
                   void *usr)
{
    uint8_t *uint8 = (uint8_t *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *uint8 = (uint8_t)strtol(tsv->ss, &endptr, 10);
    *tsv->se = tmp;
    return 0;
}

int tsv_read_int32(tsv_t *tsv,
                   bcf1_t *rec,
                   void *usr)
{
    int32_t *int32 = (int32_t *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *int32 = (int32_t)strtol(tsv->ss, &endptr, 10);
    *tsv->se = tmp;
    return 0;
}

int tsv_read_float(tsv_t *tsv,
                   bcf1_t *rec,
                   void *usr)
{
    float *single = (float *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *single = (float)strtof(tsv->ss, &endptr);
    *tsv->se = tmp;
    return 0;
}

int tsv_read_string(tsv_t *tsv,
                    bcf1_t *rec,
                    void *usr)
{
    char **str = (char **)usr;
    if (tsv->se == tsv->ss) *str = NULL;
    else
    {
        char tmp = *tsv->se;
        *tsv->se = 0;
        *str = strdup(tsv->ss);
        *tsv->se = tmp;
    }
    return 0;
}

// Petr Danecek's similar implementation in bcftools/tsv2vcf.c
int csv_parse(tsv_t *tsv,
              bcf1_t *rec,
              char *str)
{
    int status = 0;
    tsv->icol = 0;
    tsv->ss = tsv->se = str;
    while ( *tsv->ss && tsv->icol < tsv->ncols )
    {
        while ( *tsv->se && *tsv->se!=',' ) tsv->se++;
        if ( tsv->cols[tsv->icol].setter )
        {
            int ret = tsv->cols[tsv->icol].setter(tsv,rec,tsv->cols[tsv->icol].usr);
            if ( ret<0 ) return -1;
            status++;
        }
        tsv->se++;
        tsv->ss = tsv->se;
        tsv->icol++;
    }
    return status ? 0 : -1;
}

static uint8_t get_assay_type(const char *allele_a_probe_seq,
                              const char *allele_b_probe_seq,
                              const char *source_seq)
{
    if ( !allele_a_probe_seq || !source_seq ) return 0xFF;
    if ( !allele_b_probe_seq ) return 0;
    char *ptr1 = strchr(source_seq, '[');
    char *ptr2 = strchr(source_seq, ']');
    if ( !ptr1 || !ptr2 ) error("Source sequence is malformed: %s\n", source_seq);
    char trail1 = *(ptr1-1);
    char trail2 = *(ptr2+1);
    if ( (trail1 == 'a' || trail1 == 'A' || trail1 == 't' || trail1 == 'T') && (trail2 == 'a' || trail2 == 'A' || trail2 == 't' || trail2 == 'T') ) return 1;
    if ( (trail1 == 'c' || trail1 == 'C' || trail1 == 'g' || trail1 == 'G') && (trail2 == 'c' || trail2 == 'C' || trail2 == 'g' || trail2 == 'G') ) return 2;
    char last = allele_a_probe_seq[(strlen(allele_a_probe_seq)-1)];
    if ( (last == 'g' || last == 'G' || last == 't' || last == 'T') ) return 1;
    if ( (last == 'a' || last == 'A' || last == 'c' || last == 'C') ) return 2;
    error("Unable to retrieve assay type: %s %s\n", allele_a_probe_seq, source_seq);
}

static void locus_merge(LocusEntry *dest, LocusEntry *src)
{
    if ( src->version ) dest->version = src->version;
    if ( src->norm_id != 0xFF ) dest->norm_id = src->norm_id;
    if ( strcmp(dest->ilmn_id, src->ilmn_id) ) error("BPM and CSV manifests have conflicting IDs: %s and %s\n", dest->ilmn_id, src->ilmn_id);
    else { free(dest->ilmn_id); dest->ilmn_id = src->ilmn_id; }
    if ( src->name ) { free(dest->name); dest->name = src->name; }
    if ( src->index != 0 ) dest->index = src->index;
    if ( src->ilmn_strand ) { free(dest->ilmn_strand); dest->ilmn_strand = src->ilmn_strand; }
    if ( src->snp ) { free(dest->snp); dest->snp = src->snp; }
    if ( src->chrom ) { free(dest->chrom); dest->chrom = src->chrom; }
    if ( src->ploidy ) { free(dest->ploidy); dest->ploidy = src->ploidy; }
    if ( src->species ) { free(dest->species); dest->species = src->species; }
    if ( src->map_info ) { free(dest->map_info); dest->map_info = src->map_info; }
    if ( src->unknown_strand ) { free(dest->unknown_strand); dest->unknown_strand = src->unknown_strand; }
    if ( src->address_a != 0 ) dest->address_a = src->address_a;
    if ( src->allele_a_probe_seq ) { free(dest->allele_a_probe_seq); dest->allele_a_probe_seq = src->allele_a_probe_seq; }
    if ( src->address_b != 0 ) dest->address_b = src->address_b;
    if ( src->allele_b_probe_seq ) { free(dest->allele_b_probe_seq); dest->allele_b_probe_seq = src->allele_b_probe_seq; }
    if ( src->genome_build ) { free(dest->genome_build); dest->genome_build = src->genome_build; }
    if ( src->source ) { free(dest->source); dest->source = src->source; }
    if ( src->source_version ) { free(dest->source_version); dest->source_version = src->source_version; }
    if ( src->source_strand ) { free(dest->source_strand); dest->source_strand = src->source_strand; }
    if ( src->source_seq ) { free(dest->source_seq); dest->source_seq = src->source_seq; }
    if ( src->top_genomic_seq ) { free(dest->top_genomic_seq); dest->top_genomic_seq = src->top_genomic_seq; }
    if ( src->beadset_id ) dest->beadset_id = src->beadset_id;
    if ( src->exp_clusters ) dest->exp_clusters = src->exp_clusters;
    if ( src->intensity_only ) dest->intensity_only = src->intensity_only;
    if ( src->assay_type != 0xFF ) dest->assay_type = src->assay_type;
    if ( src->frac_a ) dest->frac_a = src->frac_a;
    if ( src->frac_c ) dest->frac_c = src->frac_c;
    if ( src->frac_g ) dest->frac_g = src->frac_g;
    if ( src->frac_t ) dest->frac_t = src->frac_t;
    if ( src->ref_strand ) { free(dest->ref_strand); dest->ref_strand = src->ref_strand; }
}

// this line will read a CSV file and if a BPM object is provided it will fill it rather than create a new one
static bpm_t *bpm_csv_init(const char *fn, bpm_t *bpm)
{
    int bpm_available = bpm != NULL;
    if ( !bpm_available ) bpm = (bpm_t *)calloc(1, sizeof(bpm_t));
    int bpm_prev_num_loci = bpm->num_loci;

    free(bpm->fn);
    bpm->fn = strdup(fn);
    if ( bpm->fp && hclose(bpm->fp) < 0 ) error("Error closing BPM file\n");
    bpm->fp = hopen(bpm->fn, "r");
    if ( bpm->fp == NULL ) error("Could not open %s: %s\n", bpm->fn, strerror(errno));
    if ( is_gzip(bpm->fp) ) error("File %s is gzip compressed and currently cannot be sought\n", bpm->fn);

    kstring_t str = {0, 0, NULL};
    if ( kgetline(&str, (kgets_func *)hgets, bpm->fp) < 0 ) error("Empty file: %s\n", bpm->fn);
    if ( strncmp(str.s, "Illumina, Inc.", 14) )
        error("Header of file %s is incorrect: %s\n", bpm->fn, str.s);
    char *tmp = NULL;
    size_t prev = 0;
    while ( strncmp(str.s + prev, "[Assay]", 7) )
    {
        if ( strncmp(str.s + prev, "Descriptor File Name,", 21) == 0 ) { free(bpm->manifest_name); bpm->manifest_name = strdup(str.s + prev + 21); }
        else if ( strncmp(str.s + prev, "Loci Count ,", 12) == 0 ) bpm->num_loci = (int)strtol(str.s + prev + 12, &tmp, 0);
        kputc('\n', &str);
        prev = str.l;
        if ( kgetline(&str, (kgets_func *)hgets, bpm->fp) < 0 ) error("Error reading from file: %s\n", bpm->fn);
    }
    if ( bpm->num_loci == 0 ) error("Could not understand number of loci from header of manifest file %s\n", bpm->fn);
    else if ( bpm_available && bpm_prev_num_loci != bpm->num_loci )
      error("BPM manifest file has %d loci while CSV manifest file %s has %d loci\n", bpm_prev_num_loci, bpm->fn, bpm->num_loci);
    int moff = 0, *off = NULL;
    for (int i=0; i<bpm->m_header; i++) free(bpm->header[i]);
    bpm->m_header = ksplit_core(str.s, '\n', &moff, &off);
    free(bpm->header);
    bpm->header = (char **)malloc(bpm->m_header * sizeof(char *));
    for (int i=0; i<bpm->m_header; i++) bpm->header[i] = strdup(&str.s[off[i]]);
    free(off);

    str.l = 0;
    if ( kgetline(&str, (kgets_func *)hgets, bpm->fp) < 0 ) error("Error reading from file: %s\n", bpm->fn);

    LocusEntry locus_entry;
    tsv_t *tsv = tsv_init(str.s);
    tsv_register(tsv, "Index", tsv_read_int32, &locus_entry.index);
    int norm_id = tsv_register(tsv, "NormID", tsv_read_uint8, &locus_entry.norm_id);
    tsv_register(tsv, "IlmnID", tsv_read_string, &locus_entry.ilmn_id);
    tsv_register(tsv, "Name", tsv_read_string, &locus_entry.name);
    tsv_register(tsv, "IlmnStrand", tsv_read_string, &locus_entry.ilmn_strand);
    tsv_register(tsv, "SNP", tsv_read_string, &locus_entry.snp);
    tsv_register(tsv, "AddressA_ID", tsv_read_int32, &locus_entry.address_a);
    tsv_register(tsv, "AlleleA_ProbeSeq", tsv_read_string, &locus_entry.allele_a_probe_seq);
    tsv_register(tsv, "AddressB_ID", tsv_read_int32, &locus_entry.address_b);
    tsv_register(tsv, "AlleleB_ProbeSeq", tsv_read_string, &locus_entry.allele_b_probe_seq);
    tsv_register(tsv, "GenomeBuild", tsv_read_string, &locus_entry.genome_build);
    tsv_register(tsv, "Chr", tsv_read_string, &locus_entry.chrom);
    tsv_register(tsv, "MapInfo", tsv_read_string, &locus_entry.map_info);
    tsv_register(tsv, "Ploidy", tsv_read_string, &locus_entry.ploidy);
    tsv_register(tsv, "Species", tsv_read_string, &locus_entry.species);
    tsv_register(tsv, "Source", tsv_read_string, &locus_entry.source);
    tsv_register(tsv, "SourceVersion", tsv_read_string, &locus_entry.source_version);
    tsv_register(tsv, "SourceStrand", tsv_read_string, &locus_entry.source_strand);
    tsv_register(tsv, "SourceSeq", tsv_read_string, &locus_entry.source_seq);
    tsv_register(tsv, "TopGenomicSeq", tsv_read_string, &locus_entry.top_genomic_seq);
    tsv_register(tsv, "BeadSetID", tsv_read_int32, &locus_entry.beadset_id);
    tsv_register(tsv, "Exp_Clusters", tsv_read_uint8, &locus_entry.exp_clusters);
    tsv_register(tsv, "Intensity_Only", tsv_read_uint8, &locus_entry.intensity_only);
    tsv_register(tsv, "Frac A", tsv_read_float, &locus_entry.frac_a);
    tsv_register(tsv, "Frac C", tsv_read_float, &locus_entry.frac_c);
    tsv_register(tsv, "Frac G", tsv_read_float, &locus_entry.frac_g);
    tsv_register(tsv, "Frac T", tsv_read_float, &locus_entry.frac_t);
    tsv_register(tsv, "RefStrand", tsv_read_string, &locus_entry.ref_strand);

    if ( !bpm_available ) bpm->locus_entries = (LocusEntry *)malloc(bpm->num_loci * sizeof(LocusEntry));
    for (int i=0; i<bpm->num_loci; i++)
    {
        memset(&locus_entry, 0, sizeof(LocusEntry));
        locus_entry.norm_id = 0xFF;
        str.l = 0;
        if ( kgetline(&str, (kgets_func *)hgets, bpm->fp) < 0 ) error("Error reading from file: %s\n", bpm->fn);
        if ( csv_parse(tsv, NULL, str.s) < 0 ) error("Could not parse the manifest file: %s\n", str.s);
        locus_entry.assay_type = get_assay_type(locus_entry.allele_a_probe_seq, locus_entry.allele_b_probe_seq, locus_entry.source_seq);
        int idx = locus_entry.index ? locus_entry.index - 1 : i;
        if ( idx < 0 || idx >= bpm->num_loci ) error("Locus entry index %d is out of boundaries\n", idx);
        if ( !bpm_available ) memcpy(&bpm->locus_entries[ idx ], &locus_entry, sizeof(LocusEntry));
        else locus_merge(&bpm->locus_entries[ idx ], &locus_entry);
    }
    tsv_destroy(tsv);
    free(str.s);

    if ( norm_id == 0 )
    {
        free(bpm->norm_lookups);
        bpm->norm_lookups = bpm_norm_lookups(bpm);
    }

    return bpm;
}

/****************************************
 * EGT FILE IMPLEMENTATION              *
 ****************************************/

// see ClusterFile.py from https://github.com/Illumina/BeadArrayFiles

typedef struct
{
    int32_t N; // Number of samples assigned to cluster during training
    float r_dev; // R (intensity) std deviation value
    float r_mean; // R (intensity) mean value
    float theta_dev; // Theta std devation value
    float theta_mean; // Theta mean value
}
ClusterStats;

typedef struct
{
    float cluster_separation; // A score measure the separation between genotype clusters
    float total_score; // The GenTrain score
    float original_score; // The original score before editing this cluster
    uint8_t edited; // Whether this cluster has been manually manipulated
}
ClusterScore;

typedef struct
{
    ClusterStats aa_cluster_stats; // Describes AA genotype cluster
    ClusterStats ab_cluster_stats; // Describes AB genotype cluster
    ClusterStats bb_cluster_stats; // Describes BB genotype cluster
    float intensity_threshold; // Intensity threshold for no-call
    ClusterScore cluster_score; // Various scores for cluster
    int32_t address; // Bead type identifier for probe A
}
ClusterRecord;

typedef struct
{
    char *fn;
    hFILE *fp;
    int32_t version;
    char *gencall_version; // The GenCall version
    char *cluster_version; // The clustering algorithm version
    char *call_version; // The genotyping algorithm version
    char *normalization_version; // The normalization algorithm version
    char *date_created; // The date the cluster file was created (e.g., 3/9/2017 2:18:30 PM)
    uint8_t is_wgt;
    int32_t data_block_version;
    char *manifest_name; // The manifest name used to build this cluster file
    int32_t num_records;
    ClusterRecord *cluster_records;
    char **loci_names;
}
egt_t;

static void clusterscore_read(ClusterScore *clusterscore, hFILE *fp)
{
    read_bytes(fp, (void *)&clusterscore->cluster_separation, sizeof(float));
    read_bytes(fp, (void *)&clusterscore->total_score, sizeof(float));
    read_bytes(fp, (void *)&clusterscore->original_score, sizeof(float));
    read_bytes(fp, (void *)&clusterscore->edited, sizeof(uint8_t));
}

static void clusterrecord_read(ClusterRecord *clusterrecord, hFILE *fp)
{
    read_bytes(fp, (void *)&clusterrecord->aa_cluster_stats.N, sizeof(int32_t));
    read_bytes(fp, (void *)&clusterrecord->ab_cluster_stats.N, sizeof(int32_t));
    read_bytes(fp, (void *)&clusterrecord->bb_cluster_stats.N, sizeof(int32_t));
    read_bytes(fp, (void *)&clusterrecord->aa_cluster_stats.r_dev, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->ab_cluster_stats.r_dev, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->bb_cluster_stats.r_dev, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->aa_cluster_stats.r_mean, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->ab_cluster_stats.r_mean, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->bb_cluster_stats.r_mean, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->aa_cluster_stats.theta_dev, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->ab_cluster_stats.theta_dev, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->bb_cluster_stats.theta_dev, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->aa_cluster_stats.theta_mean, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->ab_cluster_stats.theta_mean, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->bb_cluster_stats.theta_mean, sizeof(float));
    read_bytes(fp, (void *)&clusterrecord->intensity_threshold, sizeof(float));
    read_bytes(fp, NULL, 14 * sizeof(float));
}

static egt_t *egt_init(const char *fn)
{
    egt_t *egt = (egt_t *)calloc(1, sizeof(egt_t));
    egt->fn = strdup(fn);
    egt->fp = hopen(egt->fn, "rb");
    if ( egt->fp == NULL ) error("Could not open %s: %s\n", egt->fn, strerror(errno));
    if ( is_gzip(egt->fp) ) error("File %s is gzip compressed and currently cannot be sought\n", egt->fn);

    read_bytes(egt->fp, (void *)&egt->version, sizeof(int32_t));
    if ( egt->version != 3 ) error("EGT cluster file version %d not supported\n", egt->version);

    read_pfx_string(egt->fp, &egt->gencall_version, NULL);
    read_pfx_string(egt->fp, &egt->cluster_version, NULL);
    read_pfx_string(egt->fp, &egt->call_version, NULL);
    read_pfx_string(egt->fp, &egt->normalization_version, NULL);
    read_pfx_string(egt->fp, &egt->date_created, NULL);

    read_bytes(egt->fp, (void *)&egt->is_wgt, sizeof(uint8_t));
    if ( egt->is_wgt != 1 ) error("Only WGT cluster file version supported\n");

    read_pfx_string(egt->fp, &egt->manifest_name, NULL);

    read_bytes(egt->fp, (void *)&egt->data_block_version, sizeof(int32_t));
    if ( egt->data_block_version < 7 || egt->data_block_version > 9 ) error("Data block version %d in cluster file not supported\n", egt->data_block_version);
    read_pfx_string(egt->fp, NULL, NULL); // opa

    read_bytes(egt->fp, (void *)&egt->num_records, sizeof(int32_t));
    egt->cluster_records = (ClusterRecord *)malloc(egt->num_records * sizeof(ClusterRecord));
    for (int i=0; i<egt->num_records; i++) clusterrecord_read(&egt->cluster_records[i], egt->fp);
    for (int i=0; i<egt->num_records; i++) clusterscore_read(&egt->cluster_records[i].cluster_score, egt->fp);

    // toss useless strings such as aa_ab_bb/aa_ab/aa_bb/ab_bb
    for (int i=0; i<egt->num_records; i++) read_pfx_string(egt->fp, NULL, NULL);

    egt->loci_names = (char **)malloc(egt->num_records * sizeof(char *));
    for (int i=0; i<egt->num_records; i++)
    {
        read_pfx_string(egt->fp, &egt->loci_names[i], NULL);
    }
    for (int i=0; i<egt->num_records; i++) read_bytes(egt->fp, (void *)&egt->cluster_records[i].address, sizeof(int32_t));

    int32_t aa_n, ab_n, bb_n;
    for (int i=0; i<egt->num_records; i++)
    {
        read_bytes(egt->fp, (void *)&aa_n, sizeof(int32_t));
        read_bytes(egt->fp, (void *)&ab_n, sizeof(int32_t));
        read_bytes(egt->fp, (void *)&bb_n, sizeof(int32_t));
        if ( egt->cluster_records[i].aa_cluster_stats.N != aa_n ||
             egt->cluster_records[i].ab_cluster_stats.N != ab_n ||
             egt->cluster_records[i].bb_cluster_stats.N != bb_n )
             error("Cluster counts don't match with EGT cluster file %s\n", egt->fn);
    }

    if ( egt->data_block_version == 9 ) read_bytes(egt->fp, NULL, egt->num_records * sizeof(float));
    if ( !heof(egt->fp) ) error("EGT reader did not reach the end of file %s at position %ld\n", egt->fn, htell(egt->fp));

    return egt;
}

static void egt_destroy(egt_t *egt)
{
    if (!egt) return;
    free(egt->fn);
    if ( hclose(egt->fp) < 0 ) error("Error closing EGT file\n");
    free(egt->gencall_version);
    free(egt->cluster_version);
    free(egt->call_version);
    free(egt->normalization_version);
    free(egt->date_created);
    free(egt->manifest_name);
    free(egt->cluster_records);
    for (int i=0; i<egt->num_records; i++) free(egt->loci_names[i]);
    free(egt->loci_names);
    free(egt);
}

static void egt_summary(const egt_t *egt, FILE *stream)
{
    fprintf(stream, "EGT cluster file version = %d\n", egt->version);
    fprintf(stream, "GenCall version = %s\n", egt->gencall_version);
    fprintf(stream, "Clustering algorithm version = %s\n", egt->cluster_version);
    fprintf(stream, "Genotyping algorithm version = %s\n", egt->call_version);
    fprintf(stream, "Normalization algorithm version = %s\n", egt->normalization_version);
    fprintf(stream, "Date the cluster file was created = %s\n", egt->date_created);
    fprintf(stream, "Manifest name used to build this cluster file = %s\n", egt->manifest_name);
    fprintf(stream, "Data block version = %d\n", egt->data_block_version);
    fprintf(stream, "Number of records = %d\n", egt->num_records);
}

static void egt_to_csv(const egt_t *egt, FILE *stream)
{
    fprintf(stream, "Illumina, Inc.\n");
    fprintf(stream, "[Heading]\n");
    fprintf(stream, "Descriptor File Name,%s\n", strrchr(egt->fn, '/') ? strrchr(egt->fn, '/') + 1 : egt->fn);
    fprintf(stream, "GenCall version,%s\n", egt->gencall_version);
    fprintf(stream, "Clustering algorithm version,%s\n", egt->cluster_version);
    fprintf(stream, "Genotyping algorithm version,%s\n", egt->call_version);
    fprintf(stream, "Normalization algorithm version,%s\n", egt->normalization_version);
    fprintf(stream, "Date Manufactured,%s\n", egt->date_created);
    fprintf(stream, "Records Count ,%d\n", egt->num_records);
    fprintf(stream, "[Assay]\n");
    fprintf(stream, "Name,AA.N,AA.R_dev,AA.R_mean,AA.Theta_dev,AA.Theta_mean,AB.N,AB.R_dev,AB.R_mean,AB.Theta_dev,AB.Theta_mean,BB.N,BB.R_dev,BB.R_mean,BB.Theta_dev,BB.Theta_mean,Intensity Threshold,Cluster Separation,GenTrain Score,Original Score,Edited,Address\n");
    for (int i=0; i<egt->num_records; i++)
    {
        ClusterRecord *cluster_record = &egt->cluster_records[i];
        fprintf(stream, "%s,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\n",
            egt->loci_names[i],
            cluster_record->aa_cluster_stats.N,
            cluster_record->aa_cluster_stats.r_dev,
            cluster_record->aa_cluster_stats.r_mean,
            cluster_record->aa_cluster_stats.theta_dev,
            cluster_record->aa_cluster_stats.theta_mean,
            cluster_record->ab_cluster_stats.N,
            cluster_record->ab_cluster_stats.r_dev,
            cluster_record->ab_cluster_stats.r_mean,
            cluster_record->ab_cluster_stats.theta_dev,
            cluster_record->ab_cluster_stats.theta_mean,
            cluster_record->bb_cluster_stats.N,
            cluster_record->bb_cluster_stats.r_dev,
            cluster_record->bb_cluster_stats.r_mean,
            cluster_record->bb_cluster_stats.theta_dev,
            cluster_record->bb_cluster_stats.theta_mean,
            cluster_record->intensity_threshold,
            cluster_record->cluster_score.cluster_separation,
            cluster_record->cluster_score.total_score,
            cluster_record->cluster_score.original_score,
            cluster_record->cluster_score.edited,
            cluster_record->address);
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

/****************************************
 * GTC FILE IMPLEMENTATION              *
 ****************************************/

// see GenotypeCalls.py from https://github.com/Illumina/BeadArrayFiles

#define NUM_SNPS 1
#define PLOIDY 2
#define PLOIDY_TYPE 3
#define SAMPLE_NAME 10
#define SAMPLE_PLATE 11
#define SAMPLE_WELL 12
#define CLUSTER_FILE 100
#define SNP_MANIFEST 101
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
#define B_ALLELE_FREQS 1012
#define LOGR_RATIOS 1013
#define PERCENTILES_X 1014
#define PERCENTILES_Y 1015
#define SLIDE_IDENTIFIER 1016

const char *code2genotype[] = {
    "NC",
    "AA",
    "AB",
    "BB",
    "NULL",
    "A",
    "B",
    "AAA",
    "AAB",
    "ABB",
    "BBB",
    "AAAA",
    "AAAB",
    "AABB",
    "ABBB",
    "BBBB",
    "AAAAA",
    "AAAAB",
    "AAABB",
    "AABBB",
    "ABBBB",
    "BBBBB",
    "AAAAAA",
    "AAAAAB",
    "AAAABB",
    "AAABBB",
    "AABBBB",
    "ABBBBB",
    "BBBBBB",
    "AAAAAAA",
    "AAAAAAB",
    "AAAAABB",
    "AAAABBB",
    "AAABBBB",
    "AABBBBB",
    "ABBBBBB",
    "BBBBBBB",
    "AAAAAAAA",
    "AAAAAAAB",
    "AAAAAABB",
    "AAAAABBB",
    "AAAABBBB",
    "AAABBBBB",
    "AABBBBBB",
    "ABBBBBBB",
    "BBBBBBBB"
};

typedef struct
{
    int32_t version;
    float offset_x;
    float offset_y;
    float scale_x;
    float scale_y;
    float shear;
    float theta;
    float reserved[6];
}
XForm;

typedef char BaseCall[2];

typedef struct
{
    uint16_t raw_x;
    uint16_t raw_y;
    float norm_x;
    float norm_y;
    float ilmn_theta;
    float ilmn_r;
    float baf;
    float lrr;
} intensities_t;

typedef struct
{
    char *scanner_name;
    int32_t pmt_green;
    int32_t pmt_red;
    char *scanner_version;
    char *imaging_user;
} ScannerData;

typedef struct
{
    float p50gc;
    int32_t num_calls;
    int32_t num_no_calls;
    int32_t num_intensity_only;
} SampleData;

typedef uint16_t Percentiles[3];

typedef struct
{
    char *fn;
    hFILE *fp;
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

    float *sin_theta; // precomputed sine transforms
    float *cos_theta; // precomputed cosine transforms

    buffer_array_t *raw_x;
    buffer_array_t *raw_y;
    buffer_array_t *genotypes;
    buffer_array_t *base_calls;
    buffer_array_t *genotype_scores;
    buffer_array_t *b_allele_freqs;
    buffer_array_t *logr_ratios;
}
gtc_t;

static int gtc_read(gtc_t *gtc, uint16_t id)
{
    int i;
    for (i=0; i<gtc->number_toc_entries && id != gtc->id[i]; i++);
    if (i == gtc->number_toc_entries) return -1;
    if ( id != NUM_SNPS && id != PLOIDY && id != PLOIDY_TYPE )
    {
        if ( hseek(gtc->fp, gtc->toc[i], SEEK_SET) < 0 ) error("Fail to seek to position %d in GTC %s file \n", gtc->toc[i], gtc->fn);
    }

    switch (id)
    {
        case NUM_SNPS: gtc->num_snps = gtc->toc[i]; break;
        case PLOIDY: gtc->ploidy = gtc->toc[i]; break;
        case PLOIDY_TYPE: gtc->ploidy = gtc->toc[i]; break;
        case SAMPLE_NAME: read_pfx_string(gtc->fp, &gtc->sample_name, NULL); break;
        case SAMPLE_PLATE: read_pfx_string(gtc->fp, &gtc->sample_plate, NULL); break;
        case SAMPLE_WELL: read_pfx_string(gtc->fp, &gtc->sample_well, NULL); break;
        case CLUSTER_FILE: read_pfx_string(gtc->fp, &gtc->cluster_file, NULL); break;
        case SNP_MANIFEST: read_pfx_string(gtc->fp, &gtc->snp_manifest, NULL); break;
        case IMAGING_DATE: read_pfx_string(gtc->fp, &gtc->imaging_date, NULL); break;
        case AUTOCALL_DATE: read_pfx_string(gtc->fp, &gtc->autocall_date, NULL); break;
        case AUTOCALL_VERSION: read_pfx_string(gtc->fp, &gtc->autocall_version, NULL); break;
        case NORMALIZATION_TRANSFORMS: read_pfx_array(gtc->fp, (void **)&gtc->normalization_transforms, &gtc->m_normalization_transforms, sizeof(XForm)); break;
        case CONTROLS_X: read_pfx_array(gtc->fp, (void **)&gtc->controls_x, &gtc->m_controls_x, sizeof(uint16_t)); break;
        case CONTROLS_Y: read_pfx_array(gtc->fp, (void **)&gtc->controls_y, &gtc->m_controls_y, sizeof(uint16_t)); break;
        case RAW_X: gtc->raw_x = buffer_array_init(gtc->fp, 0, sizeof(uint16_t)); break;
        case RAW_Y: gtc->raw_y = buffer_array_init(gtc->fp, 0, sizeof(uint16_t)); break;
        case GENOTYPES: gtc->genotypes = buffer_array_init(gtc->fp, 0, sizeof(uint8_t)); break;
        case BASE_CALLS: gtc->base_calls = buffer_array_init(gtc->fp, 0, sizeof(BaseCall)); break;
        case GENOTYPE_SCORES: gtc->genotype_scores = buffer_array_init(gtc->fp, 0, sizeof(float)); break;
        case SCANNER_DATA:
            read_pfx_string(gtc->fp, &gtc->scanner_data.scanner_name, NULL);
            read_bytes(gtc->fp, (void *)&gtc->scanner_data.pmt_green, sizeof(float));
            read_bytes(gtc->fp, (void *)&gtc->scanner_data.pmt_red, sizeof(float));
            read_pfx_string(gtc->fp, &gtc->scanner_data.scanner_version, NULL);
            read_pfx_string(gtc->fp, &gtc->scanner_data.imaging_user, NULL);
            break;
        case CALL_RATE: read_bytes(gtc->fp, (void *)&gtc->call_rate, sizeof(float)); break;
        case GENDER: read_bytes(gtc->fp, (void *)&gtc->gender, sizeof(char)); break;
        case LOGR_DEV: read_bytes(gtc->fp, (void *)&gtc->logr_dev, sizeof(float)); break;
        case GC10: read_bytes(gtc->fp, (void *)&gtc->p10gc, sizeof(float)); break;
        case DX: read_bytes(gtc->fp, (void *)&gtc->dx, sizeof(int32_t)); break;
        case SAMPLE_DATA: read_bytes(gtc->fp, (void *)&gtc->sample_data, sizeof(SampleData)); break;
        case B_ALLELE_FREQS: gtc->b_allele_freqs = buffer_array_init(gtc->fp, 0, sizeof(float)); break;
        case LOGR_RATIOS: gtc->logr_ratios = buffer_array_init(gtc->fp, 0, sizeof(float)); break;
        case PERCENTILES_X: read_bytes(gtc->fp, (void *)&gtc->percentiles_x, sizeof(Percentiles)); break;
        case PERCENTILES_Y: read_bytes(gtc->fp, (void *)&gtc->percentiles_y, sizeof(Percentiles)); break;
        case SLIDE_IDENTIFIER: read_pfx_string(gtc->fp, &gtc->sentrix_id, NULL); break;
        default:
            error("GTC file format does not support TOC entry %d\n", id);
            break;
    }
    return 0;
}

static gtc_t *gtc_init(const char *fn)
{
    gtc_t *gtc = (gtc_t *)calloc(1, sizeof(gtc_t));
    gtc->fn = strdup(fn);
    gtc->fp = hopen(gtc->fn, "rb");
    if ( gtc->fp == NULL ) error("Could not open %s: %s\n", gtc->fn, strerror(errno));
    if ( is_gzip(gtc->fp) ) error("File %s is gzip compressed and currently cannot be seeked\n", gtc->fn);

    uint8_t buffer[4];
    if ( hread(gtc->fp, (void *)buffer, 4) < 4 ) error("Failed to read magic number from %s file\n", gtc->fn);
    if ( memcmp(buffer, "gtc", 3) != 0 ) error("GTC file %s format identifier is bad\n", gtc->fn);
    if ( buffer[3]>5 && buffer[3]<3 ) error("GTC file %s version %d is unsupported\n", gtc->fn, buffer[3]);
    gtc->version = (int32_t)buffer[3];

    read_bytes(gtc->fp, (void *)&gtc->number_toc_entries, sizeof(int32_t));
    gtc->id = (uint16_t *)malloc(gtc->number_toc_entries * sizeof(uint16_t));
    gtc->toc = (int32_t *)malloc(gtc->number_toc_entries * sizeof(int32_t));
    for (int i=0; i<gtc->number_toc_entries; i++)
    {
        read_bytes(gtc->fp, (void *)&gtc->id[i], sizeof(uint16_t));
        read_bytes(gtc->fp, (void *)&gtc->toc[i], sizeof(int32_t));
    }

    for (int i=0; i<gtc->number_toc_entries; i++) gtc_read(gtc, gtc->id[i]);

    gtc->sin_theta = (float *)malloc(gtc->m_normalization_transforms * sizeof(float));
    gtc->cos_theta = (float *)malloc(gtc->m_normalization_transforms * sizeof(float));
    for (int i=0; i<gtc->m_normalization_transforms; i++)
    {
        gtc->sin_theta[i] = sinf(gtc->normalization_transforms[i].theta);
        gtc->cos_theta[i] = cosf(gtc->normalization_transforms[i].theta);
    }

    return gtc;
}

static void gtc_destroy(gtc_t *gtc)
{
    if (!gtc) return;
    free(gtc->fn);
    if ( hclose(gtc->fp) < 0 ) error("Error closing GTC file\n");
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

    free(gtc->sin_theta);
    free(gtc->cos_theta);

    buffer_array_destroy(gtc->raw_x);
    buffer_array_destroy(gtc->raw_y);
    buffer_array_destroy(gtc->genotypes);
    buffer_array_destroy(gtc->base_calls);
    buffer_array_destroy(gtc->genotype_scores);
    buffer_array_destroy(gtc->b_allele_freqs);
    buffer_array_destroy(gtc->logr_ratios);
    free(gtc);
}

static void gtc_summary(const gtc_t *gtc, FILE *stream)
{
    fprintf(stream, "GTC genotype file version = %d\n", gtc->version);
    fprintf(stream, "Number of TOC entries = %d\n", gtc->number_toc_entries);
    fprintf(stream, "Number of SNPs = %d\n", gtc->num_snps);
    fprintf(stream, "Ploidy = %d\n", gtc->ploidy);
    fprintf(stream, "Ploidy Type = %d\n", gtc->ploidy_type);
    fprintf(stream, "Sample name = %s\n", gtc->sample_name);
    fprintf(stream, "Sample plate = %s\n", gtc->sample_plate);
    fprintf(stream, "Sample well = %s\n", gtc->sample_well);
    fprintf(stream, "Cluster file = %s\n", gtc->cluster_file);
    fprintf(stream, "SNP manifest = %s\n", gtc->snp_manifest);
    fprintf(stream, "Imaging date = %s\n", gtc->imaging_date);
    fprintf(stream, "AutoCall date = %s\n", gtc->autocall_date);
    fprintf(stream, "AutoCall version = %s\n", gtc->autocall_version);
    fprintf(stream, "Number of normalization transforms = %ld\n", gtc->m_normalization_transforms);
    fprintf(stream, "Number of controls X = %ld\n", gtc->m_controls_x);
    fprintf(stream, "Number of controls Y = %ld\n", gtc->m_controls_y);
    fprintf(stream, "Name of the scanner = %s\n", gtc->scanner_data.scanner_name);
    fprintf(stream, "Pmt Green = %d\n", gtc->scanner_data.pmt_green);
    fprintf(stream, "Pmt Red = %d\n", gtc->scanner_data.pmt_red);
    fprintf(stream, "Version of the scanner software used = %s\n", gtc->scanner_data.scanner_version);
    fprintf(stream, "Name of the scanner user = %s\n", gtc->scanner_data.imaging_user);
    fprintf(stream, "Call Rate = %f\n", gtc->call_rate);
    fprintf(stream, "Gender = %c\n", gtc->gender);
    fprintf(stream, "LogR deviation = %f\n", gtc->logr_dev);
    fprintf(stream, "GenCall score - 10th percentile = %f\n", gtc->p10gc);
    fprintf(stream, "DX = %d\n", gtc->dx);
    fprintf(stream, "GenCall score - 50th percentile = %f\n", gtc->sample_data.p50gc);
    fprintf(stream, "Number of valid calls = %d\n", gtc->sample_data.num_calls);
    fprintf(stream, "Number of invalid calls = %d\n", gtc->sample_data.num_no_calls);
    fprintf(stream, "Number of loci that are \"Intensity Only\" or \"Zeroed\" = %d\n", gtc->sample_data.num_intensity_only);
    fprintf(stream, "P05 X = %d\n", gtc->percentiles_x[0]);
    fprintf(stream, "P50 X = %d\n", gtc->percentiles_x[1]);
    fprintf(stream, "P95 X = %d\n", gtc->percentiles_x[2]);
    fprintf(stream, "P05 Y = %d\n", gtc->percentiles_y[0]);
    fprintf(stream, "P50 Y = %d\n", gtc->percentiles_y[1]);
    fprintf(stream, "P95 Y = %d\n", gtc->percentiles_y[2]);
    fprintf(stream, "Sentrix identifier for the slide = %s\n", gtc->sentrix_id);
}

static void gtc_to_csv(const gtc_t *gtc, FILE *stream)
{
    fprintf(stream, "Illumina, Inc.\n");
    fprintf(stream, "[Heading]\n");
    fprintf(stream, "Descriptor File Name,%s\n", strrchr(gtc->fn, '/') ? strrchr(gtc->fn, '/') + 1 : gtc->fn);
    fprintf(stream, "Cluster File,%s\n", gtc->cluster_file);
    fprintf(stream, "SNP Manifest,%s\n", gtc->snp_manifest);
    fprintf(stream, "Imaging Date,%s\n", gtc->imaging_date);
    fprintf(stream, "AutoCall Date,%s\n", gtc->autocall_date);
    fprintf(stream, "AutoCall Version,%s\n", gtc->autocall_version);
    fprintf(stream, "SNP Count ,%d\n", gtc->num_snps);
    fprintf(stream, "[Assay]\n");
    fprintf(stream, "Raw X,Raw Y,GType,Top Alleles,Score,B Allele Freq,Log R Ratio\n");
    for (int i=0; i<gtc->num_snps; i++)
    {
        uint16_t raw_x = 0, raw_y = 0;
        get_element(gtc->raw_x, (void *)&raw_x, i);
        get_element(gtc->raw_y, (void *)&raw_y, i);
        uint8_t genotype = 0;
        get_element(gtc->genotypes, (void *)&genotype, i);
        BaseCall base_call = {'-', '-'};
        get_element(gtc->base_calls, (void *)&base_call, i);
        float genotype_score = NAN, b_allele_freq = NAN, logr_ratio = NAN;
        get_element(gtc->genotype_scores, (void *)&genotype_score, i);
        get_element(gtc->b_allele_freqs, (void *)&b_allele_freq, i);
        get_element(gtc->logr_ratios, (void *)&logr_ratio, i);
        fprintf(stream, "%d,%d,%s,%c%c,%f,%f,%f\n", raw_x, raw_y, code2genotype[genotype], base_call[0], base_call[1],
            genotype_score, b_allele_freq, logr_ratio);
    }
    fprintf(stream, "[Normalization Transforms]\n");
    fprintf(stream, "Version,Offset X,Offset Y,Scale X,Scale Y,Shear,Theta\n");
    for (int i=0; i<gtc->m_normalization_transforms; i++)
    {
        fprintf(stream, "%d,%f,%f,%f,%f,%f,%f\n",
            gtc->normalization_transforms[i].version,
            gtc->normalization_transforms[i].offset_x,
            gtc->normalization_transforms[i].offset_y,
            gtc->normalization_transforms[i].scale_x,
            gtc->normalization_transforms[i].scale_y,
            gtc->normalization_transforms[i].shear,
            gtc->normalization_transforms[i].theta);
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

/****************************************
 * INTENSITIES COMPUTATIONS             *
 ****************************************/

// compute normalized X Y intensities
static void get_norm_xy(uint16_t raw_x,
                        uint16_t raw_y,
                        const gtc_t *gtc,
                        const bpm_t *bpm,
                        int idx,
                        float *norm_x,
                        float *norm_y)
{
    if ( bpm->norm_lookups && bpm->locus_entries[idx].norm_id != 0xFF )
    {
        int norm_id = bpm->norm_lookups[bpm->locus_entries[idx].norm_id];
        XForm *xform = gtc->normalization_transforms + norm_id;
        float temp_x = (float)raw_x - xform->offset_x;
        float temp_y = (float)raw_y - xform->offset_y;
        float temp_x2 =  gtc->cos_theta[norm_id] * temp_x + gtc->sin_theta[norm_id] * temp_y;
        float temp_y2 = -gtc->sin_theta[norm_id] * temp_x + gtc->cos_theta[norm_id] * temp_y;
        float temp_x3 = temp_x2 - xform->shear * temp_y2;
        *norm_x = temp_x3 < 0.0f ? 0.0f : temp_x3 / xform->scale_x;
        *norm_y = temp_y2 < 0.0f ? 0.0f : temp_y2 / xform->scale_y;
    }
    else
    {
        *norm_x = NAN;
        *norm_y = NAN;
    }
}

// compute Theta and R from raw intensities
static inline void get_ilmn_theta_r(float norm_x,
                                    float norm_y,
                                    float *ilmn_theta,
                                    float *ilmn_r)
{
    *ilmn_theta = atanf(norm_y / norm_x) * (float)M_2_PI;
    *ilmn_r = norm_x + norm_y;
}

// compute BAF and LRR from Theta and R
static void get_baf_lrr(const float ilmn_theta,
                        const float ilmn_r,
                        const ClusterRecord *cluster_record,
                        float *baf,
                        float *lrr)
{
    float aa_theta = cluster_record->aa_cluster_stats.theta_mean;
    float ab_theta = cluster_record->ab_cluster_stats.theta_mean;
    float bb_theta = cluster_record->bb_cluster_stats.theta_mean;
    float aa_r = cluster_record->aa_cluster_stats.r_mean;
    float ab_r = cluster_record->ab_cluster_stats.r_mean;
    float bb_r = cluster_record->bb_cluster_stats.r_mean;

    // compute LRR and BAF
    if ( ilmn_theta == ab_theta )
    {
        *lrr = logf(ilmn_r / ab_r) * (float)M_LOG2E;
        *baf = 0.5f;
    }
    else if ( ilmn_theta < ab_theta )
    {
        float slope = ( aa_r - ab_r ) / ( aa_theta - ab_theta );
        float b = aa_r - ( aa_theta * slope );
        float r_ref = ( slope * ilmn_theta ) + b;
        *lrr = logf(ilmn_r / r_ref) * (float)M_LOG2E;
        *baf = ilmn_theta < aa_theta ? 0.0f : 0.5f - (ab_theta - ilmn_theta) * 0.5f / (ab_theta - aa_theta);
    }
    else if ( ilmn_theta > ab_theta )
    {
        float slope = ( ab_r - bb_r ) / ( ab_theta - bb_theta );
        float b = ab_r - ( ab_theta * slope );
        float r_ref = ( slope * ilmn_theta ) + b;
        *lrr = logf(ilmn_r / r_ref) * (float)M_LOG2E;
        *baf = ilmn_theta >= bb_theta ? 1.0f : 1.0f - (bb_theta - ilmn_theta) * 0.5f / (bb_theta - ab_theta);
    }
    else
    {
        *lrr = NAN;
        *baf = NAN;
    }
}

static inline void get_intensities(gtc_t *gtc,
                                   const bpm_t *bpm,
                                   const egt_t *egt,
                                   int idx,
                                   intensities_t *intensities)
{
    get_element(gtc->raw_x, (void *)&intensities->raw_x, idx);
    get_element(gtc->raw_y, (void *)&intensities->raw_y, idx);
    get_norm_xy(intensities->raw_x, intensities->raw_y, gtc, bpm, idx, &intensities->norm_x, &intensities->norm_y);
    get_ilmn_theta_r(intensities->norm_x, intensities->norm_y, &intensities->ilmn_theta, &intensities->ilmn_r);

    if ( bpm->norm_lookups && egt )
    {
        get_baf_lrr(intensities->ilmn_theta, intensities->ilmn_r, &egt->cluster_records[idx], &intensities->baf, &intensities->lrr);
    }
    else if ( gtc->b_allele_freqs && gtc->logr_ratios )
    {
        get_element(gtc->b_allele_freqs, (void *)&intensities->baf, idx);
        get_element(gtc->logr_ratios, (void *)&intensities->lrr, idx);
    }
    else
    {
        intensities->baf = NAN;
        intensities->lrr = NAN;
    }
}

static void adjust_clusters(const uint8_t *gts,
                            const float *ilmn_theta,
                            const float *ilmn_r,
                            int n,
                            ClusterRecord *cluster_record)
{
    cluster_record->aa_cluster_stats.N = 0;
    cluster_record->ab_cluster_stats.N = 0;
    cluster_record->bb_cluster_stats.N = 0;
    cluster_record->aa_cluster_stats.theta_mean *= 0.2f;
    cluster_record->ab_cluster_stats.theta_mean *= 0.2f;
    cluster_record->bb_cluster_stats.theta_mean *= 0.2f;
    cluster_record->aa_cluster_stats.r_mean *= 0.2f;
    cluster_record->ab_cluster_stats.r_mean *= 0.2f;
    cluster_record->bb_cluster_stats.r_mean *= 0.2f;

    for (int i=0; i<n; i++)
    {
        switch ( gts[i] )
        {
            case 1: // AA
                cluster_record->aa_cluster_stats.N++;
                cluster_record->aa_cluster_stats.theta_mean += ilmn_theta[i];
                cluster_record->aa_cluster_stats.r_mean += ilmn_r[i];
                break;
            case 2: // AB
                cluster_record->ab_cluster_stats.N++;
                cluster_record->ab_cluster_stats.theta_mean += ilmn_theta[i];
                cluster_record->ab_cluster_stats.r_mean += ilmn_r[i];
                break;
            case 3: // BB
                cluster_record->bb_cluster_stats.N++;
                cluster_record->bb_cluster_stats.theta_mean += ilmn_theta[i];
                cluster_record->bb_cluster_stats.r_mean += ilmn_r[i];
                break;
            default: // NC
                break;
        }
    }

    cluster_record->aa_cluster_stats.theta_mean /= ((float)cluster_record->aa_cluster_stats.N + 0.2f);
    cluster_record->ab_cluster_stats.theta_mean /= ((float)cluster_record->ab_cluster_stats.N + 0.2f);
    cluster_record->bb_cluster_stats.theta_mean /= ((float)cluster_record->bb_cluster_stats.N + 0.2f);
    cluster_record->aa_cluster_stats.r_mean /= ((float)cluster_record->aa_cluster_stats.N + 0.2f);
    cluster_record->ab_cluster_stats.r_mean /= ((float)cluster_record->ab_cluster_stats.N + 0.2f);
    cluster_record->bb_cluster_stats.r_mean /= ((float)cluster_record->bb_cluster_stats.N + 0.2f);
}

/****************************************
 * INDEL IMPLEMENTATION                 *
 ****************************************/

// see BPMRecord.py from https://github.com/Illumina/GTCtoVCF

static inline int is_indel(const LocusEntry *locusentry)
{
    return locusentry->snp[1] == 'D' || locusentry->snp[3] == 'D';
}

static inline void strupper(char *str) {
    char *s = str;
    while (*s) {
        *s = toupper((unsigned char) *s);
        s++;
     }
}

static inline char complement(char c)
{
    switch(toupper(c))
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'S': return 'S';
        case 'W': return 'W';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'B': return 'V';
        case 'D': return 'H';
        case 'H': return 'D';
        case 'V': return 'B';
        default:
            error("Sequence nucleotide character %c is not in IUPAC notation\n", c);
    }
}

static void reverse_complement(char *str)
{
    size_t len = strlen(str);
    for (size_t i=0; i<len/2; i++)
    {
        char tmp = str[i];
        str[i] = complement( str[len-i-1] );
        str[len-i-1] = complement( tmp );
    }
    if (len % 2 == 1) str[len/2] = complement( str[len/2] );
}

#define NUC_A (1<<0)
#define NUC_C (1<<1)
#define NUC_G (1<<2)
#define NUC_T (1<<3)

static inline int degeneracy(char c)
{
    switch(toupper(c))
    {
        case 'A': return NUC_A;
        case 'C': return NUC_C;
        case 'G': return NUC_G;
        case 'T': return NUC_T;
        case 'N': return NUC_A | NUC_C | NUC_G | NUC_T;
        case 'R': return NUC_A || NUC_G;
        case 'Y': return NUC_C || NUC_T;
        case 'S': return NUC_C || NUC_G;
        case 'W': return NUC_A || NUC_T;
        case 'K': return NUC_G || NUC_T;
        case 'M': return NUC_A || NUC_C;
        case 'B': return NUC_C || NUC_G || NUC_T;
        case 'D': return NUC_A || NUC_G || NUC_T;
        case 'H': return NUC_A || NUC_C || NUC_T;
        case 'V': return NUC_A || NUC_C || NUC_G;
        default:
            error("Sequence nucleotide character %c is not in IUPAC notation\n", c);
    }
}

static int max_suffix_match(const char *str1, const char *str2, size_t str1_len, size_t str2_len)
{
    int ret = 0;
    for (size_t i=str1_len, j=str2_len; i > 0 && j > 0; i--, j--)
    {
        if ( degeneracy(str1[i-1]) & degeneracy(str2[j-1]) )
            ret++;
        else
        {
            return ret;
        }
    }
    return ret;
}

static int max_prefix_match(const char *str1, const char *str2, size_t str1_len, size_t str2_len)
{
    int ret = 0;
    for (size_t i=0, j=0; i < str1_len && j < str2_len; i++, j++)
    {
        if ( degeneracy(str1[i]) & degeneracy(str2[j]) )
            ret++;
        else
            return ret;
    }
    return ret;
}

typedef struct
{
    char *five_prime;
    const char *indel;
    const char *three_prime;
    const char *end;
}
indel_t;

static int inline five_prime_len(indel_t *self) { return (int)(self->indel - self->five_prime); }
static int inline indel_len(indel_t *self) { return (int)(self->three_prime - self->indel); }
static int inline three_prime_len(indel_t *self) { return (int)(self->end - self->three_prime); }

static void indel_reverse_complement(indel_t *self)
{
    reverse_complement(self->five_prime);
    int three_prime_n = three_prime_len(self);
    int indel_n = indel_len(self);
    self->indel = self->five_prime + three_prime_n;
    self->three_prime = self->indel + indel_n;
}

static void indel_left_shift(indel_t *self)
{
    while( self->indel > self->five_prime && *(self->indel-1) == *(self->three_prime-1) )
    {
        self->indel--;
        self->three_prime--;
    }
}

static indel_t indel_source_seq_init(const char *source_seq, int generate_reverse_complement)
{
    indel_t self;
    const char *left_pos = strchr(source_seq, '/');
    if ( *(left_pos-1) != '-' ) error("Not an indel: %s\n", source_seq);
    const char *right_pos = strchr(source_seq, ']');
    int five_prime_n = (int)( left_pos - source_seq - 2 );
    int indel_n = (int)( right_pos - left_pos - 1 );
    int three_prime_n = strlen( right_pos + 1 );
    self.five_prime = (char *)malloc(sizeof(char) * (five_prime_n + indel_n + three_prime_n + 1));
    sprintf(self.five_prime, "%.*s%.*s%s", five_prime_n, source_seq, indel_n, left_pos + 1, right_pos + 1);
    strupper(self.five_prime);
    self.indel = self.five_prime + five_prime_n;
    self.three_prime = self.indel + indel_n;
    self.end = self.three_prime + three_prime_n;
    if ( generate_reverse_complement ) indel_reverse_complement( &self );
    return self;
}

static indel_t indel_init(const char *five_prime,
                          int five_prime_n,
                          const char *indel,
                          int indel_n,
                          const char *three_prime,
                          int three_prime_n)
{
    indel_t self;
    self.five_prime = (char *)malloc(sizeof(char) * (five_prime_n + indel_n + three_prime_n + 1));
    sprintf(self.five_prime, "%.*s%.*s%.*s", five_prime_n, five_prime, indel_n, indel, three_prime_n, three_prime);
    strupper(self.five_prime);
    self.indel = self.five_prime + five_prime_n;
    self.three_prime = self.indel + indel_n;
    self.end = self.three_prime + three_prime_n;
    return self;
}

/****************************************
 * CONVERSION UTILITIES                 *
 ****************************************/

#define BPM_LOOKUPS     (1<<0)
#define EGT_LOADED      (1<<1)
#define ADJUST_CLUSTERS (1<<2)
#define FORMAT_IGC      (1<<3)
#define FORMAT_BAF      (1<<4)
#define FORMAT_LRR      (1<<5)
#define FORMAT_NORMX    (1<<6)
#define FORMAT_NORMY    (1<<7)
#define FORMAT_R        (1<<8)
#define FORMAT_THETA    (1<<9)
#define FORMAT_X        (1<<10)
#define FORMAT_Y        (1<<11)

static void gtcs_to_gs(gtc_t **gtc,
                       int n,
                       const bpm_t *bpm,
                       const egt_t *egt,
                       FILE *stream)
{
    // print header
    fprintf(stream, "Index\tName\tAddress\tChr\tPosition\tGenTrain Score\tFrac A\tFrac C\tFrac G\tFrac T");
    for (int i=0; i<n; i++) fprintf(stream, "\t%s.GType\t%s.Score\t%s.Theta\t%s.R\t%s.B Allele Freq\t%s.Log R Ratio\t%s.Top Alleles\t%s.Plus/Minus Alleles", gtc[i]->sample_name, gtc[i]->sample_name, gtc[i]->sample_name, gtc[i]->sample_name, gtc[i]->sample_name, gtc[i]->sample_name, gtc[i]->sample_name, gtc[i]->sample_name);
    fprintf(stream, "\n");

    // print loci
    for (int j=0; j<bpm->num_loci; j++)
    {
        fprintf(stream, "%d\t%s\t%d\t%s\t%s\t%f\t%f\t%f\t%f\t%f", bpm->indexes ? bpm->indexes[j] : j,
                                                                  bpm->locus_entries[j].name,
                                                                  bpm->locus_entries[j].address_a,
                                                                  bpm->locus_entries[j].chrom,
                                                                  bpm->locus_entries[j].map_info,
                                                                  egt ? egt->cluster_records[j].cluster_score.total_score : NAN,
                                                                  bpm->locus_entries[j].frac_a,
                                                                  bpm->locus_entries[j].frac_c,
                                                                  bpm->locus_entries[j].frac_g,
                                                                  bpm->locus_entries[j].frac_t);
        for (int i=0; i<n; i++)
        {
            uint8_t genotype;
            get_element(gtc[i]->genotypes, (void *)&genotype, j);
            float genotype_score;
            get_element(gtc[i]->genotype_scores, (void *)&genotype_score, j);
            BaseCall base_call;
            get_element(gtc[i]->base_calls, (void *)&base_call, j);
            intensities_t intensities;
            get_intensities(gtc[i], bpm, egt, j, &intensities);
            fprintf(stream, "\t%s\t%f\t%f\t%f\t%f\t%f\t%c%c\t--", code2genotype[genotype],
                                                                  genotype_score,
                                                                  intensities.ilmn_theta,
                                                                  intensities.ilmn_r,
                                                                  intensities.baf,
                                                                  intensities.lrr,
                                                                  base_call[0],
                                                                  base_call[1]);
        }
        fprintf(stream, "\n");
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

static int bcf_hdr_name2id_flexible(const bcf_hdr_t *hdr, char *chr)
{
    static kstring_t str = {0, 0, NULL};
    if (!str.s) kputs("chr", &str);
    if (!chr) { free(str.s); return 0; }
    int rid = bcf_hdr_name2id(hdr, chr);
    if ( rid >= 0 ) return rid;
    str.l = 3;
    kputs(chr, &str);
    rid = bcf_hdr_name2id(hdr, str.s);
    if ( rid >= 0 ) return rid;
    if ( strcmp(chr, "XY") == 0 || strcmp(chr, "XX") == 0 )
    {
        rid = bcf_hdr_name2id(hdr, "X");
        if ( rid >= 0 ) return rid;
        rid = bcf_hdr_name2id(hdr, "chrX");
    }
    else if ( strcmp(chr, "MT") == 0 )
    {
        rid = bcf_hdr_name2id(hdr, "chrM");
    }
    return rid;
}

static bcf_hdr_t *get_hdr(const faidx_t *fai, int flags)
{
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    int n = faidx_nseq(fai);
    for (int i=0; i<n; i++)
    {
        const char *seq = faidx_iseq(fai, i);
        int len = faidx_seq_len(fai, seq);
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq, len);
    }
    bcf_hdr_append(hdr, "##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description=\"A allele\">");
    bcf_hdr_append(hdr, "##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description=\"B allele\">");
    if ( flags & EGT_LOADED )
    {
        bcf_hdr_append(hdr, "##INFO=<ID=N_AA,Number=1,Type=Integer,Description=\"Number of AA calls in training set\">");
        bcf_hdr_append(hdr, "##INFO=<ID=N_AB,Number=1,Type=Integer,Description=\"Number of AB calls in training set\">");
        bcf_hdr_append(hdr, "##INFO=<ID=N_BB,Number=1,Type=Integer,Description=\"Number of BB calls in training set\">");
        bcf_hdr_append(hdr, "##INFO=<ID=devR_AA,Number=1,Type=Float,Description=\"Standard deviation of normalized R for AA cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=devR_AB,Number=1,Type=Float,Description=\"Standard deviation of normalized R for AB cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=devR_BB,Number=1,Type=Float,Description=\"Standard deviation of normalized R for BB cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=devTHETA_AA,Number=1,Type=Float,Description=\"Standard deviation of normalized THETA for AA cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=devTHETA_AB,Number=1,Type=Float,Description=\"Standard deviation of normalized THETA for AB cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=devTHETA_BB,Number=1,Type=Float,Description=\"Standard deviation of normalized THETA for BB cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=meanR_AA,Number=1,Type=Float,Description=\"Mean of normalized R for AA cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=meanR_AB,Number=1,Type=Float,Description=\"Mean of normalized R for AB cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=meanR_BB,Number=1,Type=Float,Description=\"Mean of normalized R for BB cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=meanTHETA_AA,Number=1,Type=Float,Description=\"Mean of normalized THETA for AA cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=meanTHETA_AB,Number=1,Type=Float,Description=\"Mean of normalized THETA for AB cluster\">");
        bcf_hdr_append(hdr, "##INFO=<ID=meanTHETA_BB,Number=1,Type=Float,Description=\"Mean of normalized THETA for BB cluster\">");
    }
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    if ( flags & FORMAT_IGC ) bcf_hdr_append(hdr, "##FORMAT=<ID=IGC,Number=1,Type=Float,Description=\"Illumina GenCall Confidence Score\">");
    if ( flags & FORMAT_BAF ) bcf_hdr_append(hdr, "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B Allele Frequency\">");
    if ( flags & FORMAT_LRR ) bcf_hdr_append(hdr, "##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R Ratio\">");
    if ( flags & BPM_LOOKUPS )
    {
        if ( flags & FORMAT_NORMX ) bcf_hdr_append(hdr, "##FORMAT=<ID=NORMX,Number=1,Type=Float,Description=\"Normalized X intensity\">");
        if ( flags & FORMAT_NORMY ) bcf_hdr_append(hdr, "##FORMAT=<ID=NORMY,Number=1,Type=Float,Description=\"Normalized Y intensity\">");
        if ( flags & FORMAT_R     ) bcf_hdr_append(hdr, "##FORMAT=<ID=R,Number=1,Type=Float,Description=\"Normalized R value\">");
        if ( flags & FORMAT_THETA ) bcf_hdr_append(hdr, "##FORMAT=<ID=THETA,Number=1,Type=Float,Description=\"Normalized Theta value\">");
    }
    if ( flags & FORMAT_X ) bcf_hdr_append(hdr, "##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"Raw X intensity\">");
    if ( flags & FORMAT_Y ) bcf_hdr_append(hdr, "##FORMAT=<ID=Y,Number=1,Type=Integer,Description=\"Raw Y intensity\">");
    return hdr;
}

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

static void gtcs_to_vcf(faidx_t *fai,
                        const bpm_t *bpm,
                        const egt_t *egt,
                        gtc_t **gtc,
                        int n,
                        htsFile *out_fh,
                        bcf_hdr_t *hdr,
                        int flags)
{
    if ( bcf_hdr_write(out_fh, hdr) < 0 ) error("Unable to write to output VCF file\n");
    bcf1_t *rec = bcf_init();
    rec->n_sample = n;
    bcf_float_set_missing(rec->qual);

    char reference_base[] = "\0\0";
    kstring_t allele_a = {0, 0, NULL};
    kstring_t allele_b = {0, 0, NULL};
    kstring_t indel_seq = {0, 0, NULL};

    uint8_t *gts = (uint8_t *) malloc(n * sizeof(uint8_t));
    int32_t *gt_arr = (int32_t *) malloc(n*2 * sizeof(int32_t));
    int32_t *gq_arr = (int32_t *) malloc(n * sizeof(int32_t));
    float *igc_arr = (float *) malloc(n * sizeof(float));
    float *baf_arr = (float *) malloc(n * sizeof(float));
    float *lrr_arr = (float *) malloc(n * sizeof(float));
    float *norm_x_arr = (float *) malloc(n * sizeof(float));
    float *norm_y_arr = (float *) malloc(n * sizeof(float));
    float *ilmn_r_arr = (float *) malloc(n * sizeof(float));
    float *ilmn_theta_arr = (float *) malloc(n * sizeof(float));
    int32_t *raw_x_arr = (int32_t *) malloc(n * sizeof(int32_t));
    int32_t *raw_y_arr = (int32_t *) malloc(n * sizeof(int32_t));

    int skipped_markers = 0;
    for (int j=0; j<bpm->num_loci; j++)
    {
        bcf_clear(rec);
        rec->rid = bcf_hdr_name2id_flexible(hdr, bpm->locus_entries[j].chrom);
        char *endptr;
        rec->pos = strtol(bpm->locus_entries[j].map_info, &endptr, 10) - 1;
        if ( bpm->locus_entries[j].map_info==endptr )
            error("Map info %s for marker %s is not understood\n",bpm->locus_entries[j].map_info, bpm->locus_entries[j].ilmn_id);
        if ( rec->rid < 0 || rec->pos < 0 )
        {
            fprintf(stderr, "Unable to determine location for marker %s with chromosome %s and position %s\n",
                bpm->locus_entries[j].ilmn_id, bpm->locus_entries[j].chrom, bpm->locus_entries[j].map_info);
            skipped_markers++;
            continue;
        }
        bcf_update_id( hdr, rec, bpm->locus_entries[j].name );

        int allele_a_idx, allele_b_idx;
        allele_a.l = allele_b.l = 0;
        if ( is_indel(&bpm->locus_entries[j]) )
        {
            if ( !bpm->locus_entries[j].source_seq )
            {
                fprintf(stderr, "Unable to determine reference allele for indel %s as source sequence is missing\n",
                    bpm->locus_entries[j].ilmn_id);
                skipped_markers++;
                continue;
            }
            // see BPMRecord.py from https://github.com/Illumina/GTCtoVCF
            char source_strand = toupper(bpm->locus_entries[j].source_strand[0]);
            char ilmn_strand = toupper(bpm->locus_entries[j].ilmn_strand[0]);
            char ref_strand = bpm->locus_entries[j].ref_strand[0];
            if ( ( source_strand != 'P' && source_strand != 'M' && source_strand != 'T' && source_strand != 'B' ) ||
                 ( ( source_strand == 'P' || source_strand == 'M' ) && ilmn_strand != 'P' && ilmn_strand != 'M' ) ||
                 ( ( source_strand == 'T' || source_strand == 'B' ) && ilmn_strand != 'T' && ilmn_strand != 'B' ) ||
                 ( ref_strand != '+' && ref_strand != '-' ) )
                error("Unable to process indel %s with source strand %s, Illumina strand %s, and reference strand %s\n",
                    bpm->locus_entries[j].ilmn_id, bpm->locus_entries[j].source_strand,
                    bpm->locus_entries[j].ilmn_strand, bpm->locus_entries[j].ref_strand);

            int reverse_complement = (source_strand == ilmn_strand) != (ref_strand == '+');
            indel_t source_seq = indel_source_seq_init(bpm->locus_entries[j].source_seq, reverse_complement);
            indel_seq.l = 0;
            ksprintf(&indel_seq, "%.*s", indel_len(&source_seq), source_seq.indel);
            indel_left_shift( &source_seq );

            int len;
            char *genomic_seq = faidx_fetch_seq(fai, bcf_seqname(hdr, rec),
                rec->pos-five_prime_len(&source_seq), rec->pos-1+indel_seq.l+three_prime_len(&source_seq), &len);
            strupper(genomic_seq);
            int indel_sequence_match = strncmp(genomic_seq + five_prime_len(&source_seq), indel_seq.s, indel_seq.l) == 0;

            indel_t genomic_del = indel_init(genomic_seq, five_prime_len(&source_seq), indel_seq.s, indel_seq.l,
                genomic_seq + five_prime_len(&source_seq) + indel_seq.l, three_prime_len(&source_seq));
            indel_left_shift( &genomic_del );
            int del_context_match_lengths_five_prime = max_suffix_match(genomic_del.five_prime,
                source_seq.five_prime, five_prime_len(&genomic_del), five_prime_len(&source_seq));
            int del_context_match_lengths_three_prime = max_prefix_match(genomic_del.three_prime,
                source_seq.three_prime, three_prime_len(&genomic_del), three_prime_len(&source_seq));
            int max_del_context = min( five_prime_len(&genomic_del), five_prime_len(&source_seq) ) +
                min( three_prime_len(&genomic_del), three_prime_len(&source_seq) ) + indel_seq.l;
            float del_context_score = (float)(del_context_match_lengths_five_prime + del_context_match_lengths_three_prime +
                ( indel_sequence_match ? indel_seq.l : 0) ) / (float)max_del_context;

            indel_t genomic_ins = indel_init(genomic_seq + 1, five_prime_len(&source_seq), indel_seq.s, indel_seq.l,
                genomic_seq + five_prime_len(&source_seq) + 1, three_prime_len(&source_seq));
            indel_left_shift( &genomic_ins );
            int ins_context_match_lengths_five_prime = max_suffix_match(genomic_ins.five_prime,
                source_seq.five_prime, five_prime_len(&genomic_ins), five_prime_len(&source_seq));
            int ins_context_match_lengths_three_prime = max_prefix_match(genomic_ins.three_prime,
                source_seq.three_prime, three_prime_len(&genomic_ins), three_prime_len(&source_seq));
            int max_ins_context = min( five_prime_len(&genomic_ins), five_prime_len(&source_seq) ) +
                min( three_prime_len(&genomic_ins), three_prime_len(&source_seq) );
            float ins_context_score = (float)(ins_context_match_lengths_five_prime + ins_context_match_lengths_three_prime) /
                (float)max_ins_context;

            int is_del = indel_sequence_match && del_context_score > ins_context_score &&
                ( del_context_match_lengths_five_prime >= 1 && del_context_match_lengths_three_prime >= 1);
            int is_ins = ins_context_score > del_context_score &&
                ( ins_context_match_lengths_five_prime >= 1 && ins_context_match_lengths_three_prime >= 1);

            if ( is_del ) rec->pos--;
            reference_base[0] = genomic_seq[five_prime_len(&source_seq) - is_del];
            free(source_seq.five_prime);
            free(genomic_del.five_prime);
            free(genomic_ins.five_prime);
            free(genomic_seq);

            if ( is_del == is_ins )
            {
                fprintf(stderr, "Unable to determine reference allele for indel %s with source sequence %s\n",
                    bpm->locus_entries[j].ilmn_id, bpm->locus_entries[j].source_seq);
                skipped_markers++;
                continue;
            }

            kputc(reference_base[0], &allele_a);
            kputc(reference_base[0], &allele_b);
            if ( bpm->locus_entries[j].snp[1] == 'I' )
            {
                ksprintf(&allele_a, "%s", indel_seq.s);
                allele_a_idx = !is_del;
                allele_b_idx = is_del;
            }
            else
            {
                ksprintf(&allele_b, "%s", indel_seq.s);
                allele_a_idx = is_del;
                allele_b_idx = !is_del;
            }

        }
        else
        {
            if ( strcmp(bpm->locus_entries[j].ref_strand, "+") == 0 )
            {
                kputc(bpm->locus_entries[j].snp[1], &allele_a);
                kputc(bpm->locus_entries[j].snp[3], &allele_b);
            }
            else if ( strcmp(bpm->locus_entries[j].ref_strand, "-") == 0 )
            {
                kputc(revnt(bpm->locus_entries[j].snp[1]), &allele_a);
                kputc(revnt(bpm->locus_entries[j].snp[3]), &allele_b);
            }
            int len;
            char *genomic_seq = faidx_fetch_seq(fai, bcf_seqname(hdr, rec), rec->pos, rec->pos, &len);
            reference_base[0] = genomic_seq[0];
            free(genomic_seq);
            if ( reference_base[0] == allele_a.s[0] )
            {
                allele_a_idx = 0;
                allele_b_idx = 1;
            }
            else if ( reference_base[0] == allele_b.s[0] )
            {
                allele_a_idx = 1;
                allele_b_idx = 0;
            }
            else
            {
                allele_a_idx = 1;
                allele_b_idx = 2;
            }
        }

        const char *alleles[3];
        int n_alleles;
        switch ( allele_b_idx )
        {
            case 0:
                alleles[0] = allele_b.s;
                alleles[1] = allele_a.s;
                n_alleles = 2;
                break;
            case 1:
                alleles[0] = allele_a.s;
                alleles[1] = allele_b.s;
                n_alleles = 2;
                break;
            case 2:
                alleles[0] = reference_base;
                alleles[1] = allele_a.s;
                alleles[2] = allele_b.s;
                n_alleles = 3;
                break;
            default:
                error("Unable to process marker %s\n", bpm->locus_entries[j].ilmn_id);
                break;
        }
        bcf_update_alleles(hdr, rec, alleles, n_alleles);
        bcf_update_info_int32(hdr, rec, "ALLELE_A", &allele_a_idx, 1);
        bcf_update_info_int32(hdr, rec, "ALLELE_B", &allele_b_idx, 1);

        for (int i=0; i<n; i++)
        {
            get_element(gtc[i]->genotypes, (void *)&gts[i], j);
            switch ( gts[i] )
            {
                case 1: // AA
                    gt_arr[2*i] = bcf_gt_unphased(allele_a_idx);
                    gt_arr[2*i+1] = bcf_gt_unphased(allele_a_idx);
                    break;
                case 2: // AB
                    gt_arr[2*i] = bcf_gt_unphased( ( allele_a_idx + allele_b_idx ) / 2 );
                    gt_arr[2*i+1] = bcf_gt_unphased( ( allele_a_idx + allele_b_idx + 1 ) / 2 );
                    break;
                case 3: // BB
                    gt_arr[2*i] = bcf_gt_unphased(allele_b_idx);
                    gt_arr[2*i+1] = bcf_gt_unphased(allele_b_idx);
                    break;
                default: // NC
                    gt_arr[2*i] = bcf_gt_missing;
                    gt_arr[2*i+1] = bcf_gt_missing;
                    break;
            }
            get_element(gtc[i]->genotype_scores, (void *)&igc_arr[i], j);
            gq_arr[i] = (int)(-10*log10(1 - igc_arr[i]) + .5);
            if ( gq_arr[i] < 0 ) gq_arr[i] = 0;
            if ( gq_arr[i] > 50 ) gq_arr[i] = 50;
            intensities_t intensities;
            get_intensities(gtc[i], bpm, egt, j, &intensities);
            baf_arr[i] = intensities.baf;
            lrr_arr[i] = intensities.lrr;
            norm_x_arr[i] = intensities.norm_x;
            norm_y_arr[i] = intensities.norm_y;
            ilmn_r_arr[i] = intensities.ilmn_r;
            ilmn_theta_arr[i] = intensities.ilmn_theta;
            raw_x_arr[i] = (int32_t)intensities.raw_x;
            raw_y_arr[i] = (int32_t)intensities.raw_y;
        }

        if ( flags & EGT_LOADED )
        {
            if ( flags & ADJUST_CLUSTERS )
            {
                adjust_clusters(gts, ilmn_theta_arr, ilmn_r_arr, n, &egt->cluster_records[j]);
                for (int i=0; i<n; i++)
                {
                    get_baf_lrr(ilmn_theta_arr[i], ilmn_r_arr[i], &egt->cluster_records[j], &baf_arr[i], &lrr_arr[i]);
                }
            }
            bcf_update_info_int32(hdr, rec, "N_AA", &egt->cluster_records[j].aa_cluster_stats.N, 1);
            bcf_update_info_int32(hdr, rec, "N_AB", &egt->cluster_records[j].ab_cluster_stats.N, 1);
            bcf_update_info_int32(hdr, rec, "N_BB", &egt->cluster_records[j].bb_cluster_stats.N, 1);
            bcf_update_info_float(hdr, rec, "devR_AA", &egt->cluster_records[j].aa_cluster_stats.r_dev, 1);
            bcf_update_info_float(hdr, rec, "devR_AB", &egt->cluster_records[j].ab_cluster_stats.r_dev, 1);
            bcf_update_info_float(hdr, rec, "devR_BB", &egt->cluster_records[j].bb_cluster_stats.r_dev, 1);
            bcf_update_info_float(hdr, rec, "devTHETA_AA", &egt->cluster_records[j].aa_cluster_stats.theta_dev, 1);
            bcf_update_info_float(hdr, rec, "devTHETA_AB", &egt->cluster_records[j].ab_cluster_stats.theta_dev, 1);
            bcf_update_info_float(hdr, rec, "devTHETA_BB", &egt->cluster_records[j].bb_cluster_stats.theta_dev, 1);
            bcf_update_info_float(hdr, rec, "meanR_AA", &egt->cluster_records[j].aa_cluster_stats.r_mean, 1);
            bcf_update_info_float(hdr, rec, "meanR_AB", &egt->cluster_records[j].ab_cluster_stats.r_mean, 1);
            bcf_update_info_float(hdr, rec, "meanR_BB", &egt->cluster_records[j].bb_cluster_stats.r_mean, 1);
            bcf_update_info_float(hdr, rec, "meanTHETA_AA", &egt->cluster_records[j].aa_cluster_stats.theta_mean, 1);
            bcf_update_info_float(hdr, rec, "meanTHETA_AB", &egt->cluster_records[j].ab_cluster_stats.theta_mean, 1);
            bcf_update_info_float(hdr, rec, "meanTHETA_BB", &egt->cluster_records[j].bb_cluster_stats.theta_mean, 1);
        }

        bcf_update_genotypes(hdr, rec, gt_arr, n*2);
        bcf_update_format_int32(hdr, rec, "GQ", gq_arr, n);
        if ( flags & FORMAT_IGC   ) bcf_update_format_float(hdr, rec, "IGC", igc_arr, n);
        if ( flags & FORMAT_BAF   ) bcf_update_format_float(hdr, rec, "BAF", baf_arr, n);
        if ( flags & FORMAT_LRR   ) bcf_update_format_float(hdr, rec, "LRR", lrr_arr, n);
        if ( flags & FORMAT_NORMX ) bcf_update_format_float(hdr, rec, "NORMX", norm_x_arr, n);
        if ( flags & FORMAT_NORMY ) bcf_update_format_float(hdr, rec, "NORMY", norm_y_arr, n);
        if ( flags & FORMAT_R     ) bcf_update_format_float(hdr, rec, "R", ilmn_r_arr, n);
        if ( flags & FORMAT_THETA ) bcf_update_format_float(hdr, rec, "THETA", ilmn_theta_arr, n);
        if ( flags & FORMAT_X     ) bcf_update_format_int32(hdr, rec, "X", raw_x_arr, n);
        if ( flags & FORMAT_Y     ) bcf_update_format_int32(hdr, rec, "Y", raw_y_arr, n);
        if ( bcf_write(out_fh, hdr, rec) < 0 ) error("Unable to write to output VCF file\n");
    }
    fprintf(stderr, "%d markers skipped in the conversion to VCF\n", skipped_markers);
    bcf_hdr_name2id_flexible(hdr, NULL);

    free(gts);
    free(gt_arr);
    free(gq_arr);
    free(igc_arr);
    free(baf_arr);
    free(lrr_arr);
    free(norm_x_arr);
    free(norm_y_arr);
    free(ilmn_r_arr);
    free(ilmn_theta_arr);
    free(raw_x_arr);
    free(raw_y_arr);

    free(allele_a.s);
    free(allele_b.s);
    free(indel_seq.s);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(out_fh);
}

#define STRAND_PLUS  0
#define STRAND_MINUS 1
#define STRAND_TOP   2
#define STRAND_BOT   3

#define GS_GT    0
#define GS_TOP   1
#define GS_IGC   2
#define GS_BAF   3
#define GS_LRR   4
#define GS_NORMX 5
#define GS_NORMY 6
#define GS_R     7
#define GS_THETA 8
#define GS_X     9
#define GS_Y     10

typedef struct
{
    int *col2sample;
    int type;
    void *ptr;
}
gs_col_t;

static int tsv_setter_gs_col(tsv_t *tsv,
                             bcf1_t *rec,
                             void *usr)
{
    gs_col_t *gs_col = (gs_col_t *)usr;
    int32_t *gts;
    char *endptr;
    switch ( gs_col->type )
    {
        case GS_GT:
            gts = (int32_t *)gs_col->ptr + 2 * gs_col->col2sample[ tsv->icol ];
            if ( tsv->ss[0]=='A' && tsv->ss[1]=='A' )
            {
                gts[0] = bcf_gt_unphased(0);
                gts[1] = bcf_gt_unphased(0);
            }
            else if ( tsv->ss[0]=='A' && tsv->ss[1]=='B' )
            {
                gts[0] = bcf_gt_unphased(0);
                gts[1] = bcf_gt_unphased(1);
            }
            else if ( tsv->ss[0]=='B' && tsv->ss[1]=='B' )
            {
                gts[0] = bcf_gt_unphased(1);
                gts[1] = bcf_gt_unphased(1);
            }
            else if ( tsv->ss[0]=='N' && tsv->ss[1]=='C' )
            {
                gts[0] = bcf_gt_missing;
                gts[1] = bcf_gt_missing;
            }
            else return -1;
            break;
        case GS_TOP:
            gts = (int32_t *)gs_col->ptr + 2 * gs_col->col2sample[ tsv->icol ];
            if ( gts[0] == bcf_gt_unphased(0) ) rec->d.allele[0][0] = tsv->ss[0];
            if ( gts[0] == bcf_gt_unphased(1) ) rec->d.allele[1][0] = tsv->ss[0];
            if ( gts[1] == bcf_gt_unphased(0) ) rec->d.allele[0][0] = tsv->ss[1];
            if ( gts[1] == bcf_gt_unphased(1) ) rec->d.allele[1][0] = tsv->ss[1];
            break;
        case GS_IGC:
        case GS_BAF:
        case GS_LRR:
        case GS_NORMX:
        case GS_NORMY:
        case GS_R:
        case GS_THETA:
            ((float *)gs_col->ptr + gs_col->col2sample[ tsv->icol ])[0] = strtof(tsv->ss, &endptr);
            if ( tsv->ss==endptr ) return -1;
            break;
        case GS_X:
        case GS_Y:
            ((int32_t *)gs_col->ptr + gs_col->col2sample[ tsv->icol ])[0] = strtod(tsv->ss, &endptr);
            if ( tsv->ss==endptr ) return -1;
            break;
        default:
            return -1;
    }
    return 0;
}

static int tsv_setter_chrom_flexible(tsv_t *tsv,
                                     bcf1_t *rec,
                                     void *usr)
{
    char tmp = *tsv->se;
    *tsv->se = 0;
    rec->rid = bcf_hdr_name2id_flexible((bcf_hdr_t*)usr, tsv->ss);
    *tsv->se = tmp;
    return rec->rid==-1 ? -1 : 0;
}

static int tsv_setter_strand(tsv_t *tsv,
                             bcf1_t *rec,
                             void *usr)
{
    if ( strncmp(tsv->ss, "PLUS", 4) == 0 ) ((int *)usr)[0] = STRAND_PLUS;
    if ( strncmp(tsv->ss, "MINUS", 5) == 0 ) ((int *)usr)[0] = STRAND_MINUS;
    else if ( strncmp(tsv->ss, "TOP", 3) == 0 ) ((int *)usr)[0] = STRAND_TOP;
    else if ( strncmp(tsv->ss, "BOT", 3) == 0 ) ((int *)usr)[0] = STRAND_BOT;
    else return -1;
    return 0;
}

static int tsv_setter_snp(tsv_t *tsv,
                          bcf1_t *rec,
                          void *usr)
{
    char **snp = (char **)usr;
    if ( strncmp(tsv->ss, "[N/A]", 5) == 0 ) *snp = NULL;
    else *snp = tsv->ss;
    return 0;
}

int tsv_register_all(tsv_t *tsv,
                     const char *id,
                     tsv_setter_t setter, void *usr)
{
    int i, n = 0;
    for (i=0; i<tsv->ncols; i++)
    {
        if ( !tsv->cols[i].name || strcasecmp(tsv->cols[i].name,id) ) continue;
        tsv->cols[i].setter = setter;
        tsv->cols[i].usr    = usr;
        n++;
    }
    return n ? 0 : -1;
}

int tsv_parse_delimiter(tsv_t *tsv,
                        bcf1_t *rec,
                        char *str,
                        int delimiter)
{
    int status = 0;
    tsv->icol = 0;
    tsv->ss = tsv->se = str;
    while ( *tsv->ss && tsv->icol < tsv->ncols )
    {
        if ( delimiter ) while ( *tsv->se && (*tsv->se)!=delimiter ) tsv->se++;
        else while ( *tsv->se && !isspace(*tsv->se) ) tsv->se++;
        if ( tsv->cols[tsv->icol].setter )
        {
            int ret = tsv->cols[tsv->icol].setter(tsv,rec,tsv->cols[tsv->icol].usr);
            if ( ret<0 ) return -1;
            status++;
        }
        if ( delimiter ) tsv->se++;
        else while ( *tsv->se && isspace(*tsv->se) ) tsv->se++;
        tsv->ss = tsv->se;
        tsv->icol++;
    }
    return status ? 0 : -1;
}

static void gs_to_vcf(htsFile *gs_fh,
                      htsFile *out_fh,
                      bcf_hdr_t *hdr,
                      int flags)
{
    // read the header of the table
    kstring_t line = {0, 0, NULL};
    if ( hts_getline(gs_fh, KS_SEP_LINE, &line) <= 0 ) error("Empty file: %s\n", gs_fh->fn);
    int moff = 0, *off = NULL, ncols = ksplit_core(line.s, '\t', &moff, &off);
    kstring_t str = {0, 0, NULL};
    int *col2sample = (int *) malloc(sizeof(int)*ncols);
    for (int i=0; i<ncols; i++)
    {
        char *ptr;
        if ( i>0 ) ksprintf(&str, ",");
        if ( ( ptr = strrchr(&line.s[off[i]], '.') ) )
        {
            *ptr++ = '\0';
            if( ( bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, &line.s[off[i]]) < 0 ) )
                bcf_hdr_add_sample(hdr, &line.s[off[i]]);
            if ( strcmp(ptr, "GType")==0 ) ksprintf(&str, "GT");
            else if ( strcmp(ptr, "Score")==0 ) ksprintf(&str, "IGC");
            else if ( strcmp(ptr, "Theta")==0 ) ksprintf(&str, "THETA");
            else if ( strcmp(ptr, "R")==0 ) ksprintf(&str, "R");
            else if ( strcmp(ptr, "X Raw")==0 ) ksprintf(&str, "X");
            else if ( strcmp(ptr, "Y Raw")==0 ) ksprintf(&str, "Y");
            else if ( strcmp(ptr, "X")==0 ) ksprintf(&str, "NORMX");
            else if ( strcmp(ptr, "Y")==0 ) ksprintf(&str, "NORMY");
            else if ( strcmp(ptr, "B Allele Freq")==0 ) ksprintf(&str, "BAF");
            else if ( strcmp(ptr, "Log R Ratio")==0 ) ksprintf(&str, "LRR");
            else if ( strcmp(ptr, "Top Alleles")==0 ) ksprintf(&str, "TOP");
            else if ( strcmp(ptr, "Plus/Minus Alleles")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Import Calls")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Concordance")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Orig Call")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "CNV Value")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "CNV Confidence")==0 ) ksprintf(&str, "-");
            else error("Could not recognize FORMAT field: %s\n", ptr);
            col2sample[ i ] = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, &line.s[off[i]]);
        }
        else
        {
            ptr = &line.s[off[i]];
            if ( strcmp(ptr, "Index")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Name")==0 ) ksprintf(&str, "ID");
            else if ( strcmp(ptr, "Address")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Chr")==0 || strcmp(ptr, "Chromosome")==0 ) ksprintf(&str, "CHROM");
            else if ( strcmp(ptr, "Manifest")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Position")==0 ) ksprintf(&str, "POS");
            else if ( strcmp(ptr, "GenTrain Score")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Frac A")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Frac C")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Frac G")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "Frac T")==0 ) ksprintf(&str, "-");
            else if ( strcmp(ptr, "IlmnStrand")==0 ) ksprintf(&str, "STRAND");
            else if ( strcmp(ptr, "SNP")==0 ) ksprintf(&str, "SNP");
            else error("Could not recognize INFO field: %s\n", ptr);
            col2sample[ i ] = -1;
        }
    }
    free(off);
    if ( bcf_hdr_sync( hdr ) < 0 ) error_errno("[%s] Failed to update header", __func__); // updates the number of samples
    int nsamples = bcf_hdr_nsamples(hdr);

    tsv_t *tsv = tsv_init(str.s);
    if ( tsv_register(tsv, "CHROM", tsv_setter_chrom_flexible, hdr) < 0 ) error("Expected CHROM column\n");
    if ( tsv_register(tsv, "POS", tsv_setter_pos, NULL) < 0 ) error("Expected POS column\n");
    tsv_register(tsv, "ID", tsv_setter_id, hdr);

    int strand;
    tsv_register(tsv, "STRAND", tsv_setter_strand, &strand);
    char *snp;
    tsv_register(tsv, "SNP", tsv_setter_snp, &snp);

    int32_t *gt_arr = (int32_t *)malloc(nsamples*2 * sizeof(int32_t));
    gs_col_t gs_gt_arr = {col2sample, GS_GT, gt_arr};
    if ( tsv_register_all(tsv, "GT", tsv_setter_gs_col, &gs_gt_arr) < 0 ) error("Expected GType column\n");
    gs_col_t gs_top = {col2sample, GS_TOP, gt_arr};
    tsv_register_all(tsv, "TOP", tsv_setter_gs_col, &gs_top);

    gs_col_t gs_igc = {col2sample, GS_IGC, NULL};
    if ( (flags & FORMAT_IGC) && tsv_register_all(tsv, "IGC", tsv_setter_gs_col, &gs_igc) == 0 ) gs_igc.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_baf = {col2sample, GS_BAF, NULL};
    if ( (flags & FORMAT_BAF) && tsv_register_all(tsv, "BAF", tsv_setter_gs_col, &gs_baf) == 0 ) gs_baf.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_lrr = {col2sample, GS_LRR, NULL};
    if ( (flags & FORMAT_LRR) && tsv_register_all(tsv, "LRR", tsv_setter_gs_col, &gs_lrr) == 0 ) gs_lrr.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_norm_x = {col2sample, GS_NORMX, NULL};
    if ( (flags & FORMAT_NORMX) && tsv_register_all(tsv, "NORMX", tsv_setter_gs_col, &gs_norm_x) == 0 ) gs_norm_x.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_norm_y = {col2sample, GS_NORMY, NULL};
    if ( (flags & FORMAT_NORMY) && tsv_register_all(tsv, "NORMY", tsv_setter_gs_col, &gs_norm_y) == 0 ) gs_norm_y.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_ilmn_r = {col2sample, GS_R, NULL};
    if ( (flags & FORMAT_R) && tsv_register_all(tsv, "R", tsv_setter_gs_col, &gs_ilmn_r) == 0 ) gs_ilmn_r.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_ilmn_theta = {col2sample, GS_THETA, NULL};
    if ( (flags & FORMAT_THETA) && tsv_register_all(tsv, "THETA", tsv_setter_gs_col, &gs_ilmn_theta) == 0 ) gs_ilmn_theta.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_raw_x = {col2sample, GS_X, NULL};
    if ( (flags & FORMAT_X) && tsv_register_all(tsv, "X", tsv_setter_gs_col, &gs_raw_x) == 0 ) gs_raw_x.ptr = malloc(nsamples * sizeof(int32_t));

    gs_col_t gs_raw_y = {col2sample, GS_Y, NULL};
    if ( (flags & FORMAT_Y) && tsv_register_all(tsv, "Y", tsv_setter_gs_col, &gs_raw_y) == 0 ) gs_raw_y.ptr = malloc(nsamples * sizeof(int32_t));

    if ( bcf_hdr_write(out_fh, hdr) < 0 ) error("Unable to write to output VCF file\n");

    bcf1_t *rec = bcf_init();
    // bcf_float_set_missing(rec->qual);
    int n_total = 0, n_skipped = 0;
    while ( hts_getline(gs_fh, KS_SEP_LINE, &line) > 0 )
    {
        if ( line.s[0]=='#' ) continue;     // skip comments
        strand = STRAND_PLUS;
        snp = NULL;
        bcf_clear(rec);
        rec->n_sample = nsamples;
        bcf_update_alleles_str(hdr, rec, "N,N");

        n_total++;
        if ( !tsv_parse_delimiter(tsv, rec, line.s, '\t') )
        {
            if ( snp )
            {
                if ( strand == STRAND_TOP )
                {
                    rec->d.allele[0][0] = snp[1];
                    rec->d.allele[1][0] = snp[3];
                }
                else if ( strand == STRAND_BOT )
                {
                    rec->d.allele[0][0] = revnt(snp[1]);
                    rec->d.allele[1][0] = revnt(snp[3]);
                }
            }
            bcf_update_genotypes(hdr, rec, gt_arr, nsamples*2);
            if ( gs_igc.ptr ) bcf_update_format_float(hdr, rec, "IGC", (float *)gs_igc.ptr, nsamples);
            if ( gs_baf.ptr ) bcf_update_format_float(hdr, rec, "BAF", (float *)gs_baf.ptr, nsamples);
            if ( gs_lrr.ptr ) bcf_update_format_float(hdr, rec, "LRR", (float *)gs_lrr.ptr, nsamples);
            if ( gs_norm_x.ptr ) bcf_update_format_float(hdr, rec, "NORMX", (float *)gs_norm_x.ptr, nsamples);
            if ( gs_norm_y.ptr ) bcf_update_format_float(hdr, rec, "NORMY", (float *)gs_norm_y.ptr, nsamples);
            if ( gs_ilmn_r.ptr ) bcf_update_format_float(hdr, rec, "R", (float *)gs_ilmn_r.ptr, nsamples);
            if ( gs_ilmn_theta.ptr ) bcf_update_format_float(hdr, rec, "THETA", (float *)gs_ilmn_theta.ptr, nsamples);
            if ( gs_raw_x.ptr ) bcf_update_format_int32(hdr, rec, "X", (int32_t *)gs_raw_x.ptr, nsamples);
            if ( gs_raw_y.ptr ) bcf_update_format_int32(hdr, rec, "Y", (int32_t *)gs_raw_y.ptr, nsamples);
            if ( bcf_write(out_fh, hdr, rec) < 0 ) error("Unable to write to output VCF file\n");
        }
        else
            n_skipped++;
    }
    if ( hts_close(gs_fh) ) error("Close failed: %s\n", gs_fh->fn);
    free(line.s);

    bcf_hdr_destroy(hdr);
    hts_close(out_fh);
    tsv_destroy(tsv);
    bcf_destroy(rec);
    free(col2sample);
    free(gt_arr);
    free(gs_igc.ptr);
    free(gs_baf.ptr);
    free(gs_lrr.ptr);
    free(gs_norm_x.ptr);
    free(gs_norm_y.ptr);
    free(gs_ilmn_r.ptr);
    free(gs_ilmn_theta.ptr);
    free(gs_raw_x.ptr);
    free(gs_raw_y.ptr);

    fprintf(stderr, "Rows total: \t%d\n", n_total);
    fprintf(stderr, "Rows skipped: \t%d\n", n_skipped);
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void)
{
    return "Convert GTC files to VCF.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: convert Illumina GTC files containing intensity data into VCF. (version "GTC2VCF_VERSION" https://github.com/freeseek/gtc2vcf)\n"
        "\n"
        "Usage: bcftools +gtc2vcf [options] <A.gtc> [...]\n"
        "\n"
        "Plugin options:\n"
        "    -l, --list-tags                    list available tags with description for VCF output\n"
        "    -t, --tags LIST                    list of output tags [IGC,BAF,LRR,NORMX,NORMY,R,THETA,X,Y]\n"
        "    -i  --idat <file>                  IDAT file\n"
        "    -b  --bpm <file>                   BPM manifest file\n"
        "    -c  --csv <file>                   CSV manifest file\n"
        "    -e  --egt <file>                   EGT cluster file\n"
        "    -f, --fasta-ref <file>             reference sequence in fasta format\n"
        "        --set-cache-size <int>         select fasta cache size in bytes\n"
        "    -g, --gtc-list <file>              read GTC file names from file\n"
        "        --adjust-clusters              adjust cluster centers in (Theta, R) space\n"
        "    -x, --sex <file>                   output GenCall gender estimate into file\n"
        "        --do-not-check-bpm             do not check whether BPM and GTC files match manifest file name\n"
        "        --genome-studio                input a genome studio final report file (in matrix format)\n"
        "        --no-version                   do not append version and command line to the header\n"
        "    -o, --output <file>                write output to a file [standard output]\n"
        "    -O, --output-type b|u|z|v|g        b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF, g GenomeStudio [v]\n"
        "        --threads <int>                number of extra output compression threads [0]\n"
        "\n";
}

static FILE *get_file_handle(const char *str)
{
    FILE *ret;
    if ( strcmp(str, "-") == 0 )
        ret = stdout;
    else
    {
        ret = fopen(str, "w");
        if ( !ret ) error("Failed to open %s: %s\n", str, strerror(errno));
    }
    return ret;
}

static int parse_tags(const char *str)
{
    int flags = 0, n;
    char **tags = hts_readlist(str, 0, &n);
    for (int i=0; i<n; i++)
    {
        if ( !strcasecmp(tags[i], "IGC") ) flags |= FORMAT_IGC;
        else if ( !strcasecmp(tags[i], "X") ) flags |= FORMAT_X;
        else if ( !strcasecmp(tags[i], "Y") ) flags |= FORMAT_Y;
        else if ( !strcasecmp(tags[i], "NORMX") ) flags |= FORMAT_NORMX;
        else if ( !strcasecmp(tags[i], "NORMY") ) flags |= FORMAT_NORMY;
        else if ( !strcasecmp(tags[i], "R") ) flags |= FORMAT_R;
        else if ( !strcasecmp(tags[i], "THETA") ) flags |= FORMAT_THETA;
        else if ( !strcasecmp(tags[i], "LRR") ) flags |= FORMAT_LRR;
        else if ( !strcasecmp(tags[i], "BAF") ) flags |= FORMAT_BAF;
        else error("Error parsing \"--tags %s\": the tag \"%s\" is not supported\n", str, tags[i]);
        free(tags[i]);
    }
    if (n) free(tags);
    return flags;
}

void list_tags(void)
{
    error(
        "FORMAT/IGC      Number:1  Type:Float    ..  Illumina GenCall Confidence Score\n"
        "FORMAT/BAF      Number:1  Type:Float    ..  B Allele Frequency\n"
        "FORMAT/LRR      Number:1  Type:Float    ..  Log R Ratio\n"
        "FORMAT/NORMX    Number:1  Type:Float    ..  Normalized X intensity\n"
        "FORMAT/NORMY    Number:1  Type:Float    ..  Normalized Y intensity\n"
        "FORMAT/R        Number:1  Type:Float    ..  Normalized R value\n"
        "FORMAT/THETA    Number:1  Type:Float    ..  Normalized Theta value\n"
        "FORMAT/X        Number:1  Type:Integer  ..  Raw X intensity\n"
        "FORMAT/Y        Number:1  Type:Integer  ..  Raw Y intensity\n"
        );
}

int run(int argc, char *argv[])
{
    char *tag_list = "IGC,BAF,LRR,NORMX,NORMY,R,THETA,X,Y";
    char *idat_fname = NULL;
    char *bpm_fname = NULL;
    char *csv_fname = NULL;
    char *egt_fname = NULL;
    char *gs_fname = NULL;
    char *output_fname = "-";
    char *ref_fname = NULL;
    char *gtc_list = NULL;
    char *sex_fname = NULL;
    int flags = 0;
    int output_type = FT_VCF;
    int cache_size = 0;
    int bpm_check = 1;
    int n_threads = 0;
    int record_cmd_line = 1;
    int binary_to_csv = 0;
    faidx_t *fai = NULL;
    htsFile *out_fh = NULL;
    FILE *out_txt = NULL;
    FILE *out_sex = NULL;

    static struct option loptions[] =
    {
        {"list-tags", no_argument, NULL, 'l'},
        {"tags", required_argument, NULL, 't'},
        {"idat", required_argument, NULL, 'i'},
        {"bpm", required_argument, NULL, 'b'},
        {"csv", required_argument, NULL, 'c'},
        {"egt", required_argument, NULL, 'e'},
        {"output", required_argument, NULL, 'o'},
        {"output-type", required_argument, NULL, 'O'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"set-cache-size", required_argument, NULL, 1},
        {"gtc-list", required_argument, NULL, 'g'},
        {"adjust-clusters", no_argument, NULL, 2},
        {"sex", required_argument, NULL, 'x'},
        {"do-not-check-bpm", no_argument, NULL, 3},
        {"genome-studio", required_argument, NULL, 4},
        {"no-version", no_argument, NULL, 8},
        {"threads", required_argument, NULL, 9},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "h?lt:o:O:i:b:c:e:f:g:x:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'l': list_tags(); break;
            case 't': tag_list = optarg; break;
            case 'i': idat_fname = optarg; break;
            case 'b': bpm_fname = optarg; break;
            case 'c': csv_fname = optarg; break;
            case 'e': egt_fname = optarg; break;
            case 'o': output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': output_type = FT_BCF_GZ; break;
                          case 'u': output_type = FT_BCF; break;
                          case 'z': output_type = FT_VCF_GZ; break;
                          case 'v': output_type = FT_VCF; break;
                          case 'g': output_type = FT_GS; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'f': ref_fname = optarg; break;
            case  1 : cache_size = strtol(optarg, NULL, 0); break;
            case 'g': gtc_list = optarg; break;
            case  2 : flags |= ADJUST_CLUSTERS; break;
            case 'x': sex_fname = optarg; break;
            case  3 : bpm_check = 0; break;
            case  4 : gs_fname = optarg; break;
            case  9 : n_threads = strtol(optarg, NULL, 0); break;
            case  8 : record_cmd_line = 0; break;
            case 'h':
            case '?':
            default: error("%s", usage_text());
        }
    }
    if ( (idat_fname!=NULL) + ( (bpm_fname!=NULL) || (csv_fname!=NULL) ) + (egt_fname!=NULL) + argc-optind == 1 ) binary_to_csv = 1;
    if ( !binary_to_csv )
    {
        if ( idat_fname ) error("IDAT file only allowed when converting to CSV\n%s", usage_text());
        if ( !bpm_fname && !csv_fname && !gs_fname ) error("Manifest file required when converting to VCF\n%s", usage_text());
        if ( !egt_fname && ( flags & ADJUST_CLUSTERS ) ) error("Cluster file required when adjusting cluster centers\n%s", usage_text());
        if ( gs_fname && (argc-optind>0 || gtc_list || output_type & FT_GS) ) error("If Genome Studio file provided, do not pass GTC files and do not output to GenomeStudio format\n%s", usage_text());
        if ( argc-optind>0 && gtc_list ) error("GTC files cannot be listed through both command interface and file list\n%s", usage_text());
        if ( !gs_fname && !(output_type & FT_GS) && sex_fname ) out_sex = get_file_handle( sex_fname );
    }
    flags |= parse_tags(tag_list);

    int nfiles;
    char **files;
    if ( gtc_list )
    {
        files = hts_readlines(gtc_list, &nfiles);
        if ( !files ) error("Failed to read from file %s\n", gtc_list);
    }
    else
    {
        nfiles = argc - optind;
        files = argv + optind;
    }
    // make sure the process is allowed to open enough files
    struct rlimit lim;
    getrlimit(RLIMIT_NOFILE, &lim);
    if (nfiles + 7 > lim.rlim_max) error("On this system you cannot open more than %ld files at once while %d is required\n", lim.rlim_max, nfiles + 7);
    if (nfiles + 7 > lim.rlim_cur)
    {
        lim.rlim_cur = nfiles + 7;
        setrlimit(RLIMIT_NOFILE, &lim);
    }
    if ( ( flags & ADJUST_CLUSTERS ) && nfiles < 100 )
        fprintf(stderr, "Warning: adjusting clusters with %d sample(s) is not recommended\n", nfiles);

    if ( binary_to_csv || output_type & FT_GS ) out_txt = get_file_handle( output_fname );
    else
    {
        out_fh = hts_open(output_fname, hts_bcf_wmode(output_type));
        if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
        if ( n_threads ) hts_set_threads(out_fh, n_threads);
        if ( !ref_fname ) error("VCF output requires the --fasta-ref option\n");
        fai = fai_load(ref_fname);
        if ( !fai ) error("Could not load the reference %s\n", ref_fname);
        if ( cache_size ) fai_set_cache_size( fai, cache_size );
    }

    if ( idat_fname )
    {
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "Reading IDAT file %s\n", idat_fname);
        idat_t *idat = idat_init(idat_fname);
        fprintf(stderr, "================================================================================\n");
        idat_summary(idat, stderr);
        if ( binary_to_csv ) idat_to_csv(idat, out_txt);
        idat_destroy(idat);
    }

    bpm_t *bpm = NULL;
    if ( bpm_fname )
    {
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "Reading BPM file %s\n", bpm_fname);
        bpm = bpm_init(bpm_fname);
        bpm_summary(bpm, stderr);
        if ( binary_to_csv && !csv_fname ) bpm_to_csv(bpm, out_txt);
    }

    if ( csv_fname )
    {
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "Reading CSV file %s\n", csv_fname);
        bpm = bpm_csv_init(csv_fname, bpm);
        bpm_summary(bpm, stderr);
        if ( binary_to_csv ) bpm_to_csv(bpm, out_txt);
    }

    if ( bpm->norm_lookups ) flags |= BPM_LOOKUPS;
    if ( ( flags & ADJUST_CLUSTERS ) && !( flags & BPM_LOOKUPS ) )
        error("Cannot adjust clusters as couldn't generate the normalization lookup table\n");

    egt_t *egt = NULL;
    if ( egt_fname )
    {
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "Reading EGT file %s\n", egt_fname);
        egt = egt_init(egt_fname);
        egt_summary(egt, stderr);
        if ( binary_to_csv ) egt_to_csv(egt, out_txt);
        else flags |= EGT_LOADED;
    }

    gtc_t **gtc = (gtc_t **)malloc(nfiles * sizeof(gtc_t *));
    for (int i=0; i<nfiles; i++)
    {
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "Reading GTC file %s\n", files[i]);
        gtc[i] = gtc_init(files[i]);
        // GenCall fills the GTC SNP manifest with the BPM file name rather than the BPM manifest name
        if ( bpm_check && bpm && strcmp( bpm->manifest_name, gtc[i]->snp_manifest ) &&
             strcmp( strrchr(bpm->fn, '/') ? strrchr(bpm->fn, '/') + 1 : bpm->fn, gtc[i]->snp_manifest ) )
            error("Manifest name %s in BPM file %s does not match manifest name %s in GTC file %s\n", bpm->manifest_name, bpm->fn, gtc[i]->snp_manifest, gtc[i]->fn);
        gtc_summary(gtc[i], stderr);
        if ( binary_to_csv ) gtc_to_csv(gtc[i], out_txt);
    }

    if ( !binary_to_csv )
    {
        if ( output_type & FT_GS ) gtcs_to_gs(gtc, nfiles, bpm, egt, out_txt);
        else
        {
            // if it couldn't generate the normalization lookup table from the BPM / CSV manifests
            bcf_hdr_t *hdr = get_hdr(fai, flags);
            if ( bpm_fname ) bcf_hdr_printf(hdr, "##BPM=%s", strrchr(bpm_fname, '/') ? strrchr(bpm_fname, '/') + 1 : bpm_fname);
            if ( csv_fname ) bcf_hdr_printf(hdr, "##CSV=%s", strrchr(csv_fname, '/') ? strrchr(csv_fname, '/') + 1 : csv_fname);
            if ( egt_fname ) bcf_hdr_printf(hdr, "##EGT=%s", strrchr(egt_fname, '/') ? strrchr(egt_fname, '/') + 1 : egt_fname);
            if ( record_cmd_line ) bcf_hdr_append_version(hdr, argc, argv, "bcftools_+gtc2vcf");
            if ( gs_fname )
            {
                htsFile *gs_fh = hts_open(gs_fname, "r");
                gs_to_vcf(gs_fh, out_fh, hdr, flags);
            }
            else
            {
                for (int i=0; i<nfiles; i++)
                {
                    char *ptr = strrchr(gtc[i]->fn, '/');
                    char *str = ptr ? ptr+1 : gtc[i]->fn;
                    ptr = strstr(str, ".gtc");
                    if (ptr) *ptr = '\0';
                    // TODO add here a check to make sure the sample names are distinct
                    bcf_hdr_add_sample(hdr, str);
                    if ( out_sex ) fprintf(out_sex, "%s\t%c\n", str, gtc[i]->gender);
                    if (ptr) *ptr = '.';
                }
                gtcs_to_vcf(fai, bpm, egt, gtc, nfiles, out_fh, hdr, flags);
            }
        }
     }

    if ( fai ) fai_destroy(fai);
    egt_destroy(egt);
    bpm_destroy(bpm);
    if ( gtc_list )
    {
        for (int i=0; i<nfiles; i++) free(files[i]);
        free(files);
    }
    for (int i=0; i<nfiles; i++) gtc_destroy(gtc[i]);
    free(gtc);
    if (out_sex && out_sex != stdout && out_sex != stderr) fclose(out_sex);
    return 0;
}
