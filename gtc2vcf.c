/* The MIT License

   Copyright (c) 2018-2021 Giulio Genovese

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
#include <sys/resource.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "tsv2vcf.h"
#include "gtc2vcf.h"

#define GTC2VCF_VERSION "2021-01-20"

#define GT_NC 0
#define GT_AA 1
#define GT_AB 2
#define GT_BB 3

#define TAG_LIST_DFLT "GT,GQ,IGC,BAF,LRR,NORMX,NORMY,R,THETA,X,Y"
#define GC_WIN_DFLT "200"
#define CAPACITY_DFLT "32768"
#define GENOME_BUILD_DFLT "GRCh38"

#define VERBOSE (1 << 0)
#define BPM_LOADED (1 << 1)
#define CSV_LOADED (1 << 2)
#define EGT_LOADED (1 << 3)
#define LOAD_IDAT (1 << 4)
#define ADJUST_CLUSTERS (1 << 5)
#define GENOME_STUDIO (1 << 6)
#define NO_INFO_GC (1 << 7)
#define FORMAT_GT (1 << 8)
#define FORMAT_GQ (1 << 9)
#define FORMAT_IGC (1 << 10)
#define FORMAT_BAF (1 << 11)
#define FORMAT_LRR (1 << 12)
#define FORMAT_NORMX (1 << 13)
#define FORMAT_NORMY (1 << 14)
#define FORMAT_R (1 << 15)
#define FORMAT_THETA (1 << 16)
#define FORMAT_X (1 << 17)
#define FORMAT_Y (1 << 18)

/****************************************
 * hFILE READING FUNCTIONS              *
 ****************************************/

// read or skip a fixed length array
static inline void read_array(hFILE *hfile, void **arr, size_t *m_arr, size_t nmemb, size_t size, size_t term) {
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
        if (hread(hfile, *arr, nmemb * size) < nmemb * size) {
            error("Failed to read %ld bytes from stream\n", nmemb * size);
        }
    } else {
        int c = 0;
        for (int i = 0; i < nmemb * size; i++) c = hgetc(hfile);
        if (c == EOF) error("Failed to reposition stream forward %ld bytes\n", nmemb * size);
    }
}

// read or skip a length-prefixed array
static inline void read_pfx_array(hFILE *hfile, void **arr, size_t *m_arr, size_t item_size) {
    int32_t n;
    if (hread(hfile, (void *)&n, 4) < 4) {
        error("Failed to read 4 bytes from stream\n");
    }
    read_array(hfile, arr, m_arr, n, item_size, 0);
}

// read or skip a length-prefixed string
// http://en.wikipedia.org/wiki/LEB128#Decode_unsigned_integer
static inline void read_pfx_string(hFILE *hfile, char **str, size_t *m_str) {
    uint8_t byte;
    size_t n = 0, shift = 0;
    while (1) {
        if (hread(hfile, (void *)&byte, 1) < 1) {
            error("Failed to read 1 byte from stream\n");
        }
        n |= (size_t)(byte & 0x7F) << shift;
        if (!(byte & 0x80)) break;
        shift += 7;
    }
    if (n || m_str) {
        read_array(hfile, (void **)str, m_str, n, 1, 1);
        if (str) (*str)[n] = '\0';
    }
}

// check whether file is compressed with gzip
static inline int is_gzip(hFILE *hfile) {
    uint8_t buffer[2];
    if (hpeek(hfile, (void *)buffer, 2) < 2) error("Failed to read 2 bytes from stream\n");
    return (buffer[0] == 0x1f && buffer[1] == 0x8b);
}

/****************************************
 * BUFFER ARRAY IMPLEMENTATION          *
 ****************************************/

typedef struct {
    hFILE *hfile;
    off_t offset;
    int32_t item_num;
    int32_t item_offset;
    size_t item_capacity;
    size_t item_size;
    char *buffer;
} buffer_array_t;

static buffer_array_t *buffer_array_init(hFILE *hfile, size_t capacity, size_t item_size) {
    buffer_array_t *arr = (buffer_array_t *)malloc(1 * sizeof(buffer_array_t));
    arr->hfile = hfile;
    read_bytes(hfile, (void *)&arr->item_num, sizeof(int32_t));
    arr->offset = htell(arr->hfile);
    arr->item_offset = 0;
    arr->item_capacity = (capacity <= 0) ? (size_t)strtol(CAPACITY_DFLT, NULL, 0) : capacity;
    arr->item_size = item_size;
    arr->buffer = (char *)malloc(arr->item_capacity * item_size);
    read_bytes(hfile, (void *)arr->buffer,
               (arr->item_num < arr->item_capacity ? arr->item_num : arr->item_capacity) * item_size);
    return arr;
}

static inline int get_element(buffer_array_t *arr, void *dst, size_t item_idx) {
    if (!arr || item_idx >= arr->item_num) {
        return -1;
    } else if (item_idx - arr->item_offset < arr->item_capacity) {
        memcpy(dst, (void *)(arr->buffer + (item_idx - arr->item_offset) * arr->item_size), arr->item_size);
        return 0;
    }
    arr->item_offset = item_idx;
    if (hseek(arr->hfile, arr->offset + item_idx * arr->item_size, SEEK_SET) < 0) {
        error("Fail to seek to position %ld in file\n", arr->offset + item_idx * arr->item_size);
    }
    read_bytes(arr->hfile, (void *)arr->buffer,
               ((arr->item_num - arr->item_offset) < arr->item_capacity ? (arr->item_num - arr->item_offset)
                                                                        : arr->item_capacity)
                   * arr->item_size);
    memcpy(dst, (void *)arr->buffer, arr->item_size);
    return 0;
}

static void buffer_array_destroy(buffer_array_t *arr) {
    if (!arr) return;
    free(arr->buffer);
    free(arr);
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
    float frac_a;
    float frac_c;
    float frac_g;
    float frac_t;
    char *ref_strand; // RefStrand annotation
} LocusEntry;

// retrieve assay type according to the following (allele_a_probe_seq, source_seq) -> assay_type
// map
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
    char trail_a_probe_seq = toupper(allele_a_probe_seq[strlen(allele_a_probe_seq) - 2]);
    if (trail_a_probe_seq == 'C' || trail_a_probe_seq == 'G') return 1;
    if (trail_a_probe_seq == 'A' || trail_a_probe_seq == 'T') return 2;
    error("Unable to retrieve assay type: %s %s\n", allele_a_probe_seq, source_seq);
}

static void locusentry_read(LocusEntry *locus_entry, hFILE *hfile) {
    locus_entry->norm_id = 0xFF;
    read_bytes(hfile, (void *)&locus_entry->version, sizeof(int32_t));
    if (locus_entry->version < 4 || locus_entry->version == 5 || locus_entry->version > 8)
        error("Locus version %d in manifest file not supported\n", locus_entry->version);
    read_pfx_string(hfile, &locus_entry->ilmn_id, NULL);
    read_pfx_string(hfile, &locus_entry->name, NULL);
    read_pfx_string(hfile, NULL, NULL);
    read_pfx_string(hfile, NULL, NULL);
    read_pfx_string(hfile, NULL, NULL);
    read_bytes(hfile, (void *)&locus_entry->index, sizeof(int32_t));
    read_pfx_string(hfile, NULL, NULL);
    read_pfx_string(hfile, &locus_entry->ilmn_strand, NULL);
    read_pfx_string(hfile, &locus_entry->snp, NULL);
    read_pfx_string(hfile, &locus_entry->chrom, NULL);
    read_pfx_string(hfile, &locus_entry->ploidy, NULL);
    read_pfx_string(hfile, &locus_entry->species, NULL);
    read_pfx_string(hfile, &locus_entry->map_info, NULL);
    read_pfx_string(hfile, &locus_entry->top_genomic_seq, NULL); // only version 4
    read_pfx_string(hfile, &locus_entry->customer_strand, NULL);
    read_bytes(hfile, (void *)&locus_entry->address_a, sizeof(int32_t));
    read_bytes(hfile, (void *)&locus_entry->address_b, sizeof(int32_t));
    read_pfx_string(hfile, &locus_entry->allele_a_probe_seq, NULL); // only version 4
    read_pfx_string(hfile, &locus_entry->allele_b_probe_seq, NULL); // only version 4
    read_pfx_string(hfile, &locus_entry->genome_build, NULL);
    read_pfx_string(hfile, &locus_entry->source, NULL);
    read_pfx_string(hfile, &locus_entry->source_version, NULL);
    read_pfx_string(hfile, &locus_entry->source_strand, NULL);
    read_pfx_string(hfile, &locus_entry->source_seq, NULL); // only version 4
    if (locus_entry->source_seq) {
        char *ptr = strchr(locus_entry->source_seq, '-');
        if (ptr && *(ptr - 1) == '/') {
            *ptr = *(ptr - 2);
            *(ptr - 2) = '-';
        }
    }

    if (locus_entry->version >= 6) {
        read_bytes(hfile, NULL, 1);
        read_bytes(hfile, (void *)&locus_entry->exp_clusters, sizeof(int8_t));
        read_bytes(hfile, (void *)&locus_entry->intensity_only, sizeof(int8_t));
        read_bytes(hfile, (void *)&locus_entry->assay_type, sizeof(uint8_t));

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
        read_bytes(hfile, &locus_entry->frac_a, sizeof(float));
        read_bytes(hfile, &locus_entry->frac_c, sizeof(float));
        read_bytes(hfile, &locus_entry->frac_t, sizeof(float));
        read_bytes(hfile, &locus_entry->frac_g, sizeof(float));
    }
    if (locus_entry->version >= 8) read_pfx_string(hfile, &locus_entry->ref_strand, NULL);
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
    uint8_t *norm_ids;
    LocusEntry *locus_entries;
    uint8_t *norm_lookups;
    char **header;
    size_t m_header;
} bpm_t;

static uint8_t *bpm_norm_lookups(bpm_t *bpm) {
    uint8_t sorted_norm_ids[256];
    for (int i = 0; i < 256; i++) sorted_norm_ids[i] = 0xFF;
    for (int i = 0; i < bpm->num_loci; i++) {
        int norm_id = bpm->locus_entries[i].norm_id;
        sorted_norm_ids[norm_id] = norm_id;
    }
    int j = 0;
    for (int i = 0; i < 256; i++)
        if (sorted_norm_ids[i] != 0xFF) sorted_norm_ids[j++] = sorted_norm_ids[i];
    uint8_t *norm_lookups = (uint8_t *)malloc(256 * sizeof(uint8_t *));
    memset((void *)norm_lookups, 0xFF, 256 * sizeof(uint8_t *));
    for (int i = 0; i < j; i++) norm_lookups[sorted_norm_ids[i]] = i;
    return norm_lookups;
}

static bpm_t *bpm_init(const char *fn) {
    bpm_t *bpm = (bpm_t *)calloc(1, sizeof(bpm_t));
    bpm->fn = strdup(fn);
    bpm->hfile = hopen(bpm->fn, "rb");
    if (bpm->hfile == NULL) error("Could not open %s: %s\n", bpm->fn, strerror(errno));
    if (is_gzip(bpm->hfile)) error("File %s is gzip compressed and currently cannot be sought\n", bpm->fn);

    uint8_t buffer[4];
    if (hread(bpm->hfile, (void *)buffer, 4) < 4) error("Failed to read magic number from %s file\n", bpm->fn);
    if (memcmp(buffer, "BPM", 3) != 0) error("BPM file %s format identifier is bad\n", bpm->fn);
    if (buffer[3] != 1) error("BPM file %s version is unknown\n", bpm->fn);

    read_bytes(bpm->hfile, (void *)&bpm->version, sizeof(int32_t));
    if (bpm->version & 0x1000) bpm->version ^= 0x1000;
    if (bpm->version > 5 || bpm->version < 3) error("BPM file %s version %d is unsupported\n", bpm->fn, bpm->version);
    read_pfx_string(bpm->hfile, &bpm->manifest_name, NULL);

    if (bpm->version > 1) read_pfx_string(bpm->hfile, &bpm->control_config, NULL);

    read_bytes(bpm->hfile, (void *)&bpm->num_loci, sizeof(int32_t));
    read_array(bpm->hfile, (void **)&bpm->indexes, NULL, bpm->num_loci, sizeof(int32_t), 0);
    bpm->names = (char **)malloc(bpm->num_loci * sizeof(char *));
    for (int i = 0; i < bpm->num_loci; i++) read_pfx_string(bpm->hfile, &bpm->names[i], NULL);
    read_array(bpm->hfile, (void **)&bpm->norm_ids, NULL, bpm->num_loci, sizeof(uint8_t), 0);

    bpm->locus_entries = (LocusEntry *)malloc(bpm->num_loci * sizeof(LocusEntry));
    LocusEntry locus_entry;
    for (int i = 0; i < bpm->num_loci; i++) {
        memset(&locus_entry, 0, sizeof(LocusEntry));
        locusentry_read(&locus_entry, bpm->hfile);
        int idx = locus_entry.index - 1;
        if (idx < 0 || idx >= bpm->num_loci) error("Locus entry index %d is out of boundaries\n", locus_entry.index);
        if (bpm->norm_ids[idx] > 100)
            error("Manifest format error: read invalid normalization ID %d\n", bpm->norm_ids[idx]);
        // To mimic the flawed byte-wrapping behavior from GenomeStudio, AutoCall, and
        // IAAP, this value is allowed to overflow beyond 255, which happens with some
        // probes in the Omni5 arrays
        locus_entry.norm_id = bpm->norm_ids[idx] + 100 * locus_entry.assay_type;
        memcpy(&bpm->locus_entries[idx], &locus_entry, sizeof(LocusEntry));
    }
    bpm->norm_lookups = bpm_norm_lookups(bpm);
    for (int i = 0; i < bpm->num_loci; i++) {
        if (i != bpm->locus_entries[i].index - 1)
            error("Manifest format error: read invalid number of assay entries\n");
    }
    if (bpm->locus_entries[0].version < 8)
        fprintf(stderr, "Warning: RefStrand annotation missing from manifest file %s\n", bpm->fn);

    read_bytes(bpm->hfile, (void *)&bpm->m_header, sizeof(int32_t));
    bpm->header = (char **)malloc(bpm->m_header * sizeof(char *));
    for (int i = 0; i < bpm->m_header; i++) read_pfx_string(bpm->hfile, &bpm->header[i], NULL);

    if (!heof(bpm->hfile))
        error("BPM reader did not reach the end of file %s at position %ld\n", bpm->fn, htell(bpm->hfile));

    return bpm;
}

static void bpm_destroy(bpm_t *bpm) {
    if (!bpm) return;
    if (bpm->hfile && hclose(bpm->hfile) < 0) error("Error closing BPM file %s\n", bpm->fn);
    free(bpm->fn);
    if (bpm->fp && hts_close(bpm->fp) < 0) error("Error closing CSV file %s\n", bpm->fp->fn);
    free(bpm->manifest_name);
    free(bpm->control_config);
    free(bpm->indexes);
    if (bpm->names) {
        for (int i = 0; i < bpm->num_loci; i++) free(bpm->names[i]);
        free(bpm->names);
    }
    free(bpm->norm_ids);
    for (int i = 0; i < bpm->num_loci; i++) {
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
    for (int i = 0; i < bpm->m_header; i++) free(bpm->header[i]);
    free(bpm->header);
    free(bpm);
}

static void bpm_to_csv(const bpm_t *bpm, FILE *stream, int flags) {
    for (int i = 0; i < bpm->m_header; i++) fprintf(stream, "%s\n", bpm->header[i]);
    if (flags & BPM_LOADED) {
        fprintf(stream,
                "Index,NormID,IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_"
                "ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,"
                "SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID,Exp_Clusters,"
                "Intensity_Only,Assay_Type,Frac A,Frac C,Frac G,Frac T,RefStrand\n");
    } else {
        fprintf(stream,
                "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_"
                "ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,"
                "SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand\n");
    }
    if (flags & VERBOSE) {
        kstring_t address_b = {0, 0, NULL};
        if (flags & BPM_LOADED) {
            for (int i = 0; i < bpm->num_loci; i++) {
                LocusEntry *locus_entry = &bpm->locus_entries[i];
                address_b.l = 0;
                ksprintf(&address_b, locus_entry->address_b ? "%010d" : "", locus_entry->address_b);
                fprintf(stream,
                        "%d,%d,%s,%s,%s,%s,%010d,%-s,%s,%-s,%s,%s,%s,%s,%s,%s,%s,%s,%-s,%-s,%d,"
                        "%d,%d,%d,%f,%f,%f,%f,%s\n",
                        locus_entry->index, locus_entry->norm_id, locus_entry->ilmn_id, locus_entry->name,
                        locus_entry->ilmn_strand, locus_entry->snp, locus_entry->address_a,
                        locus_entry->allele_a_probe_seq ? locus_entry->allele_a_probe_seq : "", address_b.s,
                        locus_entry->allele_b_probe_seq ? locus_entry->allele_b_probe_seq : "",
                        locus_entry->genome_build, locus_entry->chrom, locus_entry->map_info, locus_entry->ploidy,
                        locus_entry->species, locus_entry->source, locus_entry->source_version,
                        locus_entry->source_strand, locus_entry->source_seq ? locus_entry->source_seq : "",
                        locus_entry->top_genomic_seq ? locus_entry->top_genomic_seq : "", locus_entry->beadset_id,
                        locus_entry->exp_clusters, locus_entry->intensity_only, locus_entry->assay_type,
                        locus_entry->frac_a, locus_entry->frac_c, locus_entry->frac_g, locus_entry->frac_t,
                        locus_entry->ref_strand ? locus_entry->ref_strand : "");
            }
        } else {
            for (int i = 0; i < bpm->num_loci; i++) {
                LocusEntry *locus_entry = &bpm->locus_entries[i];
                address_b.l = 0;
                ksprintf(&address_b, locus_entry->address_b ? "%010d" : "", locus_entry->address_b);
                fprintf(stream, "%s,%s,%s,%s,%010d,%-s,%s,%-s,%s,%s,%s,%s,%s,%s,%s,%s,%-s,%-s,%d,%d,%s\n",
                        locus_entry->ilmn_id, locus_entry->name, locus_entry->ilmn_strand, locus_entry->snp,
                        locus_entry->address_a, locus_entry->allele_a_probe_seq, address_b.s,
                        locus_entry->allele_b_probe_seq ? locus_entry->allele_b_probe_seq : "",
                        locus_entry->genome_build, locus_entry->chrom, locus_entry->map_info, locus_entry->ploidy,
                        locus_entry->species, locus_entry->source, locus_entry->source_version,
                        locus_entry->source_strand, locus_entry->source_seq, locus_entry->top_genomic_seq,
                        locus_entry->beadset_id, locus_entry->exp_clusters,
                        locus_entry->ref_strand ? locus_entry->ref_strand : "");
            }
        }
        free(address_b.s);
    } else {
        fprintf(stream, "... use --verbose to visualize Assay data ...\n");
    }
    fprintf(stream, "[Controls]\n");
    fprintf(stream, "%s", bpm->control_config);
}

/****************************************
 * CSV FILE IMPLEMENTATION              *
 ****************************************/

static int tsv_read_uint8(tsv_t *tsv, bcf1_t *rec, void *usr) {
    uint8_t *uint8 = (uint8_t *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *uint8 = (uint8_t)strtol(tsv->ss, &endptr, 0);
    *tsv->se = tmp;
    return 0;
}

static int tsv_read_int32(tsv_t *tsv, bcf1_t *rec, void *usr) {
    int32_t *int32 = (int32_t *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *int32 = (int32_t)strtol(tsv->ss, &endptr, 10);
    *tsv->se = tmp;
    return 0;
}

static int tsv_read_float(tsv_t *tsv, bcf1_t *rec, void *usr) {
    float *single = (float *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *single = (float)strtof(tsv->ss, &endptr);
    *tsv->se = tmp;
    return 0;
}

static int tsv_read_string(tsv_t *tsv, bcf1_t *rec, void *usr) {
    char **str = (char **)usr;
    if (tsv->se == tsv->ss) {
        *str = NULL;
    } else {
        char tmp = *tsv->se;
        *tsv->se = 0;
        *str = strdup(tsv->ss);
        *tsv->se = tmp;
    }
    return 0;
}

// Petr Danecek's similar implementation in bcftools/tsv2vcf.c
static int csv_parse(tsv_t *tsv, bcf1_t *rec, char *str) {
    int status = 0;
    tsv->icol = 0;
    tsv->ss = tsv->se = str;
    while (*tsv->ss && tsv->icol < tsv->ncols) {
        while (*tsv->se && *tsv->se != ',') tsv->se++;
        if (tsv->cols[tsv->icol].setter) {
            int ret = tsv->cols[tsv->icol].setter(tsv, rec, tsv->cols[tsv->icol].usr);
            if (ret < 0) return -1;
            status++;
        }
        tsv->se++;
        tsv->ss = tsv->se;
        tsv->icol++;
    }
    return status ? 0 : -1;
}

static void locus_merge(LocusEntry *dest, LocusEntry *src) {
    if (src->version) dest->version = src->version;
    if (src->norm_id != 0xFF) dest->norm_id = src->norm_id;
    if (strcmp(dest->ilmn_id, src->ilmn_id)) {
        error("BPM and CSV manifests have conflicting IDs: %s and %s\n", dest->ilmn_id, src->ilmn_id);
    } else {
        free(dest->ilmn_id);
        dest->ilmn_id = src->ilmn_id;
    }
    if (src->name) {
        free(dest->name);
        dest->name = src->name;
    }
    if (src->index != 0) dest->index = src->index;
    if (src->ilmn_strand) {
        free(dest->ilmn_strand);
        dest->ilmn_strand = src->ilmn_strand;
    }
    if (src->snp) {
        free(dest->snp);
        dest->snp = src->snp;
    }
    if (src->chrom) {
        free(dest->chrom);
        dest->chrom = src->chrom;
    }
    if (src->ploidy) {
        free(dest->ploidy);
        dest->ploidy = src->ploidy;
    }
    if (src->species) {
        free(dest->species);
        dest->species = src->species;
    }
    if (src->map_info) {
        free(dest->map_info);
        dest->map_info = src->map_info;
    }
    if (src->customer_strand) {
        free(dest->customer_strand);
        dest->customer_strand = src->customer_strand;
    }
    if (src->address_a != 0) dest->address_a = src->address_a;
    if (src->allele_a_probe_seq) {
        free(dest->allele_a_probe_seq);
        dest->allele_a_probe_seq = src->allele_a_probe_seq;
    }
    if (src->address_b != 0) dest->address_b = src->address_b;
    if (src->allele_b_probe_seq) {
        free(dest->allele_b_probe_seq);
        dest->allele_b_probe_seq = src->allele_b_probe_seq;
    }
    if (src->genome_build) {
        free(dest->genome_build);
        dest->genome_build = src->genome_build;
    }
    if (src->source) {
        free(dest->source);
        dest->source = src->source;
    }
    if (src->source_version) {
        free(dest->source_version);
        dest->source_version = src->source_version;
    }
    if (src->source_strand) {
        free(dest->source_strand);
        dest->source_strand = src->source_strand;
    }
    if (src->source_seq) {
        free(dest->source_seq);
        dest->source_seq = src->source_seq;
    }
    if (src->top_genomic_seq) {
        free(dest->top_genomic_seq);
        dest->top_genomic_seq = src->top_genomic_seq;
    }
    if (src->beadset_id) dest->beadset_id = src->beadset_id;
    if (src->exp_clusters) dest->exp_clusters = src->exp_clusters;
    if (src->intensity_only) dest->intensity_only = src->intensity_only;
    if (src->assay_type != 0xFF) dest->assay_type = src->assay_type;
    if (src->frac_a) dest->frac_a = src->frac_a;
    if (src->frac_c) dest->frac_c = src->frac_c;
    if (src->frac_g) dest->frac_g = src->frac_g;
    if (src->frac_t) dest->frac_t = src->frac_t;
    if (src->ref_strand) {
        free(dest->ref_strand);
        dest->ref_strand = src->ref_strand;
    }
}

// this line will read a CSV file and if a BPM object is provided it will fill it rather than
// create a new one
static bpm_t *bpm_csv_init(const char *fn, bpm_t *bpm) {
    int bpm_available = bpm != NULL;
    if (!bpm_available) bpm = (bpm_t *)calloc(1, sizeof(bpm_t));
    int bpm_prev_num_loci = bpm->num_loci;

    bpm->fp = hts_open(fn, "r");
    if (bpm->fp == NULL) error("Could not open %s: %s\n", fn, strerror(errno));

    kstring_t str = {0, 0, NULL};
    kstring_t hdr = {0, 0, NULL};
    if (hts_getline(bpm->fp, KS_SEP_LINE, &str) <= 0) error("Empty file: %s\n", fn);
    if (strncmp(str.s, "Illumina", 8) && strncmp(str.s, "\"Illumina", 9))
        error("Header of file %s is incorrect: %s\n", fn, str.s);
    kputs(str.s, &hdr);
    kputc('\n', &hdr);

    char *tmp = NULL;
    size_t prev = 0;
    while (strncmp(str.s + prev, "[Assay]", 7)) {
        if (strncmp(str.s + prev, "Descriptor File Name,", 21) == 0) {
            free(bpm->manifest_name);
            bpm->manifest_name = strdup(str.s + prev + 21);
            char *ptr = strchr(bpm->manifest_name, ',');
            if (ptr) *ptr = '\0';
        } else if (strncmp(str.s + prev, "Loci Count ,", 12) == 0) {
            bpm->num_loci = (int)strtol(str.s + prev + 12, &tmp, 0);
        } else if (strncmp(str.s + prev, "Loci Count,", 11) == 0) {
            bpm->num_loci = (int)strtol(str.s + prev + 11, &tmp, 0);
        }
        if (hts_getline(bpm->fp, KS_SEP_LINE, &str) <= 0) error("Error reading from file: %s\n", fn);
        kputs(str.s, &hdr);
        kputc('\n', &hdr);
    }
    if (bpm->num_loci == 0)
        error("Could not understand number of loci from header of manifest file %s\n", fn);
    else if (bpm_available && bpm_prev_num_loci != bpm->num_loci)
        error("BPM manifest file has %d loci while CSV manifest file %s has %d loci\n", bpm_prev_num_loci, fn,
              bpm->num_loci);

    int moff = 0, *off = NULL;
    for (int i = 0; i < bpm->m_header; i++) free(bpm->header[i]);
    bpm->m_header = ksplit_core(hdr.s, '\n', &moff, &off);
    free(bpm->header);
    bpm->header = (char **)malloc(bpm->m_header * sizeof(char *));
    for (int i = 0; i < bpm->m_header; i++) bpm->header[i] = strdup(&hdr.s[off[i]]);
    free(off);
    free(hdr.s);

    if (hts_getline(bpm->fp, KS_SEP_LINE, &str) <= 0) error("Error reading from file: %s\n", fn);

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
    int ref_strand = tsv_register(tsv, "RefStrand", tsv_read_string, &locus_entry.ref_strand);
    if (ref_strand < 0) fprintf(stderr, "Warning: RefStrand annotation missing from manifest file %s\n", fn);

    if (!bpm_available) bpm->locus_entries = (LocusEntry *)malloc(bpm->num_loci * sizeof(LocusEntry));
    for (int i = 0; i < bpm->num_loci; i++) {
        memset(&locus_entry, 0, sizeof(LocusEntry));
        locus_entry.norm_id = 0xFF;
        if (hts_getline(bpm->fp, KS_SEP_LINE, &str) <= 0) error("Error reading from file: %s\n", fn);
        if (csv_parse(tsv, NULL, str.s) < 0) error("Could not parse the manifest file: %s\n", str.s);
        if (locus_entry.source_seq) {
            char *ptr = strchr(locus_entry.source_seq, '-');
            if (ptr && *(ptr - 1) == '/') {
                *ptr = *(ptr - 2);
                *(ptr - 2) = '-';
            }
        }
        locus_entry.assay_type =
            get_assay_type(locus_entry.allele_a_probe_seq, locus_entry.allele_b_probe_seq, locus_entry.source_seq);
        if (locus_entry.index == 0) locus_entry.index = i + 1;
        int idx = locus_entry.index - 1;
        if (idx < 0 || idx >= bpm->num_loci) error("Locus entry index %d is out of boundaries\n", idx);
        if (!bpm_available)
            memcpy(&bpm->locus_entries[idx], &locus_entry, sizeof(LocusEntry));
        else
            locus_merge(&bpm->locus_entries[idx], &locus_entry);
    }
    tsv_destroy(tsv);

    if (hts_getline(bpm->fp, KS_SEP_LINE, &str) <= 0) error("Error reading from file: %s\n", fn);
    if (strncmp(str.s, "[Controls]", 10) != 0)
        error(
            "Missing [Controls] section from manifest file: %s\n"
            "Found the following line instead: %s\n",
            fn, str.s);
    while (hts_getline(bpm->fp, KS_SEP_LINE, &str) > 0) kputc('\n', &str);
    free(bpm->control_config);
    bpm->control_config = str.s;

    if (norm_id == 0) {
        free(bpm->norm_lookups);
        bpm->norm_lookups = bpm_norm_lookups(bpm);
    }

    return bpm;
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
    char **loci_names;
} egt_t;

static void clusterscore_read(ClusterScore *clusterscore, hFILE *hfile) {
    read_bytes(hfile, (void *)&clusterscore->cluster_separation, sizeof(float));
    read_bytes(hfile, (void *)&clusterscore->total_score, sizeof(float));
    read_bytes(hfile, (void *)&clusterscore->original_score, sizeof(float));
    read_bytes(hfile, (void *)&clusterscore->edited, sizeof(uint8_t));
}

static void clusterrecord_read(ClusterRecord *clusterrecord, hFILE *hfile, int32_t data_block_version) {
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.N, sizeof(int32_t));
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.N, sizeof(int32_t));
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.N, sizeof(int32_t));
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.r_dev, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.r_dev, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.r_dev, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.r_mean, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.r_mean, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.r_mean, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.theta_dev, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.theta_dev, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.theta_dev, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->aa_cluster_stats.theta_mean, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->ab_cluster_stats.theta_mean, sizeof(float));
    read_bytes(hfile, (void *)&clusterrecord->bb_cluster_stats.theta_mean, sizeof(float));
    if (data_block_version >= 7) {
        read_bytes(hfile, (void *)&clusterrecord->intensity_threshold, sizeof(float));
        read_bytes(hfile, NULL, 14 * sizeof(float));
    } else {
        clusterrecord->intensity_threshold = NAN;
    }
}

static egt_t *egt_init(const char *fn) {
    egt_t *egt = (egt_t *)calloc(1, sizeof(egt_t));
    egt->fn = strdup(fn);
    egt->hfile = hopen(egt->fn, "rb");
    if (egt->hfile == NULL) error("Could not open %s: %s\n", egt->fn, strerror(errno));
    if (is_gzip(egt->hfile)) error("File %s is gzip compressed and currently cannot be sought\n", egt->fn);

    read_bytes(egt->hfile, (void *)&egt->version, sizeof(int32_t));
    if (egt->version != 3) error("EGT cluster file version %d not supported\n", egt->version);

    read_pfx_string(egt->hfile, &egt->gencall_version, NULL);
    read_pfx_string(egt->hfile, &egt->cluster_version, NULL);
    read_pfx_string(egt->hfile, &egt->call_version, NULL);
    read_pfx_string(egt->hfile, &egt->normalization_version, NULL);
    read_pfx_string(egt->hfile, &egt->date_created, NULL);

    read_bytes(egt->hfile, (void *)&egt->is_wgt, sizeof(uint8_t));
    if (egt->is_wgt != 1) error("Only WGT cluster file version supported\n");

    read_pfx_string(egt->hfile, &egt->manifest_name, NULL);

    read_bytes(egt->hfile, (void *)&egt->data_block_version, sizeof(int32_t));
    if (egt->data_block_version < 5 || egt->data_block_version == 6 || egt->data_block_version > 9)
        error("Data block version %d in cluster file not supported\n", egt->data_block_version);
    read_pfx_string(egt->hfile, &egt->opa, NULL);

    read_bytes(egt->hfile, (void *)&egt->num_records, sizeof(int32_t));
    egt->cluster_records = (ClusterRecord *)malloc(egt->num_records * sizeof(ClusterRecord));
    for (int i = 0; i < egt->num_records; i++)
        clusterrecord_read(&egt->cluster_records[i], egt->hfile, egt->data_block_version);
    for (int i = 0; i < egt->num_records; i++) clusterscore_read(&egt->cluster_records[i].cluster_score, egt->hfile);

    // toss useless strings such as aa_ab_bb/aa_ab/aa_bb/ab_bb
    for (int i = 0; i < egt->num_records; i++) read_pfx_string(egt->hfile, NULL, NULL);

    egt->loci_names = (char **)malloc(egt->num_records * sizeof(char *));
    for (int i = 0; i < egt->num_records; i++) {
        read_pfx_string(egt->hfile, &egt->loci_names[i], NULL);
    }
    for (int i = 0; i < egt->num_records; i++)
        read_bytes(egt->hfile, (void *)&egt->cluster_records[i].address, sizeof(int32_t));

    int32_t aa_n, ab_n, bb_n;
    for (int i = 0; i < egt->num_records; i++) {
        read_bytes(egt->hfile, (void *)&aa_n, sizeof(int32_t));
        read_bytes(egt->hfile, (void *)&ab_n, sizeof(int32_t));
        read_bytes(egt->hfile, (void *)&bb_n, sizeof(int32_t));
        if (egt->cluster_records[i].aa_cluster_stats.N != aa_n || egt->cluster_records[i].ab_cluster_stats.N != ab_n
            || egt->cluster_records[i].bb_cluster_stats.N != bb_n)
            error("Cluster counts don't match with EGT cluster file %s\n", egt->fn);
    }

    if (egt->data_block_version == 9) read_bytes(egt->hfile, NULL, egt->num_records * sizeof(float));
    if (!heof(egt->hfile))
        error("EGT reader did not reach the end of file %s at position %ld\n", egt->fn, htell(egt->hfile));

    return egt;
}

static void egt_destroy(egt_t *egt) {
    if (!egt) return;
    if (hclose(egt->hfile) < 0) error("Error closing EGT file %s\n", egt->fn);
    free(egt->fn);
    free(egt->gencall_version);
    free(egt->cluster_version);
    free(egt->call_version);
    free(egt->normalization_version);
    free(egt->date_created);
    free(egt->manifest_name);
    free(egt->cluster_records);
    for (int i = 0; i < egt->num_records; i++) free(egt->loci_names[i]);
    free(egt->loci_names);
    free(egt);
}

static void egt_to_csv(const egt_t *egt, FILE *stream, int verbose) {
    fprintf(stream, "Illumina, Inc.\n");
    fprintf(stream, "[Heading]\n");
    fprintf(stream, "Descriptor File Name,%s\n", strrchr(egt->fn, '/') ? strrchr(egt->fn, '/') + 1 : egt->fn);
    fprintf(stream, "GenCall version,%s\n", egt->gencall_version);
    fprintf(stream, "Clustering algorithm version,%s\n", egt->cluster_version);
    fprintf(stream, "Genotyping algorithm version,%s\n", egt->call_version);
    fprintf(stream, "Normalization algorithm version,%s\n", egt->normalization_version);
    fprintf(stream, "Date Manufactured,%s\n", egt->date_created);
    fprintf(stream, "Manifest name used to build this cluster file,%s\n", egt->manifest_name);
    fprintf(stream, "OPA,%s\n", egt->opa ? egt->opa : "");
    fprintf(stream, "Loci Count,%d\n", egt->num_records);
    fprintf(stream, "[Assay]\n");
    fprintf(stream,
            "Name,AA.N,AA.R_dev,AA.R_mean,AA.Theta_dev,AA.Theta_mean,AB.N,AB.R_dev,AB.R_mean,AB."
            "Theta_dev,AB.Theta_mean,BB.N,BB.R_dev,BB.R_mean,BB.Theta_dev,BB.Theta_mean,Intensity "
            "Threshold,Cluster Separation,GenTrain Score,Original Score,Edited,Address\n");
    if (verbose) {
        for (int i = 0; i < egt->num_records; i++) {
            ClusterRecord *cluster_record = &egt->cluster_records[i];
            fprintf(stream, "%s,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\n", egt->loci_names[i],
                    cluster_record->aa_cluster_stats.N, cluster_record->aa_cluster_stats.r_dev,
                    cluster_record->aa_cluster_stats.r_mean, cluster_record->aa_cluster_stats.theta_dev,
                    cluster_record->aa_cluster_stats.theta_mean, cluster_record->ab_cluster_stats.N,
                    cluster_record->ab_cluster_stats.r_dev, cluster_record->ab_cluster_stats.r_mean,
                    cluster_record->ab_cluster_stats.theta_dev, cluster_record->ab_cluster_stats.theta_mean,
                    cluster_record->bb_cluster_stats.N, cluster_record->bb_cluster_stats.r_dev,
                    cluster_record->bb_cluster_stats.r_mean, cluster_record->bb_cluster_stats.theta_dev,
                    cluster_record->bb_cluster_stats.theta_mean, cluster_record->intensity_threshold,
                    cluster_record->cluster_score.cluster_separation, cluster_record->cluster_score.total_score,
                    cluster_record->cluster_score.original_score, cluster_record->cluster_score.edited,
                    cluster_record->address);
        }
    } else {
        fprintf(stream, "... use --verbose to visualize Assay data ...\n");
    }
}

/****************************************
 * IDAT FILE IMPLEMENTATION             *
 ****************************************/

// http://github.com/snewhouse/glu-genetics/blob/master/glu/lib/illumina.py
// http://github.com/HenrikBengtsson/illuminaio/blob/master/R/readIDAT.R

#define NUM_SNPS_READ 1000
#define ILLUMINA_ID 102
#define SD 103
#define MEAN 104
#define NBEADS 107
#define MID_BLOCK 200
#define RED_GREEN 400
#define IDAT_SNP_MANIFEST 401
#define SENTRIX_BARCODE 402
#define CHIP_TYPE 403
#define SENTRIX_POSITION 404
#define OPA 405
#define IDAT_SAMPLE_NAME 406
#define DESCRIPTION 407
#define IDAT_SAMPLE_PLATE 408
#define IDAT_SAMPLE_WELL 409
#define UNKNOWN_1 410
#define UNKNOWN_2 510
#define RUN_INFO 300

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
    {"1-95um_multi-swath_for_8x2-5M", 2550870, 2550870, "HumanOmni2.5-8v1"},
    {"1-95um_multi-swath_for_8x2-5M", 2575219, 2575219, "HumanOmni2.5-8v1"},
    {"1-95um_multi-swath_for_8x2-5M", 2563064, 2563064, "HumanOmni25M-8v1-1"},
    {"BeadChip 12x1", 55300, 55300, "humanmethylation27_270596_v1-2 ???"},
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
    {"BeadChip 1x12", 661182, 49163, "HumanHap650Yv3"},
    {"BeadChip 1x40", 1129736, 57373, "Human1Mv1"},
    {"BeadChip 1x40 66", 1078890, 52497, "Human1Mv1"},
    {"BeadChip 24x1x4", 306776, 306776, "InfiniumCore-24v1-2"},
    {"BeadChip 24x1x4", 527136, 527136, "OncoArray-500K"},
    {"BeadChip 24x1x4", 577781, 577781, "HumanCoreExome-24v1-0"},
    {"BeadChip 24x1x4", 623302, 623302, "PsychChip_15048346"},
    {"BeadChip 24x1x4", 623513, 623513, "InfiniumPsychArray-24v1-1"},
    {"BeadChip 24x1x4", 638714, 638714, "PsychChip_v1-1_15073391"},
    {"BeadChip 24x1x4", 663209, 663209, "GSA-24v1-0"},
    {"BeadChip 24x1x4", 710576, 710576, "GSAMD-24v1-0_20011747"},
    {"BeadChip 24x1x4", 710606, 710606, "GSAMD-24v1-0_20011747"},
    {"BeadChip 24x1x4", 710608, 710608, "GSAMD-24v1-0_20011747"},
    {"BeadChip 24x1x4", 716279, 716279, "InfiniumOmniExpress-24v1-2"},
    {"BeadChip 24x1x4", 718963, 718963, "HumanOmniExpress-24-v1-0"},
    {"BeadChip 24x1x4", 719234, 719234, "HumanOmniExpress-24-v1-0"},
    {"BeadChip 24x1x4", 751614, 751614, "GSA-24v3-0"},
    {"BeadChip 24x1x4", 780509, 780509, "GSAMD-24v2-0_20024620"},
    {"BeadChip 24x1x4", 818205, 818205, "GSA-24v2-0"},
    {"BeadChip 2x10", 321354, 37161, "HumanHap300v2"},
    {"BeadChip 2x12", 381079, 29275, "HumanCNV370v1"},
    {"BeadChip 2x20", 561686, 54936, "HumanHap550v3"},
    {"BeadChip 2x6Q", 1224000, 180026, "Human1M-Duov3"},
    {"BeadChip 2x6Q", 1224629, 180026, "Human1M-Duov3"},
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
    {"BeadChip 8x5", 996003, 996003, "HumanOmniExpressExome-8-v1-2"},
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
    size_t capacity;
    buffer_array_t *ilmn_id;
    buffer_array_t *sd;
    buffer_array_t *mean;
    buffer_array_t *nbeads;
    buffer_array_t *mid_block;
    uint8_t red_green[4];
    char *snp_manifest;
    char *sentrix_barcode;
    char *chip_type;
    char *sentrix_position;
    char *opa;
    char *sample_name;
    char *description;
    char *sample_plate;
    char *sample_well;
    char *unknown1;
    char *unknown2;
    RunInfo *run_infos;
    int32_t m_run_infos;
    const char *chip_type_guess;
    const char *imaging_date;
    const char *scanner_data;
} idat_t;

static int idat_read(idat_t *idat, uint16_t id) {
    int i;
    for (i = 0; i < idat->number_toc_entries && id != idat->id[i]; i++)
        ;
    if (i == idat->number_toc_entries) return -1;
    if (hseek(idat->hfile, idat->toc[i], SEEK_SET) < 0)
        error("Fail to seek to position %ld in IDAT %s file\n", idat->toc[i], idat->fn);

    switch (id) {
    case NUM_SNPS_READ:
        read_bytes(idat->hfile, (void *)&idat->num_snps, sizeof(int32_t));
        break;
    case ILLUMINA_ID:
        idat->ilmn_id = buffer_array_init(idat->hfile, idat->capacity, sizeof(int32_t));
        break;
    case SD:
        idat->sd = buffer_array_init(idat->hfile, idat->capacity, sizeof(uint16_t));
        break;
    case MEAN:
        idat->mean = buffer_array_init(idat->hfile, idat->capacity, sizeof(uint16_t));
        break;
    case NBEADS:
        idat->nbeads = buffer_array_init(idat->hfile, idat->capacity, sizeof(uint8_t));
        break;
    case MID_BLOCK:
        idat->mid_block = buffer_array_init(idat->hfile, idat->capacity, sizeof(uint8_t));
        break;
    case RED_GREEN:
        read_bytes(idat->hfile, (void *)&idat->red_green, 4 * sizeof(uint8_t));
        break;
    case IDAT_SNP_MANIFEST:
        read_pfx_string(idat->hfile, &idat->snp_manifest, NULL);
        break;
    case SENTRIX_BARCODE:
        read_pfx_string(idat->hfile, &idat->sentrix_barcode, NULL);
        break;
    case CHIP_TYPE:
        read_pfx_string(idat->hfile, &idat->chip_type, NULL);
        break;
    case SENTRIX_POSITION:
        read_pfx_string(idat->hfile, &idat->sentrix_position, NULL);
        break;
    case OPA:
        read_pfx_string(idat->hfile, &idat->opa, NULL);
        break;
    case IDAT_SAMPLE_NAME:
        read_pfx_string(idat->hfile, &idat->sample_name, NULL);
        break;
    case DESCRIPTION:
        read_pfx_string(idat->hfile, &idat->description, NULL);
        break;
    case IDAT_SAMPLE_PLATE:
        read_pfx_string(idat->hfile, &idat->sample_plate, NULL);
        break;
    case IDAT_SAMPLE_WELL:
        read_pfx_string(idat->hfile, &idat->sample_well, NULL);
        break;
    case UNKNOWN_1:
        read_pfx_string(idat->hfile, &idat->unknown1, NULL);
        break;
    case UNKNOWN_2:
        read_pfx_string(idat->hfile, &idat->unknown2, NULL);
        break;
    case RUN_INFO:
        read_bytes(idat->hfile, (void *)&idat->m_run_infos, sizeof(int32_t));
        idat->run_infos = (RunInfo *)malloc(idat->m_run_infos * sizeof(RunInfo));
        for (int i = 0; i < idat->m_run_infos; i++) {
            read_pfx_string(idat->hfile, &idat->run_infos[i].run_time, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].block_type, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].block_pars, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].block_code, NULL);
            read_pfx_string(idat->hfile, &idat->run_infos[i].code_version, NULL);
        }
        break;
    default:
        error("IDAT file format does not support TOC entry %d\n", id);
        break;
    }
    return 0;
}

static idat_t *idat_init(const char *fn, size_t capacity) {
    idat_t *idat = (idat_t *)calloc(1, sizeof(idat_t));
    idat->fn = strdup(fn);
    idat->hfile = hopen(idat->fn, "rb");
    if (idat->hfile == NULL) error("Could not open %s: %s\n", idat->fn, strerror(errno));
    if (is_gzip(idat->hfile)) error("File %s is gzip compressed and currently cannot be sought\n", idat->fn);

    uint8_t buffer[4];
    if (hread(idat->hfile, (void *)buffer, 4) < 4) error("Failed to read magic number from %s file\n", idat->fn);
    if (memcmp(buffer, "IDAT", 4) != 0) error("IDAT file %s format identifier is bad\n", idat->fn);

    read_bytes(idat->hfile, (void *)&idat->version, sizeof(int64_t));
    if (idat->version < 3)
        error("Cannot read IDAT file %s. Unsupported IDAT file format version: %ld\n", idat->fn, idat->version);

    read_bytes(idat->hfile, (void *)&idat->number_toc_entries, sizeof(int32_t));
    idat->id = (uint16_t *)malloc(idat->number_toc_entries * sizeof(uint16_t));
    idat->toc = (int64_t *)malloc(idat->number_toc_entries * sizeof(int64_t));
    for (int i = 0; i < idat->number_toc_entries; i++) {
        read_bytes(idat->hfile, (void *)&idat->id[i], sizeof(uint16_t));
        read_bytes(idat->hfile, (void *)&idat->toc[i], sizeof(int64_t));
    }

    idat->capacity = capacity;
    for (int i = 0; i < idat->number_toc_entries; i++) idat_read(idat, idat->id[i]);

    for (const chip_type_t *ptr = chip_types; ptr->chip_type; ptr++) {
        if (strcmp(idat->chip_type, ptr->chip_type) == 0 && ptr->num_snps == idat->num_snps
            && ptr->num_mid_blocks == idat->mid_block->item_num)
            idat->chip_type_guess = ptr->chip_type_guess;
    }

    for (int i = 0; i < idat->m_run_infos; i++) {
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
    free(idat->opa);
    free(idat->sample_name);
    free(idat->description);
    free(idat->sample_plate);
    free(idat->sample_well);
    free(idat->unknown1);
    free(idat->unknown2);
    for (int i = 0; i < idat->m_run_infos; i++) {
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

static void idat_to_csv(const idat_t *idat, FILE *stream, int verbose) {

    fprintf(stream, "Illumina, Inc.\n");
    fprintf(stream, "[Heading]\n");
    fprintf(stream, "Descriptor File Name,%s\n", strrchr(idat->fn, '/') ? strrchr(idat->fn, '/') + 1 : idat->fn);
    fprintf(stream, "IDAT file version,%ld\n", idat->version);
    fprintf(stream, "Number of TOC entries,%d\n", idat->number_toc_entries);
    fprintf(stream, "Probes Count,%d\n", idat->num_snps);
    fprintf(stream, "Mid Blocks Count,%d\n", idat->mid_block->item_num);
    fprintf(stream, "Red Green,%02x %02x %02x %02x\n", idat->red_green[0], idat->red_green[1], idat->red_green[2],
            idat->red_green[3]);
    fprintf(stream, "SNP Manifest,%s\n", idat->snp_manifest ? idat->snp_manifest : "");
    fprintf(stream, "Sentrix Barcode,%s\n", idat->sentrix_barcode);
    fprintf(stream, "Chip Type,%s\n", idat->chip_type);
    fprintf(stream, "Sentrix Position,%s\n", idat->sentrix_position);
    fprintf(stream, "OPA,%s\n", idat->opa ? idat->opa : "");
    fprintf(stream, "Sample Name,%s\n", idat->sample_name ? idat->sample_name : "");
    fprintf(stream, "Description,%s\n", idat->description ? idat->description : "");
    fprintf(stream, "Sample Plate,%s\n", idat->sample_plate ? idat->sample_plate : "");
    fprintf(stream, "Sample Well,%s\n", idat->sample_well ? idat->sample_well : "");
    fprintf(stream, "Unknown 1,%s\n", idat->unknown1 ? idat->unknown1 : "");
    fprintf(stream, "Unknown 2,%s\n", idat->unknown2 ? idat->unknown2 : "");
    fprintf(stream, "Chip Prefix (Guess),%s\n", idat->chip_type_guess ? idat->chip_type_guess : "Unknown");
    fprintf(stream, "[Assay]\n");
    fprintf(stream, "IlmnID,Sd,Mean,Nbeads\n");
    if (verbose) {
        for (int i = 0; i < idat->num_snps; i++) {
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
        for (int i = 0; i < idat->mid_block->item_num; i++) {
            int8_t mid_block;
            get_element(idat->mid_block, (void *)&mid_block, i);
            fprintf(stream, "%d\n", mid_block);
        }
    } else {
        fprintf(stream, "... use --verbose to visualize Assay data ...\n");
        fprintf(stream, "[Mid Blocks]\n");
        fprintf(stream, "... use --verbose to visualize Mid Blocks data ...\n");
    }
    fprintf(stream, "[Run Infos]\n");
    for (int i = 0; i < idat->m_run_infos; i++) {
        fprintf(stream, "%s\t%s\t%s\t%s\t%s\n", idat->run_infos[i].run_time, idat->run_infos[i].block_type,
                idat->run_infos[i].block_pars, idat->run_infos[i].block_code, idat->run_infos[i].code_version);
    }
}

static void idats_to_tsv(idat_t **idats, int n, FILE *stream) {
    fprintf(stream,
            "idat\tnumber_probes\tnumber_mid_blocks\tred_green\tmanifest_file\tsentrix_"
            "barcode\tchip_type\t"
            "sentrix_position\topa\tsample_name\tdescription\tsample_plate\tsample_"
            "well\tunknown1\tunknown2\t"
            "chip_type_guess\tscan_date\tscanner_data\n");
    for (int i = 0; i < n; i++) {
        idat_t *idat = idats[i];
        fprintf(stream,
                "%s\t%d\t%d\t%02x %02x %02x "
                "%02x\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                strrchr(idat->fn, '/') ? strrchr(idat->fn, '/') + 1 : idat->fn, idat->num_snps,
                idat->mid_block->item_num, idat->red_green[0], idat->red_green[1], idat->red_green[2],
                idat->red_green[3], idat->snp_manifest ? idat->snp_manifest : "", idat->sentrix_barcode,
                idat->chip_type, idat->sentrix_position, idat->opa ? idat->opa : "",
                idat->sample_name ? idat->sample_name : "", idat->description ? idat->description : "",
                idat->sample_plate ? idat->sample_plate : "", idat->sample_well ? idat->sample_well : "",
                idat->unknown1 ? idat->unknown1 : "", idat->unknown2 ? idat->unknown2 : "",
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
#define PLOIDY 2
#define PLOIDY_TYPE 3
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
#define B_ALLELE_FREQS 1012
#define LOGR_RATIOS 1013
#define PERCENTILES_X 1014
#define PERCENTILES_Y 1015
#define SLIDE_IDENTIFIER 1016

static const char *code2genotype[] = {
    "NC",       "AA",       "AB",       "BB",       "NULL",     "A",        "B",        "AAA",
    "AAB",      "ABB",      "BBB",      "AAAA",     "AAAB",     "AABB",     "ABBB",     "BBBB",
    "AAAAA",    "AAAAB",    "AAABB",    "AABBB",    "ABBBB",    "BBBBB",    "AAAAAA",   "AAAAAB",
    "AAAABB",   "AAABBB",   "AABBBB",   "ABBBBB",   "BBBBBB",   "AAAAAAA",  "AAAAAAB",  "AAAAABB",
    "AAAABBB",  "AAABBBB",  "AABBBBB",  "ABBBBBB",  "BBBBBBB",  "AAAAAAAA", "AAAAAAAB", "AAAAAABB",
    "AAAAABBB", "AAAABBBB", "AAABBBBB", "AABBBBBB", "ABBBBBBB", "BBBBBBBB"};

typedef struct {
    int32_t version;
    float offset_x;
    float offset_y;
    float scale_x;
    float scale_y;
    float shear;
    float theta;
    float reserved[6];
} XForm;

typedef char BaseCall[2];

typedef struct {
    uint16_t raw_x;
    uint16_t raw_y;
    float norm_x;
    float norm_y;
    float ilmn_theta;
    float ilmn_r;
    float baf;
    float lrr;
} intensities_t;

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

    size_t capacity;
    buffer_array_t *raw_x;
    buffer_array_t *raw_y;
    buffer_array_t *genotypes;
    buffer_array_t *base_calls;
    buffer_array_t *genotype_scores;
    buffer_array_t *b_allele_freqs;
    buffer_array_t *logr_ratios;
} gtc_t;

static int gtc_read(gtc_t *gtc, uint16_t id) {
    int i;
    for (i = 0; i < gtc->number_toc_entries && id != gtc->id[i]; i++)
        ;
    if (i == gtc->number_toc_entries) return -1;
    if (id != NUM_SNPS && id != PLOIDY && id != PLOIDY_TYPE) {
        if (hseek(gtc->hfile, gtc->toc[i], SEEK_SET) < 0)
            error("Fail to seek to position %d in GTC %s file \n", gtc->toc[i], gtc->fn);
    }

    switch (id) {
    case NUM_SNPS:
        gtc->num_snps = gtc->toc[i];
        break;
    case PLOIDY:
        gtc->ploidy = gtc->toc[i];
        break;
    case PLOIDY_TYPE:
        gtc->ploidy = gtc->toc[i];
        break;
    case GTC_SAMPLE_NAME:
        read_pfx_string(gtc->hfile, &gtc->sample_name, NULL);
        break;
    case GTC_SAMPLE_PLATE:
        read_pfx_string(gtc->hfile, &gtc->sample_plate, NULL);
        break;
    case GTC_SAMPLE_WELL:
        read_pfx_string(gtc->hfile, &gtc->sample_well, NULL);
        break;
    case CLUSTER_FILE:
        read_pfx_string(gtc->hfile, &gtc->cluster_file, NULL);
        break;
    case GTC_SNP_MANIFEST:
        read_pfx_string(gtc->hfile, &gtc->snp_manifest, NULL);
        break;
    case IMAGING_DATE:
        read_pfx_string(gtc->hfile, &gtc->imaging_date, NULL);
        break;
    case AUTOCALL_DATE:
        read_pfx_string(gtc->hfile, &gtc->autocall_date, NULL);
        break;
    case AUTOCALL_VERSION:
        read_pfx_string(gtc->hfile, &gtc->autocall_version, NULL);
        break;
    case NORMALIZATION_TRANSFORMS:
        read_pfx_array(gtc->hfile, (void **)&gtc->normalization_transforms, &gtc->m_normalization_transforms,
                       sizeof(XForm));
        break;
    case CONTROLS_X:
        read_pfx_array(gtc->hfile, (void **)&gtc->controls_x, &gtc->m_controls_x, sizeof(uint16_t));
        break;
    case CONTROLS_Y:
        read_pfx_array(gtc->hfile, (void **)&gtc->controls_y, &gtc->m_controls_y, sizeof(uint16_t));
        break;
    case RAW_X:
        gtc->raw_x = buffer_array_init(gtc->hfile, gtc->capacity, sizeof(uint16_t));
        break;
    case RAW_Y:
        gtc->raw_y = buffer_array_init(gtc->hfile, gtc->capacity, sizeof(uint16_t));
        break;
    case GENOTYPES:
        gtc->genotypes = buffer_array_init(gtc->hfile, gtc->capacity, sizeof(uint8_t));
        break;
    case BASE_CALLS:
        gtc->base_calls = buffer_array_init(gtc->hfile, gtc->capacity, sizeof(BaseCall));
        break;
    case GENOTYPE_SCORES:
        gtc->genotype_scores = buffer_array_init(gtc->hfile, gtc->capacity, sizeof(float));
        break;
    case SCANNER_DATA:
        read_pfx_string(gtc->hfile, &gtc->scanner_data.scanner_name, NULL);
        read_bytes(gtc->hfile, (void *)&gtc->scanner_data.pmt_green, sizeof(float));
        read_bytes(gtc->hfile, (void *)&gtc->scanner_data.pmt_red, sizeof(float));
        read_pfx_string(gtc->hfile, &gtc->scanner_data.scanner_version, NULL);
        read_pfx_string(gtc->hfile, &gtc->scanner_data.imaging_user, NULL);
        break;
    case CALL_RATE:
        read_bytes(gtc->hfile, (void *)&gtc->call_rate, sizeof(float));
        break;
    case GENDER:
        read_bytes(gtc->hfile, (void *)&gtc->gender, sizeof(char));
        break;
    case LOGR_DEV:
        read_bytes(gtc->hfile, (void *)&gtc->logr_dev, sizeof(float));
        break;
    case GC10:
        read_bytes(gtc->hfile, (void *)&gtc->p10gc, sizeof(float));
        break;
    case DX:
        read_bytes(gtc->hfile, (void *)&gtc->dx, sizeof(int32_t));
        break;
    case SAMPLE_DATA:
        read_bytes(gtc->hfile, (void *)&gtc->sample_data, sizeof(SampleData));
        break;
    case B_ALLELE_FREQS:
        gtc->b_allele_freqs = buffer_array_init(gtc->hfile, gtc->capacity, sizeof(float));
        break;
    case LOGR_RATIOS:
        gtc->logr_ratios = buffer_array_init(gtc->hfile, gtc->capacity, sizeof(float));
        break;
    case PERCENTILES_X:
        read_bytes(gtc->hfile, (void *)&gtc->percentiles_x, sizeof(Percentiles));
        break;
    case PERCENTILES_Y:
        read_bytes(gtc->hfile, (void *)&gtc->percentiles_y, sizeof(Percentiles));
        break;
    case SLIDE_IDENTIFIER:
        read_pfx_string(gtc->hfile, &gtc->sentrix_id, NULL);
        break;
    default:
        error("GTC file format does not support TOC entry %d\n", id);
        break;
    }
    return 0;
}

static gtc_t *gtc_init(const char *fn, size_t capacity) {
    gtc_t *gtc = (gtc_t *)calloc(1, sizeof(gtc_t));
    gtc->fn = strdup(fn);
    gtc->hfile = hopen(gtc->fn, "rb");
    if (gtc->hfile == NULL) error("Could not open %s: %s\n", gtc->fn, strerror(errno));
    if (is_gzip(gtc->hfile)) error("File %s is gzip compressed and currently cannot be sought\n", gtc->fn);

    uint8_t buffer[4];
    if (hread(gtc->hfile, (void *)buffer, 4) < 4) error("Failed to read magic number from %s file\n", gtc->fn);
    if (memcmp(buffer, "gtc", 3) != 0) error("GTC file %s format identifier is bad\n", gtc->fn);
    if (buffer[3] > 5 && buffer[3] < 3) error("GTC file %s version %d is unsupported\n", gtc->fn, buffer[3]);
    gtc->version = (int32_t)buffer[3];

    read_bytes(gtc->hfile, (void *)&gtc->number_toc_entries, sizeof(int32_t));
    gtc->id = (uint16_t *)malloc(gtc->number_toc_entries * sizeof(uint16_t));
    gtc->toc = (int32_t *)malloc(gtc->number_toc_entries * sizeof(int32_t));
    for (int i = 0; i < gtc->number_toc_entries; i++) {
        read_bytes(gtc->hfile, (void *)&gtc->id[i], sizeof(uint16_t));
        read_bytes(gtc->hfile, (void *)&gtc->toc[i], sizeof(int32_t));
    }

    gtc->capacity = capacity;
    for (int i = 0; i < gtc->number_toc_entries; i++) gtc_read(gtc, gtc->id[i]);

    const char *ptr = strrchr(gtc->fn, '/') ? strrchr(gtc->fn, '/') + 1 : gtc->fn;
    gtc->display_name = strndup(ptr, strlen(ptr) - 4);

    gtc->sin_theta = (float *)malloc(gtc->m_normalization_transforms * sizeof(float));
    gtc->cos_theta = (float *)malloc(gtc->m_normalization_transforms * sizeof(float));
    for (int i = 0; i < gtc->m_normalization_transforms; i++) {
        gtc->sin_theta[i] = sinf(gtc->normalization_transforms[i].theta);
        gtc->cos_theta[i] = cosf(gtc->normalization_transforms[i].theta);
    }

    return gtc;
}

static void gtc_destroy(gtc_t *gtc) {
    if (!gtc) return;
    if (hclose(gtc->hfile) < 0) error("Error closing GTC file %s\n", gtc->fn);
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

    buffer_array_destroy(gtc->raw_x);
    buffer_array_destroy(gtc->raw_y);
    buffer_array_destroy(gtc->genotypes);
    buffer_array_destroy(gtc->base_calls);
    buffer_array_destroy(gtc->genotype_scores);
    buffer_array_destroy(gtc->b_allele_freqs);
    buffer_array_destroy(gtc->logr_ratios);
    free(gtc);
}

static void gtc_to_csv(const gtc_t *gtc, FILE *stream, int verbose) {
    fprintf(stream, "Illumina, Inc.\n");
    fprintf(stream, "[Heading]\n");
    fprintf(stream, "Descriptor File Name,%s\n", strrchr(gtc->fn, '/') ? strrchr(gtc->fn, '/') + 1 : gtc->fn);
    fprintf(stream, "GTC genotype file version,%d\n", gtc->version);
    fprintf(stream, "Number of TOC entries,%d\n", gtc->number_toc_entries);
    fprintf(stream, "Number of SNPs,%d\n", gtc->num_snps);
    fprintf(stream, "Ploidy,%d\n", gtc->ploidy);
    fprintf(stream, "Ploidy Type,%d\n", gtc->ploidy_type);
    fprintf(stream, "Sample name,%s\n", gtc->sample_name ? gtc->sample_name : "");
    fprintf(stream, "Sample plate,%s\n", gtc->sample_plate ? gtc->sample_plate : "");
    fprintf(stream, "Sample well,%s\n", gtc->sample_well ? gtc->sample_well : "");
    fprintf(stream, "Cluster file,%s\n", gtc->cluster_file);
    fprintf(stream, "SNP manifest,%s\n", gtc->snp_manifest);
    fprintf(stream, "Imaging date,%s\n", gtc->imaging_date);
    fprintf(stream, "AutoCall date,%s\n", gtc->autocall_date);
    fprintf(stream, "AutoCall version,%s\n", gtc->autocall_version);
    fprintf(stream, "Number of normalization transforms,%ld\n", gtc->m_normalization_transforms);
    fprintf(stream, "Number of controls X,%ld\n", gtc->m_controls_x);
    fprintf(stream, "Number of controls Y,%ld\n", gtc->m_controls_y);
    fprintf(stream, "Name of the scanner,%s\n", gtc->scanner_data.scanner_name);
    fprintf(stream, "Pmt Green,%d\n", gtc->scanner_data.pmt_green);
    fprintf(stream, "Pmt Red,%d\n", gtc->scanner_data.pmt_red);
    fprintf(stream, "Version of the scanner software used,%s\n", gtc->scanner_data.scanner_version);
    fprintf(stream, "Name of the scanner user,%s\n", gtc->scanner_data.imaging_user);
    fprintf(stream, "Call Rate,%f\n", gtc->call_rate);
    fprintf(stream, "Computed Gender,%c\n", gtc->gender);
    fprintf(stream, "LogR deviation,%f\n", gtc->logr_dev);
    fprintf(stream, "GenCall score - 10th percentile,%f\n", gtc->p10gc);
    fprintf(stream, "DX,%d\n", gtc->dx);
    fprintf(stream, "GenCall score - 50th percentile,%f\n", gtc->sample_data.p50gc);
    fprintf(stream, "Number of valid calls,%d\n", gtc->sample_data.num_calls);
    fprintf(stream, "Number of invalid calls,%d\n", gtc->sample_data.num_no_calls);
    fprintf(stream, "Number of loci that are \"Intensity Only\" or \"Zeroed\",%d\n",
            gtc->sample_data.num_intensity_only);
    fprintf(stream, "P05 X,%d\n", gtc->percentiles_x[0]);
    fprintf(stream, "P50 X,%d\n", gtc->percentiles_x[1]);
    fprintf(stream, "P95 X,%d\n", gtc->percentiles_x[2]);
    fprintf(stream, "P05 Y,%d\n", gtc->percentiles_y[0]);
    fprintf(stream, "P50 Y,%d\n", gtc->percentiles_y[1]);
    fprintf(stream, "P95 Y,%d\n", gtc->percentiles_y[2]);
    fprintf(stream, "Sentrix identifier for the slide,%s\n", gtc->sentrix_id ? gtc->sentrix_id : "");
    fprintf(stream, "[Assay]\n");
    fprintf(stream, "Raw X,Raw Y,GType,Top Alleles,Score,B Allele Freq,Log R Ratio\n");
    if (verbose) {
        for (int i = 0; i < gtc->num_snps; i++) {
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
            fprintf(stream, "%d,%d,%s,%c%c,%f,%f,%f\n", raw_x, raw_y, code2genotype[genotype], base_call[0],
                    base_call[1], genotype_score, b_allele_freq, logr_ratio);
        }
    } else {
        fprintf(stream, "... use --verbose to visualize assay data ...\n");
    }
    fprintf(stream, "[Normalization Transforms]\n");
    fprintf(stream, "Version,Offset X,Offset Y,Scale X,Scale Y,Shear,Theta\n");
    if (verbose) {
        for (int i = 0; i < gtc->m_normalization_transforms; i++)
            fprintf(stream, "%d,%f,%f,%f,%f,%f,%f\n", gtc->normalization_transforms[i].version,
                    gtc->normalization_transforms[i].offset_x, gtc->normalization_transforms[i].offset_y,
                    gtc->normalization_transforms[i].scale_x, gtc->normalization_transforms[i].scale_y,
                    gtc->normalization_transforms[i].shear, gtc->normalization_transforms[i].theta);
    } else {
        fprintf(stream, "... use --verbose to visualize assay data ...\n");
    }
}

static void gtcs_to_tsv(gtc_t **gtcs, int n, FILE *stream) {
    fprintf(stream,
            "gtc\tnumber_snps\tploidy\tploidy_type\tsample_name\tsample_plate\tsample_"
            "well\tcluster_file\tsnp_manifest\t"
            "scan_date\tautocall_date\tautocall_version\tnumber_normalization_"
            "transforms\tnumber_x_controls\t"
            "number_y_controls\tscanner_name\tpmt_green\tpmt_red\tscanner_software_"
            "version\tscanner_username\tcall_rate\t"
            "computed_gender\tlogr_deviation\tgencall_score_10_percentile\tdx\tgencall_score_"
            "50_percentile\t"
            "number_valid_calls\tnumber_invalid_calls\tnumber_intensity_only_or_zeroed_"
            "loci\tp05_x\tp50_x\tp95_x\tp05_y\t"
            "p50_y\tp95_y\tsentrix_barcode\n");
    for (int i = 0; i < n; i++) {
        gtc_t *gtc = gtcs[i];
        fprintf(stream,
                "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%ld\t%ld\t%ld\t%s\t%d\t%d\t%"
                "s\t%s\t%f\t%c\t%f\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
                strrchr(gtc->fn, '/') ? strrchr(gtc->fn, '/') + 1 : gtc->fn, gtc->num_snps, gtc->ploidy,
                gtc->ploidy_type, gtc->sample_name ? gtc->sample_name : "", gtc->sample_plate ? gtc->sample_plate : "",
                gtc->sample_well ? gtc->sample_well : "", gtc->cluster_file, gtc->snp_manifest, gtc->imaging_date,
                gtc->autocall_date, gtc->autocall_version, gtc->m_normalization_transforms, gtc->m_controls_x,
                gtc->m_controls_y, gtc->scanner_data.scanner_name, gtc->scanner_data.pmt_green,
                gtc->scanner_data.pmt_red, gtc->scanner_data.scanner_version, gtc->scanner_data.imaging_user,
                gtc->call_rate, gtc->gender, gtc->logr_dev, gtc->p10gc, gtc->dx, gtc->sample_data.p50gc,
                gtc->sample_data.num_calls, gtc->sample_data.num_no_calls, gtc->sample_data.num_intensity_only,
                gtc->percentiles_x[0], gtc->percentiles_x[1], gtc->percentiles_x[2], gtc->percentiles_y[0],
                gtc->percentiles_y[1], gtc->percentiles_y[2], gtc->sentrix_id ? gtc->sentrix_id : "");
    }
}

/****************************************
 * SAM FILE IMPLEMENTATION              *
 ****************************************/

static bpm_t *sam_csv_init(const char *fn, bpm_t *bpm, const char *genome_build, int flags) {
    htsFile *hts = hts_open(fn, "r");
    if (hts == NULL || hts_get_format(hts)->category != sequence_data)
        error("File %s does not contain sequence data\n", fn);
    sam_hdr_t *sam_hdr = sam_hdr_read(hts);
    if (sam_hdr == NULL) error("Reading header from \"%s\" failed", fn);
    bam1_t *b = bam_init1();
    if (b == NULL) error("Cannot create SAM record\n");

    kstring_t str = {0, 0, NULL};
    const char *chromosome = NULL;
    int strand = -1, position = 0, n_unmapped = 0;
    for (int i = 0; i < bpm->num_loci; i++) {
        LocusEntry *locus_entry = &bpm->locus_entries[i];
        int idx = get_position(hts, sam_hdr, b, locus_entry->ilmn_id, locus_entry->source_seq, 1, &chromosome,
                               &position, &strand);
        if (idx < 0) {
            error("Reading from %s failed", fn);
        } else if (idx == 0) {
            if (flags & VERBOSE) fprintf(stderr, "Unable to determine position for marker %s\n", locus_entry->ilmn_id);
            n_unmapped++;
        }
        free(locus_entry->genome_build);
        locus_entry->genome_build = strdup(genome_build);
        free(locus_entry->chrom);
        locus_entry->chrom = strdup(chromosome ? chromosome : "0");
        free(locus_entry->map_info);
        str.l = 0;
        kputw(position, &str);
        locus_entry->map_info = strdup(str.s);
        free(locus_entry->ref_strand);
        locus_entry->ref_strand =
            ((strand < 0) || ((strcasecmp(locus_entry->ilmn_strand, locus_entry->source_strand) != 0) == strand))
                ? strdup("+")
                : strdup("-");
    }
    fprintf(stderr, "Lines   total/unmapped:\t%d/%d\n", bpm->num_loci, n_unmapped);
    free(str.s);

    bam_destroy1(b);
    sam_hdr_destroy(sam_hdr);
    if (hts_close(hts) < 0) error("closing \"%s\" failed", fn);
    return bpm;
}

/****************************************
 * INTENSITIES COMPUTATIONS             *
 ****************************************/

// compute normalized X Y intensities
static void get_norm_xy(uint16_t raw_x, uint16_t raw_y, const gtc_t *gtc, const bpm_t *bpm, int idx, float *norm_x,
                        float *norm_y) {
    if (bpm->norm_lookups && bpm->locus_entries[idx].norm_id != 0xFF) {
        int norm_id = bpm->norm_lookups[bpm->locus_entries[idx].norm_id];
        XForm *xform = gtc->normalization_transforms + norm_id;
        float temp_x = (float)raw_x - xform->offset_x;
        float temp_y = (float)raw_y - xform->offset_y;
        float temp_x2 = gtc->cos_theta[norm_id] * temp_x + gtc->sin_theta[norm_id] * temp_y;
        float temp_y2 = -gtc->sin_theta[norm_id] * temp_x + gtc->cos_theta[norm_id] * temp_y;
        float temp_x3 = temp_x2 - xform->shear * temp_y2;
        *norm_x = temp_x3 < 0.0f ? 0.0f : temp_x3 / xform->scale_x;
        *norm_y = temp_y2 < 0.0f ? 0.0f : temp_y2 / xform->scale_y;
    } else {
        *norm_x = NAN;
        *norm_y = NAN;
    }
}

// compute Theta and R from raw intensities
static inline void get_ilmn_theta_r(float norm_x, float norm_y, float *ilmn_theta, float *ilmn_r) {
    *ilmn_theta = atanf(norm_y / norm_x) * (float)M_2_PI;
    *ilmn_r = norm_x + norm_y;
}

static inline void get_intensities(gtc_t *gtc, const bpm_t *bpm, const egt_t *egt, int idx,
                                   intensities_t *intensities) {
    get_element(gtc->raw_x, (void *)&intensities->raw_x, idx);
    get_element(gtc->raw_y, (void *)&intensities->raw_y, idx);
    get_norm_xy(intensities->raw_x, intensities->raw_y, gtc, bpm, idx, &intensities->norm_x, &intensities->norm_y);
    get_ilmn_theta_r(intensities->norm_x, intensities->norm_y, &intensities->ilmn_theta, &intensities->ilmn_r);

    if (bpm->norm_lookups && egt) {
        get_baf_lrr(intensities->ilmn_theta, intensities->ilmn_r, egt->cluster_records[idx].aa_cluster_stats.theta_mean,
                    egt->cluster_records[idx].ab_cluster_stats.theta_mean,
                    egt->cluster_records[idx].bb_cluster_stats.theta_mean,
                    egt->cluster_records[idx].aa_cluster_stats.r_mean,
                    egt->cluster_records[idx].ab_cluster_stats.r_mean,
                    egt->cluster_records[idx].bb_cluster_stats.r_mean, &intensities->baf, &intensities->lrr);
    } else if (gtc->b_allele_freqs && gtc->logr_ratios) {
        get_element(gtc->b_allele_freqs, (void *)&intensities->baf, idx);
        get_element(gtc->logr_ratios, (void *)&intensities->lrr, idx);
    } else {
        intensities->baf = NAN;
        intensities->lrr = NAN;
    }
}

static void adjust_clusters(const uint8_t *gts, const float *ilmn_theta, const float *ilmn_r, int n,
                            ClusterRecord *cluster_record) {
    cluster_record->aa_cluster_stats.N = 0;
    cluster_record->ab_cluster_stats.N = 0;
    cluster_record->bb_cluster_stats.N = 0;
    cluster_record->aa_cluster_stats.theta_mean *= 0.2f;
    cluster_record->ab_cluster_stats.theta_mean *= 0.2f;
    cluster_record->bb_cluster_stats.theta_mean *= 0.2f;
    cluster_record->aa_cluster_stats.r_mean *= 0.2f;
    cluster_record->ab_cluster_stats.r_mean *= 0.2f;
    cluster_record->bb_cluster_stats.r_mean *= 0.2f;

    for (int i = 0; i < n; i++) {
        switch (gts[i]) {
        case GT_AA:
            cluster_record->aa_cluster_stats.N++;
            cluster_record->aa_cluster_stats.theta_mean += ilmn_theta[i];
            cluster_record->aa_cluster_stats.r_mean += ilmn_r[i];
            break;
        case GT_AB:
            cluster_record->ab_cluster_stats.N++;
            cluster_record->ab_cluster_stats.theta_mean += ilmn_theta[i];
            cluster_record->ab_cluster_stats.r_mean += ilmn_r[i];
            break;
        case GT_BB:
            cluster_record->bb_cluster_stats.N++;
            cluster_record->bb_cluster_stats.theta_mean += ilmn_theta[i];
            cluster_record->bb_cluster_stats.r_mean += ilmn_r[i];
            break;
        default:
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
 * CONVERSION UTILITIES                 *
 ****************************************/

static inline char rev_allele(char allele) {
    static const char allele_complement[128] = {
        0, 0,   0, 0,   0,   0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0,   0, 0,   0,   0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 'T', 0, 'G', 'I', 0, 0, 'C', 0, 'D', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    if (allele > 95) return 0;
    return allele_complement[(int)allele];
}

static void gtcs_to_gs(gtc_t **gtc, int n, const bpm_t *bpm, const egt_t *egt, FILE *stream, int flags) {
    // print header
    fputs("Index\tName\tAddress\tChr\tPosition", stream);
    if (flags & EGT_LOADED) fputs("\tGenTrain Score", stream);
    if (flags & BPM_LOADED) fputs("\tFrac A\tFrac C\tFrac G\tFrac T", stream);
    for (int i = 0; i < n; i++) {
        if (flags & FORMAT_GT) fprintf(stream, "\t%s.GType", gtc[i]->display_name);
        if (flags & FORMAT_IGC) fprintf(stream, "\t%s.Score", gtc[i]->display_name);
        if ((flags & BPM_LOADED) && (flags & FORMAT_THETA)) fprintf(stream, "\t%s.Theta", gtc[i]->display_name);
        if ((flags & BPM_LOADED) && (flags & FORMAT_R)) fprintf(stream, "\t%s.R", gtc[i]->display_name);
        if (flags & FORMAT_BAF) fprintf(stream, "\t%s.B Allele Freq", gtc[i]->display_name);
        if (flags & FORMAT_LRR) fprintf(stream, "\t%s.Log R Ratio", gtc[i]->display_name);
        if (flags & FORMAT_X) fprintf(stream, "\t%s.X Raw", gtc[i]->display_name);
        if (flags & FORMAT_Y) fprintf(stream, "\t%s.Y Raw", gtc[i]->display_name);
        if ((flags & BPM_LOADED) && (flags & FORMAT_NORMX)) fprintf(stream, "\t%s.X", gtc[i]->display_name);
        if ((flags & BPM_LOADED) && (flags & FORMAT_NORMY)) fprintf(stream, "\t%s.Y", gtc[i]->display_name);
        if (flags & FORMAT_GT)
            fprintf(stream, "\t%s.Top Alleles\t%s.Plus/Minus Alleles", gtc[i]->display_name, gtc[i]->display_name);
    }
    fputc('\n', stream);

    // print loci
    for (int j = 0; j < bpm->num_loci; j++) {
        LocusEntry *locus_entry = &bpm->locus_entries[j];
        int strand = !locus_entry->ref_strand ? -1
                                              : (strcmp(locus_entry->ref_strand, "+") == 0
                                                     ? 0
                                                     : (strcmp(locus_entry->ref_strand, "-") == 0 ? 1 : -1));
        if (strand < 0) error("Unable to process reference strand %s\n", locus_entry->ref_strand);
        fprintf(stream, "%d\t%s\t%d\t%s\t%s", bpm->indexes ? bpm->indexes[j] : j, locus_entry->name,
                locus_entry->address_a, locus_entry->chrom, locus_entry->map_info);
        if (flags & EGT_LOADED) fprintf(stream, "\t%f", egt->cluster_records[j].cluster_score.total_score);
        if (flags & BPM_LOADED)
            fprintf(stream, "\t%f\t%f\t%f\t%f", locus_entry->frac_a, locus_entry->frac_c, locus_entry->frac_g,
                    locus_entry->frac_t);
        for (int i = 0; i < n; i++) {
            uint8_t genotype;
            get_element(gtc[i]->genotypes, (void *)&genotype, j);
            float genotype_score;
            get_element(gtc[i]->genotype_scores, (void *)&genotype_score, j);
            BaseCall base_call;
            get_element(gtc[i]->base_calls, (void *)&base_call, j);
            intensities_t intensities;
            get_intensities(gtc[i], bpm, egt, j, &intensities);
            char allele_a = strand ? rev_allele(locus_entry->snp[1]) : locus_entry->snp[1];
            char allele_b = strand ? rev_allele(locus_entry->snp[3]) : locus_entry->snp[3];
            BaseCall ref_call;
            switch (genotype) {
            case GT_NC:
                ref_call[0] = '-';
                ref_call[1] = '-';
                break;
            case GT_AA:
                ref_call[0] = allele_a;
                ref_call[1] = allele_a;
                break;
            case GT_AB:
                ref_call[0] = allele_a;
                ref_call[1] = allele_b;
                break;
            case GT_BB:
                ref_call[0] = allele_b;
                ref_call[1] = allele_b;
                break;
            default:
                error("Unable to process marker %s\n", locus_entry->name);
                break;
            }
            if (flags & FORMAT_GT) fprintf(stream, "\t%s", code2genotype[genotype]);
            if (flags & FORMAT_IGC) fprintf(stream, "\t%f", genotype_score);
            if ((flags & BPM_LOADED) && (flags & FORMAT_THETA)) fprintf(stream, "\t%f", intensities.ilmn_theta);
            if ((flags & BPM_LOADED) && (flags & FORMAT_R)) fprintf(stream, "\t%f", intensities.ilmn_r);
            if (flags & FORMAT_BAF) fprintf(stream, "\t%f", intensities.baf);
            if (flags & FORMAT_LRR) fprintf(stream, "\t%f", intensities.lrr);
            if (flags & FORMAT_X) fprintf(stream, "\t%d", intensities.raw_x);
            if (flags & FORMAT_Y) fprintf(stream, "\t%d", intensities.raw_y);
            if ((flags & BPM_LOADED) && (flags & FORMAT_NORMX)) fprintf(stream, "\t%f", intensities.norm_x);
            if ((flags & BPM_LOADED) && (flags & FORMAT_NORMY)) fprintf(stream, "\t%f", intensities.norm_y);
            if (flags & FORMAT_GT)
                fprintf(stream, "\t%c%c\t%c%c", base_call[0], base_call[1], ref_call[0], ref_call[1]);
        }
        fputc('\n', stream);
    }
}

static bcf_hdr_t *hdr_init(const faidx_t *fai, int flags) {
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    int n = faidx_nseq(fai);
    for (int i = 0; i < n; i++) {
        const char *seq = faidx_iseq(fai, i);
        int len = faidx_seq_len(fai, seq);
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq, len);
    }

    bcf_hdr_append(hdr, "##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description=\"A allele\">");
    bcf_hdr_append(hdr, "##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description=\"B allele\">");
    if (flags & BPM_LOADED) {
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_A,Number=1,Type=Float,Description=\"Fraction of the A "
                       "nucleotide in the top genomic sequence\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_C,Number=1,Type=Float,Description=\"Fraction of the C "
                       "nucleotide in the top genomic sequence\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_G,Number=1,Type=Float,Description=\"Fraction of the G "
                       "nucleotide in the top genomic sequence\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_T,Number=1,Type=Float,Description=\"Fraction of the T "
                       "nucleotide in the top genomic sequence\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=NORM_ID,Number=1,Type=Integer,Description=\"Normalization "
                       "lookups from manifest\">");
    }
    if (flags & CSV_LOADED) {
        bcf_hdr_append(hdr,
                       "##INFO=<ID=BEADSET_ID,Number=1,Type=Integer,Description=\"Bead set ID "
                       "for normalization\">");
    }
    if ((flags & BPM_LOADED) | (flags & CSV_LOADED)) {
        bcf_hdr_append(hdr,
                       "##INFO=<ID=ASSAY_TYPE,Number=1,Type=Integer,Description=\"Identifies type of "
                       "assay (0 - Infinium II, 1 - Infinium I (A/T), 2 - Infinium I (G/C)\">");
    }
    if (flags & EGT_LOADED) {
        bcf_hdr_append(hdr,
                       "##INFO=<ID=GenTrain_Score,Number=1,Type=Float,Description=\"The SNP "
                       "cluster quality from the GenTrain clustering algorithm\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=Orig_Score,Number=1,Type=Float,Description=\"The original "
                       "GenTrain score for the SNP before edits\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=Edited,Number=0,Type=Flag,Description=\"The SNP was edited "
                       "after identifying clustering positions\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=Cluster_Sep,Number=1,Type=Float,Description=\"The cluster "
                       "separation measurement for the SNP that ranges between 0 and 1\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=N_AA,Number=1,Type=Integer,Description=\"Number of AA calls "
                       "in training set\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=N_AB,Number=1,Type=Integer,Description=\"Number of AB calls "
                       "in training set\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=N_BB,Number=1,Type=Integer,Description=\"Number of BB calls "
                       "in training set\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=devR_AA,Number=1,Type=Float,Description=\"Standard "
                       "deviation of normalized R for AA cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=devR_AB,Number=1,Type=Float,Description=\"Standard "
                       "deviation of normalized R for AB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=devR_BB,Number=1,Type=Float,Description=\"Standard "
                       "deviation of normalized R for BB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=devTHETA_AA,Number=1,Type=Float,Description=\"Standard "
                       "deviation of normalized THETA for AA cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=devTHETA_AB,Number=1,Type=Float,Description=\"Standard "
                       "deviation of normalized THETA for AB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=devTHETA_BB,Number=1,Type=Float,Description=\"Standard "
                       "deviation of normalized THETA for BB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=meanR_AA,Number=1,Type=Float,Description=\"Mean of "
                       "normalized R for AA cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=meanR_AB,Number=1,Type=Float,Description=\"Mean of "
                       "normalized R for AB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=meanR_BB,Number=1,Type=Float,Description=\"Mean of "
                       "normalized R for BB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=meanTHETA_AA,Number=1,Type=Float,Description=\"Mean of "
                       "normalized THETA for AA cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=meanTHETA_AB,Number=1,Type=Float,Description=\"Mean of "
                       "normalized THETA for AB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=meanTHETA_BB,Number=1,Type=Float,Description=\"Mean of "
                       "normalized THETA for BB cluster\">");
        bcf_hdr_append(hdr,
                       "##INFO=<ID=Intensity_Threshold,Number=1,Type=Float,Description=\"The "
                       "intensity threshold value\">");
    }
    if (!(flags & NO_INFO_GC))
        bcf_hdr_append(hdr,
                       "##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC ratio content "
                       "around the variant\">");

    if (flags & FORMAT_GT) bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if (flags & FORMAT_GQ)
        bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    if (flags & FORMAT_IGC)
        bcf_hdr_append(hdr,
                       "##FORMAT=<ID=IGC,Number=1,Type=Float,Description=\"Illumina GenCall "
                       "Confidence Score\">");
    if (flags & FORMAT_BAF)
        bcf_hdr_append(hdr, "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B Allele Frequency\">");
    if (flags & FORMAT_LRR) bcf_hdr_append(hdr, "##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R Ratio\">");
    if ((flags & BPM_LOADED) | (flags & GENOME_STUDIO)) {
        if (flags & FORMAT_NORMX)
            bcf_hdr_append(hdr,
                           "##FORMAT=<ID=NORMX,Number=1,Type=Float,Description=\"Normalized X "
                           "intensity\">");
        if (flags & FORMAT_NORMY)
            bcf_hdr_append(hdr,
                           "##FORMAT=<ID=NORMY,Number=1,Type=Float,Description=\"Normalized Y "
                           "intensity\">");
        if (flags & FORMAT_R)
            bcf_hdr_append(hdr, "##FORMAT=<ID=R,Number=1,Type=Float,Description=\"Normalized R value\">");
        if (flags & FORMAT_THETA)
            bcf_hdr_append(hdr,
                           "##FORMAT=<ID=THETA,Number=1,Type=Float,Description=\"Normalized "
                           "Theta value\">");
    }
    if (flags & FORMAT_X) bcf_hdr_append(hdr, "##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"Raw X intensity\">");
    if (flags & FORMAT_Y) bcf_hdr_append(hdr, "##FORMAT=<ID=Y,Number=1,Type=Integer,Description=\"Raw Y intensity\">");
    return hdr;
}

static int gts_to_gt_arr(int32_t *gt_arr, const uint8_t *gts, int n, int allele_a_idx, int allele_b_idx) {
    for (int i = 0; i < n; i++) {
        switch (gts[i]) {
        case GT_NC:
            gt_arr[2 * i] = bcf_gt_missing;
            gt_arr[2 * i + 1] = bcf_gt_missing;
            break;
        case GT_AA:
            gt_arr[2 * i] = bcf_gt_unphased(allele_a_idx);
            gt_arr[2 * i + 1] = bcf_gt_unphased(allele_a_idx);
            break;
        case GT_AB:
            gt_arr[2 * i] = bcf_gt_unphased(min(allele_a_idx, allele_b_idx));
            gt_arr[2 * i + 1] = bcf_gt_unphased(max(allele_a_idx, allele_b_idx));
            break;
        case GT_BB:
            gt_arr[2 * i] = bcf_gt_unphased(allele_b_idx);
            gt_arr[2 * i + 1] = bcf_gt_unphased(allele_b_idx);
            break;
        default:
            return -1;
        }
    }
    return 0;
}

static void gtcs_to_vcf(faidx_t *fai, const bpm_t *bpm, const egt_t *egt, gtc_t **gtc, int n, htsFile *out_fh,
                        bcf_hdr_t *hdr, int flags, int gc_win) {
    if (bcf_hdr_write(out_fh, hdr) < 0) error("Unable to write to output VCF file\n");
    bcf1_t *rec = bcf_init();
    char ref_base[] = {'\0', '\0'};
    kstring_t allele_a = {0, 0, NULL};
    kstring_t allele_b = {0, 0, NULL};
    kstring_t flank = {0, 0, NULL};

    uint8_t *gts = (uint8_t *)malloc(n * sizeof(uint8_t));
    int32_t *gt_arr = (int32_t *)malloc(n * 2 * sizeof(int32_t));
    int32_t *gq_arr = (int32_t *)malloc(n * sizeof(int32_t));
    float *igc_arr = (float *)malloc(n * sizeof(float));
    float *baf_arr = (float *)malloc(n * sizeof(float));
    float *lrr_arr = (float *)malloc(n * sizeof(float));
    float *norm_x_arr = (float *)malloc(n * sizeof(float));
    float *norm_y_arr = (float *)malloc(n * sizeof(float));
    float *ilmn_r_arr = (float *)malloc(n * sizeof(float));
    float *ilmn_theta_arr = (float *)malloc(n * sizeof(float));
    int32_t *raw_x_arr = (int32_t *)malloc(n * sizeof(int32_t));
    int32_t *raw_y_arr = (int32_t *)malloc(n * sizeof(int32_t));

    int n_missing = 0, n_skipped = 0;
    for (int j = 0; j < bpm->num_loci; j++) {
        LocusEntry *locus_entry = &bpm->locus_entries[j];
        bcf_clear(rec);
        rec->n_sample = n;
        rec->rid = bcf_hdr_name2id_flexible(hdr, locus_entry->chrom);
        char *endptr;
        rec->pos = strtol(locus_entry->map_info, &endptr, 0) - 1;
        if (locus_entry->map_info == endptr)
            error("Map info %s for marker %s is not understood\n", locus_entry->map_info, locus_entry->ilmn_id);
        int strand = !locus_entry->ref_strand ? -1
                                              : (strcmp(locus_entry->ref_strand, "+") == 0
                                                     ? 0
                                                     : (strcmp(locus_entry->ref_strand, "-") == 0 ? 1 : -1));
        if (rec->rid < 0 || rec->pos < 0) {
            if (flags & VERBOSE) fprintf(stderr, "Skipping unlocalized marker %s\n", locus_entry->ilmn_id);
            n_skipped++;
            continue;
        }
        bcf_update_id(hdr, rec, locus_entry->name);

        int len, win = min(max(100, locus_entry->source_seq ? max(gc_win, strlen(locus_entry->source_seq)) : gc_win),
                           rec->pos);
        char *ref = faidx_fetch_seq(fai, bcf_seqname(hdr, rec), rec->pos - win, rec->pos + win, &len);
        if (!ref || len == 1)
            error("faidx_fetch_seq failed at %s:%" PRId64 " (are you using the correct reference genome?)\n",
                  bcf_seqname(hdr, rec), rec->pos + 1);
        strupper(ref);
        if (!(flags & NO_INFO_GC)) {
            float gc_ratio = get_gc_ratio(&ref[max(win - gc_win, 0)], &ref[min(win + gc_win, len)]);
            bcf_update_info_float(hdr, rec, "GC", &gc_ratio, 1);
        }
        ref_base[0] = ref[win];
        int32_t allele_b_idx;
        allele_a.l = allele_b.l = 0;
        kputc(locus_entry->snp[1], &allele_a);
        kputc(locus_entry->snp[3], &allele_b);
        int is_indel = allele_a.s[0] == 'D' || allele_a.s[0] == 'I' || allele_b.s[0] == 'D' || allele_b.s[0] == 'I';
        int ref_is_del = -1;
        if (is_indel && strand >= 0 && locus_entry->source_seq && strchr(locus_entry->source_seq, '-')) {
            flank.l = 0;
            kputs(locus_entry->source_seq, &flank);
            strupper(flank.s);
            if ((strcasecmp(locus_entry->ilmn_strand, locus_entry->source_strand) != 0) != strand)
                flank_reverse_complement(flank.s);
            flank_left_shift(flank.s);

            ref_is_del = get_indel_alleles(&allele_a, &allele_b, flank.s, ref, win, len);
            if (ref_is_del == 0) {
                rec->pos--;
                ref_base[0] = ref[win - 1];
            }
            allele_b_idx = ref_is_del < 0 ? 1 : ref_is_del ^ (locus_entry->snp[3] == 'D');
        } else {
            if (allele_a.s[0] == 'N' && allele_b.s[0] == 'A') {
                allele_a.s[0] = '.';
                allele_b.s[0] = '.';
            } else if (is_indel) {
                ref_base[0] = allele_a.s[0];
            } else {
                if (strand < 0) {
                    if (strcmp(locus_entry->ilmn_strand, "BOT") == 0 || strcmp(locus_entry->ilmn_strand, "Bot") == 0) {
                        allele_a.s[0] = rev_nt(allele_a.s[0]);
                        allele_b.s[0] = rev_nt(allele_b.s[0]);
                    }
                    strand = get_strand_from_top_alleles(allele_a.s, allele_b.s, ref, win, len);
                    if (strand < 0) {
                        if (flags & VERBOSE)
                            fprintf(stderr, "Unable to determine reference strand for SNP %s\n", locus_entry->ilmn_id);
                        allele_a.s[0] = '.';
                        allele_b.s[0] = '.';
                    }
                }
                if (strand == 1) {
                    allele_a.s[0] = rev_nt(allele_a.s[0]);
                    allele_b.s[0] = rev_nt(allele_b.s[0]);
                }
            }
            allele_b_idx = get_allele_b_idx(ref_base[0], allele_a.s, allele_b.s);
        }
        if (is_indel && ref_is_del < 0) {
            if (flags & VERBOSE) fprintf(stderr, "Unable to determine alleles for indel %s\n", locus_entry->ilmn_id);
            n_missing++;
        }
        free(ref);

        int32_t allele_a_idx = get_allele_a_idx(allele_b_idx);
        const char *alleles[3];
        int nals = alleles_ab_to_vcf(alleles, ref_base, allele_a.s, allele_b.s, allele_b_idx);
        if (nals < 0) error("Unable to process marker %s\n", locus_entry->ilmn_id);
        bcf_update_alleles(hdr, rec, alleles, nals);
        bcf_update_info_int32(hdr, rec, "ALLELE_A", &allele_a_idx, 1);
        bcf_update_info_int32(hdr, rec, "ALLELE_B", &allele_b_idx, 1);

        for (int i = 0; i < n; i++) {
            get_element(gtc[i]->genotypes, (void *)&gts[i], j);
            get_element(gtc[i]->genotype_scores, (void *)&igc_arr[i], j);
            gq_arr[i] = (int)(-10 * log10(1 - igc_arr[i]) + .5);
            if (gq_arr[i] < 0) gq_arr[i] = 0;
            if (gq_arr[i] > 50) gq_arr[i] = 50;
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

        if (flags & BPM_LOADED) {
            bcf_update_info_float(hdr, rec, "FRAC_A", &locus_entry->frac_a, 1);
            bcf_update_info_float(hdr, rec, "FRAC_C", &locus_entry->frac_c, 1);
            bcf_update_info_float(hdr, rec, "FRAC_G", &locus_entry->frac_g, 1);
            bcf_update_info_float(hdr, rec, "FRAC_T", &locus_entry->frac_t, 1);
            bcf_update_info_int32(hdr, rec, "NORM_ID", &locus_entry->norm_id, 1);
        }
        if (flags & CSV_LOADED) {
            bcf_update_info_int32(hdr, rec, "BEADSET_ID", &locus_entry->beadset_id, 1);
        }
        if ((flags & BPM_LOADED) | (flags & CSV_LOADED)) {
            int32_t assay_type = (int32_t)locus_entry->assay_type;
            bcf_update_info_int32(hdr, rec, "ASSAY_TYPE", &assay_type, 1);
        }
        if (flags & EGT_LOADED) {
            if (flags & ADJUST_CLUSTERS) {
                adjust_clusters(gts, ilmn_theta_arr, ilmn_r_arr, n, &egt->cluster_records[j]);
                for (int i = 0; i < n; i++) {
                    get_baf_lrr(ilmn_theta_arr[i], ilmn_r_arr[i], egt->cluster_records[j].aa_cluster_stats.theta_mean,
                                egt->cluster_records[j].ab_cluster_stats.theta_mean,
                                egt->cluster_records[j].bb_cluster_stats.theta_mean,
                                egt->cluster_records[j].aa_cluster_stats.r_mean,
                                egt->cluster_records[j].ab_cluster_stats.r_mean,
                                egt->cluster_records[j].bb_cluster_stats.r_mean, &baf_arr[i], &lrr_arr[i]);
                }
            }
            bcf_update_info_float(hdr, rec, "GenTrain_Score", &egt->cluster_records[j].cluster_score.total_score, 1);
            bcf_update_info_float(hdr, rec, "Orig_Score", &egt->cluster_records[j].cluster_score.original_score, 1);
            if (egt->cluster_records[j].cluster_score.edited) bcf_update_info_flag(hdr, rec, "Edited", NULL, 1);
            bcf_update_info_float(hdr, rec, "Cluster_Sep", &egt->cluster_records[j].cluster_score.cluster_separation,
                                  1);
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
            bcf_update_info_float(hdr, rec, "Intensity_Threshold", &egt->cluster_records[j].intensity_threshold, 1);
        }

        gts_to_gt_arr(gt_arr, gts, n, allele_a_idx, allele_b_idx);
        bcf_update_genotypes(hdr, rec, gt_arr, n * 2);
        bcf_update_format_int32(hdr, rec, "GQ", gq_arr, n);
        bcf_update_format_float(hdr, rec, "IGC", igc_arr, n);
        bcf_update_format_float(hdr, rec, "BAF", baf_arr, n);
        bcf_update_format_float(hdr, rec, "LRR", lrr_arr, n);
        bcf_update_format_float(hdr, rec, "NORMX", norm_x_arr, n);
        bcf_update_format_float(hdr, rec, "NORMY", norm_y_arr, n);
        bcf_update_format_float(hdr, rec, "R", ilmn_r_arr, n);
        bcf_update_format_float(hdr, rec, "THETA", ilmn_theta_arr, n);
        bcf_update_format_int32(hdr, rec, "X", raw_x_arr, n);
        bcf_update_format_int32(hdr, rec, "Y", raw_y_arr, n);
        if (bcf_write(out_fh, hdr, rec) < 0) error("Unable to write to output VCF file\n");
    }
    fprintf(stderr, "Lines   total/missing-reference/skipped:\t%d/%d/%d\n", bpm->num_loci, n_missing, n_skipped);

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
    free(flank.s);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    if (hts_close(out_fh) < 0) error("Close failed: %s\n", out_fh->fn);
}

#define GS_GT 0
#define GS_TOP_STRAND 1
#define GS_REF_STRAND 2
#define GS_IGC 3
#define GS_BAF 4
#define GS_LRR 5
#define GS_NORMX 6
#define GS_NORMY 7
#define GS_R 8
#define GS_THETA 9
#define GS_X 10
#define GS_Y 11

typedef struct {
    int *col2sample;
    int type;
    void *ptr;
} gs_col_t;

static int tsv_setter_gs_col(tsv_t *tsv, bcf1_t *rec, void *usr) {
    gs_col_t *gs_col = (gs_col_t *)usr;
    uint8_t *gts;
    char *strand_alleles, *endptr;
    switch (gs_col->type) {
    case GS_GT:
        gts = (uint8_t *)gs_col->ptr + gs_col->col2sample[tsv->icol];
        if ((tsv->ss[0] == 'N' && tsv->ss[1] == 'C') || (tsv->ss[0] == '-' && tsv->ss[1] == '-'))
            *gts = GT_NC;
        else if (tsv->ss[0] == 'A' && tsv->ss[1] == 'A')
            *gts = GT_AA;
        else if (tsv->ss[0] == 'A' && tsv->ss[1] == 'B')
            *gts = GT_AB;
        else if (tsv->ss[0] == 'B' && tsv->ss[1] == 'B')
            *gts = GT_BB;
        else
            return -1;
        break;
    case GS_TOP_STRAND:
    case GS_REF_STRAND:
        strand_alleles = (char *)gs_col->ptr + 2 * gs_col->col2sample[tsv->icol];
        strand_alleles[0] = tsv->ss[0];
        strand_alleles[1] = tsv->ss[1];
        break;
    case GS_IGC:
    case GS_BAF:
    case GS_LRR:
    case GS_NORMX:
    case GS_NORMY:
    case GS_R:
    case GS_THETA:
        ((float *)gs_col->ptr + gs_col->col2sample[tsv->icol])[0] = strtof(tsv->ss, &endptr);
        if (tsv->ss == endptr) return -1;
        break;
    case GS_X:
    case GS_Y:
        ((int32_t *)gs_col->ptr + gs_col->col2sample[tsv->icol])[0] = strtol(tsv->ss, &endptr, 0);
        if (tsv->ss == endptr) return -1;
        break;
    default:
        return -1;
    }
    return 0;
}

static int tsv_setter_chrom_flexible(tsv_t *tsv, bcf1_t *rec, void *usr) {
    char tmp = *tsv->se;
    *tsv->se = 0;
    rec->rid = bcf_hdr_name2id_flexible((bcf_hdr_t *)usr, tsv->ss);
    *tsv->se = tmp;
    return rec->rid == -1 ? -1 : 0;
}

static int tsv_register_all(tsv_t *tsv, const char *id, tsv_setter_t setter, void *usr) {
    int i, n = 0;
    for (i = 0; i < tsv->ncols; i++) {
        if (!tsv->cols[i].name || strcasecmp(tsv->cols[i].name, id)) continue;
        tsv->cols[i].setter = setter;
        tsv->cols[i].usr = usr;
        n++;
    }
    return n ? 0 : -1;
}

static int tsv_parse_delimiter(tsv_t *tsv, bcf1_t *rec, char *str, int delimiter) {
    int status = 0;
    tsv->icol = 0;
    tsv->ss = tsv->se = str;
    while (*tsv->ss && tsv->icol < tsv->ncols) {
        if (delimiter)
            while (*tsv->se && (*tsv->se) != delimiter) tsv->se++;
        else
            while (*tsv->se && !isspace(*tsv->se)) tsv->se++;
        if (tsv->cols[tsv->icol].setter) {
            int ret = tsv->cols[tsv->icol].setter(tsv, rec, tsv->cols[tsv->icol].usr);
            if (ret < 0) return -1;
            status++;
        }
        if (delimiter)
            tsv->se++;
        else
            while (*tsv->se && isspace(*tsv->se)) tsv->se++;
        tsv->ss = tsv->se;
        tsv->icol++;
    }
    return status ? 0 : -1;
}

static void gs_to_vcf(faidx_t *fai, htsFile *gs_fh, htsFile *out_fh, bcf_hdr_t *hdr, int flags, int gc_win) {
    // read the header of the table
    kstring_t line = {0, 0, NULL};
    if (hts_getline(gs_fh, KS_SEP_LINE, &line) <= 0) error("Empty file: %s\n", gs_fh->fn);
    int moff = 0, *off = NULL, ncols = ksplit_core(line.s, '\t', &moff, &off);
    kstring_t str = {0, 0, NULL};
    int *col2sample = (int *)malloc(sizeof(int) * ncols);
    for (int i = 0; i < ncols; i++) {
        char *ptr;
        if (i > 0) kputc(',', &str);
        if ((ptr = strrchr(&line.s[off[i]], '.'))) {
            *ptr++ = '\0';
            if ((bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, &line.s[off[i]]) < 0)) bcf_hdr_add_sample(hdr, &line.s[off[i]]);
            if (strcmp(ptr, "GType") == 0)
                kputs("GT", &str);
            else if (strcmp(ptr, "Score") == 0 || strcmp(ptr, "GC Score") == 0)
                kputs("IGC", &str);
            else if (strcmp(ptr, "Theta") == 0 || strcmp(ptr, "Theta Illumina") == 0)
                kputs("THETA", &str);
            else if (strcmp(ptr, "R") == 0 || strcmp(ptr, "R Illumina") == 0)
                kputc('R', &str);
            else if (strcmp(ptr, "X Raw") == 0 || strcmp(ptr, "Raw X") == 0)
                kputc('X', &str);
            else if (strcmp(ptr, "Y Raw") == 0 || strcmp(ptr, "Raw Y") == 0)
                kputc('Y', &str);
            else if (strcmp(ptr, "X") == 0)
                kputs("NORMX", &str);
            else if (strcmp(ptr, "Y") == 0)
                kputs("NORMY", &str);
            else if (strcmp(ptr, "B Allele Freq") == 0)
                kputs("BAF", &str);
            else if (strcmp(ptr, "Log R Ratio") == 0)
                kputs("LRR", &str);
            else if (strcmp(ptr, "Top Alleles") == 0)
                kputs("TOP_STRAND", &str);
            else if (strcmp(ptr, "Plus/Minus Alleles") == 0)
                kputs("REF_STRAND", &str);
            else if (strcmp(ptr, "Import Calls") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "Concordance") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "Orig Call") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "CNV Value") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "CNV Confidence") == 0)
                kputc('-', &str);
            else
                error("Could not recognize FORMAT field: %s\n", ptr);
            col2sample[i] = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, &line.s[off[i]]);
        } else {
            ptr = &line.s[off[i]];
            if (strcmp(ptr, "Index") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "Name") == 0 || strcmp(ptr, "SNP Name") == 0)
                kputs("ID", &str);
            else if (strcmp(ptr, "Address") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "Chr") == 0 || strcmp(ptr, "Chromosome") == 0)
                kputs("CHROM", &str);
            else if (strcmp(ptr, "Manifest") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "Position") == 0)
                kputs("POS", &str);
            else if (strcmp(ptr, "GenTrain Score") == 0)
                kputs("GENTRAIN_SCORE", &str);
            else if (strcmp(ptr, "Frac A") == 0)
                kputs("FRAC_A", &str);
            else if (strcmp(ptr, "Frac C") == 0)
                kputs("FRAC_C", &str);
            else if (strcmp(ptr, "Frac G") == 0)
                kputs("FRAC_G", &str);
            else if (strcmp(ptr, "Frac T") == 0)
                kputs("FRAC_T", &str);
            else if (strcmp(ptr, "IlmnStrand") == 0)
                kputc('-', &str);
            else if (strcmp(ptr, "SNP") == 0)
                kputc('-', &str);
            else
                error("Could not recognize INFO field: %s\n", ptr);
            col2sample[i] = -1;
        }
    }
    free(off);
    if (bcf_hdr_sync(hdr) < 0) error_errno("[%s] Failed to update header",
                                           __func__); // updates the number of samples
    int nsamples = bcf_hdr_nsamples(hdr);

    tsv_t *tsv = tsv_init(str.s);
    if (tsv_register(tsv, "CHROM", tsv_setter_chrom_flexible, hdr) < 0) error("Expected CHROM column\n");
    if (tsv_register(tsv, "POS", tsv_setter_pos, NULL) < 0) error("Expected POS column\n");
    tsv_register(tsv, "ID", tsv_setter_id, hdr);

    float total_score;
    int gentrain_score = tsv_register(tsv, "GENTRAIN_SCORE", tsv_read_float, &total_score);
    if (gentrain_score)
        bcf_hdr_append(hdr,
                       "##INFO=<ID=GenTrain_Score,Number=1,Type=Float,Description=\"The SNP "
                       "cluster quality from the GenTrain clustering algorithm\">");
    float frac[4];
    int frac_a = tsv_register(tsv, "FRAC_A", tsv_read_float, &frac[0]);
    if (frac_a == 0)
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_A,Number=1,Type=Float,Description=\"Fraction of the A "
                       "nucleotide in the top genomic sequence\">");
    int frac_c = tsv_register(tsv, "FRAC_C", tsv_read_float, &frac[1]);
    if (frac_c == 0)
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_C,Number=1,Type=Float,Description=\"Fraction of the C "
                       "nucleotide in the top genomic sequence\">");
    int frac_g = tsv_register(tsv, "FRAC_G", tsv_read_float, &frac[2]);
    if (frac_g == 0)
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_G,Number=1,Type=Float,Description=\"Fraction of the G "
                       "nucleotide in the top genomic sequence\">");
    int frac_t = tsv_register(tsv, "FRAC_T", tsv_read_float, &frac[3]);
    if (frac_t == 0)
        bcf_hdr_append(hdr,
                       "##INFO=<ID=FRAC_T,Number=1,Type=Float,Description=\"Fraction of the T "
                       "nucleotide in the top genomic sequence\">");

    uint8_t *gts = (uint8_t *)malloc(nsamples * sizeof(uint8_t));
    int32_t *gt_arr = (int32_t *)malloc(nsamples * 2 * sizeof(int32_t));
    gs_col_t gs_gts = {col2sample, GS_GT, gts};
    if (tsv_register_all(tsv, "GT", tsv_setter_gs_col, &gs_gts) < 0) error("Expected GType column\n");
    int32_t *gq_arr = (int32_t *)malloc(nsamples * sizeof(int32_t));

    char *top_strand_alleles = (char *)malloc(nsamples * 2 * sizeof(char));
    gs_col_t gs_top_strand = {col2sample, GS_TOP_STRAND, top_strand_alleles};
    tsv_register_all(tsv, "TOP_STRAND", tsv_setter_gs_col, &gs_top_strand);
    const char *strand_alleles = top_strand_alleles;

    char *ref_strand_alleles = (char *)malloc(nsamples * 2 * sizeof(char));
    gs_col_t gs_ref_strand = {col2sample, GS_REF_STRAND, ref_strand_alleles};
    if (tsv_register_all(tsv, "REF_STRAND", tsv_setter_gs_col, &gs_ref_strand) >= 0)
        strand_alleles = ref_strand_alleles;

    gs_col_t gs_igc = {col2sample, GS_IGC, NULL};
    if ((flags & FORMAT_IGC) && tsv_register_all(tsv, "IGC", tsv_setter_gs_col, &gs_igc) == 0)
        gs_igc.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_baf = {col2sample, GS_BAF, NULL};
    if ((flags & FORMAT_BAF) && tsv_register_all(tsv, "BAF", tsv_setter_gs_col, &gs_baf) == 0)
        gs_baf.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_lrr = {col2sample, GS_LRR, NULL};
    if ((flags & FORMAT_LRR) && tsv_register_all(tsv, "LRR", tsv_setter_gs_col, &gs_lrr) == 0)
        gs_lrr.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_norm_x = {col2sample, GS_NORMX, NULL};
    if ((flags & FORMAT_NORMX) && tsv_register_all(tsv, "NORMX", tsv_setter_gs_col, &gs_norm_x) == 0)
        gs_norm_x.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_norm_y = {col2sample, GS_NORMY, NULL};
    if ((flags & FORMAT_NORMY) && tsv_register_all(tsv, "NORMY", tsv_setter_gs_col, &gs_norm_y) == 0)
        gs_norm_y.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_ilmn_r = {col2sample, GS_R, NULL};
    if ((flags & FORMAT_R) && tsv_register_all(tsv, "R", tsv_setter_gs_col, &gs_ilmn_r) == 0)
        gs_ilmn_r.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_ilmn_theta = {col2sample, GS_THETA, NULL};
    if ((flags & FORMAT_THETA) && tsv_register_all(tsv, "THETA", tsv_setter_gs_col, &gs_ilmn_theta) == 0)
        gs_ilmn_theta.ptr = malloc(nsamples * sizeof(float));

    gs_col_t gs_raw_x = {col2sample, GS_X, NULL};
    if ((flags & FORMAT_X) && tsv_register_all(tsv, "X", tsv_setter_gs_col, &gs_raw_x) == 0)
        gs_raw_x.ptr = malloc(nsamples * sizeof(int32_t));

    gs_col_t gs_raw_y = {col2sample, GS_Y, NULL};
    if ((flags & FORMAT_Y) && tsv_register_all(tsv, "Y", tsv_setter_gs_col, &gs_raw_y) == 0)
        gs_raw_y.ptr = malloc(nsamples * sizeof(int32_t));

    if (bcf_hdr_write(out_fh, hdr) < 0) error("Unable to write to output VCF file\n");

    bcf1_t *rec = bcf_init();
    char ref_base[] = {'\0', '\0'};
    char allele_a[] = {'\0', '\0'};
    char allele_b[] = {'\0', '\0'};
    int32_t allele_a_idx, allele_b_idx;
    int n_total = 0, n_missing = 0, n_skipped = 0;
    while (hts_getline(gs_fh, KS_SEP_LINE, &line) > 0) {
        if (line.s[0] == '#') continue; // skip comments
        bcf_clear(rec);
        rec->n_sample = nsamples;

        n_total++;
        if (!tsv_parse_delimiter(tsv, rec, line.s, '\t')) {
            if (rec->rid < 0 || rec->pos < 0) {
                if (flags & VERBOSE) fprintf(stderr, "Skipping unlocalized marker %s\n", rec->d.id);
                n_skipped++;
                continue;
            }
            // determine A and B alleles
            allele_a[0] = '.';
            allele_b[0] = '.';
            for (int i = 0; i < nsamples; i++) {
                switch (gts[i]) {
                case GT_NC:
                    break;
                case GT_AA:
                    allele_a[0] = strand_alleles[2 * i];
                    break;
                case GT_AB:
                    allele_a[0] = strand_alleles[2 * i];
                    allele_b[0] = strand_alleles[2 * i + 1];
                    break;
                case GT_BB:
                    allele_b[0] = strand_alleles[2 * i];
                    break;
                default:
                    error("Unable to process marker %s\n", rec->d.id);
                    break;
                }
            }

            int len, win = min(max(100, gc_win), rec->pos);
            char *ref = faidx_fetch_seq(fai, bcf_seqname(hdr, rec), rec->pos - win, rec->pos + win, &len);
            if (!ref || len == 1)
                error("faidx_fetch_seq failed at %s:%" PRId64 " (are you using the correct reference genome?)\n",
                      bcf_seqname(hdr, rec), rec->pos + 1);
            strupper(ref);
            if (!(flags & NO_INFO_GC)) {
                float gc_ratio = get_gc_ratio(&ref[max(win - gc_win, 0)], &ref[min(win + gc_win, len)]);
                bcf_update_info_float(hdr, rec, "GC", &gc_ratio, 1);
            }
            ref_base[0] = ref[win];
            int is_indel = allele_a[0] == 'D' || allele_a[0] == 'I' || allele_b[0] == 'D' || allele_b[0] == 'I';
            if (is_indel) {
                if (allele_a[0] == '.') {
                    allele_a[0] = allele_b[0] == 'D' ? 'I' : 'D';
                }
                if (allele_b[0] == '.') {
                    allele_b[0] = allele_a[0] == 'D' ? 'I' : 'D';
                }
                ref_base[0] = allele_a[0];
                n_missing++;
            } else if (strand_alleles == top_strand_alleles) {
                if (allele_a[0] == '.' || allele_b[0] == '.') {
                    allele_a[0] = '.';
                    allele_b[0] = '.';
                } else {
                    int strand = get_strand_from_top_alleles(allele_a, allele_b, ref, win, len);
                    if (strand < 0) {
                        if (flags & VERBOSE)
                            fprintf(stderr, "Unable to determine reference strand for SNP %s\n", rec->d.id);
                        allele_a[0] = '.';
                        allele_b[0] = '.';
                    } else if (strand == 1) {
                        allele_a[0] = rev_nt(allele_a[0]);
                        allele_b[0] = rev_nt(allele_b[0]);
                    }
                }
            }
            free(ref);

            allele_b_idx = get_allele_b_idx(ref_base[0], allele_a, allele_b);
            allele_a_idx = get_allele_a_idx(allele_b_idx);
            const char *alleles[3];
            int nals = alleles_ab_to_vcf(alleles, ref_base, allele_a, allele_b, allele_b_idx);
            if (nals < 0) error("Unable to process marker %s\n", rec->d.id);
            bcf_update_alleles(hdr, rec, alleles, nals);
            bcf_update_info_int32(hdr, rec, "ALLELE_A", &allele_a_idx, 1);
            bcf_update_info_int32(hdr, rec, "ALLELE_B", &allele_b_idx, 1);

            if (allele_a_idx >= 0 && allele_b_idx >= 0) {
                gts_to_gt_arr(gt_arr, gts, nsamples, allele_a_idx, allele_b_idx);
            } else {
                for (int i = 0; i < nsamples; i++) {
                    gt_arr[2 * i] = bcf_gt_missing;
                    gt_arr[2 * i + 1] = bcf_gt_missing;
                }
            }
            if (gentrain_score == 0) bcf_update_info_float(hdr, rec, "GenTrain_Score", &total_score, 1);
            if (frac_a == 0) bcf_update_info_float(hdr, rec, "FRAC_A", &frac[0], 1);
            if (frac_c == 0) bcf_update_info_float(hdr, rec, "FRAC_C", &frac[1], 1);
            if (frac_g == 0) bcf_update_info_float(hdr, rec, "FRAC_G", &frac[2], 1);
            if (frac_t == 0) bcf_update_info_float(hdr, rec, "FRAC_T", &frac[3], 1);

            bcf_update_genotypes(hdr, rec, gt_arr, nsamples * 2);

            if (gs_igc.ptr) {
                for (int i = 0; i < nsamples; i++) {
                    gq_arr[i] = (int)(-10 * log10(1 - ((float *)gs_igc.ptr)[i]) + .5);
                    if (gq_arr[i] < 0) gq_arr[i] = 0;
                    if (gq_arr[i] > 50) gq_arr[i] = 50;
                }
                bcf_update_format_int32(hdr, rec, "GQ", gq_arr, nsamples);
                bcf_update_format_float(hdr, rec, "IGC", (float *)gs_igc.ptr, nsamples);
            }
            if (gs_baf.ptr) bcf_update_format_float(hdr, rec, "BAF", (float *)gs_baf.ptr, nsamples);
            if (gs_lrr.ptr) bcf_update_format_float(hdr, rec, "LRR", (float *)gs_lrr.ptr, nsamples);
            if (gs_norm_x.ptr) bcf_update_format_float(hdr, rec, "NORMX", (float *)gs_norm_x.ptr, nsamples);
            if (gs_norm_y.ptr) bcf_update_format_float(hdr, rec, "NORMY", (float *)gs_norm_y.ptr, nsamples);
            if (gs_ilmn_r.ptr) bcf_update_format_float(hdr, rec, "R", (float *)gs_ilmn_r.ptr, nsamples);
            if (gs_ilmn_theta.ptr) bcf_update_format_float(hdr, rec, "THETA", (float *)gs_ilmn_theta.ptr, nsamples);
            if (gs_raw_x.ptr) bcf_update_format_int32(hdr, rec, "X", (int32_t *)gs_raw_x.ptr, nsamples);
            if (gs_raw_y.ptr) bcf_update_format_int32(hdr, rec, "Y", (int32_t *)gs_raw_y.ptr, nsamples);
            if (bcf_write(out_fh, hdr, rec) < 0) error("Unable to write to output VCF file\n");
        } else {
            if (flags & VERBOSE) fprintf(stderr, "Failed to process marker %s\n", rec->d.id);
            n_skipped++;
        }
    }
    fprintf(stderr, "Lines   total/missing-reference/skipped:\t%d/%d/%d\n", n_total, n_missing, n_skipped);
    free(line.s);

    free(col2sample);
    free(gts);
    free(gt_arr);
    free(gq_arr);
    free(top_strand_alleles);
    free(ref_strand_alleles);
    free(gs_igc.ptr);
    free(gs_baf.ptr);
    free(gs_lrr.ptr);
    free(gs_norm_x.ptr);
    free(gs_norm_y.ptr);
    free(gs_ilmn_r.ptr);
    free(gs_ilmn_theta.ptr);
    free(gs_raw_x.ptr);
    free(gs_raw_y.ptr);
    tsv_destroy(tsv);
    free(str.s);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    if (hts_close(out_fh) < 0) error("Close failed: %s\n", out_fh->fn);
    if (hts_close(gs_fh) < 0) error("Close failed: %s\n", gs_fh->fn);
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Convert Illumina GTC files to VCF.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: convert Illumina GTC files containing intensity data into VCF. "
           "(version " GTC2VCF_VERSION
           " https://github.com/freeseek/gtc2vcf)\n"
           "Usage: bcftools +gtc2vcf [options] [<A.gtc> ...]\n"
           "\n"
           "Plugin options:\n"
           "    -l, --list-tags                 list available FORMAT tags with description for VCF output\n"
           "    -t, --tags LIST                 list of output FORMAT tags [" TAG_LIST_DFLT
           "]\n"
           "    -b, --bpm <file>                BPM manifest file\n"
           "    -c, --csv <file>                CSV manifest file (can be gzip compressed)\n"
           "    -e, --egt <file>                EGT cluster file\n"
           "    -f, --fasta-ref <file>          reference sequence in fasta format\n"
           "        --set-cache-size <int>      select fasta cache size in bytes\n"
           "        --gc-window-size <int>      window size in bp used to compute the GC content (-1 for no estimate) "
           "[" GC_WIN_DFLT
           "]\n"
           "    -g, --gtcs <dir|file>           GTC genotype files from directory or list from file\n"
           "    -i, --idat                      input IDAT files rather than GTC files\n"
           "        --capacity <int>            number of variants to read from intensity files per I/O operation "
           "[" CAPACITY_DFLT
           "]\n"
           "        --adjust-clusters           adjust cluster centers in (Theta, R) space (requires --bpm and --egt)\n"
           "        --use-gtc-sample-names      use sample name in GTC files rather than GTC file name\n"
           "        --do-not-check-bpm          do not check whether BPM and GTC files match manifest file name\n"
           "        --genome-studio <file>      input a GenomeStudio final report file (in matrix format)\n"
           "        --no-version                do not append version and command line to the header\n"
           "    -o, --output <file>             write output to a file [standard output]\n"
           "    -O, --output-type <b|u|z|v|t>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF\n"
           "                                    v: uncompressed VCF, t: GenomeStudio tab-delimited text output [v]\n"
           "        --threads <int>             number of extra output compression threads [0]\n"
           "    -x, --extra <file>              write GTC metadata to a file\n"
           "    -v, --verbose                   print verbose information\n"
           "\n"
           "Manifest options:\n"
           "        --beadset-order             output BeadSetID normalization order (requires --bpm and --csv)\n"
           "        --fasta-flank               output flank sequence in FASTA format (requires --csv)\n"
           "    -s, --sam-flank <file>          input flank sequence alignment in SAM/BAM format (requires --csv)\n"
           "        --genome-build <assembly>   genome build ID used to update the manifest file [" GENOME_BUILD_DFLT
           "]\n"
           "\n"
           "Examples:\n"
           "    bcftools +gtc2vcf -i 5434246082_R03C01_Grn.idat\n"
           "    bcftools +gtc2vcf 5434246082_R03C01.gtc\n"
           "    bcftools +gtc2vcf -b HumanOmni2.5-4v1_H.bpm -c HumanOmni2.5-4v1_H.csv\n"
           "    bcftools +gtc2vcf -e HumanOmni2.5-4v1_H.egt\n"
           "    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv -e GSA-24v3-0_A1_ClusterFile.egt -f human_g1k_v37.fasta -o "
           "GSA-24v3-0_A1.vcf\n"
           "    bcftools +gtc2vcf -c HumanOmni2.5-4v1_H.csv -f human_g1k_v37.fasta 5434246082_R03C01.gtc -o "
           "5434246082_R03C01.vcf\n"
           "    bcftools +gtc2vcf -f human_g1k_v37.fasta --genome-studio GenotypeReport.txt -o GenotypeReport.vcf\n"
           "\n"
           "Examples of manifest file options:\n"
           "    bcftools +gtc2vcf -b GSA-24v3-0_A1.bpm -c GSA-24v3-0_A1.csv --beadset-order\n"
           "    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv --fasta-source-seq -o GSA-24v3-0_A1.fasta\n"
           "    bwa mem -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GSA-24v3-0_A1.fasta -o GSA-24v3-0_A1.sam\n"
           "    bcftools +gtc2vcf -c GSA-24v3-0_A1.csv --sam-source-seq GSA-24v3-0_A1.sam -o GSA-24v3-0_A1.GRCh38.csv\n"
           "\n";
}

static int parse_tags(const char *str) {
    int flags = 0, n;
    char **tags = hts_readlist(str, 0, &n);
    for (int i = 0; i < n; i++) {
        if (!strcasecmp(tags[i], "GT"))
            flags |= FORMAT_GT;
        else if (!strcasecmp(tags[i], "GQ"))
            flags |= FORMAT_GQ;
        else if (!strcasecmp(tags[i], "IGC"))
            flags |= FORMAT_IGC;
        else if (!strcasecmp(tags[i], "X"))
            flags |= FORMAT_X;
        else if (!strcasecmp(tags[i], "Y"))
            flags |= FORMAT_Y;
        else if (!strcasecmp(tags[i], "NORMX"))
            flags |= FORMAT_NORMX;
        else if (!strcasecmp(tags[i], "NORMY"))
            flags |= FORMAT_NORMY;
        else if (!strcasecmp(tags[i], "R"))
            flags |= FORMAT_R;
        else if (!strcasecmp(tags[i], "THETA"))
            flags |= FORMAT_THETA;
        else if (!strcasecmp(tags[i], "LRR"))
            flags |= FORMAT_LRR;
        else if (!strcasecmp(tags[i], "BAF"))
            flags |= FORMAT_BAF;
        else
            error("Error parsing \"--tags %s\": the tag \"%s\" is not supported\n", str, tags[i]);
        free(tags[i]);
    }
    if (n) free(tags);
    return flags;
}

static void list_tags(void) {
    error(
        "FORMAT/GT       Number:1  Type:String   ..  Genotype\n"
        "FORMAT/GQ       Number:1  Type:Integer  ..  Genotype Quality\n"
        "FORMAT/IGC      Number:1  Type:Float    ..  Illumina GenCall Confidence Score\n"
        "FORMAT/BAF      Number:1  Type:Float    ..  B Allele Frequency\n"
        "FORMAT/LRR      Number:1  Type:Float    ..  Log R Ratio\n"
        "FORMAT/NORMX    Number:1  Type:Float    ..  Normalized X intensity\n"
        "FORMAT/NORMY    Number:1  Type:Float    ..  Normalized Y intensity\n"
        "FORMAT/R        Number:1  Type:Float    ..  Normalized R value\n"
        "FORMAT/THETA    Number:1  Type:Float    ..  Normalized Theta value\n"
        "FORMAT/X        Number:1  Type:Integer  ..  Raw X intensity\n"
        "FORMAT/Y        Number:1  Type:Integer  ..  Raw Y intensity\n");
}

int run(int argc, char *argv[]) {
    const char *tag_list = TAG_LIST_DFLT;
    const char *bpm_fname = NULL;
    const char *csv_fname = NULL;
    const char *egt_fname = NULL;
    const char *gs_fname = NULL;
    const char *output_fname = "-";
    const char *ref_fname = NULL;
    const char *pathname = NULL;
    const char *extra_fname = NULL;
    const char *sam_fname = NULL;
    const char *genome_build = GENOME_BUILD_DFLT;
    char *tmp;
    int flags = 0;
    int output_type = FT_VCF;
    size_t capacity = 0;
    int cache_size = 0;
    int gc_win = (int)strtol(GC_WIN_DFLT, NULL, 0);
    int gtc_sample_names = 0;
    int bpm_check = 1;
    int n_threads = 0;
    int record_cmd_line = 1;
    int binary_to_csv = 0;
    int beadset_order = 0;
    int fasta_flank = 0;
    faidx_t *fai = NULL;
    htsFile *out_fh = NULL;
    FILE *out_txt = NULL;

    static struct option loptions[] = {{"list-tags", no_argument, NULL, 'l'},
                                       {"tags", required_argument, NULL, 't'},
                                       {"bpm", required_argument, NULL, 'b'},
                                       {"csv", required_argument, NULL, 'c'},
                                       {"egt", required_argument, NULL, 'e'},
                                       {"fasta-ref", required_argument, NULL, 'f'},
                                       {"set-cache-size", required_argument, NULL, 1},
                                       {"gc-window-size", required_argument, NULL, 2},
                                       {"gtcs", required_argument, NULL, 'g'},
                                       {"idat", no_argument, NULL, 'i'},
                                       {"capacity", required_argument, NULL, 3},
                                       {"adjust-clusters", no_argument, NULL, 4},
                                       {"use-gtc-sample-names", no_argument, NULL, 5},
                                       {"do-not-check-bpm", no_argument, NULL, 6},
                                       {"genome-studio", required_argument, NULL, 7},
                                       {"no-version", no_argument, NULL, 8},
                                       {"output", required_argument, NULL, 'o'},
                                       {"output-type", required_argument, NULL, 'O'},
                                       {"threads", required_argument, NULL, 9},
                                       {"extra", required_argument, NULL, 'x'},
                                       {"verbose", no_argument, NULL, 'v'},
                                       {"beadset-order", no_argument, NULL, 10},
                                       {"fasta-flank", no_argument, NULL, 11},
                                       {"sam-flank", required_argument, NULL, 's'},
                                       {"genome-build", required_argument, NULL, 12},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?lt:b:c:e:f:g:io:O:x:vs:", loptions, NULL)) >= 0) {
        switch (c) {
        case 'l':
            list_tags();
            break;
        case 't':
            tag_list = optarg;
            break;
        case 'b':
            bpm_fname = optarg;
            break;
        case 'c':
            csv_fname = optarg;
            break;
        case 'e':
            egt_fname = optarg;
            break;
        case 'f':
            ref_fname = optarg;
            break;
        case 1:
            cache_size = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --set-cache-size %s\n", optarg);
            break;
        case 2:
            gc_win = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --gc-window-size %s\n", optarg);
            if (gc_win <= 0) flags |= NO_INFO_GC;
            break;
        case 'g':
            pathname = optarg;
            break;
        case 'i':
            flags |= LOAD_IDAT;
            break;
        case 3:
            capacity = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --capacity %s\n", optarg);
            break;
        case 4:
            flags |= ADJUST_CLUSTERS;
            break;
        case 5:
            gtc_sample_names = 1;
            break;
        case 6:
            bpm_check = 0;
            break;
        case 7:
            gs_fname = optarg;
            break;
        case 8:
            record_cmd_line = 0;
            break;
        case 'o':
            output_fname = optarg;
            break;
        case 'O':
            switch (optarg[0]) {
            case 'b':
                output_type = FT_BCF_GZ;
                break;
            case 'u':
                output_type = FT_BCF;
                break;
            case 'z':
                output_type = FT_VCF_GZ;
                break;
            case 'v':
                output_type = FT_VCF;
                break;
            case 't':
                output_type = FT_TAB_TEXT;
                break;
            default:
                error("The output type \"%s\" not recognised\n", optarg);
            }
            break;
        case 9:
            n_threads = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse argument: --threads %s\n", optarg);
            break;
        case 'x':
            extra_fname = optarg;
            break;
        case 'v':
            flags |= VERBOSE;
            break;
        case 10:
            beadset_order = 1;
            break;
        case 11:
            fasta_flank = 1;
            break;
        case 's':
            sam_fname = optarg;
            break;
        case 12:
            genome_build = optarg;
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
        }
    }
    if (((bpm_fname != NULL) || (csv_fname != NULL)) + (egt_fname != NULL) + (ref_fname != NULL) + (gs_fname != NULL)
            + (argc - optind > 0) + (pathname != NULL)
        == 1)
        binary_to_csv = 1;
    if (sam_fname && (csv_fname == NULL)) error("The --sam-flank option requires the --csv option\n%s", usage_text());
    if (binary_to_csv) {
        if (beadset_order && (bpm_fname == NULL || csv_fname == NULL))
            error("The --beadset-order option requires both the --bpm and the --csv options\n%s", usage_text());
        if (fasta_flank && (csv_fname == NULL))
            error("The --fasta-flank option requires the --csv option\n%s", usage_text());
        if (beadset_order + fasta_flank + (sam_fname != NULL) > 1)
            error(
                "Only one of --beadset-order or --fasta-flank or --sam-flank options can be "
                "used at once\n%s",
                usage_text());
    } else {
        if (beadset_order)
            error("The --beadset-order option can only be used with options --bpm and --csv\n%s", usage_text());
        if (fasta_flank)
            error("The --fasta-flank option can only be used with options --bpm and --csv\n%s", usage_text());
        if (!bpm_fname && !csv_fname && !gs_fname)
            error("Manifest file required when converting to VCF\n%s", usage_text());
        if (!egt_fname && (flags & ADJUST_CLUSTERS))
            error("Cluster file required when adjusting cluster centers\n%s", usage_text());
        if (gs_fname && (bpm_fname || csv_fname || egt_fname || sam_fname))
            error(
                "If a GenomeStudio final report file is provided, you cannot use --bpm/--csv/--egt/--sam-flank "
                "options\n%s",
                usage_text());
        if (gs_fname && (argc - optind > 0 || pathname))
            error("If a GenomeStudio final report file is provided, do not pass GTC files\n%s", usage_text());
        if (gs_fname && output_type == FT_TAB_TEXT)
            error("If a GenomeStudio final report file is provided, you cannot output in GenomeStudio format\n%s",
                  usage_text());
        if (argc - optind > 0 && pathname)
            error("GTC files cannot be listed through both command interface and file list\n%s", usage_text());
        if (!gs_fname && output_type != FT_TAB_TEXT && extra_fname) out_txt = get_file_handle(extra_fname);
    }
    flags |= parse_tags(tag_list);

    // beginning of plugin run
    fprintf(stderr, "gtc2vcf " GTC2VCF_VERSION " https://github.com/freeseek/gtc2vcf\n");

    int nfiles = 0;
    char **filenames = NULL;
    if (pathname) {
        filenames = get_file_list(pathname, flags & LOAD_IDAT ? "idat" : "gtc", &nfiles);
    } else {
        nfiles = argc - optind;
        filenames = argv + optind;
    }
    void **files = (void **)malloc(nfiles * sizeof(void *));

    // make sure the process is allowed to open enough files
    struct rlimit lim;
    getrlimit(RLIMIT_NOFILE, &lim);
    if (nfiles + 10 > lim.rlim_max)
        error("On this system you cannot open more than %ld files at once while %d is required\n", lim.rlim_max,
              nfiles + 10);
    if (nfiles + 10 > lim.rlim_cur) {
        lim.rlim_cur = nfiles + 10;
        fprintf(stderr, "Adjusting the limit of how many files can be open at once to %ld\n", lim.rlim_cur);
        setrlimit(RLIMIT_NOFILE, &lim);
    }

    if ((flags & ADJUST_CLUSTERS) && nfiles < 100)
        fprintf(stderr, "Warning: adjusting clusters with %d sample(s) is not recommended\n", nfiles);

    if (binary_to_csv || output_type == FT_TAB_TEXT) {
        out_txt = get_file_handle(output_fname);
    } else {
        out_fh = hts_open(output_fname, hts_bcf_wmode(output_type));
        if (out_fh == NULL) error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
        if (n_threads) hts_set_threads(out_fh, n_threads);
        if (!ref_fname) error("VCF output requires the --fasta-ref option\n");
        fai = fai_load(ref_fname);
        if (!fai) error("Could not load the reference %s\n", ref_fname);
        if (cache_size) fai_set_cache_size(fai, cache_size);
        if (extra_fname) out_txt = get_file_handle(extra_fname);
    }

    bpm_t *bpm = NULL;
    if (bpm_fname) {
        fprintf(stderr, "Reading BPM file %s\n", bpm_fname);
        bpm = bpm_init(bpm_fname);
        flags |= BPM_LOADED;
        if (binary_to_csv && !csv_fname) bpm_to_csv(bpm, out_txt, flags);
    }

    if (csv_fname) {
        fprintf(stderr, "Reading CSV file %s\n", csv_fname);
        bpm = bpm_csv_init(csv_fname, bpm);
        flags |= CSV_LOADED;
        if (binary_to_csv && !sam_fname && !beadset_order && !fasta_flank) bpm_to_csv(bpm, out_txt, flags);
    }

    // output source sequences in FASTA format to be realigned by bwa mem
    if (fasta_flank) {
        for (int i = 0; i < bpm->num_loci; i++)
            flank2fasta(bpm->locus_entries[i].ilmn_id, bpm->locus_entries[i].source_seq, out_txt);
    }

    // input source sequence alignments in SAM format to generate new coordinates for the
    // CSV manifest file
    if (sam_fname) {
        fprintf(stderr, "Reading SAM file %s\n", sam_fname);
        bpm = sam_csv_init(sam_fname, bpm, genome_build, flags);
        if (binary_to_csv) bpm_to_csv(bpm, out_txt, flags);
    }

    // the BeadSet normalization order is the only information in the BPM manifest file
    // missing from the CSV manifest file
    kstring_t str = {0, 0, NULL};
    if ((flags & BPM_LOADED) && (flags & CSV_LOADED)) {
        int32_t norm_id_to_beadset_id[100] = {0};
        for (int i = 0; i < bpm->num_loci; i++) {
            uint8_t norm_id = bpm->norm_ids[i];
            if (norm_id_to_beadset_id[norm_id] != 0
                && norm_id_to_beadset_id[norm_id] != bpm->locus_entries[i].beadset_id)
                error("Normalization ID %d corresponds to multiple BeadSet IDs %d and %d\n", norm_id,
                      norm_id_to_beadset_id[norm_id], bpm->locus_entries[i].beadset_id);
            else
                norm_id_to_beadset_id[norm_id] = bpm->locus_entries[i].beadset_id;
        }
        for (int i = 0, j = 0; i < 100; i++) {
            if (norm_id_to_beadset_id[i] == 0) continue;
            if (i != j) error("Normalization ID %d not corresponding to any BeadSet ID", j);
            if (i > 0) kputc(',', &str);
            kputw(norm_id_to_beadset_id[i], &str);
            j++;
        }
        if (beadset_order && out_txt) fprintf(out_txt, "%s,%s\n", bpm->manifest_name, str.s);
    }
    if ((flags & ADJUST_CLUSTERS) && !(flags & BPM_LOADED))
        error("Cannot adjust clusters as couldn't generate the normalization lookup table\n");

    egt_t *egt = NULL;
    if (egt_fname) {
        fprintf(stderr, "Reading EGT file %s\n", egt_fname);
        egt = egt_init(egt_fname);
        if (binary_to_csv)
            egt_to_csv(egt, out_txt, flags & VERBOSE);
        else
            flags |= EGT_LOADED;
    }

    if (gs_fname) flags |= GENOME_STUDIO;

    for (int i = 0; i < nfiles; i++) {
        if (flags & LOAD_IDAT) {
            fprintf(stderr, "Reading IDAT file %s\n", filenames[i]);
            idat_t *idat = idat_init(filenames[i], capacity);
            files[i] = (void *)idat;
        } else {
            fprintf(stderr, "Reading GTC file %s\n", filenames[i]);
            gtc_t *gtc = gtc_init(filenames[i], capacity);
            // GenCall fills the GTC SNP manifest with the BPM file name rather than
            // the BPM manifest name
            if (bpm_check && bpm && strcmp(bpm->manifest_name, gtc->snp_manifest)
                && strcmp(strrchr(bpm->fn, '/') ? strrchr(bpm->fn, '/') + 1 : bpm->fn, gtc->snp_manifest))
                error(
                    "Manifest name %s in BPM file %s does not match manifest name %s in GTC "
                    "file %s\nUse --do-not-check-bpm to suppress this check\n",
                    bpm->manifest_name, bpm->fn, gtc->snp_manifest, gtc->fn);
            files[i] = (void *)gtc;
        }
    }

    if (binary_to_csv && nfiles > 0) {
        if (flags & LOAD_IDAT) {
            if (nfiles == 1)
                idat_to_csv((idat_t *)files[0], out_txt, flags & VERBOSE);
            else
                idats_to_tsv((idat_t **)files, nfiles, out_txt);
        } else {
            if (nfiles == 1)
                gtc_to_csv((gtc_t *)files[0], out_txt, flags & VERBOSE);
            else
                gtcs_to_tsv((gtc_t **)files, nfiles, out_txt);
        }
    }

    if (!binary_to_csv) {
        if (nfiles == 1) fprintf(stderr, "Warning: it is recommended to convert multiple GTC files at once\n");
        if (output_type == FT_TAB_TEXT) {
            fprintf(stderr, "Writing GenomeStudio final report file\n");
            gtcs_to_gs((gtc_t **)files, nfiles, bpm, egt, out_txt, flags);
        } else {
            fprintf(stderr, "Writing VCF file\n");
            bcf_hdr_t *hdr = hdr_init(fai, flags);
            if (bpm_fname)
                bcf_hdr_printf(hdr, "##BPM=%s", strrchr(bpm_fname, '/') ? strrchr(bpm_fname, '/') + 1 : bpm_fname);
            if (csv_fname)
                bcf_hdr_printf(hdr, "##CSV=%s", strrchr(csv_fname, '/') ? strrchr(csv_fname, '/') + 1 : csv_fname);
            if (egt_fname)
                bcf_hdr_printf(hdr, "##EGT=%s", strrchr(egt_fname, '/') ? strrchr(egt_fname, '/') + 1 : egt_fname);
            if (sam_fname)
                bcf_hdr_printf(hdr, "##SAM=%s", strrchr(sam_fname, '/') ? strrchr(sam_fname, '/') + 1 : sam_fname);
            if ((flags & BPM_LOADED) && (flags & CSV_LOADED)) bcf_hdr_printf(hdr, "##BeadSet_Order=%s", str.s);
            if (record_cmd_line) bcf_hdr_append_version(hdr, argc, argv, "bcftools_+gtc2vcf");
            if (gs_fname) {
                htsFile *gs_fh = hts_open(gs_fname, "r");
                bcf_hdr_printf(hdr, "##GenomeStudio=%s",
                               strrchr(gs_fname, '/') ? strrchr(gs_fname, '/') + 1 : gs_fname);
                gs_to_vcf(fai, gs_fh, out_fh, hdr, flags, gc_win);
            } else {
                if (extra_fname) gtcs_to_tsv((gtc_t **)files, nfiles, out_txt);
                for (int i = 0; i < nfiles; i++) {
                    gtc_t *gtc = (gtc_t *)files[i];
                    const char *sample_name =
                        (gtc_sample_names && gtc->sample_name) ? gtc->sample_name : gtc->display_name;
                    if (bcf_hdr_add_sample(hdr, sample_name) < 0)
                        error("GTC files must correspond to different samples\n");
                }
                gtcs_to_vcf(fai, bpm, egt, (gtc_t **)files, nfiles, out_fh, hdr, flags, gc_win);
            }
        }
    }

    free(str.s);
    fai_destroy(fai);
    egt_destroy(egt);
    bpm_destroy(bpm);
    if (pathname) {
        for (int i = 0; i < nfiles; i++) free(filenames[i]);
        free(filenames);
    }
    for (int i = 0; i < nfiles; i++) {
        if (flags & LOAD_IDAT)
            idat_destroy((idat_t *)files[i]);
        else
            gtc_destroy((gtc_t *)files[i]);
    }
    free(files);
    if (out_txt && out_txt != stdout && out_txt != stderr) fclose(out_txt);
    return 0;
}
