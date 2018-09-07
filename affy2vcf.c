/* The MIT License

   Copyright (c) 2018 Giulio Genovese

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
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "htslib/khash_str2int.h"

static inline char revnt(char nt)
{
    if ( nt=='A' ) return 'T';
    if ( nt=='C' ) return 'G';
    if ( nt=='G' ) return 'C';
    if ( nt=='T' ) return 'A';
    return -1;
}

static htsFile *unheader(const char *fn, kstring_t *str)
{
    htsFile *fp = hts_open(fn, "r");
    if ( !fp ) error("Could not read: %s\n", fn);

    if ( hts_getline(fp, KS_SEP_LINE, str) <= 0 ) error("Empty file: %s\n", fn);

    // skip header
    while (str->s[0]=='#') hts_getline(fp, KS_SEP_LINE, str);
    return fp;
}

/****************************************
 * SNP-POSTERIORS FILE IMPLEMENTATION   *
 ****************************************/

typedef struct
{
    float aa_delta;
    float aa_size;
    float ab_delta;
    float ab_size;
    float bb_delta;
    float bb_size;
}
cluster_t;

typedef struct
{
    cluster_t *clusters;
    int n_clusters, m_clusters;
}
snp_posteriors_t;

static snp_posteriors_t *snp_posteriors_init(const char *fn)
{
    kstring_t str = {0, 0, NULL};
    htsFile *fp = unheader(fn, &str);

    int moff = 0, *off = NULL;
    int ncols = ksplit_core(str.s, '\t', &moff, &off);
    if ( strcmp(&str.s[off[0]], "id") ) error("Malformed header in posterior models file: %s\n", fn);

    snp_posteriors_t *snp_posteriors = (snp_posteriors_t *)calloc(1, sizeof(snp_posteriors_t));
    cluster_t *cluster = NULL;
    int moff2 = 0, *off2 = NULL, ncols2;
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        ncols = ksplit_core(str.s, '\t', &moff, &off);
        if ( ncols < 4 ) error("Missing information in posterior models file: %s\n", fn);
        snp_posteriors->n_clusters++;
        hts_expand(cluster_t, snp_posteriors->n_clusters, snp_posteriors->m_clusters, snp_posteriors->clusters);
        cluster = &snp_posteriors->clusters[snp_posteriors->n_clusters - 1];
        cluster->aa_delta = NAN;
        cluster->aa_size = NAN;
        cluster->ab_delta = NAN;
        cluster->ab_size = NAN;
        cluster->bb_delta = NAN;
        cluster->bb_size = NAN;

        char *tmp;
        ncols2 = ksplit_core(&str.s[off[1]], ',', &moff2, &off2);
        if ( ncols2 < 5 ) error("Missing information for cluster BB for probeset %s in file: %s\n", &str.s[off[0]], fn);
        cluster->bb_delta = strtof( &str.s[off[1] + off2[0]], &tmp );
        cluster->bb_size = strtof( &str.s[off[1] + off2[4]], &tmp );
        ncols2 = ksplit_core(&str.s[off[2]], ',', &moff2, &off2);
        if ( ncols2 < 5 ) error("Missing information for cluster AB for probeset %s in file: %s\n", &str.s[off[0]], fn);
        cluster->ab_delta = strtof( &str.s[off[2] + off2[0]], &tmp );
        cluster->ab_size = strtof( &str.s[off[2] + off2[4]], &tmp );
        ncols2 = ksplit_core(&str.s[off[3]], ',', &moff2, &off2);
        if ( ncols2 < 5 ) error("Missing information for cluster AA for probeset %s in file: %s\n", &str.s[off[0]], fn);
        cluster->aa_delta = strtof( &str.s[off[3] + off2[0]], &tmp );
        cluster->aa_size = strtof( &str.s[off[3] + off2[4]], &tmp );
    }

    free(off2);
    free(off);
    free(str.s);
    hts_close(fp);
    return snp_posteriors;
}

static void snp_posteriors_destroy(snp_posteriors_t *snp_posteriors)
{
    free(snp_posteriors->clusters);
}

/****************************************
 * ANNOT.CSV FILE IMPLEMENTATION        *
 ****************************************/

typedef struct
{
    char *probe_set;
    char *dbsnp;
    int chrom;
    int pos;
    char strand;
    char allele_a;
    char allele_b;
}
record_t;

typedef struct
{
    void *probeset_id;
    record_t *records;
    int n_records, m_records;
}
annot_t;

static char *unquote(char *str)
{
    char *ptr = strrchr(str, '"');
    if ( ptr ) *ptr = '\0';
    return str + 1;
}

static int bcf_hdr_name2id_flexible(const bcf_hdr_t *hdr, const char *chr)
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

static annot_t *annot_init(const char *fn, const bcf_hdr_t *hdr)
{
    kstring_t str = {0, 0, NULL};
    htsFile *fp = unheader(fn, &str);

    int probeset_id = -1, dbsnp_id = -1, chrom_id = -1, position_id = -1, strand_id = -1, allele_a_id = -1, allele_b_id = -1;

    int moff = 0, *off = NULL;
    int ncols = ksplit_core(str.s, ',', &moff, &off);
    for (int i=0; i<ncols; i++)
    {
        if ( strcmp(&str.s[off[i]], "\"Probe Set ID\"")==0 ) probeset_id = i;
        else if ( strcmp(&str.s[off[i]], "\"dbSNP RS ID\"")==0 ) dbsnp_id = i;
        else if ( strcmp(&str.s[off[i]], "\"Chromosome\"")==0 ) chrom_id = i;
        else if ( strcmp(&str.s[off[i]], "\"Physical Position\"")==0 ) position_id = i;
        else if ( strcmp(&str.s[off[i]], "\"Strand\"")==0 ) strand_id = i;
        else if ( strcmp(&str.s[off[i]], "\"Allele A\"")==0 ) allele_a_id = i;
        else if ( strcmp(&str.s[off[i]], "\"Allele B\"")==0 ) allele_b_id = i;
    }
    if (probeset_id == -1) error("Probe Set ID missing from file: %s\n", fn);
    if (dbsnp_id == -1) error("dbSNP RS ID missing from file: %s\n", fn);
    if (chrom_id == -1) error("Chromosome missing from file: %s\n", fn);
    if (position_id == -1) error("Physical Position missing from file: %s\n", fn);
    if (strand_id == -1) error("Strand missing from file: %s\n", fn);
    if (allele_a_id == -1) error("Allele A missing from file: %s\n", fn);
    if (allele_b_id == -1) error("Allele B missing from file: %s\n", fn);

    annot_t *annot = (annot_t *)calloc(1, sizeof(annot_t));
    annot->probeset_id = khash_str2int_init();
    char *tmp;
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        hts_expand0(record_t, annot->n_records + 1, annot->m_records, annot->records);
        ncols = ksplit_core(str.s, ',', &moff, &off);
        char *ptr = unquote(&str.s[off[probeset_id]]);
        if ( strcmp(ptr, "---") )
        {
            annot->records[annot->n_records].probe_set = strdup(ptr);
            khash_str2int_inc(annot->probeset_id, annot->records[annot->n_records].probe_set);
        }
        ptr = unquote(&str.s[off[dbsnp_id]]);
        if ( strcmp(ptr, "---") ) annot->records[annot->n_records].dbsnp = strdup(ptr);
        ptr = unquote(&str.s[off[chrom_id]]);
        if ( strcmp(ptr, "---") ) annot->records[annot->n_records].chrom = bcf_hdr_name2id_flexible(hdr, ptr);
        ptr = unquote(&str.s[off[position_id]]);
        if ( strcmp(ptr, "---") ) annot->records[annot->n_records].pos = strtol(ptr, &tmp, 10);
        ptr = unquote(&str.s[off[strand_id]]);
        if ( strcmp(ptr, "---") ) annot->records[annot->n_records].strand = ptr[0];
        ptr = unquote(&str.s[off[allele_a_id]]);
        annot->records[annot->n_records].allele_a = ptr[0];
        ptr = unquote(&str.s[off[allele_b_id]]);
        annot->records[annot->n_records].allele_b = ptr[0];
        annot->n_records++;
    }

    free(off);
    free(str.s);
    hts_close(fp);
    return annot;
}

static void annot_destroy(annot_t *annot)
{
    khash_str2int_destroy(annot->probeset_id);
    for (int i=0; i<annot->n_records; i++)
    {
        free(annot->records[i].probe_set);
        free(annot->records[i].dbsnp);
    }
    free(annot->records);
}

/****************************************
 * REPORT.TXT FILE IMPLEMENTATION       *
 ****************************************/

typedef struct
{
    char **cel_files;
    int8_t *genders;
    int n_samples, m_cel_file, m_genders;
}
report_t;

static report_t *report_init(const char *fn)
{
    kstring_t str = {0, 0, NULL};
    htsFile *fp = unheader(fn, &str);
    int moff = 0, *off = NULL, ncols;
    ncols = ksplit_core(str.s, '\t', &moff, &off);
    if ( ncols < 2 ) error("Missing information in report file: %s\n", fn);
    if ( strcmp(&str.s[off[1]], "computed_gender") ) error("Second column not genders in file: %s\n", fn);

    report_t *report = (report_t *)calloc(1, sizeof(report_t));
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        ncols = ksplit_core(str.s, '\t', &moff, &off);
        if ( ncols < 2 ) error("Missing information in report file: %s\n", fn);
        hts_expand(report->cel_files, report->n_samples + 1, report->m_cel_file, report->cel_files);
        hts_expand0(int8_t, report->n_samples + 1, report->m_genders, report->genders);
        report->cel_files[report->n_samples] = strdup(&str.s[off[0]]);
        if ( strcmp(&str.s[off[1]], "male")==0 ) report->genders[report->n_samples] = 1;
        else if ( strcmp(&str.s[off[1]], "female")==0 ) report->genders[report->n_samples] = 2;
        report->n_samples++;
    }

    free(off);
    free(str.s);
    hts_close(fp);
    return report;
}

static void report_destroy(report_t *report)
{
    for (int i=0; i<report->n_samples; i++) free(report->cel_files[i]);
    free(report->genders);
}

/****************************************
 * OUTPUT FUNCTIONS                     *
 ****************************************/

static bcf_hdr_t *hdr_init(const faidx_t *fai)
{
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    int n = faidx_nseq(fai);
    for (int i=0; i<n; i++)
    {
        const char *seq = faidx_iseq(fai, i);
        int len = faidx_seq_len(fai, seq);
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq, len);
    }
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=CONF,Number=1,Type=Float,Description=\"Genotype confidences\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=NORMX,Number=1,Type=Float,Description=\"Normalized X intensity\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=NORMY,Number=1,Type=Float,Description=\"Normalized Y intensity\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DELTA,Number=1,Type=Float,Description=\"Normalized contrast value\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=SIZE,Number=1,Type=Float,Description=\"Normalized size value\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B Allele Frequency\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R Ratio\">");
    return hdr;
}

static void get_lrr_baf(const float *norm_x,
                        const float *norm_y,
                        int n,
                        const cluster_t *cluster,
                        float *delta,
                        float *size,
                        float *baf,
                        float *lrr)
{
    for (int i=0; i<n; i++) // compute THETA and R
    {
        float log2x = logf(norm_x[i]) * (float)M_LOG2E;
        float log2y = logf(norm_y[i]) * (float)M_LOG2E;
        delta[i] = log2x - log2y;
        size[i] = ( log2x + log2y ) * 0.5f;
    }

    for (int i=0; i<n; i++) // compute LRR and BAF
    {
        if ( delta[i] == cluster->ab_delta )
        {
            lrr[i] = size[i] - cluster->ab_size;
            baf[i] = 0.5f;
        }
        else if ( ( delta[i] < cluster->ab_delta && cluster->aa_delta < cluster->ab_delta ) ||
                  ( delta[i] > cluster->ab_delta && cluster->aa_delta > cluster->ab_delta ) )
        {
            float slope = ( cluster->aa_size - cluster->ab_size ) / ( cluster->aa_delta - cluster->ab_delta );
            float b = cluster->aa_size - ( cluster->aa_delta * slope );
            float size_ref = ( slope * delta[i] ) + b;
            lrr[i] = size[i] - size_ref;
            baf[i] = 0.5f - (cluster->ab_delta - delta[i]) * 0.5f / (cluster->ab_delta - cluster->aa_delta);
        }
        else if ( ( delta[i] > cluster->ab_delta && cluster->bb_delta > cluster->ab_delta ) ||
                  ( delta[i] < cluster->ab_delta && cluster->bb_delta < cluster->ab_delta ) )
        {
            float slope = ( cluster->ab_size - cluster->bb_size ) / ( cluster->ab_delta - cluster->bb_delta );
            float b = cluster->ab_size - ( cluster->ab_delta * slope );
            float size_ref = ( slope * delta[i] ) + b;
            lrr[i] = size[i] - size_ref;
            baf[i] = 1.0f - (cluster->bb_delta - delta[i]) * 0.5f / (cluster->bb_delta - cluster->ab_delta);
        }
        else
        {
            lrr[i] = NAN;
            baf[i] = NAN;
        }
        if ( baf[i] < 0.0f ) baf[i] = 0.0f;
        if ( baf[i] > 1.0f ) baf[i] = 1.0f;
    }
}

static void process(htsFile *out_fh,
                    bcf_hdr_t *hdr,
                    const annot_t *annot,
                    const snp_posteriors_t *snp_posteriors,
                    const char *summary_fn,
                    const char *calls_fn,
                    const char *confidences_fn)
{
    kstring_t str = {0, 0, NULL};
    int moff = 0, *off = NULL, ncols;

    htsFile *summary_fp = unheader(summary_fn, &str);
    ncols = ksplit_core(str.s, '\t', &moff, &off);
    if ( strcmp(&str.s[off[0]], "probeset_id") ) error("Malformed first line from summary file: %s\n%s\n", summary_fn, str.s);

    for (int i=1; i<ncols; i++) bcf_hdr_add_sample(hdr, &str.s[off[i]]);
    if ( bcf_hdr_write(out_fh, hdr) < 0 ) error("Unable to write to output VCF file\n");
    bcf_hdr_sync( hdr ); // updates the number of samples

    htsFile *confidences_fp = unheader(confidences_fn, &str);
    ncols = ksplit_core(str.s, '\t', &moff, &off);
    if ( strcmp(&str.s[off[0]], "probeset_id") ) error("Malformed first line from confidences file: %s\n%s\n", confidences_fn, str.s);

    htsFile *calls_fp = unheader(calls_fn, &str);
    ncols = ksplit_core(str.s, '\t', &moff, &off);
    if ( strcmp(&str.s[off[0]], "probeset_id") ) error("Malformed first line from calls file: %s\n%s\n", calls_fn, str.s);

    bcf1_t *rec = bcf_init();
    int32_t *gts = (int32_t *)malloc(bcf_hdr_nsamples(hdr)*2 * sizeof(int32_t));
    float *conf_arr = (float *)malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
    float *norm_x_arr = (float *)malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
    float *norm_y_arr = (float *)malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
    float *delta_arr = (float *)malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
    float *size_arr = (float *)malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
    float *baf_arr = (float *)malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
    float *lrr_arr = (float *)malloc(bcf_hdr_nsamples(hdr) * sizeof(float));
    cluster_t *cluster = snp_posteriors->clusters;
    while ( hts_getline(calls_fp, KS_SEP_LINE, &str) > 0 )
    {
        int idx;
        ncols = ksplit_core(str.s, '\t', &moff, &off);
        if ( ncols != bcf_hdr_nsamples(hdr) + 1 ) error("Expected %d columns but %d columns found in the call file\n", bcf_hdr_nsamples(hdr) + 1, ncols);
        int ret = khash_str2int_get( annot->probeset_id, &str.s[off[0]], &idx );
        if ( ret < 0 ) error("Probeset ID %s not found in manifest file\n", &str.s[off[0]]);
        record_t *record = &annot->records[idx];
        bcf_clear(rec);
        rec->n_sample = bcf_hdr_nsamples(hdr);
        bcf_update_alleles_str(hdr, rec, "N,N");
        rec->rid = record->chrom;
        rec->pos = record->pos - 1;
        if (record->dbsnp) bcf_update_id(hdr, rec, record->dbsnp);
        else bcf_update_id(hdr, rec, &str.s[off[0]]);
        if ( record->strand == '+' )
        {
            rec->d.allele[0][0] = record->allele_a;
            rec->d.allele[1][0] = record->allele_b;
        }
        else if ( record->strand == '-' )
        {
            rec->d.allele[0][0] = revnt(record->allele_a);
            rec->d.allele[1][0] = revnt(record->allele_b);
        }
        char *tmp;
        // read genotypes
        for (int i=1; i<ncols; i++)
        {
            int gt = strtol( &str.s[off[i]], &tmp, 10 );
            switch ( gt )
            {
                case -1:
                    gts[2*(i-1)] = bcf_gt_missing;
                    gts[2*(i-1)+1] = bcf_gt_missing;
                    break;
                case 0:
                    gts[2*(i-1)] = bcf_gt_unphased(0);
                    gts[2*(i-1)+1] = bcf_gt_unphased(0);
                    break;
                case 1:
                    gts[2*(i-1)] = bcf_gt_unphased(0);
                    gts[2*(i-1)+1] = bcf_gt_unphased(1);
                    break;
                case 2:
                    gts[2*(i-1)] = bcf_gt_unphased(1);
                    gts[2*(i-1)+1] = bcf_gt_unphased(1);
                    break;
                default:
                    error("Genotype for probeset ID %s is malformed: %s\n", &str.s[off[0]], &str.s[off[i]]);
                    break;
            }
        }
        // read confidences
        if ( hts_getline(confidences_fp, KS_SEP_LINE, &str) <=0 ) error("Confidences file ended prematurely\n");
        ncols = ksplit_core(str.s, '\t', &moff, &off);
        if ( ncols != bcf_hdr_nsamples(hdr) + 1 ) error("Expected %d columns but %d columns found in the confidences file\n", bcf_hdr_nsamples(hdr) + 1, ncols);
        for (int i=1; i<ncols; i++) conf_arr[i-1] = strtof(&str.s[off[i]], &tmp);

        // read intensities
        if ( hts_getline(summary_fp, KS_SEP_LINE, &str) <=0 ) error("Summary file ended prematurely\n");
        ncols = ksplit_core(str.s, '\t', &moff, &off);
        if ( ncols != bcf_hdr_nsamples(hdr) + 1 ) error("Expected %d columns but %d columns found in the confidences file\n", bcf_hdr_nsamples(hdr) + 1, ncols);
        for (int i=1; i<ncols; i++) norm_x_arr[i-1] = strtof(&str.s[off[i]], &tmp);
        if ( hts_getline(summary_fp, KS_SEP_LINE, &str) <=0 ) error("Summary file ended prematurely\n");
        ncols = ksplit_core(str.s, '\t', &moff, &off);
        if ( ncols != bcf_hdr_nsamples(hdr) + 1 ) error("Expected %d columns but %d columns found in the confidences file\n", bcf_hdr_nsamples(hdr) + 1, ncols);
        for (int i=1; i<ncols; i++) norm_y_arr[i-1] = strtof(&str.s[off[i]], &tmp);

        // compute LRR and BAF
        get_lrr_baf(norm_x_arr, norm_y_arr, bcf_hdr_nsamples(hdr), cluster, delta_arr, size_arr, baf_arr, lrr_arr);

        bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr)*2);
        bcf_update_format_float(hdr, rec, "CONF", conf_arr, bcf_hdr_nsamples(hdr));
        bcf_update_format_float(hdr, rec, "NORMX", norm_x_arr, bcf_hdr_nsamples(hdr));
        bcf_update_format_float(hdr, rec, "NORMY", norm_y_arr, bcf_hdr_nsamples(hdr));
        bcf_update_format_float(hdr, rec, "DELTA", delta_arr, bcf_hdr_nsamples(hdr));
        bcf_update_format_float(hdr, rec, "SIZE", size_arr, bcf_hdr_nsamples(hdr));
        bcf_update_format_float(hdr, rec, "BAF", baf_arr, bcf_hdr_nsamples(hdr));
        bcf_update_format_float(hdr, rec, "LRR", lrr_arr, bcf_hdr_nsamples(hdr));
        if ( bcf_write(out_fh, hdr, rec) < 0 ) error("Unable to write to output VCF file\n");
        cluster++;
    }
    free(conf_arr);
    free(norm_x_arr);
    free(norm_y_arr);
    free(delta_arr);
    free(size_arr);
    free(baf_arr);
    free(lrr_arr);
    free(gts);
    bcf_destroy(rec);
    free(off);
    free(str.s);
    hts_close(summary_fp);
    hts_close(confidences_fp);
    hts_close(calls_fp);
    return;
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void)
{
    return "convert Affymetrix files to VCF\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: convert Affymetrix apt-probeset-genotype output files to VCF (2018-09-07)\n"
        "\n"
        "Usage: bcftools +affy2vcf [options] --fasta-ref <fasta> --annot <file> --snp-posteriors <file>\n"
        "                                         --summary <file> --calls <file> --confidences <file> \n"
        "\n"
        "Plugin options:\n"
        "    -f, --fasta-ref <file>                     reference sequence in fasta format\n"
        "        --annot <file>                         probeset annotation file\n"
        "        --snp-posteriors <snp-posteriors.txt>  apt-probeset-genotype snp-posteriors output\n"
        "        --summary <summary.txt>                apt-probeset-genotype summary output\n"
        "        --report <report.txt>                  apt-probeset-genotype report output\n"
        "        --calls <calls.txt>                    apt-probeset-genotype calls output\n"
        "        --confidences <confidences.txt>        apt-probeset-genotype confidences output\n"
        "    -x, --sex <file>                           output apt-probeset-genotype gender estimate into file\n"
        "        --no-version                           do not append version and command line to the header\n"
        "    -o, --output <file>                        write output to a file [standard output]\n"
        "    -O, --output-type b|u|z|v                  b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF\n"
        "        --threads <int>                        number of extra output compression threads [0]\n"
        "\n";
}

int run(int argc, char *argv[])
{
    char *ref_fname = NULL;
    char *sex_fname = NULL;
    char *annot_fname = NULL;
    char *snp_posteriors_fname = NULL;
    char *summary_fname = NULL;
    char *report_fname = NULL;
    char *calls_fname = NULL;
    char *confidences_fname = NULL;
    char *output_fname = "-";
    int output_type = FT_VCF;
    int n_threads = 0;
    int record_cmd_line = 1;
    faidx_t *fai = NULL;

    static struct option loptions[] =
    {
        {"annot", required_argument, NULL, 1},
        {"snp-posteriors", required_argument, NULL, 2},
        {"summary", required_argument, NULL, 3},
        {"report", required_argument, NULL, 4},
        {"calls", required_argument, NULL, 5},
        {"confidences", required_argument, NULL, 6},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"sex", required_argument, NULL, 'x'},
        {"output", required_argument, NULL, 'o'},
        {"output-type", required_argument, NULL, 'O'},
        {"no-version", no_argument, NULL, 8},
        {"threads", required_argument, NULL, 9},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "h?lt:o:O:i:b:c:e:f:g:x:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'o': output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': output_type = FT_BCF_GZ; break;
                          case 'u': output_type = FT_BCF; break;
                          case 'z': output_type = FT_VCF_GZ; break;
                          case 'v': output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'f': ref_fname = optarg; break;
            case 'x': sex_fname = optarg; break;
            case  1 : annot_fname = optarg; break;
            case  2 : snp_posteriors_fname = optarg; break;
            case  3 : summary_fname = optarg; break;
            case  4 : report_fname = optarg; break;
            case  5 : calls_fname = optarg; break;
            case  6 : confidences_fname = optarg; break;
            case  9 : n_threads = strtol(optarg, NULL, 0); break;
            case  8 : record_cmd_line = 0; break;
            case 'h':
            case '?':
            default: error("%s", usage_text());
        }
    }
    if ( !ref_fname ) error("Expected --fasta-ref option\n");
    if ( !annot_fname ) error("Expected --annot option\n");
    if ( !snp_posteriors_fname ) error("Expected --snp-posteriors option\n");
    if ( !summary_fname ) error("Expected --summary option\n");
    if ( sex_fname && !report_fname ) error("Expected --report option with --sex option\n");
    if ( !calls_fname ) error("Expected --calls option\n");
    if ( !confidences_fname ) error("Expected --confidences option\n");

    fai = fai_load(ref_fname);
    if ( !fai ) error("Could not load the reference %s\n", ref_fname);
    bcf_hdr_t *hdr = hdr_init(fai);
    fai_destroy(fai);

    if (record_cmd_line) bcf_hdr_append_version(hdr, argc, argv, "bcftools_+affy2vcf");
    htsFile *out_fh = hts_open(output_fname, hts_bcf_wmode(output_type));
    if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
    if ( n_threads ) hts_set_threads(out_fh, n_threads);

    annot_t *annot = annot_init(annot_fname, hdr);

    snp_posteriors_t *snp_posteriors = snp_posteriors_init(snp_posteriors_fname);

    if ( sex_fname )
    {
        report_t *report = report_init(report_fname);
        FILE *sex_fh = fopen(sex_fname, "w");
        if ( !sex_fh ) error("Failed to open %s: %s\n", sex_fname, strerror(errno));
        for (int i=0; i<report->n_samples; i++) fprintf(sex_fh, "%s\t%d\n", report->cel_files[i], report->genders[i]);
        fclose(sex_fh);
        report_destroy(report);
    }

    process(out_fh, hdr, annot, snp_posteriors, summary_fname , calls_fname, confidences_fname);

    snp_posteriors_destroy(snp_posteriors);
    annot_destroy(annot);
    bcf_hdr_destroy(hdr);
    hts_close(out_fh);
    return 0;
}