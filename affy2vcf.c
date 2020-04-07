/* The MIT License

   Copyright (c) 2018-2020 Giulio Genovese

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
#include <htslib/hfile.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include "bcftools.h"
#include "htslib/khash_str2int.h"
#include "gtc2vcf.h"

#define AFFY2VCF_VERSION "2020-04-07"

#define GT_NC -1
#define GT_AA 0
#define GT_AB 1
#define GT_BB 2

#define VERBOSE (1 << 0)
#define CALLS_LOADED (1 << 1)
#define CONFIDENCES_LOADED (1 << 2)
#define SUMMARY_LOADED (1 << 3)
#define SNP_POSTERIORS_LOADED (1 << 4)
#define ADJUST_CLUSTERS (1 << 5)

/****************************************
 * htsFILE READING FUNCTIONS            *
 ****************************************/

static htsFile *unheader(const char *fn, kstring_t *str)
{
	htsFile *fp = hts_open(fn, "r");
	if (!fp)
		error("Could not read: %s\n", fn);

	if (hts_getline(fp, KS_SEP_LINE, str) <= 0)
		error("Empty file: %s\n", fn);

	// skip header
	while (str->s[0] == '#')
		hts_getline(fp, KS_SEP_LINE, str);

	return fp;
}

/****************************************
 * ANNOT.CSV FILE IMPLEMENTATION        *
 ****************************************/

typedef struct {
	char *probe_set_id;
	char *affy_snp_id;
	char *dbsnp_rs_id;
	char *chromosome;
	int position;
	int strand;
	char *flank;
} record_t;

typedef struct {
	void *probe_set_id;
	record_t *records;
	int n_records, m_records;
} annot_t;

static inline char *unquote(char *str)
{
	if (strcmp(str, "\"---\"") == 0)
		return NULL;
	char *ptr = strrchr(str, '"');
	if (ptr)
		*ptr = '\0';
	return str + 1;
}

static annot_t *annot_init(const char *fn, const char *sam_fn, const char *out_fn, int flags)
{
	annot_t *annot = NULL;
	FILE *out_txt = get_file_handle(out_fn);
	htsFile *hts = NULL;
	sam_hdr_t *sam_hdr = NULL;
	bam1_t *b = NULL;
	if (sam_fn) {
		hts = hts_open(sam_fn, "r");
		if (hts == NULL || hts_get_format(hts)->category != sequence_data)
			error("File %s does not contain sequence data\n", sam_fn);
		sam_hdr = sam_hdr_read(hts);
		if (sam_hdr == NULL)
			error("Reading header from \"%s\" failed", sam_fn);
		b = bam_init1();
		if (b == NULL)
			error("Cannot create SAM record\n");
	}
	kstring_t str = {0, 0, NULL};

	htsFile *fp = hts_open(fn, "r");
	if (!fp)
		error("Could not read: %s\n", fn);
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0)
		error("Empty file: %s\n", fn);
	const char *null_strand = "---";
	while (str.s[0] == '#') {
		if (strcmp(str.s, "#%netaffx-annotation-tabular-format-version=1.0") == 0)
			null_strand = "---";
		if (strcmp(str.s, "#%netaffx-annotation-tabular-format-version=1.5") == 0)
			null_strand = "+";
		if (hts && out_txt)
			fprintf(out_txt, "%s\n", str.s);
		hts_getline(fp, KS_SEP_LINE, &str);
	}

	if (hts && out_txt)
		fprintf(out_txt, "%s\n", str.s);

	int probe_set_id_idx = -1;
	int affy_snp_id_idx = -1;
	int dbsnp_rs_id_idx = -1;
	int chromosome_idx = -1;
	int position_idx = -1;
	int position_end_idx = -1;
	int strand_idx = -1;
	int flank_idx = -1;
	int allele_a_idx = -1;
	int allele_b_idx = -1;

	int moff = 0, *off = NULL;
	int ncols = ksplit_core(str.s, ',', &moff, &off);
	for (int i = 0; i < ncols; i++) {
		if (strcmp(&str.s[off[i]], "\"Probe Set ID\"") == 0)
			probe_set_id_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Affy SNP ID\"") == 0)
			affy_snp_id_idx = i;
		else if (strcmp(&str.s[off[i]], "\"dbSNP RS ID\"") == 0)
			dbsnp_rs_id_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Chromosome\"") == 0)
			chromosome_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Physical Position\"") == 0)
			position_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Position End\"") == 0)
			position_end_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Strand\"") == 0)
			strand_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Flank\"") == 0)
			flank_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Allele A\"") == 0)
			allele_a_idx = i;
		else if (strcmp(&str.s[off[i]], "\"Allele B\"") == 0)
			allele_b_idx = i;
	}
	if (probe_set_id_idx != 0)
		error("Probe Set ID not the first column in file: %s\n", fn);
	if (flank_idx == -1)
		error("Flank missing from file: %s\n", fn);
	if (allele_a_idx == -1)
		error("Allele A missing from file: %s\n", fn);
	if (allele_b_idx == -1)
		error("Allele B missing from file: %s\n", fn);
	const char *probe_set_id, *flank, *allele_a, *allele_b;

	if (!hts && out_txt) {

		while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
			ncols = ksplit_core(str.s, ',', &moff, &off);
			probe_set_id = unquote(&str.s[off[probe_set_id_idx]]);
			flank = unquote(&str.s[off[flank_idx]]);
			if (flank)
				flank2fasta(probe_set_id, flank, out_txt);
		}
	} else {
		if (dbsnp_rs_id_idx == -1)
			error("dbSNP RS ID missing from file: %s\n", fn);
		if (chromosome_idx == -1)
			error("Chromosome missing from file: %s\n", fn);
		if (position_idx == -1)
			error("Physical Position missing from file: %s\n", fn);
		if (strand_idx == -1)
			error("Strand missing from file: %s\n", fn);

		if (!out_txt) {
			annot = (annot_t *)calloc(1, sizeof(annot_t));
			annot->probe_set_id = khash_str2int_init();
		}

		char *tmp;
		int n_total = 0, n_unmapped = 0;
		while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
			ncols = ksplit_core(str.s, ',', &moff, &off);
			probe_set_id = unquote(&str.s[off[probe_set_id_idx]]);
			flank = unquote(&str.s[off[flank_idx]]);
			allele_a = unquote(&str.s[off[allele_a_idx]]);
			allele_b = unquote(&str.s[off[allele_b_idx]]);
			const char *chromosome = NULL;
			int strand = -1, position = 0, idx = -1;
			if (hts) {
				if (!flank) {
					if (flags & VERBOSE)
						fprintf(stderr,
							"Missing flank sequence for marker %s\n",
							probe_set_id);
					n_unmapped++;
				} else {
					idx = get_position(hts, sam_hdr, b, probe_set_id, flank,
							   0, &chromosome, &position, &strand);
					if (idx < 0)
						error("Reading from %s failed", sam_fn);
					else if (idx == 0) {
						if (flags & VERBOSE)
							fprintf(stderr,
								"Unable to determine position for marker %s\n",
								probe_set_id);
						n_unmapped++;
					}
				}
				n_total++;
			} else {
				chromosome = unquote(&str.s[off[chromosome_idx]]);
				const char *ptr = unquote(&str.s[off[position_idx]]);
				position = ptr ? strtol(ptr, &tmp, 10) : 0;
				ptr = unquote(&str.s[off[strand_idx]]);
				if (!ptr)
					strand = -1;
				else if (strcmp(ptr, "+") == 0)
					strand = 0;
				else if (strcmp(ptr, "-") == 0)
					strand = 1;
				else
					strand = -1;
			}

			if (out_txt) {
				// "Ref Allele" and "Alt Allele" will not be updated
				fprintf(out_txt, "\"%s\"", probe_set_id);
				for (int i = 1; i < ncols; i++) {
					if (i == flank_idx)
						fprintf(out_txt, ",\"%s\"", flank);
					if (i == allele_a_idx)
						fprintf(out_txt, ",\"%s\"", allele_a);
					if (i == allele_b_idx) {
						fprintf(out_txt, ",\"%s\"", allele_b);
					} else if (i == chromosome_idx) {
						if (chromosome)
							fprintf(out_txt, ",\"%s\"", chromosome);
						else
							fprintf(out_txt, ",\"---\"");
					} else if (i == position_idx) {
						if (position)
							fprintf(out_txt, ",\"%d\"", position);
						else
							fprintf(out_txt, ",\"---\"");
					} else if (i == position_end_idx) {
						if (flank && position && idx > 0) {
							const char *left = strchr(flank, '[');
							const char *middle = strchr(flank, '/');
							const char *right = strchr(flank, ']');
							if (!left || !middle || !right)
								error("Flank sequence is malformed: %s\n",
								      flank);

							fprintf(out_txt, ",\"%d\"",
								position
									+ (int)(idx > 1 ? right - middle
											: middle - left
												  + (*(left
												       + 1)
												     == '-'))
									- 2);
						} else {
							fprintf(out_txt, ",\"---\"");
						}
					} else if (i == strand_idx) {
						fprintf(out_txt, ",\"%s\"",
							strand == 0
								? "+"
								: (strand == 1 ? "-"
									       : null_strand));
					} else {
						fprintf(out_txt, ",%s", &str.s[off[i]]);
					}
				}
				fprintf(out_txt, "\n");
			} else {
				hts_expand0(record_t, annot->n_records + 1, annot->m_records,
					    annot->records);
				annot->records[annot->n_records].probe_set_id =
					strdup(probe_set_id);
				khash_str2int_inc(
					annot->probe_set_id,
					annot->records[annot->n_records].probe_set_id);
				const char *dbsnp_rs_id = unquote(&str.s[off[dbsnp_rs_id_idx]]);
				if (dbsnp_rs_id)
					annot->records[annot->n_records].dbsnp_rs_id =
						strdup(dbsnp_rs_id);
				if (affy_snp_id_idx >= 0) {
					const char *affy_snp_id =
						unquote(&str.s[off[affy_snp_id_idx]]);
					if (affy_snp_id)
						annot->records[annot->n_records].affy_snp_id =
							strdup(affy_snp_id);
				}
				if (chromosome)
					annot->records[annot->n_records].chromosome =
						strdup(chromosome);
				annot->records[annot->n_records].position = position;
				if (flank) {
					annot->records[annot->n_records].flank = strdup(flank);
					// check whether alleles A and B need to be flipped in
					// the flank sequence (happens with T/C and T/G SNPs
					// only)
					char *left = strchr(
						annot->records[annot->n_records].flank, '[');
					char *middle = strchr(
						annot->records[annot->n_records].flank, '/');
					char *right = strchr(
						annot->records[annot->n_records].flank, ']');
					if (strncmp(left + 1, allele_b, middle - left - 1) == 0
					    && strncmp(middle + 1, allele_a, right - middle - 1)
						       == 0) {
						memcpy(left + 1, allele_a, right - middle - 1);
						*(left + (right - middle)) = '/';
						memcpy(left + (right - middle) + 1, allele_b,
						       middle - left - 1);
					}
				}
				annot->records[annot->n_records].strand = strand;
				annot->n_records++;
			}
		}
		if (hts)
			fprintf(stderr, "Lines   total/unmapped:\t%d/%d\n", n_total,
				n_unmapped);

		bam_destroy1(b);
		sam_hdr_destroy(sam_hdr);
		if (hts && hts_close(hts) < 0)
			error("closing \"%s\" failed", fn);
	}

	free(off);
	free(str.s);
	hts_close(fp);

	if (out_txt && out_txt != stdout && out_txt != stderr)
		fclose(out_txt);
	return annot;
}

static void annot_destroy(annot_t *annot)
{
	khash_str2int_destroy(annot->probe_set_id);
	for (int i = 0; i < annot->n_records; i++) {
		free(annot->records[i].probe_set_id);
		free(annot->records[i].affy_snp_id);
		free(annot->records[i].dbsnp_rs_id);
		free(annot->records[i].chromosome);
		free(annot->records[i].flank);
	}
	free(annot->records);
	free(annot);
}

/****************************************
 * REPORT.TXT FILE IMPLEMENTATION       *
 ****************************************/

typedef struct {
	char **cel_files;
	int8_t *genders;
	int n_samples, m_cel_file, m_genders;
} report_t;

static report_t *report_init(const char *fn)
{
	kstring_t str = {0, 0, NULL};
	htsFile *fp = unheader(fn, &str);
	int moff = 0, *off = NULL, ncols;
	ncols = ksplit_core(str.s, '\t', &moff, &off);
	if (ncols < 2)
		error("Missing information in report file: %s\n", fn);
	if (strcmp(&str.s[off[1]], "computed_gender"))
		error("Second column not genders in file: %s\n", fn);

	report_t *report = (report_t *)calloc(1, sizeof(report_t));
	while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
		ncols = ksplit_core(str.s, '\t', &moff, &off);
		if (ncols < 2)
			error("Missing information in report file: %s\n", fn);
		hts_expand(report->cel_files, report->n_samples + 1, report->m_cel_file,
			   report->cel_files);
		hts_expand0(int8_t, report->n_samples + 1, report->m_genders, report->genders);
		report->cel_files[report->n_samples] = strdup(&str.s[off[0]]);
		if (strcmp(&str.s[off[1]], "male") == 0)
			report->genders[report->n_samples] = 1;
		else if (strcmp(&str.s[off[1]], "female") == 0)
			report->genders[report->n_samples] = 2;
		report->n_samples++;
	}

	free(off);
	free(str.s);
	hts_close(fp);
	return report;
}

static void report_destroy(report_t *report)
{
	for (int i = 0; i < report->n_samples; i++)
		free(report->cel_files[i]);
	free(report->cel_files);
	free(report->genders);
	free(report);
}

/****************************************
 * OUTPUT FUNCTIONS                     *
 ****************************************/

static bcf_hdr_t *hdr_init(const faidx_t *fai, int flags)
{
	bcf_hdr_t *hdr = bcf_hdr_init("w");
	int n = faidx_nseq(fai);
	for (int i = 0; i < n; i++) {
		const char *seq = faidx_iseq(fai, i);
		int len = faidx_seq_len(fai, seq);
		bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq, len);
	}
	bcf_hdr_append(hdr,
		       "##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description=\"A allele\">");
	bcf_hdr_append(hdr,
		       "##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description=\"B allele\">");
	bcf_hdr_append(
		hdr,
		"##INFO=<ID=PROBE_SET_ID,Number=1,Type=String,Description=\"Affymetrix Probe Set ID\">");
	bcf_hdr_append(
		hdr,
		"##INFO=<ID=AFFY_SNP_ID,Number=1,Type=String,Description=\"Affymetrix SNP ID\">");
	if (flags & SNP_POSTERIORS_LOADED) {
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanDELTA_AA,Number=1,Type=Float,Description=\"Mean of normalized DELTA for AA cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanDELTA_AB,Number=1,Type=Float,Description=\"Mean of normalized DELTA for AB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanDELTA_BB,Number=1,Type=Float,Description=\"Mean of normalized DELTA for BB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devDELTA_AA,Number=1,Type=Float,Description=\"Standard deviation of normalized DELTA for AA cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devDELTA_AB,Number=1,Type=Float,Description=\"Standard deviation of normalized DELTA for AB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devDELTA_BB,Number=1,Type=Float,Description=\"Standard deviation of normalized DELTA for BB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanN_AA,Number=1,Type=Float,Description=\"Number of AA calls in training set for mean\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanN_AB,Number=1,Type=Float,Description=\"Number of AB calls in training set for mean\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanN_BB,Number=1,Type=Float,Description=\"Number of BB calls in training set for mean\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devN_AA,Number=1,Type=Float,Description=\"Number of AA calls in training set for standard deviation\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devN_AB,Number=1,Type=Float,Description=\"Number of AB calls in training set for standard deviation\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devN_BB,Number=1,Type=Float,Description=\"Number of BB calls in training set for standard deviation\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanSIZE_AA,Number=1,Type=Float,Description=\"Mean of normalized SIZE for AA cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanSIZE_AB,Number=1,Type=Float,Description=\"Mean of normalized SIZE for AB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=meanSIZE_BB,Number=1,Type=Float,Description=\"Mean of normalized SIZE for BB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devSIZE_AA,Number=1,Type=Float,Description=\"Standard deviation of normalized SIZE for AA cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devSIZE_AB,Number=1,Type=Float,Description=\"Standard deviation of normalized SIZE for AB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=devSIZE_BB,Number=1,Type=Float,Description=\"Standard deviation of normalized SIZE for BB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=covar_AA,Number=1,Type=Float,Description=\"Covariance for AA cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=covar_AB,Number=1,Type=Float,Description=\"Covariance for AB cluster\">");
		bcf_hdr_append(
			hdr,
			"##INFO=<ID=covar_BB,Number=1,Type=Float,Description=\"Covariance for BB cluster\">");
	}
	if (flags & CALLS_LOADED)
		bcf_hdr_append(
			hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	if (flags & CONFIDENCES_LOADED)
		bcf_hdr_append(
			hdr,
			"##FORMAT=<ID=CONF,Number=1,Type=Float,Description=\"Genotype confidences\">");
	if (flags & SUMMARY_LOADED) {
		bcf_hdr_append(
			hdr,
			"##FORMAT=<ID=NORMX,Number=1,Type=Float,Description=\"Normalized X intensity\">");
		bcf_hdr_append(
			hdr,
			"##FORMAT=<ID=NORMY,Number=1,Type=Float,Description=\"Normalized Y intensity\">");
		bcf_hdr_append(
			hdr,
			"##FORMAT=<ID=DELTA,Number=1,Type=Float,Description=\"Normalized contrast value\">");
		bcf_hdr_append(
			hdr,
			"##FORMAT=<ID=SIZE,Number=1,Type=Float,Description=\"Normalized size value\">");
	}
	if ((flags & SUMMARY_LOADED) && (flags & SNP_POSTERIORS_LOADED)) {
		bcf_hdr_append(
			hdr,
			"##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B Allele Frequency\">");
		bcf_hdr_append(
			hdr,
			"##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R Ratio\">");
	}
	return hdr;
}

// compute DELTA (contrast) and SIZE
static void get_delta_size(const float *norm_x, const float *norm_y, int n, float *delta,
			   float *size)
{
	for (int i = 0; i < n; i++) {
		float log2x = logf(norm_x[i]) * (float)M_LOG2E;
		float log2y = logf(norm_y[i]) * (float)M_LOG2E;
		delta[i] = log2x - log2y;
		size[i] = (log2x + log2y) * 0.5f;
	}
}

typedef struct {
	float aa_delta_mean;
	float aa_delta_dev;
	float aa_n_mean;
	float aa_n_dev;
	float aa_size_mean;
	float aa_size_dev;
	float aa_cov;

	float ab_delta_mean;
	float ab_delta_dev;
	float ab_n_mean;
	float ab_n_dev;
	float ab_size_mean;
	float ab_size_dev;
	float ab_cov;

	float bb_delta_mean;
	float bb_delta_dev;
	float bb_n_mean;
	float bb_n_dev;
	float bb_size_mean;
	float bb_size_dev;
	float bb_cov;
} cluster_t;

// adjust cluster centers (using apt-probeset-genotype posteriors as priors)
// similar to
// https://github.com/WGLab/PennCNV/blob/master/affy/bin/generate_affy_geno_cluster.pl
static void adjust_clusters(const int *gts, const float *delta, const float *size, int n,
			    cluster_t *cluster)
{
	cluster->aa_delta_mean *= 0.2f;
	cluster->ab_delta_mean *= 0.2f;
	cluster->bb_delta_mean *= 0.2f;
	cluster->aa_size_mean *= 0.2f;
	cluster->ab_size_mean *= 0.2f;
	cluster->bb_size_mean *= 0.2f;
	cluster->aa_n_mean = 0.2f;
	cluster->ab_n_mean = 0.2f;
	cluster->bb_n_mean = 0.2f;

	for (int i = 0; i < n; i++) {
		switch (gts[i]) {
		case GT_AA:
			cluster->aa_n_mean++;
			cluster->aa_delta_mean += delta[i];
			cluster->aa_size_mean += size[i];
			break;
		case GT_AB:
			cluster->ab_n_mean++;
			cluster->ab_delta_mean += delta[i];
			cluster->ab_size_mean += size[i];
			break;
		case GT_BB:
			cluster->bb_n_mean++;
			cluster->bb_delta_mean += delta[i];
			cluster->bb_size_mean += size[i];
			break;
		default:
			break;
		}
	}

	cluster->aa_delta_mean /= cluster->aa_n_mean;
	cluster->ab_delta_mean /= cluster->ab_n_mean;
	cluster->bb_delta_mean /= cluster->bb_n_mean;
	cluster->aa_size_mean /= cluster->aa_n_mean;
	cluster->ab_size_mean /= cluster->ab_n_mean;
	cluster->bb_size_mean /= cluster->bb_n_mean;
}

// compute LRR and BAF
// similar to
// https://github.com/WGLab/PennCNV/blob/master/affy/bin/normalize_affy_geno_cluster.pl
static void get_baf_lrr(const float *delta, const float *size, int n, const cluster_t *cluster,
			float *baf, float *lrr)
{
	for (int i = 0; i < n; i++) {
		if (delta[i] == cluster->ab_delta_mean) {
			lrr[i] = size[i] - cluster->ab_size_mean;
			baf[i] = 0.5f;
		} else if ((delta[i] < cluster->ab_delta_mean
			    && cluster->aa_delta_mean < cluster->ab_delta_mean)
			   || (delta[i] > cluster->ab_delta_mean
			       && cluster->aa_delta_mean > cluster->ab_delta_mean)) {
			float slope = (cluster->aa_size_mean - cluster->ab_size_mean)
				      / (cluster->aa_delta_mean - cluster->ab_delta_mean);
			float b = cluster->aa_size_mean - (cluster->aa_delta_mean * slope);
			float size_ref = (slope * delta[i]) + b;
			lrr[i] = size[i] - size_ref;
			baf[i] = 0.5f
				 - (cluster->ab_delta_mean - delta[i]) * 0.5f
					   / (cluster->ab_delta_mean - cluster->aa_delta_mean);
		} else if ((delta[i] > cluster->ab_delta_mean
			    && cluster->bb_delta_mean > cluster->ab_delta_mean)
			   || (delta[i] < cluster->ab_delta_mean
			       && cluster->bb_delta_mean < cluster->ab_delta_mean)) {
			float slope = (cluster->ab_size_mean - cluster->bb_size_mean)
				      / (cluster->ab_delta_mean - cluster->bb_delta_mean);
			float b = cluster->ab_size_mean - (cluster->ab_delta_mean * slope);
			float size_ref = (slope * delta[i]) + b;
			lrr[i] = size[i] - size_ref;
			baf[i] = 1.0f
				 - (cluster->bb_delta_mean - delta[i]) * 0.5f
					   / (cluster->bb_delta_mean - cluster->ab_delta_mean);
		} else {
			lrr[i] = NAN;
			baf[i] = NAN;
		}
		if (baf[i] < 0.0f)
			baf[i] = 0.0f;
		if (baf[i] > 1.0f)
			baf[i] = 1.0f;
	}
}

#define MAX_LENGTH_PROBE_SET_ID 24
static void process(faidx_t *fai, const annot_t *annot, const char *calls_fn,
		    const char *confidences_fn, const char *summary_fn,
		    const char *snp_posteriors_fn, htsFile *out_fh, bcf_hdr_t *hdr, int flags)
{
	kstring_t str = {0, 0, NULL};
	int moff = 0, *off = NULL, ncols;
	int moff2 = 0, *off2 = NULL, ncols2;

	htsFile *calls_fp = NULL;
	if (calls_fn) {
		calls_fp = unheader(calls_fn, &str);
		ncols = ksplit_core(str.s, '\t', &moff, &off);
		if (strcmp(&str.s[off[0]], "probeset_id"))
			error("Malformed first line from calls file: %s\n%s\n", calls_fn,
			      str.s);
		for (int i = 1; i < ncols; i++)
			bcf_hdr_add_sample(hdr, &str.s[off[i]]);
	}

	htsFile *confidences_fp = NULL;
	if (confidences_fn) {
		confidences_fp = unheader(confidences_fn, &str);
		ncols = ksplit_core(str.s, '\t', &moff, &off);
		if (strcmp(&str.s[off[0]], "probeset_id"))
			error("Malformed first line from confidences file: %s\n%s\n",
			      confidences_fn, str.s);
		if (!calls_fp)
			for (int i = 1; i < ncols; i++)
				bcf_hdr_add_sample(hdr, &str.s[off[i]]);
	}

	htsFile *summary_fp = NULL;
	if (summary_fn) {
		summary_fp = unheader(summary_fn, &str);
		ncols = ksplit_core(str.s, '\t', &moff, &off);
		if (strcmp(&str.s[off[0]], "probeset_id"))
			error("Malformed first line from summary file: %s\n%s\n", summary_fn,
			      str.s);
		if (!calls_fp && !confidences_fp)
			for (int i = 1; i < ncols; i++)
				bcf_hdr_add_sample(hdr, &str.s[off[i]]);
	}

	htsFile *snp_posteriors_fp = NULL;
	if (snp_posteriors_fn) {
		snp_posteriors_fp = unheader(snp_posteriors_fn, &str);
		ncols = ksplit_core(str.s, '\t', &moff, &off);
		if (strcmp(&str.s[off[0]], "id"))
			error("Malformed first line from posterior models file: %s\n%s\n",
			      snp_posteriors_fn, str.s);
	}

	if (bcf_hdr_write(out_fh, hdr) < 0)
		error("Unable to write to output VCF file\n");
	if (bcf_hdr_sync(hdr) < 0)
		error_errno("[%s] Failed to update header",
			    __func__); // updates the number of samples
	int nsmpl = bcf_hdr_nsamples(hdr);
	if ((flags & ADJUST_CLUSTERS) && (nsmpl < 100))
		fprintf(stderr,
			"Warning: adjusting clusters with %d sample(s) is not recommended\n",
			nsmpl);

	bcf1_t *rec = bcf_init();
	char ref_base[] = {'\0', '\0'};
	kstring_t allele_a = {0, 0, NULL};
	kstring_t allele_b = {0, 0, NULL};
	kstring_t flank = {0, 0, NULL};

	int *gts = (int *)malloc(nsmpl * sizeof(int));
	int32_t *gt_arr = (int32_t *)malloc(nsmpl * 2 * sizeof(int32_t));
	float *conf_arr = (float *)malloc(nsmpl * sizeof(float));
	float *norm_x_arr = (float *)malloc(nsmpl * sizeof(float));
	float *norm_y_arr = (float *)malloc(nsmpl * sizeof(float));
	float *delta_arr = (float *)malloc(nsmpl * sizeof(float));
	float *size_arr = (float *)malloc(nsmpl * sizeof(float));
	float *baf_arr = (float *)malloc(nsmpl * sizeof(float));
	float *lrr_arr = (float *)malloc(nsmpl * sizeof(float));

	kstring_t probe_set_id = {0, 0, NULL};
	int i = 0, n_missing = 0, n_skipped = 0;
	for (i = 0; i < annot->n_records; i++) {
		char *tmp, buf[MAX_LENGTH_PROBE_SET_ID];
		probe_set_id.l = 0;

		// read genotypes
		if (calls_fp) {
			if (hts_getline(calls_fp, KS_SEP_LINE, &str) < 0)
				break;
			ncols = ksplit_core(str.s, '\t', &moff, &off);
			if (ncols != 1 + nsmpl)
				error("Expected %d columns but %d columns found in the calls file\n",
				      1 + nsmpl, ncols);
			kputs(&str.s[off[0]], &probe_set_id);
			for (int i = 1; i < 1 + nsmpl; i++)
				gts[i - 1] = strtol(&str.s[off[i]], &tmp, 10);
		}

		// read confidences
		if (confidences_fp) {
			if (hts_getline(confidences_fp, KS_SEP_LINE, &str) < 0)
				break;
			ncols = ksplit_core(str.s, '\t', &moff, &off);
			if (ncols != 1 + nsmpl)
				error("Expected %d columns but %d columns found in the confidences file\n",
				      1 + nsmpl, ncols);
			if (probe_set_id.l == 0)
				kputs(&str.s[off[0]], &probe_set_id);
			else if (strcmp(&str.s[off[0]], probe_set_id.s) != 0)
				error("Found Probe Set ID %s in the confidences file while %s expected\n",
				      &str.s[off[0]], probe_set_id.s);
			for (int i = 1; i < 1 + nsmpl; i++)
				conf_arr[i - 1] = strtof(&str.s[off[i]], &tmp);
			bcf_update_format_float(hdr, rec, "CONF", conf_arr, nsmpl);
		}

		// read intensities
		if (summary_fp) {
			int ret, len;
			do {
				if ((ret = hts_getline(summary_fp, KS_SEP_LINE, &str)) < 0)
					break;
				ncols = ksplit_core(str.s, '\t', &moff, &off);
				if (ncols != 1 + nsmpl)
					error("Expected %d columns but %d columns found in the summary file\n",
					      1 + nsmpl, ncols);
				len = strlen(&str.s[off[0]]);
				if (str.s[off[0] + len - 2] != '-'
				    || str.s[off[0] + len - 1] != 'A')
					error("Found Probe Set ID %s while a -A was expected\n",
					      &str.s[off[0]]);
				str.s[off[0] + len - 2] = '\0';
				// check whether the next line contains the expected -B
				// probeset_id
				if (len > MAX_LENGTH_PROBE_SET_ID)
					error("Cannot read Probe Set %s intensities\n",
					      &str.s[off[0]]);
				ret = hpeek(summary_fp->fp.hfile, buf, len);
			} while (ret < len || strncmp(&str.s[off[0]], buf, len - 2) != 0
				 || buf[len - 2] != '-' || buf[len - 1] != 'B');
			if (ret < 0)
				break;

			if (probe_set_id.l == 0)
				kputs(&str.s[off[0]], &probe_set_id);
			else if (strcmp(&str.s[off[0]], probe_set_id.s) != 0)
				error("Found Probe Set ID %s in the summary file while %s expected\n",
				      &str.s[off[0]], probe_set_id.s);
			for (int i = 1; i < 1 + nsmpl; i++)
				norm_x_arr[i - 1] = strtof(&str.s[off[i]], &tmp);
			if (hts_getline(summary_fp, KS_SEP_LINE, &str) <= 0)
				error("Summary file ended prematurely\n");
			ncols = ksplit_core(str.s, '\t', &moff, &off);
			if (ncols != 1 + nsmpl)
				error("Expected %d columns but %d columns found in the summary file\n",
				      1 + nsmpl, ncols);
			len = strlen(&str.s[off[0]]);
			str.s[off[0] + len - 2] = '\0';
			for (int i = 1; i < 1 + nsmpl; i++)
				norm_y_arr[i - 1] = strtof(&str.s[off[i]], &tmp);
			get_delta_size(norm_x_arr, norm_y_arr, nsmpl, delta_arr, size_arr);
		}

		cluster_t cluster;
		// read posteriors
		if (snp_posteriors_fp) {
			if (hts_getline(snp_posteriors_fp, KS_SEP_LINE, &str) < 0)
				break;
			ncols = ksplit_core(str.s, '\t', &moff, &off);
			if (ncols < 4)
				error("Expected %d columns but %d columns found in the SNP posteriors file\n",
				      4, ncols);
			int len = strlen(&str.s[off[0]]);
			if (str.s[off[0] + len - 2] == ':' && str.s[off[0] + len - 1] == '1')
				str.s[off[0] + len - 2] = '\0';
			if (probe_set_id.l == 0)
				kputs(&str.s[off[0]], &probe_set_id);
			else if (strcmp(&str.s[off[0]], probe_set_id.s) != 0)
				error("Found Probe Set ID %s in the SNP posteriors file while %s expected\n",
				      &str.s[off[0]], probe_set_id.s);

			ncols2 = ksplit_core(&str.s[off[1]], ',', &moff2, &off2);
			if (ncols2 < 7)
				error("Missing information for cluster BB for Probe Set %s in the SNP posteriors file\n",
				      &str.s[off[0]]);
			cluster.bb_delta_mean = strtof(&str.s[off[1] + off2[0]], &tmp);
			cluster.bb_delta_dev = strtof(&str.s[off[1] + off2[1]], &tmp);
			cluster.bb_n_mean = strtof(&str.s[off[1] + off2[2]], &tmp);
			cluster.bb_n_dev = strtof(&str.s[off[1] + off2[3]], &tmp);
			cluster.bb_size_mean = strtof(&str.s[off[1] + off2[4]], &tmp);
			cluster.bb_size_dev = strtof(&str.s[off[1] + off2[5]], &tmp);
			cluster.bb_cov = strtof(&str.s[off[1] + off2[6]], &tmp);

			ncols2 = ksplit_core(&str.s[off[2]], ',', &moff2, &off2);
			if (ncols2 < 7)
				error("Missing information for cluster AB for Probe Set %s in the SNP posteriors file\n",
				      &str.s[off[0]]);
			cluster.ab_delta_mean = strtof(&str.s[off[2] + off2[0]], &tmp);
			cluster.ab_delta_dev = strtof(&str.s[off[2] + off2[1]], &tmp);
			cluster.ab_n_mean = strtof(&str.s[off[2] + off2[2]], &tmp);
			cluster.ab_n_dev = strtof(&str.s[off[2] + off2[3]], &tmp);
			cluster.ab_size_mean = strtof(&str.s[off[2] + off2[4]], &tmp);
			cluster.ab_size_dev = strtof(&str.s[off[2] + off2[5]], &tmp);
			cluster.ab_cov = strtof(&str.s[off[2] + off2[6]], &tmp);

			ncols2 = ksplit_core(&str.s[off[3]], ',', &moff2, &off2);
			if (ncols2 < 7)
				error("Missing information for cluster AA for Probe Set %s in the SNP posteriors file\n",
				      &str.s[off[0]]);
			cluster.aa_delta_mean = strtof(&str.s[off[3] + off2[0]], &tmp);
			cluster.aa_delta_dev = strtof(&str.s[off[3] + off2[1]], &tmp);
			cluster.aa_n_mean = strtof(&str.s[off[3] + off2[2]], &tmp);
			cluster.aa_n_dev = strtof(&str.s[off[3] + off2[3]], &tmp);
			cluster.aa_size_mean = strtof(&str.s[off[3] + off2[4]], &tmp);
			cluster.aa_size_dev = strtof(&str.s[off[3] + off2[5]], &tmp);
			cluster.aa_cov = strtof(&str.s[off[3] + off2[6]], &tmp);

			// check whether the next line should be skipped
			if (probe_set_id.l + 2 > MAX_LENGTH_PROBE_SET_ID)
				error("Cannot read Probe Set %s SNP posteriors\n",
				      &str.s[off[0]]);
			int ret = hpeek(snp_posteriors_fp->fp.hfile, buf, probe_set_id.l + 2);
			if (ret == probe_set_id.l + 2
			    && strncmp(probe_set_id.s, buf, probe_set_id.l) == 0
			    && buf[probe_set_id.l] == ':' && buf[probe_set_id.l + 1] == '1')
				hts_getline(snp_posteriors_fp, KS_SEP_LINE, &str);
		}

		int idx;
		if (probe_set_id.l > 0) {
			int ret = khash_str2int_get(annot->probe_set_id, probe_set_id.s, &idx);
			if (ret < 0)
				error("Probe Set %s not found in manifest file\n",
				      probe_set_id.s);
		} else {
			idx = i;
		}
		record_t *record = &annot->records[idx];
		bcf_clear(rec);
		rec->n_sample = nsmpl;
		rec->rid = bcf_hdr_name2id_flexible(hdr, record->chromosome);
		rec->pos = record->position - 1;
		if (rec->rid < 0 || rec->pos < 0 || record->strand < 0 || !record->flank) {
			if (flags & VERBOSE)
				fprintf(stderr, "Skipping unlocalized marker %s\n",
					record->probe_set_id);
			n_skipped++;
			continue;
		}

		if (record->dbsnp_rs_id)
			bcf_update_id(hdr, rec, record->dbsnp_rs_id);
		else if (record->affy_snp_id)
			bcf_update_id(hdr, rec, record->affy_snp_id);
		else
			bcf_update_id(hdr, rec, record->probe_set_id);

		flank.l = 0;
		kputs(record->flank, &flank);
		strupper(flank.s);
		if (record->strand)
			flank_reverse_complement(flank.s);

		int32_t allele_b_idx;
		allele_a.l = allele_b.l = 0;
		if (strchr(flank.s, '-')) {
			int ref_is_del =
				get_indel_alleles(flank.s, fai, bcf_seqname(hdr, rec), rec->pos,
						  0, ref_base, &allele_a, &allele_b);
			if (ref_is_del < 0) {
				if (flags & VERBOSE)
					fprintf(stderr,
						"Unable to determine alleles for indel %s\n",
						record->probe_set_id);
				n_missing++;
			}
			if (ref_is_del == 0)
				rec->pos--;
			allele_b_idx = ref_is_del < 0 ? 1 : ref_is_del;
		} else {
			const char *left = strchr(flank.s, '[');
			const char *middle = strchr(flank.s, '/');
			const char *right = strchr(flank.s, ']');
			if (!left || !middle || !right)
				error("Flank sequence is malformed: %s\n", flank.s);

			kputsn(left + 1, middle - left - 1, &allele_a);
			kputsn(middle + 1, right - middle - 1, &allele_b);
			ref_base[0] = get_ref_base(fai, hdr, rec);
			allele_b_idx = get_allele_b_idx(ref_base[0], allele_a.s, allele_b.s);
		}
		int32_t allele_a_idx = get_allele_a_idx(allele_b_idx);
		const char *alleles[3];
		int nals = alleles_ab_to_vcf(alleles, ref_base, allele_a.s, allele_b.s,
					     allele_b_idx);
		if (nals < 0)
			error("Unable to process Probe Set %s\n", record->probe_set_id);
		bcf_update_alleles(hdr, rec, alleles, nals);
		bcf_update_info_int32(hdr, rec, "ALLELE_A", &allele_a_idx, 1);
		bcf_update_info_int32(hdr, rec, "ALLELE_B", &allele_b_idx, 1);
		bcf_update_info_string(hdr, rec, "PROBE_SET_ID", record->probe_set_id);
		if (record->affy_snp_id)
			bcf_update_info_string(hdr, rec, "AFFY_SNP_ID", record->affy_snp_id);

		if (calls_fp) {
			for (int i = 0; i < nsmpl; i++) {
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
					gt_arr[2 * i] = bcf_gt_unphased(
						min(allele_a_idx, allele_b_idx));
					gt_arr[2 * i + 1] = bcf_gt_unphased(
						max(allele_a_idx, allele_b_idx));
					break;
				case GT_BB:
					gt_arr[2 * i] = bcf_gt_unphased(allele_b_idx);
					gt_arr[2 * i + 1] = bcf_gt_unphased(allele_b_idx);
					break;
				default:
					error("Genotype for Probe Set ID %s is malformed: %d\n",
					      record->probe_set_id, gts[i]);
					break;
				}
			}
			bcf_update_genotypes(hdr, rec, gt_arr, nsmpl * 2);
		}

		if (confidences_fp)
			bcf_update_format_float(hdr, rec, "CONF", conf_arr, nsmpl);

		if (summary_fp) {
			bcf_update_format_float(hdr, rec, "NORMX", norm_x_arr, nsmpl);
			bcf_update_format_float(hdr, rec, "NORMY", norm_y_arr, nsmpl);
			bcf_update_format_float(hdr, rec, "DELTA", delta_arr, nsmpl);
			bcf_update_format_float(hdr, rec, "SIZE", size_arr, nsmpl);
		}

		if (flags & ADJUST_CLUSTERS)
			adjust_clusters(gts, delta_arr, size_arr, nsmpl, &cluster);

		if (snp_posteriors_fp) {
			bcf_update_info_float(hdr, rec, "meanDELTA_AA", &cluster.aa_delta_mean,
					      1);
			bcf_update_info_float(hdr, rec, "meanDELTA_AB", &cluster.ab_delta_mean,
					      1);
			bcf_update_info_float(hdr, rec, "meanDELTA_BB", &cluster.bb_delta_mean,
					      1);
			bcf_update_info_float(hdr, rec, "devDELTA_AA", &cluster.aa_delta_dev,
					      1);
			bcf_update_info_float(hdr, rec, "devDELTA_AB", &cluster.ab_delta_dev,
					      1);
			bcf_update_info_float(hdr, rec, "devDELTA_BB", &cluster.bb_delta_dev,
					      1);
			bcf_update_info_float(hdr, rec, "meanN_AA", &cluster.aa_n_mean, 1);
			bcf_update_info_float(hdr, rec, "meanN_AB", &cluster.ab_n_mean, 1);
			bcf_update_info_float(hdr, rec, "meanN_BB", &cluster.bb_n_mean, 1);
			bcf_update_info_float(hdr, rec, "devN_AA", &cluster.aa_n_dev, 1);
			bcf_update_info_float(hdr, rec, "devN_AB", &cluster.ab_n_dev, 1);
			bcf_update_info_float(hdr, rec, "devN_BB", &cluster.bb_n_dev, 1);
			bcf_update_info_float(hdr, rec, "meanSIZE_AA", &cluster.aa_size_mean,
					      1);
			bcf_update_info_float(hdr, rec, "meanSIZE_AB", &cluster.ab_size_mean,
					      1);
			bcf_update_info_float(hdr, rec, "meanSIZE_BB", &cluster.bb_size_mean,
					      1);
			bcf_update_info_float(hdr, rec, "devSIZE_AA", &cluster.aa_size_dev, 1);
			bcf_update_info_float(hdr, rec, "devSIZE_AB", &cluster.ab_size_dev, 1);
			bcf_update_info_float(hdr, rec, "devSIZE_BB", &cluster.bb_size_dev, 1);
			bcf_update_info_float(hdr, rec, "covar_AA", &cluster.aa_cov, 1);
			bcf_update_info_float(hdr, rec, "covar_AB", &cluster.ab_cov, 1);
			bcf_update_info_float(hdr, rec, "covar_BB", &cluster.bb_cov, 1);
		}

		if (summary_fp && snp_posteriors_fp) {
			get_baf_lrr(delta_arr, size_arr, nsmpl, &cluster, baf_arr, lrr_arr);
			bcf_update_format_float(hdr, rec, "BAF", baf_arr, nsmpl);
			bcf_update_format_float(hdr, rec, "LRR", lrr_arr, nsmpl);
		}

		if (bcf_write(out_fh, hdr, rec) < 0)
			error("Unable to write to output VCF file\n");
	}
	fprintf(stderr, "Lines   total/missing-reference/skipped:\t%d/%d/%d\n", i, n_missing,
		n_skipped);

	free(gts);
	free(gt_arr);
	free(conf_arr);
	free(norm_x_arr);
	free(norm_y_arr);
	free(delta_arr);
	free(size_arr);
	free(baf_arr);
	free(lrr_arr);

	free(probe_set_id.s);
	free(allele_a.s);
	free(allele_b.s);
	free(flank.s);

	bcf_destroy(rec);
	free(off);
	free(off2);
	free(str.s);
	if (calls_fp) {
		if (hgetc(calls_fp->fp.hfile) != EOF)
			fprintf(stderr, "Warning: End of calls file %s was not reached\n",
				calls_fn);
		hts_close(calls_fp);
	}
	if (confidences_fp) {
		if (hgetc(confidences_fp->fp.hfile) != EOF)
			fprintf(stderr, "Warning: End of confidences file %s was not reached\n",
				confidences_fn);
		hts_close(confidences_fp);
	}
	if (summary_fp) {
		if (hgetc(summary_fp->fp.hfile) != EOF)
			fprintf(stderr, "Warning: End of summary file %s was not reached\n",
				summary_fn);
		hts_close(summary_fp);
	}
	if (snp_posteriors_fp) {
		if (hgetc(snp_posteriors_fp->fp.hfile) != EOF)
			fprintf(stderr,
				"Warning: End of SNP posteriors file %s was not reached\n",
				snp_posteriors_fn);
		hts_close(snp_posteriors_fp);
	}
	return;
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void)
{
	return "convert Affymetrix files to VCF.\n";
}

static const char *usage_text(void)
{
	return "\n"
	       "About: convert Affymetrix apt-probeset-genotype output files to VCF. (version " AFFY2VCF_VERSION
	       " https://github.com/freeseek/gtc2vcf)\n"
	       "Usage: bcftools +affy2vcf [options] --csv <file> --fasta-ref <file> --calls <file>\n"
	       "                                    --confidences <file> --summary <file> --snp-posteriors <file>\n"
	       "\n"
	       "Plugin options:\n"
	       "    -c, --csv <file>              CSV manifest file\n"
	       "    -f, --fasta-ref <file>        reference sequence in fasta format\n"
	       "        --set-cache-size <int>    select fasta cache size in bytes\n"
	       "        --calls <file>            apt-probeset-genotype calls output\n"
	       "        --confidences <file>      apt-probeset-genotype confidences output\n"
	       "        --summary <file>          apt-probeset-genotype summary output\n"
	       "        --snp-posteriors <file>   apt-probeset-genotype snp-posteriors output\n"
	       "        --report <file>           apt-probeset-genotype report output\n"
	       "        --adjust-clusters         adjust cluster centers in (Contrast, Size) space (requires --summary and --snp-posteriors)\n"
	       "    -x, --sex <file>              output apt-probeset-genotype gender estimate into file (requires --report)\n"
	       "        --no-version              do not append version and command line to the header\n"
	       "    -o, --output <file>           write output to a file [standard output]\n"
	       "    -O, --output-type <b|u|z|v>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
	       "        --threads <int>           number of extra output compression threads [0]\n"
	       "    -v, --verbose                 print verbose information\n"
	       "\n"
	       "Manifest options:\n"
	       "        --fasta-flank             output flank sequence in FASTA format (requires --csv)\n"
	       "    -s, --sam-flank <file>        input source sequence alignment in SAM/BAM format (requires --csv)\n"
	       "\n"
	       "Example:\n"
	       "    bcftools +affy2vcf \\\n"
	       "        --csv GenomeWideSNP_6.na35.annot.csv \\\n"
	       "        --fasta-ref human_g1k_v37.fasta \\\n"
	       "        --calls AxiomGT1.calls.txt \\\n"
	       "        --confidences AxiomGT1.confidences.txt \\\n"
	       "        --summary AxiomGT1.summary.txt \\\n"
	       "        --snp-posteriors AxiomGT1.snp-posteriors.txt \\\n"
	       "        --output AxiomGT1.vcf\n"
	       "\n"
	       "Examples of manifest file options:\n"
	       "    bcftools +affy2vcf -c GenomeWideSNP_6.na35.annot.csv --fasta-flank -o GenomeWideSNP_6.fasta\n"
	       "    bwa mem -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GenomeWideSNP_6.fasta -o GenomeWideSNP_6.sam\n"
	       "    bcftools +affy2vcf -c GenomeWideSNP_6.na35.annot.csv -s GenomeWideSNP_6.sam -o GenomeWideSNP_6.na35.annot.GRCh38.csv\n"
	       "\n";
}

int run(int argc, char *argv[])
{
	const char *ref_fname = NULL;
	const char *sex_fname = NULL;
	const char *csv_fname = NULL;
	const char *snp_posteriors_fname = NULL;
	const char *summary_fname = NULL;
	const char *report_fname = NULL;
	const char *calls_fname = NULL;
	const char *confidences_fname = NULL;
	const char *output_fname = "-";
	const char *sam_fname = NULL;
	int flags = 0;
	int output_type = FT_VCF;
	int cache_size = 0;
	int n_threads = 0;
	int record_cmd_line = 1;
	int fasta_flank = 0;
	faidx_t *fai = NULL;

	static struct option loptions[] = {{"csv", required_argument, NULL, 'c'},
					   {"fasta-ref", required_argument, NULL, 'f'},
					   {"set-cache-size", required_argument, NULL, 1},
					   {"calls", required_argument, NULL, 2},
					   {"confidences", required_argument, NULL, 3},
					   {"summary", required_argument, NULL, 4},
					   {"snp-posteriors", required_argument, NULL, 5},
					   {"report", required_argument, NULL, 6},
					   {"adjust-clusters", no_argument, NULL, 7},
					   {"sex", required_argument, NULL, 'x'},
					   {"no-version", no_argument, NULL, 8},
					   {"output", required_argument, NULL, 'o'},
					   {"output-type", required_argument, NULL, 'O'},
					   {"threads", required_argument, NULL, 9},
					   {"verbose", no_argument, NULL, 'v'},
					   {"fasta-flank", no_argument, NULL, 10},
					   {"sam-flank", required_argument, NULL, 's'},
					   {NULL, 0, NULL, 0}};
	int c;
	while ((c = getopt_long(argc, argv, "h?c:f:x:o:O:vs:", loptions, NULL)) >= 0) {
		switch (c) {
		case 'c':
			csv_fname = optarg;
			break;
		case 'f':
			ref_fname = optarg;
			break;
		case 1:
			cache_size = strtol(optarg, NULL, 0);
			break;
		case 2:
			calls_fname = optarg;
			flags |= CALLS_LOADED;
			break;
		case 3:
			confidences_fname = optarg;
			flags |= CONFIDENCES_LOADED;
			break;
		case 4:
			summary_fname = optarg;
			flags |= SUMMARY_LOADED;
			break;
		case 5:
			snp_posteriors_fname = optarg;
			flags |= SNP_POSTERIORS_LOADED;
			break;
		case 6:
			report_fname = optarg;
			break;
		case 7:
			flags |= ADJUST_CLUSTERS;
			break;
		case 'x':
			sex_fname = optarg;
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
			default:
				error("The output type \"%s\" not recognised\n", optarg);
			}
			break;
		case 9:
			n_threads = strtol(optarg, NULL, 0);
			break;
		case 'v':
			flags |= VERBOSE;
			break;
		case 10:
			fasta_flank = 1;
			break;
		case 's':
			sam_fname = optarg;
			break;
		case 'h':
		case '?':
		default:
			error("%s", usage_text());
		}
	}
	if (!csv_fname)
		error("Expected --csv option\n%s", usage_text());
	if (fasta_flank && sam_fname)
		error("Only one of --fasta-flank or --sam-flank options can be used at once\n%s",
		      usage_text());
	if (!fasta_flank && !sam_fname && !ref_fname)
		error("Expected one of --fasta-flank or --sam-flank or --fasta-ref options\n%s",
		      usage_text());
	if ((flags & ADJUST_CLUSTERS) && (!summary_fname || !snp_posteriors_fname))
		error("Expected --summary and --snp-posteriors options with --adjust-clusters option\n%s",
		      usage_text());
	if (sex_fname && !report_fname)
		error("Expected --report option with --sex option\n%s", usage_text());

	fprintf(stderr,
		"================================================================================\n");
	fprintf(stderr, "Reading CSV file %s\n", csv_fname);
	annot_t *annot = annot_init(
		csv_fname, sam_fname,
		((sam_fname && !ref_fname) || fasta_flank) ? output_fname : NULL, flags);

	if (annot) {
		fai = fai_load(ref_fname);
		if (!fai)
			error("Could not load the reference %s\n", ref_fname);
		if (cache_size)
			fai_set_cache_size(fai, cache_size);
		fprintf(stderr,
			"================================================================================\n");
		fprintf(stderr, "Writing VCF file\n");
		bcf_hdr_t *hdr = hdr_init(fai, flags);
		bcf_hdr_printf(hdr, "##CSV=%s",
			       strrchr(csv_fname, '/') ? strrchr(csv_fname, '/') + 1
						       : csv_fname);
		if (sam_fname)
			bcf_hdr_printf(hdr, "##SAM=%s",
				       strrchr(sam_fname, '/') ? strrchr(sam_fname, '/') + 1
							       : sam_fname);
		if (record_cmd_line)
			bcf_hdr_append_version(hdr, argc, argv, "bcftools_+affy2vcf");
		htsFile *out_fh = hts_open(output_fname, hts_bcf_wmode(output_type));
		if (out_fh == NULL)
			error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
		if (n_threads)
			hts_set_threads(out_fh, n_threads);
		process(fai, annot, calls_fname, confidences_fname, summary_fname,
			snp_posteriors_fname, out_fh, hdr, flags);
		fai_destroy(fai);
		bcf_hdr_destroy(hdr);
		hts_close(out_fh);
		annot_destroy(annot);
	}

	if (sex_fname) {
		fprintf(stderr, "Reading apt-probeset-genotype report file %s\n", report_fname);
		report_t *report = report_init(report_fname);
		FILE *sex_fh = fopen(sex_fname, "w");
		if (!sex_fh)
			error("Failed to open %s: %s\n", sex_fname, strerror(errno));
		for (int i = 0; i < report->n_samples; i++)
			fprintf(sex_fh, "%s\t%d\n", report->cel_files[i], report->genders[i]);
		fclose(sex_fh);
		report_destroy(report);
	}

	return 0;
}
