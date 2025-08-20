/* The MIT License

   Copyright (C) 2024-2025 Giulio Genovese

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

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/ksort.h>
#include "bcftools.h"

#define BAFREGRESS_VERSION "2025-08-19"

#define GT_NC 0
#define GT_AA 1
#define GT_AB 2
#define GT_BB 3

KSORT_INIT_GENERIC(float)

/******************************************
 * PLUGIN                                 *
 ******************************************/

inline static double sqr(double x) { return x * x; }

const char *about(void) { return "Detects and estimates sample contamination using BAF intensity data.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: Detects and estimates sample contamination. (version " BAFREGRESS_VERSION
           " http://github.com/freeseek/gtc2vcf)\n"
           "[ Jun, G. et al. Detecting and Estimating Contamination of Human DNA Samples in Sequencing\n"
           "and Array-Based Genotype Data. AJHG 91, 839-848 (2012) http://doi.org/10.1016/j.ajhg.2012.09.004 ]\n"
           "\n"
           "Usage: bcftools +BAFregress [options] <in.vcf.gz>\n"
           "\n"
           "Plugin options:\n"
           "        --threshold <float>         minimum allele frequency for BAF regression [0.1]\n"
           "    -a, --af <file>                 file with allele frequency information\n"
           "        --tag <string>              allele frequency INFO tag [AC/AN]\n"
           "        --adjust-BAF                minimum number of genotypes for a cluster to median adjust BAF (-1 for "
           "no adjustment) [5]\n"
           "        --truncate-BAF              truncates BAF values between 0 and 1 and turns off adjustment to "
           "recover original behavior\n"
           "        --use-MAF                   uses minor allele frequency rather than A/B allele frequency to "
           "recover original behavior\n"
           "    -e, --estimates <file>          write BAF regression estimates to a file [standard output]\n"
           "    -o, --output <file>             write VCF output to a file\n"
           "    -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "    -r, --regions <region>          restrict to comma-separated list of regions\n"
           "    -R, --regions-file <file>       restrict to regions listed in a file\n"
           "        --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant "
           "overlaps (2) [1]\n"
           "    -t, --targets [^]<region>       similar to -r but streams rather than index-jumps. Exclude regions "
           "with \"^\" prefix\n"
           "    -T, --targets-file [^]<file>    similar to -R but streams rather than index-jumps. Exclude regions "
           "with \"^\" prefix\n"
           "        --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant "
           "overlaps (2) [0]\n"
           "        --threads <int>             number of extra output compression threads [0]\n"
           "    -s, --samples [^]<list>         comma separated list of samples to include (or exclude with \"^\" "
           "prefix)\n"
           "    -S, --samples-file [^]<file>    file of samples to include (or exclude with \"^\" prefix)\n"
           "        --force-samples             only warn about unknown subset samples\n"
           "    -W, --write-index[=FMT]         Automatically index the output files [off]\n"
           "\n"
           "Example:\n"
           "    bcftools +BAFregress file.bcf\n"
           "    bcftools +BAFregress --tag AF file.bcf\n"
           "    bcftools +BAFregress --af 1kGP_high_coverage_Illumina.sites.bcf file.bcf\n"
           "    bcftools +BAFregress --af 1kGP_high_coverage_Illumina.sites.bcf --truncate-BAF --use-MAF file.bcf\n"
           "\n";
}

int run(int argc, char **argv) {
    float af_threshold = 0.1;
    char *af_fname = NULL;
    char *af_tag = NULL;
    int adj_baf = 5;
    int truncate_baf = 0;
    int use_maf = 0;
    char *estimate_fname = "-";
    char *output_fname = NULL;
    int output_type = FT_VCF;
    int clevel = -1;
    int regions_overlap = 1;
    int targets_overlap = 0;
    int n_threads = 0;
    char *targets_list = NULL;
    int targets_is_file = 0;
    char *regions_list = NULL;
    int regions_is_file = 0;
    char *sample_names = NULL;
    int sample_is_file = 0;
    int force_samples = 0;
    int write_index = 0;
    char *index_fname;
    htsFile *out_fh = NULL;

    static struct option loptions[] = {{"threshold", required_argument, NULL, 1},
                                       {"af", required_argument, NULL, 'a'},
                                       {"tag", required_argument, NULL, 2},
                                       {"adjust-BAF", required_argument, NULL, 3},
                                       {"truncate-BAF", no_argument, NULL, 4},
                                       {"use-MAF", no_argument, NULL, 5},
                                       {"estimates", required_argument, NULL, 'e'},
                                       {"output", required_argument, NULL, 'o'},
                                       {"output-type", required_argument, NULL, 'O'},
                                       {"threads", required_argument, NULL, 6},
                                       {"regions", required_argument, NULL, 'r'},
                                       {"regions-file", required_argument, NULL, 'R'},
                                       {"regions-overlap", required_argument, NULL, 7},
                                       {"targets", required_argument, NULL, 't'},
                                       {"targets-file", required_argument, NULL, 'T'},
                                       {"targets-overlap", required_argument, NULL, 8},
                                       {"samples", required_argument, NULL, 's'},
                                       {"samples-file", required_argument, NULL, 'S'},
                                       {"force-samples", no_argument, NULL, 9},
                                       {"write-index", optional_argument, NULL, 'W'},
                                       {0, 0, 0, 0}};
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "h?a:e:o:O:r:R:t:T:s:S:", loptions, NULL)) >= 0) {
        switch (c) {
        case 1:
            af_threshold = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --threshold %s\n", optarg);
            if (af_threshold <= 0.0 || af_threshold >= 1.0) error("--threshold must input a value between 0 and 1\n");
            break;
        case 'a':
            af_fname = optarg;
            break;
        case 2:
            af_tag = optarg;
            break;
        case 3:
            adj_baf = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --adjust-BAF %s\n", optarg);
            break;
        case 4:
            truncate_baf = 1;
            break;
        case 5:
            use_maf = 1;
            break;
        case 'e':
            estimate_fname = optarg;
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
            default: {
                clevel = strtol(optarg, &tmp, 10);
                if (*tmp || clevel < 0 || clevel > 9) error("The output type \"%s\" not recognised\n", optarg);
            }
            };
            if (optarg[1]) {
                clevel = strtol(optarg + 1, &tmp, 10);
                if (*tmp || clevel < 0 || clevel > 9)
                    error("Could not parse argument: --compression-level %s\n", optarg + 1);
            }
            break;
        case 6:
            n_threads = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse argument: --threads %s\n", optarg);
            break;
        case 'r':
            regions_list = optarg;
            break;
        case 'R':
            regions_list = optarg;
            regions_is_file = 1;
            break;
        case 7:
            if (!strcasecmp(optarg, "0"))
                regions_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                regions_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                regions_overlap = 2;
            else
                error("Could not parse: --regions-overlap %s\n", optarg);
            break;
        case 't':
            targets_list = optarg;
            break;
        case 'T':
            targets_list = optarg;
            targets_is_file = 1;
            break;
        case 8:
            if (!strcasecmp(optarg, "0"))
                targets_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                targets_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                targets_overlap = 2;
            else
                error("Could not parse: --targets-overlap %s\n", optarg);
            break;
        case 's':
            sample_names = optarg;
            break;
        case 'S':
            sample_names = optarg;
            sample_is_file = 1;
            break;
        case 9:
            force_samples = 1;
            break;
        case 'W':
            if (!(write_index = write_index_parse(optarg))) error("Unsupported index format '%s'\n", optarg);
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
            break;
        }
    }

    if (truncate_baf) adj_baf = -1;

    char *input_fname = NULL;
    if (optind == argc) {
        if (!isatty(fileno((FILE *)stdin))) {
            input_fname = "-"; // reading from stdin
        } else {
            error("%s", usage_text());
        }
    } else if (optind + 1 != argc) {
        error("%s", usage_text());
    } else {
        input_fname = argv[optind];
    }

    bcf_srs_t *srs = bcf_sr_init();
    if (af_fname) {
        bcf_sr_set_opt(srs, BCF_SR_REQUIRE_IDX);
        bcf_sr_set_opt(srs, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);
    }

    if (regions_list) {
        bcf_sr_set_opt(srs, BCF_SR_REGIONS_OVERLAP, regions_overlap);
        if (bcf_sr_set_regions(srs, regions_list, regions_is_file) < 0)
            error("Failed to read the regions: %s\n", regions_list);
    }
    if (targets_list) {
        bcf_sr_set_opt(srs, BCF_SR_TARGETS_OVERLAP, targets_overlap);
        if (bcf_sr_set_targets(srs, targets_list, targets_is_file, 0) < 0)
            error("Failed to read the targets: %s\n", targets_list);
    }
    if (bcf_sr_set_threads(srs, n_threads) < 0) error("Failed to create threads\n");
    if (!bcf_sr_add_reader(srs, input_fname))
        error("Failed to open %s: %s\n", input_fname, bcf_sr_strerror(srs->errnum));
    if (af_fname && !bcf_sr_add_reader(srs, af_fname))
        error("Failed to open %s: %s\n", af_fname, bcf_sr_strerror(srs->errnum));

    bcf_hdr_t *hdr = bcf_sr_get_header(srs, 0);
    bcf_hdr_t *af_hdr = af_fname ? bcf_sr_get_header(srs, 1) : NULL;

    if (sample_names) {
        int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
        if (ret < 0)
            error("Error parsing the list of samples: %s\n", sample_names);
        else if (force_samples && ret > 0)
            error("Sample name mismatch: sample #%d not found in the header\n", ret);
    }

    // get IDs for all VCF formats
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    if (gt_id < 0) error("Format GT was not found in the input header\n");
    int baf_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "BAF");
    if (baf_id < 0) error("Format BAF was not found in the input header\n");
    int allele_a_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ALLELE_A");
    if (allele_a_id < 0) error("Format ALLELE_A was not found in the input header\n");
    int allele_b_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ALLELE_B");
    if (allele_b_id < 0) error("Format ALLELE_B was not found in the input header\n");
    int af_id = -1;
    if (af_tag) {
        af_id = bcf_hdr_id2int(af_hdr ? af_hdr : hdr, BCF_DT_ID, af_tag);
        if (af_id < 0) error("Format %s was not found in the allele frequency header\n", af_tag);
    }

    FILE *est_fh = strcmp("-", estimate_fname) ? fopen(estimate_fname, "w") : stdout;
    if (!est_fh) error("Error: cannot write to %s\n", estimate_fname);

    // output VCF
    if (output_fname) {
        char wmode[8];
        set_wmode(wmode, output_type, output_fname, clevel);
        out_fh = hts_open(output_fname, wmode);
        if (out_fh == NULL) error("[%s] Error: cannot write to \"%s\": %s\n", __func__, output_fname, strerror(errno));
        if (n_threads) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, srs->p);
        if (bcf_hdr_write(out_fh, hdr) < 0) error("Unable to write to output VCF file\n");
        if (init_index2(out_fh, hdr, output_fname, &index_fname, write_index) < 0)
            error("Error: failed to initialise index for %s\n", output_fname);
    }

    int n_smpls = bcf_hdr_nsamples(hdr);
    if (!af_hdr && !af_tag && n_smpls < 30)
        fprintf(
            stderr,
            "Input VCF only includes %d samples. We recommend using a separate VCF to infer marker allele frequency\n",
            n_smpls);

    int *arr = NULL;
    int marr = 0;
    float *baf_arr = NULL;
    int nbaf_arr = 0;
    int8_t *gts = (int8_t *)calloc(n_smpls, sizeof(int8_t));
    float *tmp_arr = (float *)calloc(n_smpls, sizeof(float));
    float *sumx2 = (float *)calloc(n_smpls, sizeof(float));
    float *sumxy = (float *)calloc(n_smpls, sizeof(float));
    float *sumx = (float *)calloc(n_smpls, sizeof(float));
    float *sumy = (float *)calloc(n_smpls, sizeof(float));
    int *n = (int *)calloc(n_smpls, sizeof(int));

    // run through each record present in both VCFs
    int i, j;
    while (bcf_sr_next_line(srs)) {
        bcf1_t *line = bcf_sr_get_line(srs, 0);
        if (!line) continue;
        if (out_fh && bcf_write1(out_fh, hdr, line) != 0)
            error("[%s] Error: cannot write to %s\n", __func__, output_fname);

        bcf1_t *af_line = af_hdr ? bcf_sr_get_line(srs, 1) : line;
        if (line->n_allele != 2 || !af_line || af_line->n_allele != 2) continue;

        // skip lines where the allele frequency is less than 0.01 (or greater than 0.99)
        double af;
        if (af_tag) {
            bcf_info_t *af_info = bcf_get_info_id(af_line, af_id);
            af = af_info ? (double)af_info->v1.f : NAN;
        } else {
            hts_expand(int, af_line->n_allele, marr, arr);
            int ret = bcf_calc_ac(af_hdr ? af_hdr : hdr, af_line, arr, BCF_UN_INFO | BCF_UN_FMT);
            if (ret <= 0) continue;
            int an = 0;
            for (i = 0; i < af_line->n_allele; i++) an += arr[i];
            af = (double)arr[1] / (double)an;
        }
        if (isnan(af) || af < af_threshold || af > 1.0 - af_threshold) continue;
        if (use_maf && af > 0.5) af = 1.0 - af; // uses MAF instead of AF to avoid problems with flipped Illumina probes

        // skip lines where ALLELE_A and ALLELE_B refer to alleles missing from the record (it should not happen)
        bcf_info_t *allele_a_info = bcf_get_info_id(line, allele_a_id);
        int8_t allele_a = allele_a_info ? (int8_t)allele_a_info->v1.i : bcf_int8_missing;
        bcf_info_t *allele_b_info = bcf_get_info_id(line, allele_b_id);
        int8_t allele_b = allele_b_info ? (int8_t)allele_b_info->v1.i : bcf_int8_missing;
        if (allele_a < 0 || allele_a >= line->n_allele || allele_b < 0 || allele_b >= line->n_allele) continue;
        if (allele_b == 0) af = 1.0 - af; // flip the allele frequency if ALLELE_B is the reference

        // skip lines missing genotypes (e.g. intensity only sites) or with ploidy other than 2
        int n_aa = 0, n_ab = 0, n_bb = 0;
        bcf_fmt_t *gt_fmt = bcf_get_fmt_id(line, gt_id);
        if (!gt_fmt || gt_fmt->n != 2) continue;
#define BRANCH(type_t, bcf_type_vector_end)                                                                            \
    {                                                                                                                  \
        type_t *p = (type_t *)gt_fmt->p;                                                                               \
        for (i = 0; i < n_smpls; i++, p += 2) {                                                                        \
            gts[i] = GT_NC;                                                                                            \
            if (p[0] == bcf_type_vector_end || bcf_gt_is_missing(p[0]) || p[1] == bcf_type_vector_end                  \
                || bcf_gt_is_missing(p[1]))                                                                            \
                continue;                                                                                              \
            type_t allele_0 = bcf_gt_allele(p[0]);                                                                     \
            type_t allele_1 = bcf_gt_allele(p[1]);                                                                     \
            if (allele_0 == allele_a && allele_1 == allele_a) {                                                        \
                gts[i] = GT_AA;                                                                                        \
                n_aa++;                                                                                                \
            } else if ((allele_0 == allele_a && allele_1 == allele_b)                                                  \
                       || (allele_0 == allele_b && allele_1 == allele_a)) {                                            \
                gts[i] = GT_AB;                                                                                        \
                n_ab++;                                                                                                \
            } else if (allele_0 == allele_b && allele_1 == allele_b) {                                                 \
                gts[i] = GT_BB;                                                                                        \
                n_bb++;                                                                                                \
            }                                                                                                          \
        }                                                                                                              \
    }
        switch (gt_fmt->type) {
        case BCF_BT_INT8:
            BRANCH(int8_t, bcf_int8_vector_end);
            break;
        case BCF_BT_INT16:
            BRANCH(int16_t, bcf_int16_vector_end);
            break;
        case BCF_BT_INT32:
            BRANCH(int32_t, bcf_int32_vector_end);
            break;
        default:
            error("Unexpected type %d\n", gt_fmt->type);
        }
#undef BRANCH

        int nbaf = bcf_get_format_float(hdr, line, "BAF", &baf_arr, &nbaf_arr);
        if (nbaf != n_smpls) continue; // wrong number of BAF values

        // adjust BAF
        float adj_baf_aa = 0.0;
        float adj_baf_bb = 0.0;
        if (adj_baf != -1) {
            j = 0;
            if (n_aa >= adj_baf) {
                for (i = 0; i < n_smpls; i++)
                    if (gts[i] == GT_AA) tmp_arr[j++] = baf_arr[i];
                adj_baf_aa = ks_ksmall_float((size_t)j, tmp_arr, (size_t)j / 2);
                if (j % 2 == 0) adj_baf_aa = (adj_baf_aa + tmp_arr[j / 2 - 1]) * 0.5f;
            }
            j = 0;
            if (n_bb >= adj_baf) {
                for (i = 0; i < n_smpls; i++)
                    if (gts[i] == GT_BB) tmp_arr[j++] = baf_arr[i];
                adj_baf_bb = ks_ksmall_float((size_t)j, tmp_arr, (size_t)j / 2);
                if (j % 2 == 0) adj_baf_bb = (adj_baf_bb + tmp_arr[j / 2 - 1]) * 0.5f;
                adj_baf_bb -= 1.0;
            }
        } else if (truncate_baf) { // truncates the BAF between 0.0 and 1.0 like Illumina does
            for (i = 0; i < n_smpls; i++) {
                if (baf_arr[i] < 0.0)
                    baf_arr[i] = 0.0;
                else if (baf_arr[i] > 1.0)
                    baf_arr[i] = 1.0;
            }
        }

        for (i = 0; i < n_smpls; i++) {
            double baf;
            if (gts[i] == GT_AA) {
                baf = (double)(baf_arr[i] - adj_baf_aa);
                sumx2[i] += sqr(af);
                sumxy[i] += af * baf;
                sumx[i] += af;
                sumy[i] += baf;
            } else if (gts[i] == GT_BB) {
                baf = (double)(baf_arr[i] - adj_baf_bb);
                sumx2[i] += sqr(1.0 - af);
                sumxy[i] += (1.0 - af) * (1.0 - baf);
                sumx[i] += 1.0 - af;
                sumy[i] += 1.0 - baf;
            } else
                continue;
            n[i]++;
        }
    }

    fprintf(est_fh, "sample_id\tbaf_regress\tNhom\n");
    for (i = 0; i < n_smpls; i++) {
        double denom = (double)n[i] * sumx2[i] - sqr(sumx[i]);
        double m = denom ? (n[i] * sumxy[i] - sumx[i] * sumy[i]) / denom : NAN;
        // double b = denom ? (sumy[i] * sumx2[i] - sumx[i] * sumxy[i]) / denom : NAN;
        fprintf(est_fh, "%s\t%.4f\t%d\n", hdr->samples[i], m, n[i]);
    }

    if (est_fh != stdout && est_fh != stderr) fclose(est_fh);

    // close output VCF
    if (output_fname) {
        if (write_index) {
            if (bcf_idx_save(out_fh) < 0) {
                if (hts_close(out_fh) != 0)
                    error("Close failed %s\n", strcmp(output_fname, "-") ? output_fname : "stdout");
                error("Error: cannot write to index %s\n", index_fname);
            }
            free(index_fname);
        }
        hts_close(out_fh);
    }

    free(arr);
    free(baf_arr);
    free(gts);
    free(tmp_arr);
    free(sumx2);
    free(sumxy);
    free(sumx);
    free(sumy);
    free(n);
    bcf_sr_destroy(srs);

    return 0;
}
