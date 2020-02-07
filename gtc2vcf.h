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

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

static inline FILE *get_file_handle(const char *str)
{
    if ( !str ) return NULL;
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

static inline void flank2fasta(const char *name,
                               const char *flank,
                               FILE *stream)
{
    if ( !flank ) return;
    const char *left = strchr(flank, '[');
    const char *middle = strchr(flank, '/');
    const char *right = strchr(flank, ']');
    fprintf(stream, "@%s:1\n", name);
    if ( !left && !middle && !right )
    {
        fprintf(stream, "%s\n", flank);
        return;
    }
    if ( !left || !middle || !right ) error("Flank sequence is malformed: %s\n", flank);
    if ( *(middle-1) == '-' ) fprintf(stream, "%.*s%s\n", (int)(left - flank), flank, right + 1);
    else fprintf(stream, "%.*s%.*s%s\n", (int)(left - flank), flank, (int)(middle - left) - 1, left + 1, right + 1);
    fprintf(stream, "@%s:2\n", name);
    if ( *(middle-1) == '-' ) fprintf(stream, "%.*s%.*s%s\n", (int)(left - flank), flank, (int)(right - middle) - 1, middle + 1, right + 1);
    else fprintf(stream, "%.*s%.*s%s\n", (int)(left - flank), flank, (int)(right - middle) - 1, middle + 1, right + 1);
}

static inline int bcf_hdr_name2id_flexible(const bcf_hdr_t *hdr, char *chr)
{
    if (!chr) return -1;
    char buf[] = {'c', 'h', 'r', '\0', '\0', '\0'};
    int rid = bcf_hdr_name2id(hdr, chr);
    if ( rid >= 0 ) return rid;
    if ( strlen(chr) > 2 ) return -1;
    strcpy(buf + 3, chr);
    rid = bcf_hdr_name2id(hdr, buf);
    if ( rid >= 0 ) return rid;
    else if ( strcmp(chr, "X") == 0 || strcmp(chr, "XY") == 0 || strcmp(chr, "XX") == 0 )
    {
        rid = bcf_hdr_name2id(hdr, "X");
        if ( rid >= 0 ) return rid;
        rid = bcf_hdr_name2id(hdr, "chrX");
    }
    else if ( strcmp(chr, "Y") == 0 )
    {
        rid = bcf_hdr_name2id(hdr, "chrY");
    }
    else if ( strcmp(chr, "MT") == 0 )
    {
        rid = bcf_hdr_name2id(hdr, "chrM");
    }
    return rid;
}

static inline char revnt(char iupac)
{
    static const char iupac_complement[128] = {
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,'/',
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
          0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,']',  0,'[',  0,  0,
          0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
          0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
    };
    if ( iupac > 127 ) return 0;
    return iupac_complement[(int)iupac];
}

#define MAX_LENGTH_LEFT_ALLELE 8
static inline void flank_reverse_complement(char *flank)
{
    // swap alleles, but only if first allele is one base pair long
    char *left = strchr(flank, '[');
    char *middle = strchr(flank, '/');
    char *right = strchr(flank, ']');
    if ( !left || !middle || !right ) error("Flank sequence is malformed: %s\n", flank);

    char buf[MAX_LENGTH_LEFT_ALLELE];
    if ( middle - left - 1 > MAX_LENGTH_LEFT_ALLELE ) error("Cannot swap alleles in flank sequence %s\n", flank);
    memmove((void *)buf, left + 1, middle - left - 1);
    memmove((void *)left + 1, middle + 1, right - middle - 1);
    *(left + (right - middle)) = '/';
    memmove(left + (right - middle) + 1, (void *)buf, middle - left - 1);

    size_t len = strlen(flank);
    for (size_t i=0; i<len/2; i++)
    {
        char tmp = flank[i];
        flank[i] = revnt( flank[len-i-1] );
        flank[len-i-1] = revnt( tmp );
    }
    if (len % 2 == 1) flank[len/2] = revnt( flank[len/2] );
}

// this is the weird way Illumina left shifts indels (see https://github.com/Illumina/GTCtoVCF/blob/develop/BPMRecord.py)
static inline void flank_left_shift(char *flank)
{
    char *left = strchr(flank, '[');
    char *middle = strchr(flank, '/');
    char *right = strchr(flank, ']');
    if ( !left || !middle || !right ) error("Flank sequence is malformed: %s\n", flank);

    int len = (int)(right - middle) - 1;
    while( ( left - flank >= len ) && ( strncmp( left - len, middle + 1, len ) == 0 ) )
    {
        memmove(left - len, left, right - left + 1);
        left -= len;
        middle -= len;
        right -= len;
        memmove(right + 1, middle + 1, len);
    }

    char nt = *(middle + 1);
    for (const char *ptr = middle + 2; ptr < right; ptr++) if (*ptr != nt) nt = -1;
    while ( nt > 0 && *(left - 1) == nt )
    {
        memmove(left - 1, left, right - left + 1);
        *right = nt;
        left--;
        middle--;
        right--;
    }
}

// returns 1 if the first sequence is the best alignment, and 2 if the second sequence is
// if neither sequence is better or neither provides an alignment, it returns 0
// if it fails to read from the hts file, it returns -1
static inline int get_position(htsFile *hts,
                               sam_hdr_t *sam_hdr,
                               bam1_t *b,
                               const char *name,
                               const char *flank,
                               int left_shift,
                               const char **chromosome,
                               int *position,
                               int *strand)
{
    const char *left = strchr(flank, '[');
    const char *middle = strchr(flank, '/');
    const char *right = strchr(flank, ']');
    int cnv = !left && !middle && !right;
    if ( !cnv && ( !left || !middle || !right ) ) error("Flank sequence is malformed: %s\n", flank);
    const char *chromosome_pair[2];
    int position_pair[2], strand_pair[2];
    int64_t aln_score_pair[2];
    int idx = -1, ret;
    while ( idx < 1 - cnv && ( ret = sam_read1(hts, sam_hdr, b) ) >= 0 )
    {
        const char *qname = bam_get_qname(b);
        if ( b->core.flag & BAM_FSECONDARY || b->core.flag & BAM_FSUPPLEMENTARY ) continue;
        int qname_l = strlen(qname);
        if ( strncmp( qname, name, qname_l - 2 ) != 0 )
            error("Query ID %.*s found in SAM file but %s expected\n", qname_l - 2, qname, name );
        idx = qname[qname_l - 1] == '1' ? 0 : ( qname[qname_l - 1] == '2' ? 1 : -1 );
        if ( idx < 0 ) error("Query ID %s found in SAM file does not end with :1 or :2\n", qname);

        chromosome_pair[idx] = sam_hdr_tid2name( sam_hdr, b->core.tid );
        position_pair[idx] = 0;
        strand_pair[idx] = -1;
        if ( !(b->core.flag & BAM_FUNMAP) )
        {
            strand_pair[idx] = bam_is_rev(b);
            int n_cigar = b->core.n_cigar;
            const uint32_t *cigar = bam_get_cigar(b);
            position_pair[idx] = b->core.pos;

            int qlen = cnv ? (strlen(flank) + 1) / 2 : ( bam_is_rev(b) ? strlen(flank) - (right - flank) : left - flank + 1 );
            if ( strchr( flank, '-' ) )
            {
                if ( left_shift )
                {
                    int len = (int)(right - middle) - 1;
                    char nt = toupper(*(middle + 1));
                    const char *ptr;
                    for (ptr = middle + 2; ptr < right; ptr++) if (*ptr != nt) nt = -1;
                    if ( bam_is_rev(b) )
                    {
                        ptr = right + 1;
                        while( strncasecmp( middle + 1, ptr, len ) == 0 ) { qlen -= len; ptr += len; }
                        while ( nt > 0 && toupper(*ptr) == nt ) { qlen--; ptr++; }
                    }
                    else
                    {
                        ptr = left - len;
                        while( ptr >= flank && ( strncasecmp( ptr, middle + 1, len ) == 0 ) ) { qlen -= len; ptr -= len; }
                        ptr += len - 1;
                        while ( nt > 0 && toupper(*ptr) == nt ) { qlen--; ptr--; }
                    }
                }
                if ( idx == 0 ) qlen--;
            }

            for (int k = 0; k < n_cigar && qlen > 1; k++)
            {
                int type = bam_cigar_type(bam_cigar_op(cigar[k]));
                int len = bam_cigar_oplen(cigar[k]);
                if ( ( type & 1 ) && ( type & 2 ) ) // consume reference sequence ( case M )
                {
                    position_pair[idx] += min(len, qlen);
                    qlen -= len;
                }
                else if ( type & 1 ) // consume query sequence ( case I )
                {
                    qlen -= len;
                    if ( qlen <= 0 ) // we skipped the base pair that needed to be localized
                    {
                        position_pair[idx] = 0;
                    }
                }
                else if ( type & 2 )
                {
                    position_pair[idx] += len; // consume reference sequence ( case D )
                }
            }
            if ( qlen == 1 ) position_pair[idx]++;
        }
        uint8_t *as = bam_aux_get( b, "AS" );
        aln_score_pair[idx] = bam_aux2i(as);
    }
    if ( ret < -1 ) return -1;

    if ( !cnv && ( ( aln_score_pair[0] == aln_score_pair[1] && position_pair[0] != position_pair[1] ) || ( position_pair[0] == 0 && position_pair[1] == 0 ) ) )
    {
        idx = -1;
        *chromosome = NULL;
        *position = 0;
        *strand = -1;
    }
    else
    {
        idx = cnv ? 0 : ( aln_score_pair[1] > aln_score_pair[0] );
        *chromosome = chromosome_pair[idx];
        *position = position_pair[idx];
        *strand = strand_pair[idx];
    }
    return idx + 1;
}

static inline char get_ref_base(faidx_t *fai,
                                const bcf_hdr_t *hdr,
                                bcf1_t *rec)
{
    int len;
    char *ref = faidx_fetch_seq(fai, bcf_seqname(hdr, rec), rec->pos, rec->pos, &len);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%"PRId64"\n", bcf_seqname(hdr, rec), rec->pos + 1);
    char ref_base = ref[0];
    free(ref);
    return ref_base;
}

static inline void strupper(char *str) {
    char *s = str;
    while (*s) {
        *s = toupper((unsigned char) *s);
        s++;
     }
}

static inline int len_common_suffix(const char *s1,
                                    const char *s2,
                                    size_t n)
{
    int ret = 0;
    while (ret < n && *s1 == *s2)
    {
        s1--;
        s2--;
        ret++;
    }
    return ret;
}

static inline int len_common_prefix(const char *s1,
                                    const char *s2,
                                    size_t n)
{
    int ret = 0;
    while (ret < n && *s1 == *s2)
    {
        s1++;
        s2++;
        ret++;
    }
    return ret;
}

// see BPMRecord.py from https://github.com/Illumina/GTCtoVCF
// For an insertion relative to the reference, the position of the base immediately 5' to the insertion (on the plus strand) is given.
// For a deletion relative to the reference, the position of the most 5' deleted based (on the plus strand) is given.
static inline int get_indel_alleles(const char *flank,
                                    faidx_t *fai,
                                    const char *seqname,
                                    hts_pos_t pos,
                                    int allele_b_is_del,
                                    char *ref_base,
                                    kstring_t *allele_a,
                                    kstring_t *allele_b)
{
    const char *left = strchr(flank, '[');
    const char *middle = strchr(flank, '/');
    const char *right = strchr(flank, ']');
    if ( !left || !middle || !right ) error("Flank sequence is malformed: %s\n", flank);

    int len;
    char *ref = faidx_fetch_seq(fai, seqname,
         pos - (int)(left - flank), pos - 1 + (int)(right - middle) - 1 + strlen(right + 1), &len);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%"PRId64"\n", seqname, pos + 1);
    strupper(ref);
    int del_left  = len_common_suffix(left - 1  , ref + (left - flank)                       , left - flank);
    int del_right = len_common_prefix(right + 1 , ref + (left - flank) + 1                   , strlen(right + 1));
    int ins_match =      0 == strncmp(middle + 1, ref + (left - flank)                       , right - middle - 1);
    int ins_left  = len_common_suffix(left - 1  , ref + (left - flank) - 1                   , left - flank);
    int ins_right = len_common_prefix(right + 1 , ref + (left - flank) + (right - middle) - 1, strlen(right + 1));
    int ref_is_del = ( del_left >= ins_left ) && ( del_right >= ins_right );
    if ( ( ref_is_del && del_left * del_right == 0 ) || ( !ref_is_del && ( !ins_match || ins_left * ins_right == 0 ) ) )
    {
        ref_is_del = -1;
        ref_base[0] = ref[left - flank];
        kputc(allele_b_is_del ? 'I' : 'D', allele_a);
        kputc(allele_b_is_del ? 'D' : 'I', allele_b);
    }
    else
    {
        ref_base[0] = ref[left - flank - 1 + ref_is_del];
        kputc(ref_base[0], allele_a);
        kputc(ref_base[0], allele_b);
        kputsn(ref_is_del ? middle + 1 : ref + (left - flank), right - middle - 1, allele_b_is_del ? allele_a : allele_b);
    }
    free(ref);
    return ref_is_del;
}

static inline int get_allele_b_idx(char ref_base,
                                   char *allele_a,
                                   char *allele_b)
{
    if ( *allele_a == 'D' || *allele_a == 'I' || *allele_b == 'D' || *allele_b == 'I' ) return 1;

    if ( *allele_a == ref_base ) return 1;
    else if ( *allele_b == ref_base ) return 0;
    else if ( *allele_a == '.' )
    {
        *allele_a = ref_base;
        return 1;
    }
    else if ( *allele_b == '.' )
    {
        *allele_b = ref_base;
        return 0;
    }
    return 2;
}

static inline int get_allele_a_idx(int allele_b_idx)
{
    switch ( allele_b_idx )
    {
        case 0:
            return 1;
        case 1:
            return 0;
        case 2:
            return 1;
        default:
            return -1;
    }
}

static inline int alleles_ab_to_vcf(const char **alleles,
                                    const char *ref_base,
                                    const char *allele_a,
                                    const char *allele_b,
                                    int allele_b_idx)
{
    switch ( allele_b_idx )
    {
        case -1:
            alleles[0] = ref_base;
            return 1;
        case 0:
            alleles[0] = allele_b;
            if ( allele_a[0] == '.' ) return 1;
            alleles[1] = allele_a;
            return 2;
        case 1:
            alleles[0] = allele_a;
            if ( allele_b[0] == '.' ) return 1;
            alleles[1] = allele_b;
            return 2;
        case 2:
            alleles[0] = ref_base;
            alleles[1] = allele_a;
            alleles[2] = allele_b;
            return 3;
        default:
            return -1;
    }
}
