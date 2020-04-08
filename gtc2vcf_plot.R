#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2019-2020 Giulio Genovese
#
#  Author: Giulio Genovese <giulio.genovese@gmail.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

library(optparse)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
options(bitmapType = 'cairo')

parser <- OptionParser('usage: gtc2vcf_plot.R [options] --illumina|--affymetrix --vcf <file.vcf> --chrom <string> --pos <integer> --pdf|--png <file>')
parser <- add_option(parser, c('--vcf'), type = 'character', help = 'input VCF file', metavar = '<file.vcf>')
parser <- add_option(parser, c('--illumina'), action = 'store_true', default = FALSE, help = 'whether the input VCF file contains Illumina data')
parser <- add_option(parser, c('--affymetrix'), action = 'store_true', default = FALSE, help = 'whether the input VCF file contains Affymetrix data')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--png'), type = 'character', help = 'output PNG file', metavar = '<file.png>')
parser <- add_option(parser, c('--width'), type = 'integer', default = 7, help = 'inches width of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--height'), type = 'integer', default = 7, help = 'inches height of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 12, help = 'font size [12]', metavar = '<integer>')
parser <- add_option(parser, c('--chrom'), type = 'character', help = 'chromosome', metavar = '<string>')
parser <- add_option(parser, c('--pos'), type = 'integer', help = 'chromosome position', metavar = '<integer>')
parser <- add_option(parser, c('--id'), type = 'character', help = 'variant ID', metavar = '<string>')
parser <- add_option(parser, c('--samples'), type = 'character', help = 'comma-separated list of samples to include', metavar = '<list>')
parser <- add_option(parser, c('--samples-file'), type = 'character', help = 'file with list of samples to include', metavar = '<file>')
parser <- add_option(parser, c('--ellipses'), action = 'store_true', default = FALSE, help = 'plot ellipses around genotype clusters')
parser <- add_option(parser, c('--minimal'), action = 'store_true', default = FALSE, help = 'only plot NORMX/NORMY and BAF/LRR plots')
parser <- add_option(parser, c('--zcall'), action = 'store_true', default = FALSE, help = 'plot ZCall thresholds')
args <- parse_args(parser, commandArgs(trailingOnly = TRUE))

# make sure VCF is passed
if (is.null(args$vcf)) {print_help(parser); stop('option --vcf is required')}
if (is.null(args$chrom)) {print_help(parser); stop('option --chrom is required')}
if (is.null(args$pos)) {print_help(parser); stop('option --pos is required')}
if (!args$illumina && !args$affymetrix) {print_help(parser); stop('either --illumina or --affymetrix is required')}
if (args$illumina && args$affymetrix) {print_help(parser); stop('cannot use --illumina and --affymetrix at the same time')}
if (is.null(args$pdf) && is.null(args$png)) {print_help(parser); stop('either --pdf or --png is required')}
if (!is.null(args$pdf) && !is.null(args$png)) {print_help(parser); stop('cannot use --pdf and --png at the same time')}
if (!is.null(args$png) && !capabilities('png')) {print_help(parser); stop('unable to start device PNG: no png support in this version of R\nyou need to reinstall R with support for PNG to use the --png option\n')}
if (!is.null(args$samples) && !is.null(args$samples_file)) {print_help(parser); stop('cannot use --samples and --samples-file at the same time')}

base <- c('CHROM', 'POS', 'ID')
if (args$illumina)
{
  info <- c('meanR_AA', 'meanR_AB', 'meanR_BB', 'meanTHETA_AA', 'meanTHETA_AB', 'meanTHETA_BB')
  format <- c('GT', 'X', 'Y', 'NORMX', 'NORMY', 'R', 'THETA', 'BAF', 'LRR')
  x <- 'THETA'
  y <- 'R'
}
if (args$affymetrix)
{
  info <- c('meanDELTA_AA', 'meanDELTA_AB', 'meanDELTA_BB', 'meanSIZE_AA', 'meanSIZE_AB', 'meanSIZE_BB')
  format <- c('GT', 'NORMX', 'NORMY', 'DELTA', 'SIZE', 'BAF', 'LRR')
  x <- 'DELTA'
  y <- 'SIZE'
}

if (args$zcall) {
  info <- c(info, c('zthresh_X', 'zthresh_Y'))
  format <- c(format, c('GTA', 'GTZ'))
  gt_color <- 'GTA'
  gt_shape <- 'GTZ'
} else {
  gt_color <- 'GT'
  gt_shape <- 'GT'
}

fmt <- paste0('"[%', paste(base, collapse = '\\t%'), paste(c('', info), collapse = '\\t%INFO/'), paste(c('', format), collapse = '\\t%'), '\\n]"')
names <- c(base, info, format)
cmd <- paste0('bcftools query --format ', fmt, ' ', args$vcf, ' -r ', args$chrom, ':', args$pos, '-', args$pos)
if (!is.null(args$samples)) cmd <- paste(cmd, '--samples', args$samples)
if (!is.null(args$samples_file)) cmd <- paste(cmd, '--samples-file', args$samples_file)
write(paste('Command:', cmd), stderr())
if (packageVersion("data.table") < '1.11.6') {
  df <- setNames(fread(cmd, sep = '\t', header = FALSE, na.strings = '.', data.table = FALSE), names)
} else {
  df <- setNames(fread(cmd = cmd, sep = '\t', header = FALSE, na.strings = '.', data.table = FALSE), names)
}

if (!is.null(args$id)) {
  if (!(args$id %in% unique(df$ID))) stop('Specified ID not present at specified location')
  df <- df[df$ID == args$id,]
} else {
  if ( length(unique(df$ID)) > 1 ) stop('More than one variant at the specified position, use --id to specify which variant to plot')  
}

if (args$illumina)
{
  p1 <- ggplot(df, aes_string(x = 'Y', y = 'X', color = gt_color, shape = gt_shape)) +
    geom_point() +
    theme_bw(base_size = args$fontsize) +
    theme(legend.position = 'none')
  normx <- 'NORMY'
  normy <- 'NORMX'
} else {
  normx <- 'NORMX'
  normy <- 'NORMY'
}
p2 <- ggplot(df, aes_string(x = normx, y = normy, color = gt_color, shape = gt_shape)) +
  geom_point() +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'none')
if (args$zcall) {
  zthresh_X <- unique(df$zthresh_X)
  zthresh_Y <- unique(df$zthresh_Y)
  p2 <- p2 + geom_vline(xintercept = zthresh_Y, color = 'gray') +
             geom_hline(yintercept = zthresh_X, color = 'gray')
}
df_centers = setNames(data.frame(x = sapply(df[,paste0('mean', x, '_', c('AA', 'AB', 'BB'))], unique),
                                 y = sapply(df[,paste0('mean', y, '_', c('AA', 'AB', 'BB'))], unique)), c(x, y))
p3 <- ggplot(df, aes_string(x = x, y = y, color = gt_color, shape = gt_shape)) +
  geom_point() +
  geom_point(data = df_centers, color = 'black', size = 5, shape = 8) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'none')
if (args$illumina) p3 <- p3 + coord_cartesian(xlim = c(0,1))
p4 <- ggplot(df, aes_string(x = 'BAF', y = 'LRR', color = gt_color, shape = gt_shape)) +
  geom_point() +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal') +
  coord_cartesian(xlim = c(0,1))
if (args$ellipses) {
  p1 <- p1 + stat_ellipse(alpha = 1/2)
  p2 <- p2 + stat_ellipse(alpha = 1/2)
  p3 <- p3 + stat_ellipse(alpha = 1/2)
}

if (!is.null(args$pdf)) {
  pdf(args$pdf, width = args$width, height = args$height)
} else {
  png(args$png, width = args$width, height = args$height, units = 'in', res = 150)
}
if (args$minimal) {
  grid.arrange(p2, p4, nrow = 2, ncol = 1, heights = c(3, 4), top = unique(df$ID))
} else {
  if (args$illumina) grid.arrange(p1, p2, p3, p4, nrow = 4, ncol = 1, heights = c(3, 3, 3, 4), top = unique(df$ID))
  if (args$affymetrix) grid.arrange(p2, p3, p4, nrow = 3, ncol = 1, heights = c(3, 3, 4), top = unique(df$ID))
}
invisible(dev.off())
