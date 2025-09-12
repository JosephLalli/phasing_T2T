library(karyoploteR)
library(arrow)
library(PlotTools)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

get_script_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
}

named_group_split <- function(.tbl, ...) {
	grouped <- group_by(.tbl, ...)
	names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

	grouped %>%
		group_split() %>%
		rlang::set_names(names)
}

clamp <- function(x, min, max) {
	case_when(
		x < min ~ min,
		x > max ~ max,
		.default = x
	)
}

clip_data <- function(gr1, gr2){
	o = findOverlaps(gr1, gr2)
	if (length(o) == 0) {
		return(gr1[FALSE])  # Return empty GRanges with same structure
	}
	grl1 = split(gr1[queryHits(o)], 1:length(o))
	grl2 = split(gr2[subjectHits(o)], 1:length(o))
	foo = function(x, y) {
	   rv = x
	   start(rv) = max(start(x), start(y))
	   end(rv) = min(end(x), end(y))
	   return(rv)
	}
	unlist(mendoapply(foo, grl1, y=grl2))
}

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
	scale = (length(lut)-1)/(max-min)

	dev.new(width=1.75, height=5)
	plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
	axis(2, ticks, las=1)
	for (i in 1:(length(lut)-1)) {
		y = (i-1)/scale + min
		rect(0,y,10,y+1/scale, col=lut[i], border=NA)
	}
}


setwd(get_script_path())
set.seed(42)

if (!dir.exists("../figures/figure6_ideograms")) {
  dir.create("../figures/figure6_ideograms", recursive = TRUE)
}

rolling_data='../intermediate_data/rolling_stats_500k_window.parquet'
t2t_cytobands='../resources/chm13v2.0_cytobands_allchrs.w_header.bed'
grch38_cytobands='../resources/grch38_cytobands_allchrs.w_header.bed'
rolling_data = read_parquet(rolling_data)

rolling_data$switch_error_rate = replace_na(rolling_data$n_switch_errors/rolling_data$n_checked, 0) *100
rolling_data$gt_error_rate = replace_na(rolling_data$n_gt_errors/rolling_data$n_gt_checked, 0) *100



GRCh38_rolling_data = rolling_data[rolling_data$genome=='GRCh38',]
T2T_rolling_data = rolling_data[rolling_data$genome=='CHM13v2.0',]


grch38_data <- GRCh38_rolling_data %>% within(method <- paste(method_of_phasing,ground_truth_data_source,sep='_')) %>%
									   named_group_split(method) %>%
                                       map(makeGRangesFromDataFrame,
                                       	   seqnames.field='chrom',
                                       	   start.field='start',
                                       	   end.field='end',
                                       	   seqinfo=Seqinfo(genome="hg38"),
                                       	   keep.extra.columns=TRUE, ignore.strand=TRUE, starts.in.df.are.0based=TRUE)

chm13_data <- T2T_rolling_data %>% within(method <- paste(method_of_phasing,ground_truth_data_source,sep='_')) %>%
			  							  named_group_split(method) %>%
										  map(makeGRangesFromDataFrame,
											  seqnames.field='chrom',
											  start.field='start',
											  end.field='end',
											  seqinfo=SeqinfoForUCSCGenome("hs1"),
											  keep.extra.columns=TRUE, ignore.strand=TRUE, starts.in.df.are.0based=TRUE)

### make cnv_grch38, cnv_chm13
cnv=read_tsv('../resources/decipher_syndromes.txt')
cnv[cnv$start_grch38>cnv$end_grch38, c('start_grch38','end_grch38')] = cnv[cnv$start_grch38>cnv$end_grch38, c('end_grch38','start_grch38')]
cnv[cnv$start_chm13>cnv$end_chm13, c('start_chm13','end_chm13')] = cnv[cnv$start_chm13>cnv$end_chm13, c('end_chm13','start_chm13')]

cnv_chm13 = makeGRangesFromDataFrame(cnv,keep.extra.columns = TRUE,seqinfo=SeqinfoForUCSCGenome('hs01'),start.field='start_chm13',end.field='end_chm13',strand.field='strand_chm13',starts.in.df.are.0based=TRUE)
cnv_grch38 = makeGRangesFromDataFrame(cnv,keep.extra.columns = TRUE,seqinfo=SeqinfoForUCSCGenome('hg38'),start.field='start_grch38',end.field='end_grch38',strand.field='strand_grch38',starts.in.df.are.0based=TRUE)

pdf.options(width=3.5,
			height=3.5,
			pointsize=6)

chm13_3202_vs_HPRC_data=chm13_data[['1kgp_variation_phased_with_reference_panel_HPRC_samples']]
grch38_3202_vs_HPRC_data=grch38_data[['1kgp_variation_phased_with_reference_panel_HPRC_samples']]

ser_ymax=10
gt_ymax=10
num_ticks=6
var_density_palette= c("white", "#F2AF29","#EB6424")

default_pts = 6
linewidth_pts=.5
linewidth = linewidth_pts/.75 # Default is 1, which pdf makes 0.75 pt. So this converts pts to R linewidth units.
percent_size_pts = 5
chrom_size_pts = 8
label_size_pts = 6
base_number_size_pts=5

variants_band_size=.1
data_band_size=0.375
gap_size=(1-(data_band_size*2+variants_band_size))/2
var_band_start = 0
var_band_stop=variants_band_size
data_1_start = variants_band_size+gap_size
data_1_end = variants_band_size + gap_size + data_band_size
data_2_start = variants_band_size+gap_size+data_band_size+gap_size
data_2_end = 1

ser_color="#EB6424"
gt_color='#EEB902'
chm13_ser_color="#593F62"
chm13_gt_color='#7B6D8D'
grch38_ser_color="#8499B1"
grch38_gt_color='#A5C4D4'
max.var.density=round(max(max(grch38_3202_vs_HPRC_data$variants_per_1k), max(chm13_3202_vs_HPRC_data$variants_per_1k)))

plot.params <- getDefaultPlotParams(plot.type=2)
plot.params$leftmargin <- .2
plot.params$data1height <- 200
plot.params$data2height <- 200
plot.params$data1inmargin <- plot.params$data1height*gap_size
plot.params$data2inmargin <- plot.params$data2height*gap_size
plot.params$ideogramheight <- 30

slop=1.5e6

chm13_22q11_zoom <- toGRanges(data.frame("chr22", 19397653-slop, 23803198+slop))
chm13_anglemans_zoom <- toGRanges(data.frame("chr15", 20807475-slop, 25935855+slop))
chm13_16p12.1_zoom <- toGRanges(data.frame("chr16", 21398763-slop, 30473116+slop))
chm13_1q21.1_zoom <- toGRanges(data.frame("chr1", 145165377-slop, 147746604+slop))
chm13_8p23.1_zoom <- toGRanges(data.frame("chr8", 7833176-slop, 11501098+slop))
chm13_cytobands <- makeGRangesFromDataFrame(read_tsv("../resources/chm13v2.0_cytobands_allchrs.w_header.bed"),keep.extra.columns=TRUE)
chm13_regions = list(chm13_22q11_zoom, chm13_anglemans_zoom, chm13_1q21.1_zoom, chm13_16p12.1_zoom, chm13_8p23.1_zoom)
chm13_region_names=list('chm13_22q11', 'chm13_anglemans', 'chm13_1q21.1', 'chm13_16p12.1', "chm13_8p23.1")
ser_color = chm13_ser_color
gt_color = chm13_gt_color
for (i in 1:5){
	pdf(file=paste0("../figures/figure6_ideograms/",
					chm13_region_names[i],'.pdf'))
	region=chm13_regions[[i]]
	region_data=clip_data(chm13_3202_vs_HPRC_data, region)
	
	if (length(region_data) == 0) {
		cat("Warning: No data found for region", chm13_region_names[[i]], "\n")
		dev.off()
		next
	}
	
	kp <- plotKaryotype(genome = "hs1", zoom=region, plot.type=2, plot.params=plot.params, cytobands=clip_data(chm13_cytobands, region), lwd=linewidth, cex=(chrom_size_pts/default_pts))
	kp <- kpRect(kp, data = clip_data(cnv_chm13, region), y0=0, y1=1, col="#FFDDDD", border=NA, r0=data_1_start, r1=data_1_end, data.panel=1, clipping=FALSE, lwd=linewidth)
	kp <- kpRect(kp, data = clip_data(cnv_chm13, region), y0=0, y1=1, col="#FFDDDD", border=NA, r0=data_2_start, r1=1, data.panel=1, clipping=FALSE, lwd=linewidth)
	kp <- kpHeatmap(kp, data=region_data, r0=0, r1=var_band_stop, y=clamp(region_data$variants_per_1k, 0, max.var.density), ymax=max.var.density, colors=var_density_palette, data.panel = 1, clipping=FALSE, lwd=linewidth)
	kp <- kpArea(kp, data=region_data, data.panel = 1, r0=data_1_start, r1=data_1_end, y=clamp(region_data$switch_error_rate,0,ser_ymax), ymax=ser_ymax, col=ser_color, clipping=FALSE, lwd=linewidth)
	kp <- kpArea(kp, data=region_data, data.panel = 1, r0=data_2_start, r1=1, y=clamp(region_data$gt_error_rate,0,gt_ymax), ymax=gt_ymax, col=gt_color, clipping=FALSE, lwd=linewidth)
	kpAxis(kp, numticks = num_ticks, data.panel = 1, r0=data_1_start, r1=data_1_end, ymin=0, ymax=ser_ymax, cex=(percent_size_pts/default_pts), labels=paste0(seq(0,ser_ymax,length.out=num_ticks), "%"), lwd=linewidth)
	kpAxis(kp, numticks = num_ticks, data.panel = 1, r0=data_2_start, r1=1, ymin=0, ymax=gt_ymax, cex=(percent_size_pts/default_pts), labels=paste0(seq(0,gt_ymax,length.out=num_ticks), "%"), lwd=linewidth)
	kpAddLabels(kp, labels="# Variants/1Kbp", r0=0, r1=var_band_stop, data.panel = 1, cex=(label_size_pts/default_pts), label.margin=0.05, lwd=linewidth)
	kpAddLabels(kp, labels="Switch\nError\nRate", r0=data_1_start, r1=data_1_end, data.panel = 1, cex=(label_size_pts/default_pts), label.margin=0.09, lwd=linewidth)
	kpAddLabels(kp, labels="Genotype\nError\nRate", r0=data_2_start, r1=1, data.panel = 1, cex=(label_size_pts/default_pts), label.margin=0.09, lwd=linewidth)
	kpAddBaseNumbers(kp, tick.dist = 1000000, minor.tick.dist=200000, cex=(base_number_size_pts/default_pts), lwd=linewidth)
	dev.off()
}

grch38_22q11_zoom <- toGRanges(data.frame("chr22", 19022279-slop, 23380258+slop))
grch38_anglemans_zoom <- toGRanges(data.frame("chr15", 23123712-slop, 28193120+slop))
grch38_1q21.1_zoom <- toGRanges(data.frame("chr1", 145686995-slop, 148411223+slop))
grch38_16p12.1_zoom <- toGRanges(data.frame("chr16", 21398763-slop, 30188533+slop))
grch38_8p23.1_zoom <- toGRanges(data.frame("chr8", 8242533-slop, 11907120+slop))
grch38_cytobands <- makeGRangesFromDataFrame(read_tsv("../resources/grch38_cytobands_allchrs.w_header.bed"),keep.extra.columns=TRUE)
grch38_regions = list(grch38_22q11_zoom, grch38_anglemans_zoom, grch38_1q21.1_zoom, grch38_16p12.1_zoom, grch38_8p23.1_zoom)
grch38_region_names=list('grch38_22q11', 'grch38_anglemans', 'grch38_1q21.1', 'grch38_16p12.1', 'grch38_8p23.1')
ser_color = grch38_ser_color
gt_color = grch38_gt_color
for (i in 1:5){
	pdf(file=paste0("../figures/figure6_ideograms/",grch38_region_names[i],'.pdf'))
	region=grch38_regions[[i]]
	region_data=clip_data(grch38_3202_vs_HPRC_data, region)
	
	if (length(region_data) == 0) {
		cat("Warning: No data found for region", grch38_region_names[[i]], "\n")
		dev.off()
		next
	}
	
	kp <- plotKaryotype(genome = "hg38", zoom=region, plot.type=2, plot.params=plot.params, cex=(8/6), cytobands=clip_data(grch38_cytobands, region), lwd=linewidth)
	kp <- kpRect(kp, data = clip_data(cnv_grch38, region), y0=0, y1=1, col="#FFDDDD", border=NA, r0=data_1_start, r1=data_1_end, data.panel=2, lwd=linewidth)
	kp <- kpRect(kp, data = clip_data(cnv_grch38, region), y0=0, y1=1, col="#FFDDDD", border=NA, r0=data_2_start, r1=1, data.panel=2, lwd=linewidth)
	kp <- kpHeatmap(kp, data=region_data, r0=0, r1=var_band_stop, y=clamp(region_data$variants_per_1k, 0, max.var.density), ymax=max.var.density, colors=var_density_palette, data.panel = 2, clipping=FALSE, lwd=linewidth)
	kp <- kpArea(kp, data=region_data, data.panel = 2, r0=data_1_start, r1=data_1_end, y=clamp(region_data$switch_error_rate,0,ser_ymax), ymax=ser_ymax, col=ser_color, clipping=FALSE, lwd=linewidth)
	kp <- kpArea(kp, data=region_data, data.panel = 2, r0=data_2_start, r1=1, y=clamp(region_data$gt_error_rate,0,gt_ymax), ymax=gt_ymax, col=gt_color, clipping=FALSE, lwd=linewidth)
	kpAxis(kp, numticks = num_ticks, data.panel = 2, r0=data_1_start, r1=data_1_end, ymin=0, ymax=ser_ymax, cex=1, labels=paste0(seq(0,ser_ymax,length.out=num_ticks), "%"))
	kpAxis(kp, numticks = num_ticks, data.panel = 2, r0=data_2_start, r1=1, ymin=0, ymax=gt_ymax, cex=1, labels=paste0(seq(0,gt_ymax,length.out=num_ticks), "%"))
	kpAddLabels(kp, labels="# Variants/1Kbp", r0=0, r1=var_band_stop, data.panel = 2, cex=1, label.margin=0.05, lw=1)
	kpAddLabels(kp, labels="Switch\nError\nRate", r0=data_1_start, r1=data_1_end, data.panel = 2, cex=1, label.margin=0.09)
	kpAddLabels(kp, labels="Genotype\nError\nRate", r0=data_2_start, r1=1, data.panel = 2, cex=1, label.margin=0.09)
	kpAddBaseNumbers(kp, tick.dist = 1000000, minor.tick.dist=200000)
	PlotTools::SpectrumLegend(
		      x="top",                             # Legend position
		      palette = colorRampPalette(var_density_palette),                     # Display our chosen palette
		      legend = c(seq(0,max.var.density, length.out=6)),  # Annotate positions on legend
		      title = "Variants/kb",
		      bty = "n", # Don't frame with box
		      horiz=TRUE,
		      title.adj = 0
		      )
	dev.off()
}


kp <- plotKaryotype(genome = "hs1", plot.type=2, plot.params=plot.params, cytobands=chm13_cytobands)
kp <- kpRect(kp, data = cnv_chm13, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0.2, r1=1, data.panel=1)
kp <- kpHeatmap(kp, data=chm13_3202_vs_HPRC_data, r0=0, r1=0.2, y=clamp(chm13_3202_vs_HPRC_data$variants_per_1k, 0, max.var.density), ymax=max.var.density, colors=var_density_palette, data.panel = 1, clipping=FALSE)
kp <- kpArea(kp, data=chm13_3202_vs_HPRC_data, data.panel = 1, r0=0.3, r1=1, y=clamp(chm13_3202_vs_HPRC_data$switch_error_rate,0,ser_ymax), ymax=ser_ymax, col=ser_color, clipping=FALSE)


pdf(file=paste0("../figures/figure6_ideograms/whole_genome_SER_chm13.pdf"))
plot.params <- getDefaultPlotParams(plot.type=2)
plot.params$leftmargin <- .2
plot.params$data2outmargin = 110
kp <- plotKaryotype(genome = "hs1", plot.type=2, plot.params=plot.params, cytobands=chm13_cytobands)
kp <- kpRect(kp, data = cnv_chm13, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0.2, r1=1, data.panel=1)
kp <- kpHeatmap(kp, data=chm13_3202_vs_HPRC_data, r0=0, r1=0.2, y=clamp(chm13_3202_vs_HPRC_data$variants_per_1k, 0, max.var.density), ymax=max.var.density, colors=var_density_palette, data.panel = 1, clipping=TRUE)
kp <- kpArea(kp, data=chm13_3202_vs_HPRC_data, data.panel = 1, r0=0.3, r1=1, y=clamp(chm13_3202_vs_HPRC_data$switch_error_rate,0,ser_ymax), ymax=ser_ymax, col=chm13_ser_color, clipping=TRUE)

plot.params <- getDefaultPlotParams(plot.type=2)
plot.params$topmargin = plot.params$topmargin = 240
kp <- kpAddCytobands(kp, genome = "hg38", plot.type=2, plot.params=plot.params, labels.plotter = NULL)
kp <- kpRect(kp, data = cnv_grch38, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0.2, r1=1, data.panel=2)
kp <- kpHeatmap(kp, data=grch38_3202_vs_HPRC_data, r0=0, r1=0.2, y=clamp(grch38_3202_vs_HPRC_data$variants_per_1k, 0, max.var.density), ymax=max.var.density, colors=var_density_palette, data.panel = 2, clipping=TRUE)
kp <- kpArea(kp, data=grch38_3202_vs_HPRC_data, data.panel = 2, r0=0.3, r1=1, y=clamp(grch38_3202_vs_HPRC_data$switch_error_rate,0,ser_ymax), ymax=ser_ymax, col=grch38_ser_color, clipping=TRUE)


## Genotype
plot.params <- getDefaultPlotParams(plot.type=2)
kp <- plotKaryotype(genome = "hs1", plot.type=2, plot.params=plot.params, cytobands=chm13_cytobands)
kp <- kpRect(kp, data = cnv_chm13, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0.2, r1=1, data.panel=1)
kp <- kpHeatmap(kp, data=chm13_3202_vs_HPRC_data, r0=0, r1=0.2, y=clamp(chm13_3202_vs_HPRC_data$variants_per_1k, 0, max.var.density), ymax=max.var.density, colors=var_density_palette, data.panel = 1, clipping=TRUE)
kp <- kpArea(kp, data=chm13_3202_vs_HPRC_data, data.panel = 1, r0=0.3, r1=1, y=clamp(chm13_3202_vs_HPRC_data$gt_error_rate,0,gt_ymax), ymax=gt_ymax, col=chm13_gt_color, clipping=TRUE)

plot.params <- getDefaultPlotParams(plot.type=2)
plot.params$leftmargin <- .2
kp <- plotKaryotype(genome = "hg38", plot.type=2, plot.params=plot.params)
kp <- kpRect(kp, data = cnv_grch38, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0.2, r1=1, data.panel=2)
kp <- kpHeatmap(kp, data=grch38_3202_vs_HPRC_data, r0=0, r1=0.2, y=clamp(grch38_3202_vs_HPRC_data$variants_per_1k, 0, max.var.density), ymax=max.var.density, colors=var_density_palette, data.panel = 2, clipping=TRUE)
kp <- kpArea(kp, data=grch38_3202_vs_HPRC_data, data.panel = 2, r0=0.3, r1=1, y=clamp(grch38_3202_vs_HPRC_data$gt_error_rate,0,gt_ymax), ymax=gt_ymax, col=grch38_gt_color, clipping=TRUE)
