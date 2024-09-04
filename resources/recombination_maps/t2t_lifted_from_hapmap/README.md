to generate t2t chrom maps, hg38 chrom maps were:
	1) decompressed
	2) loaded as pandas dataframes and converted to bedfiles:
		map = pd.read_csv(mapfile, sep='\t')
		map['chr'] = 'chr' + map['chr'].astype(str)
		map['pos1'] = map['pos'] + 1
		map[['chr','pos','pos1','cM']].to_csv(mapfile+'.bed', sep='\t', index=False)
	3) Bedfiles were lifted over to t2t using UCSC's LiftOver tool
		liftOver -bedPlus=3 $contig.b38.bed hg38-chm13v2.over.chain $contig.t2t.bed $contig.t2t.bed.unmapped
	4) python used to rearrange columns in bedfiles to [['pos','chr','cM']].
	5) Saved with header as "$contig.t2t.gmap'
	6) compressed
