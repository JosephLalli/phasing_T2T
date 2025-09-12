import argparse

def main():
    parser = argparse.ArgumentParser(prog='generate_regions_for_rare_phasing.py', description=__doc__)
    parser.add_argument('-f', '--fasta_fai', required=True, help='ref fasta index')
    parser.add_argument('-c', '--chunk_size', help='size of each chunk (Mb).', type=int, default=5)
    parser.add_argument('-s', '--overlap_size', help='overlap size (Mb)', type=int, default=1)
    parser.add_argument('-o', '--output', help='output file name', default='./regions.txt')

    args = parser.parse_args()
    overlap = args.overlap_size*int(1e6)
    chunksize = args.chunk_size*int(1e6)

    contig_coords = dict()
    with open(args.fasta_fai, 'r') as fai:
        for line in fai:
            region_start=0
            line = line.split('\t')
            contig, contig_end = line[0], int(line[1])
            contig_coords[contig] = list()
            for end in range(chunksize-overlap, contig_end+chunksize-overlap, chunksize-overlap):
                end = min(end, contig_end-overlap)
                contig_coords[contig].append((int(region_start), int(end+overlap)))
                region_start = end
    
    with open(args.output, 'w') as out:
        for contig, coords in contig_coords.items():
            for start, end in coords:
                out.write(f"{contig}:{start}-{end}\n")
                



if __name__ == '__main__':
    main()
