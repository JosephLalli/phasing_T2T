import sys
import os
import gzip
import argparse
from cyvcf2 import VCF
from cyvcf2 import Writer as vcfOpen
import intervaltree as it
import numpy as np
from tqdm import tqdm
from pyliftover import LiftOver
from collections import defaultdict
from Bio import SeqIO
import subprocess


class Unliftable(Exception):
    pass


def parse_args():
    parser = argparse.ArgumentParser(
        description="Liftover variants between genome references in an indel-aware manner.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--input-vcf", required=True,
        help="VCF/BCF file to lift over")
    parser.add_argument("--ref-diffs-vcf", required=True,
        help="VCF/BCF of assembly differences (must be in target assembly coordinates)")
    parser.add_argument("--output-vcf", required=True,
        help="Output VCF/BCF path. Use /dev/stdout or - for stdout")
    parser.add_argument("--chain", required=True,
        help="Chain file for coordinate liftover")
    parser.add_argument("--target-fasta", required=True,
        help="Target reference FASTA (may be gzipped)")

    parser.add_argument("--chrom", nargs="+", default=None,
        help="Restrict liftover to these contigs (e.g. --chrom chr1 chr22)")
    parser.add_argument("--region", default=None,
        help="Restrict input liftover to source-assembly region(s); comma-separated regions are allowed")

    realign = parser.add_argument_group("haplotype realignment")
    realign.add_argument("--no-realign", action="store_true", default=False,
        help="Disable haplotype realignment near reference differences")
    realign.add_argument("--realign-distance", type=int, default=50,
        help="Max distance (bp) to search for nearby ref diffs (default: 50)")
    realign.add_argument("--realign-flank", type=int, default=20,
        help="Flanking bases added to each side of the realignment window (default: 20)")
    realign.add_argument("--realign-max-window", type=int, default=200,
        help="Maximum total realignment window size in bp (default: 200)")

    parser.add_argument("--threads", type=int, default=2,
        help="Threads for VCF/BCF reading via cyvcf2 (default: 2)")

    parser.add_argument("--debug", action="store_true", default=False,
        help="Enable verbose debug logging to stderr")
    parser.add_argument("--quiet", action="store_true", default=False,
        help="Suppress progress bars")

    return parser.parse_args()


DEBUG = False
REALIGN_ENABLED = True
REALIGN_DISTANCE = 50
REALIGN_FLANK = 20
REALIGN_MAX_WINDOW = 200
THREADS = 2
QUIET = False


def debug(msg):
    if DEBUG:
        sys.stderr.write(f"{msg}\n")



def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in seq)

def check_var_ref(var):
    if var.CHROM not in target_genome_seq:
        raise ValueError
    target_ref = target_genome_seq[var.CHROM][var.start:var.end]
    if any(base not in ("A", "C", "G", "T") for base in target_ref.upper()):
        raise Unliftable("degenerate base in target reference")
    if var.REF == target_ref:
        return var
    raise ValueError

def parse_regions(region_arg):
    if region_arg is None:
        return None
    regions = [region.strip() for region in region_arg.split(",") if region.strip()]
    return regions or None


def iter_variants(vcf, regions=None):
    if regions is None:
        yield from vcf
        return
    for region in regions:
        yield from vcf(region)


def get_num_vars(vcf_file, regions=None):
    if regions is None:
        out = subprocess.run(f'bcftools index -n {vcf_file}'.split(' '), capture_output = True, text=True)
        if out.returncode == 0 and out.stdout.strip().isdigit():
            return int(out.stdout.strip())
    vcf = VCF(vcf_file, threads=THREADS)
    return sum(1 for _ in iter_variants(vcf, regions))

def add_original_info_tags(var):
    var.INFO["SRC_CHROM"] = var.CHROM
    var.INFO["SRC_POS"] = var.POS
    var.INFO["Original_REF"] = var.REF
    var.INFO["Original_ALT"] = var.ALT[0]
    var.INFO["Original_ID"] = var.ID if var.ID is not None else "."
    var.INFO["SRC_REF_ALT"] = f"{var.REF},{var.ALT[0]}"
    return var

def perform_clean_liftover(var, liftover_obj, return_coordinates=False):
    var = add_original_info_tags(var)
    new_start = liftover_obj.convert_coordinate(var.CHROM, var.start)
    new_end = liftover_obj.convert_coordinate(var.CHROM, var.end)
    new_ref = var.REF
    new_alt = var.ALT[0]
    if len(new_start) != 1 or len(new_end) != 1:
        return [], []
    new_start = new_start[0][:3]
    new_end = new_end[0][:3]
    if new_start[2] == '-':
        new_start, new_end = new_end, new_start
        if not var.is_snp:
            new_base_nucleotide=target_genome_seq[new_start[0]][new_start[1]]
            new_ref = new_base_nucleotide + rev_comp(var.REF[1:])
            new_alt = new_base_nucleotide + rev_comp(var.ALT[0][1:])
        else:
            new_start = (new_start[0], new_start[1]+1)
            new_end = (new_start[0], new_end[1]+1)
            new_ref = rev_comp(var.REF)
            new_alt = rev_comp(var.ALT[0])
        var.REF = new_ref
        var.ALT = new_alt
    var.CHROM = new_start[0]
    var.set_pos(new_start[1])
    if not return_coordinates:
        return var
    return new_start[0:2], new_end[0:2]

def revert_variant(var):
    original_contig = var.INFO['SRC_CHROM']
    original_pos = var.INFO['SRC_POS']
    original_ref = var.INFO['Original_REF']
    original_alt = var.INFO['Original_ALT']
    var.CHROM=original_contig
    var.set_pos(original_pos-1)
    var.REF=original_ref
    var.ALT=original_alt
    return var

def flip_variant(var, all_variants_at_site):
    alt_alleles=var.genotype.array()
    for multivar in all_variants_at_site:
        missing = np.isnan(multivar.genotype.array())
        if missing.any(): # if any missing, fix at homozygous ref
            multivar.gt_types[missing]=0
        if multivar!=var:
            alt_alleles[:,:2] += multivar.genotype.array()[:,:2]
    alt_alleles[:,:2] = (alt_alleles[:,:2]==0).astype(int)
    ref_var = var.REF
    alt_var = var.ALT[0]
    var.REF=alt_var
    var.ALT=[ref_var]
    var.genotypes = alt_alleles
    var.genotypes=var.genotypes ##necessary per cyVCF documentation
    var.INFO["Flipped_during_liftover"] = 'Flipped'
    return var


def variant_edit_positions(var_ref, var_alt):
    if len(var_ref) == len(var_alt):
        return {i for i, (ref_base, alt_base) in enumerate(zip(var_ref, var_alt)) if ref_base != alt_base}
    return set(range(1, len(var_ref)))


def trim_common(ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return ref, alt


def left_shift_variant(pos, ref, alt, target_seq, window_start):
    if len(ref) == len(alt):
        return pos, ref, alt
    while pos > window_start:
        prev_base = target_seq[pos - window_start - 1]
        if len(ref) > len(alt):
            if ref[-1] != prev_base:
                break
            ref = prev_base + ref[:-1]
            alt = prev_base + alt
        else:
            if alt[-1] != prev_base:
                break
            ref = prev_base + ref
            alt = prev_base + alt[:-1]
        pos -= 1
    return pos, ref, alt


def normalize_variant(pos, ref, alt, target_seq, window_start):
    if not ref or not alt:
        if pos == window_start:
            anchor = target_seq[pos - window_start]
        else:
            anchor = target_seq[pos - window_start - 1]
            pos -= 1
        if not ref:
            ref = anchor
            alt = anchor + alt
        else:
            alt = anchor
            ref = anchor + ref
    ref, alt = trim_common(ref, alt)
    pos, ref, alt = left_shift_variant(pos, ref, alt, target_seq, window_start)
    ref, alt = trim_common(ref, alt)
    return pos, ref, alt


def build_source_reference(target_seq, window_start, window_end, diffs):
    source_parts = []
    target_to_source = {}
    source_index = 0
    cursor = window_start
    for diff in diffs:
        if diff.start < cursor:
            raise Unliftable("overlapping ref diffs in window")
        for pos in range(cursor, diff.start):
            source_parts.append(target_seq[pos - window_start])
            target_to_source[pos] = source_index
            source_index += 1
        alt = diff.ALT[0]
        ref_len = len(diff.REF)
        for i, base in enumerate(alt):
            source_parts.append(base)
            if i < ref_len:
                target_to_source[diff.start + i] = source_index
            source_index += 1
        cursor = diff.end
    for pos in range(cursor, window_end):
        source_parts.append(target_seq[pos - window_start])
        target_to_source[pos] = source_index
        source_index += 1
    return "".join(source_parts), target_to_source


def global_align(ref_seq, alt_seq, tie_break="left"):
    n = len(ref_seq)
    m = len(alt_seq)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = i
        trace[i][0] = 1
    for j in range(1, m + 1):
        dp[0][j] = j
        trace[0][j] = 2
    for i in range(1, n + 1):
        ref_base = ref_seq[i - 1]
        for j in range(1, m + 1):
            alt_base = alt_seq[j - 1]
            cost = 0 if ref_base == alt_base else 1
            diag = dp[i - 1][j - 1] + cost
            up = dp[i - 1][j] + 1
            left = dp[i][j - 1] + 1
            best = min(diag, up, left)
            dp[i][j] = best
            choices = []
            if diag == best:
                choices.append(0)
            if up == best:
                choices.append(1)
            if left == best:
                choices.append(2)
            if tie_break == "right":
                for move in (0, 2, 1):
                    if move in choices:
                        trace[i][j] = move
                        break
            else:
                for move in (0, 1, 2):
                    if move in choices:
                        trace[i][j] = move
                        break
    i = n
    j = m
    aln_ref = []
    aln_alt = []
    while i > 0 or j > 0:
        move = trace[i][j]
        if move == 0:
            aln_ref.append(ref_seq[i - 1])
            aln_alt.append(alt_seq[j - 1])
            i -= 1
            j -= 1
        elif move == 1:
            aln_ref.append(ref_seq[i - 1])
            aln_alt.append("-")
            i -= 1
        else:
            aln_ref.append("-")
            aln_alt.append(alt_seq[j - 1])
            j -= 1
    return "".join(reversed(aln_ref)), "".join(reversed(aln_alt))


def variant_from_alignment(aln_ref, aln_alt, window_start, target_seq):
    ref_idx = 0
    first_diff = None
    ref_allele = []
    alt_allele = []
    for ref_base, alt_base in zip(aln_ref, aln_alt):
        if ref_base == alt_base:
            if ref_base != "-":
                ref_idx += 1
            continue
        if first_diff is None:
            first_diff = ref_idx
        if ref_base != "-":
            ref_allele.append(ref_base)
        if alt_base != "-":
            alt_allele.append(alt_base)
        if ref_base != "-":
            ref_idx += 1
    if first_diff is None:
        return None
    pos = window_start + first_diff
    ref = "".join(ref_allele)
    alt = "".join(alt_allele)
    pos, ref, alt = normalize_variant(pos, ref, alt, target_seq, window_start)
    return pos, ref, alt


def attempt_haplotype_realignment(var, chrom, span_start, span_end, ref_diffs):
    if not REALIGN_ENABLED:
        return None
    if len(var.REF) == len(var.ALT[0]):
        return None
    if chrom not in target_genome_seq:
        return None
    search_start = max(0, span_start - REALIGN_DISTANCE)
    search_end = min(len(target_genome_seq[chrom]), span_end + REALIGN_DISTANCE)
    nearby_all = [x.data for x in ref_diffs[chrom][search_start:search_end]]
    indels = [x for x in nearby_all if len(x.REF) != len(x.ALT[0])]
    if not indels:
        return None

    def distance_to_span(diff):
        if diff.end <= span_start:
            return span_start - diff.end
        if diff.start >= span_end:
            return diff.start - span_end
        return 0

    diff = min(indels, key=distance_to_span)
    if distance_to_span(diff) > REALIGN_DISTANCE:
        return None
    flank = REALIGN_FLANK
    target_len = len(target_genome_seq[chrom])
    window_start = None
    window_end = None
    window_diffs = None
    while flank >= 0:
        window_start = max(0, min(span_start, diff.start) - flank)
        window_end = min(target_len, max(span_end, diff.end) + flank)
        if window_end - window_start > REALIGN_MAX_WINDOW:
            flank -= 1
            continue
        if not (window_start <= span_start < window_end and window_start < span_end <= window_end):
            flank -= 1
            continue
        if not (window_start <= diff.start and diff.end <= window_end):
            flank -= 1
            continue
        window_diffs = [x.data for x in ref_diffs[chrom][window_start:window_end]]
        other_diffs = [x for x in window_diffs if not (x.start == diff.start and x.end == diff.end and x.REF == diff.REF and x.ALT[0] == diff.ALT[0])]
        if not other_diffs:
            break
        flank -= 1
    if flank < 0 or window_diffs is None:
        return None
    target_seq = target_genome_seq[chrom][window_start:window_end]
    source_ref, target_to_source = build_source_reference(target_seq, window_start, window_end, [diff])
    if span_start not in target_to_source:
        debug(f"realign: no target_to_source for {chrom}:{span_start} in window {window_start}-{window_end}")
        return None
    source_idx = target_to_source[span_start]
    if source_ref[source_idx:source_idx + len(var.REF)] != var.REF:
        debug(f"realign: source ref mismatch at {chrom}:{span_start} ({source_ref[source_idx:source_idx + len(var.REF)]} != {var.REF})")
        return None
    source_alt = source_ref[:source_idx] + var.ALT[0] + source_ref[source_idx + len(var.REF):]
    best = None
    for tie_break in ("left", "right"):
        aln_ref, aln_alt = global_align(target_seq, source_alt, tie_break=tie_break)
        realigned = variant_from_alignment(aln_ref, aln_alt, window_start, target_seq)
        if realigned is None:
            continue
        pos, ref, alt = realigned
        if ref == alt:
            continue
        if target_seq[pos - window_start:pos - window_start + len(ref)] != ref:
            continue
        var_start = pos
        var_end = pos + len(ref)
        dist = 0
        if var_end <= diff.start:
            dist = diff.start - var_end
        elif diff.end <= var_start:
            dist = var_start - diff.end
        offset = abs(pos - diff.start)
        candidate = (dist, offset, pos, ref, alt, tie_break)
        if best is None or candidate < best:
            best = candidate

    ins_bases = var.ALT[0][len(var.REF):]
    extra = diff.REF[len(diff.ALT[0]):]
    if not ins_bases or len(set(ins_bases)) == 1:
        if ins_bases:
            base = ins_bases[0]
        else:
            base = None
        k = len(ins_bases)
        if base and len(var.REF) < len(var.ALT[0]) and extra.startswith(base) and 0 < k <= len(extra):
            run_start = span_start + 1
            run_end = diff.start
            if run_end <= run_start or all(b == base for b in target_genome_seq[chrom][run_start:run_end]):
                pos = diff.start + k
                ref = extra[k - 1:]
                alt = base
                if ref != alt and target_seq[pos - window_start:pos - window_start + len(ref)] == ref:
                    var_start = pos
                    var_end = pos + len(ref)
                    dist = 0
                    if var_end <= diff.start:
                        dist = diff.start - var_end
                    elif diff.end <= var_start:
                        dist = var_start - diff.end
                    offset = abs(pos - diff.start)
                    candidate = (dist, offset, pos, ref, alt, "shift")
                    if best is None or candidate < best:
                        best = candidate
    if best is None:
        debug(f"realign: no variant from alignment for {chrom}:{span_start}")
        return None
    _, _, pos, ref, alt, tie_break = best
    debug(f"realign({tie_break}): {chrom}:{span_start} -> {pos} {ref}/{alt} window {window_start}-{window_end}")
    return pos, ref, alt


def compute_adjusted_ref_alt(var_ref, var_alt, ref, lifted_start, lifted_end):
    start = min(lifted_start, lifted_end)
    end = max(lifted_start, lifted_end)
    overlap_start = max(start, ref.start)
    overlap_end = min(end, ref.end)
    if overlap_start >= overlap_end:
        return var_ref, var_alt

    ref_len = len(ref.REF)
    alt_len = len(ref.ALT[0])
    if ref.start < lifted_start < ref.end and ref_len < alt_len:
        raise Unliftable("start inside deletion ref diff")

    target_sub = ref.REF[overlap_start - ref.start:overlap_end - ref.start]
    if ref_len == alt_len:
        source_sub = ref.ALT[0][overlap_start - ref.start:overlap_end - ref.start]
    elif ref_len > alt_len:
        src_start = max(0, overlap_start - ref.start)
        src_end = min(alt_len, overlap_end - ref.start)
        source_sub = ref.ALT[0][src_start:src_end] if src_start < src_end else ""
    else:
        source_sub = ref.ALT[0][overlap_start - ref.start:overlap_end - ref.start]

    offset = overlap_start - start
    if ref_len != alt_len and offset == 0 and var_ref.startswith(ref.ALT[0]):
        prefix = ref.ALT[0]
        suffix = var_ref[len(prefix):]
        extra = ref.REF[len(prefix):]
        new_ref = ref.REF + suffix
        if len(var_ref) == len(var_alt):
            if len(var_alt) < len(prefix):
                raise Unliftable("variant alt shorter than ref diff prefix")
            new_alt = var_alt[:len(prefix)] + extra + var_alt[len(prefix):]
        else:
            new_alt = var_alt
        while len(new_ref) > 1 and len(new_alt) > 1 and new_ref[-1] == new_alt[-1]:
            new_ref = new_ref[:-1]
            new_alt = new_alt[:-1]
        return new_ref, new_alt
    if source_sub:
        if offset < 0 or offset + len(source_sub) > len(var_ref):
            raise Unliftable("ref diff out of bounds")
        if var_ref[offset:offset + len(source_sub)] != source_sub:
            raise Unliftable("source ref mismatch")

    if ref_len != alt_len:
        new_ref = var_ref[:offset] + target_sub + var_ref[offset + len(source_sub):]
        if len(var_ref) != len(var_alt):
            diff_region = set(range(offset, offset + len(source_sub)))
            if variant_edit_positions(var_ref, var_alt) & diff_region:
                raise Unliftable("variant edits overlap indel ref diff")
            if source_sub and var_alt[offset:offset + len(source_sub)] != source_sub:
                raise Unliftable("source alt mismatch in ref diff")
            new_alt = var_alt[:offset] + target_sub + var_alt[offset + len(source_sub):]
        else:
            if source_sub:
                replacement = list(target_sub)
                for i in range(min(len(source_sub), len(replacement))):
                    src_i = offset + i
                    if src_i >= len(var_alt) or src_i >= len(var_ref):
                        raise Unliftable("ref diff out of bounds")
                    if var_alt[src_i] != var_ref[src_i]:
                        replacement[i] = var_alt[src_i]
                new_alt = var_alt[:offset] + "".join(replacement) + var_alt[offset + len(source_sub):]
            else:
                new_alt = var_alt[:offset] + target_sub + var_alt[offset:]
    else:
        new_ref = var_ref[:offset] + target_sub + var_ref[offset + len(source_sub):]
        alt_list = list(var_alt)
        for i, tgt_base in enumerate(target_sub):
            src_i = offset + i
            if src_i >= len(var_ref):
                raise Unliftable("ref diff out of bounds")
            if src_i < len(var_alt) and var_alt[src_i] == var_ref[src_i]:
                alt_list[src_i] = tgt_base
        new_alt = "".join(alt_list)

    while len(new_ref) > 1 and len(new_alt) > 1 and new_ref[-1] == new_alt[-1]:
        new_ref = new_ref[:-1]
        new_alt = new_alt[:-1]
    return new_ref, new_alt


def convert_VCF_to_intervaltree(vcf_file):
    tree_dict = dict()
    vcf = VCF(vcf_file, gts012=True, threads=THREADS)
    for chrom in vcf.seqnames:
        if chrom_filter is not None and chrom not in chrom_filter:
            continue
        tree_dict[chrom] = list()
    for var in tqdm(vcf, total=get_num_vars(vcf_file), disable=QUIET):
        if chrom_filter is not None and var.CHROM not in chrom_filter:
            continue
        tree_dict[var.CHROM].append(var)
    sys.stderr.write('Organizing vcf containing variation between builds...\n')
    for chrom in tqdm(tree_dict.keys(), total=len(tree_dict), disable=QUIET):
        tree_dict[chrom] = it.IntervalTree.from_tuples([(var.start, var.end, var) for var in tree_dict[chrom]])
    return tree_dict


def iterate_over_positions(variants):
    iterator = iter(variants)
    try:
        var = next(iterator)
    except StopIteration:
        return
    current_position = var.POS-1
    current_chrom = var.CHROM
    poslist = [var]
    for var in iterator:
        pos = var.POS-1
        if pos != current_position:
            assert (pos > current_position) or (var.CHROM != current_chrom) , f'Variant at position {pos} is after a variant at {current_position}. Input variant file must be sorted before liftover.'
            try:
                if '*' in var.ALT[0]:
                    continue
                else:
                    yield ((current_chrom, current_position), poslist)
                    poslist = list()
                    current_position = pos
                    current_chrom=var.CHROM
                    poslist=[var]
            except IndexError:
                sys.stderr.write (f'{var.CHROM} {var.POS} {var.REF} {var.ALT}' +'\n')
                sys.stderr.write (str(var)+'\n')
        else:
            poslist.append(var)
    yield (current_chrom, current_position), poslist


#-------------------------------------
# start of main program

def load_target_genome(target_fasta_path, chroms):
    seqs = {}
    if target_fasta_path.endswith(".gz"):
        handle = gzip.open(target_fasta_path, "rt")
    else:
        handle = open(target_fasta_path, "r")
    with handle as tfasta:
        for record in SeqIO.parse(tfasta, "fasta"):
            if chroms is not None and record.id not in chroms:
                continue
            seqs[record.id] = str(record.seq).upper()
            if chroms is not None and len(seqs) >= len(chroms):
                break
    return seqs




args = parse_args()

DEBUG = args.debug
REALIGN_ENABLED = not args.no_realign
REALIGN_DISTANCE = args.realign_distance
REALIGN_FLANK = args.realign_flank
REALIGN_MAX_WINDOW = args.realign_max_window
THREADS = args.threads
QUIET = args.quiet

unlifted_vcf = args.input_vcf
ref_alignment = args.ref_diffs_vcf
out_vcf = args.output_vcf
chain = args.chain
target_fasta = args.target_fasta

if out_vcf == '/dev/stdout':
    out_vcf = '-'
    mode='wbu'
else:
    mode=None # infer


def get_num_unique_positions(vcf_file: str) -> int:
    vcf = VCF(vcf_file, threads=THREADS)
    positions = {(x.CHROM, x.POS) for x in iter_variants(vcf, source_regions)}
    return len(positions)

degenerate_nucleotides='UWSMKRYBDHV*'
source_regions = parse_regions(args.region)

chrom_filter = set()
if args.chrom is not None:
    chrom_filter.update(args.chrom)
if source_regions is not None:
    chrom_filter.update(region.split(":", 1)[0] for region in source_regions)
if not chrom_filter:
    chrom_filter = None

if chrom_filter is not None:
    debug(f"restricting to contigs: {sorted(chrom_filter)}")
if source_regions is not None:
    debug(f"restricting source input to regions: {source_regions}")

sys.stderr.write ('Loading target reference genome...\n')

target_genome_seq = load_target_genome(target_fasta, chrom_filter)
for chrom in target_genome_seq.keys():
    for nucleotide in degenerate_nucleotides:
        if nucleotide in target_genome_seq[chrom]:
            target_genome_seq[chrom] = target_genome_seq[chrom].replace(nucleotide, 'N')

sys.stderr.write ('Loading chainfile...\n')
ls = LiftOver(chain)

sys.stderr.write('Loading vcf containing variation between builds...\n')
ref_diffs = convert_VCF_to_intervaltree(ref_alignment)

invcf = VCF(unlifted_vcf, gts012=True, threads=THREADS)

invcf.add_info_to_header({'ID': 'SRC_CHROM', #'Original_Contig', 
                          'Description': 'Original contig of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'SRC_POS', #'Original_POS', 
                          'Description': 'Original position of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Original_REF',
                          'Description': 'Original reference sequence of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Original_ALT',
                          'Description': 'Original alt sequence of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'SRC_REF_ALT',
                          'Description': 'Original ref,alt of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Original_ID',
                          'Description': 'Original variant id of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Flipped_during_liftover',
                          'Description': 'REF/ALT were flipped during liftover. GTs were altered accordingly.', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Realigned_during_liftover',
                          'Description': 'ALT haplotype was realigned to target reference near build differences.', 
                          'Type':'Character', 'Number': '1'})

if os.path.exists(out_vcf):
    os.remove(out_vcf)

out = vcfOpen(out_vcf, invcf, mode=mode)


sys.stderr.write ('Lifting...\n')
unliftable = defaultdict(list)
multiple_overlaps = defaultdict(list)
ref_seq_problem_liftovers = defaultdict(list)
for pos, variants in tqdm(iterate_over_positions(iter_variants(invcf, source_regions)), total=get_num_unique_positions(unlifted_vcf), disable=QUIET):
    lifted_records = []
    lifted_variants = []
    for var in variants:
        try:
            lifted_start, lifted_end = perform_clean_liftover(var, ls, return_coordinates=True)
        except Exception as exc:
            debug(f"unliftable: liftover exception for {var.CHROM}:{var.POS} ({exc})")
            unliftable[pos].append(revert_variant(var))
            continue
        if len(lifted_start) == 0 or len(lifted_end) == 0:
            debug(f"unliftable: missing coords for {var.CHROM}:{var.POS}")
            unliftable[pos].append(revert_variant(var))
            continue
        lifted_records.append((var, lifted_start, lifted_end))

    if not lifted_records:
        continue

    already_flipped = False
    for var, lifted_start, lifted_end in lifted_records:
        chrom = lifted_start[0]
        span_start = min(lifted_start[1], lifted_end[1])
        span_end = max(lifted_start[1], lifted_end[1])
        overlap = list(ref_diffs[chrom][span_start:span_end])
        debug(f"processing {var.CHROM}:{var.POS} span {span_start}-{span_end} overlaps {len(overlap)}")

        if len(overlap) > 1:
            debug(f"multiple overlaps for {var.CHROM}:{var.POS}")
            multiple_overlaps[pos].append(revert_variant(var))
            continue
        if any('N' in x.data.REF for x in overlap):
            debug(f"unliftable: N in ref diff for {var.CHROM}:{var.POS}")
            unliftable[pos].append(revert_variant(var))
            continue
        try:
            if len(overlap) == 0:
                realigned = attempt_haplotype_realignment(var, chrom, span_start, span_end, ref_diffs)
                if realigned is not None:
                    new_pos, new_ref, new_alt = realigned
                    var.set_pos(new_pos)
                    if new_ref == new_alt:
                        assert not already_flipped, f"double flip at {var.CHROM}:{var.POS}"
                        var = flip_variant(var, variants)
                        already_flipped = True
                    else:
                        var.REF = new_ref
                        var.ALT = [new_alt]
                    var.INFO["Realigned_during_liftover"] = 'Realigned'
                check_var_ref(var)
            else:
                new_ref, new_alt = compute_adjusted_ref_alt(
                    var.REF,
                    var.ALT[0],
                    overlap[0].data,
                    lifted_start[1],
                    lifted_end[1],
                )
                if new_ref == new_alt:
                    assert not already_flipped, f"double flip at {var.CHROM}:{var.POS}"
                    var = flip_variant(var, variants)
                    already_flipped = True
                else:
                    var.REF = new_ref
                    var.ALT = [new_alt]
                assert len(var.REF) == len(var.ALT[0]) or var.REF[0] == var.ALT[0][0]
                check_var_ref(var)
            lifted_variants.append(var)
        except Unliftable:
            debug(f"unliftable: overlap rule for {var.CHROM}:{var.POS}")
            unliftable[pos].append(revert_variant(var))
            continue
        except (ValueError, AssertionError):
            debug(f"mismatch: ref/alt for {var.CHROM}:{var.POS}")
            ref_seq_problem_liftovers[pos].append(revert_variant(var))
            continue

    for var in lifted_variants:
        out.write_record(var)

out.close()

sys.stderr.write (f'{len(unliftable)} variants could not be lifted over because they did not have a start and/or end coordinate in the target assembly.\n')
sys.stderr.write (f'{len(multiple_overlaps)} variants could not be lifted over because they had multiple potential liftover spots.\n')
sys.stderr.write (f'{len(ref_seq_problem_liftovers)} variants had an unspecified error in liftover.\n')
sys.stderr.write (f'Writing these variants to disk...\n')

base_out='.'.join(out_vcf.replace('.gz','').split('.')[:-1])
out = vcfOpen(base_out+'.unliftable.bcf', invcf, mode=mode)
for pos, variants in unliftable.items():
    for var in variants:
        out.write_record(var)
out.close()

out = vcfOpen(base_out+'.multiple_overlaps.bcf', invcf, mode=mode)
for pos, variants in multiple_overlaps.items():
    for var in variants:
        out.write_record(var)
out.close()

out = vcfOpen(base_out+'.ref_seq_mismatches.bcf', invcf, mode=mode)
for pos, variants in ref_seq_problem_liftovers.items():
    for var in variants:
        out.write_record(var)
out.close()

sys.stderr.write ('\nDone!\n')
sys.stderr.write ('Note: lifted vcf file still requires indel normalizing, sorting, and recalculation of INFO fields. Format fields besides GT are no longer reliable.\n')
