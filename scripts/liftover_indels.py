import sys
import os
from cyvcf2 import VCF
from cyvcf2 import Writer as vcfOpen
import intervaltree as it
import numpy as np
from tqdm import tqdm
from pyliftover import LiftOver
from collections import defaultdict
from Bio import SeqIO
import subprocess



def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in seq)

def check_var_ref(var):
    if var.REF == target_genome_seq[var.CHROM][var.start:var.end]:
        return var
    else:
        raise ValueError

def get_num_vars(vcf_file):
    out = subprocess.run(f'bcftools index -n {vcf_file}'.split(' '), capture_output = True, text=True)
    return int(out.stdout.strip())

def add_original_info_tags(var):
    var.INFO["Original_Contig"] = var.CHROM
    var.INFO["Original_POS"] = var.POS
    var.INFO["Original_REF"] = var.REF
    var.INFO["Original_ALT"] = var.ALT[0]
    return var

def perform_clean_liftover(var, liftover_obj, return_coordinates=False):
    var = add_original_info_tags(var)
    new_start = liftover_obj.convert_coordinate(var.CHROM, var.start)
    new_end = liftover_obj.convert_coordinate(var.CHROM, var.end)
    new_ref = var.REF
    new_alt = var.ALT[0]
    if len(new_start) != 1 and len(new_end) != 1:
        return [], []
    elif len(new_start) != 1:
        new_end = new_end[0][:3]
        new_start = (new_end[0], new_end[1]-len(var.REF), new_end[2])
    elif len(new_end) != 1:
        new_start = new_start[0][:3]
        new_end = (new_start[0], new_start[1] + len(var.REF), new_start[2])
    else:
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
    original_contig = var.INFO['Original_Contig']
    original_pos = var.INFO['Original_POS']
    original_ref = var.INFO['Original_REF']
    original_alt = var.INFO['Original_ALT']
    var.CHROM=original_contig
    var.set_pos(original_pos-1)
    var.REF=original_ref
    var.ALT=original_alt
    return var

def flip_variant(var, all_variants_at_site):
    alt_alleles=var.genotype.array()
    # alt_alleles[:,:2]=0
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


def liftover_variants_adjusting_for_overlapping_ref_diff(variants, ref):
    already_flipped=False
    for var in variants:
        if var.start <= ref.start:
            new_ref, new_alt = precision_edit_lifted_ref(var, ref)
        else:
            new_ref, new_alt = get_adjusted_ref_alt(var, ref)
        if new_ref == new_alt:
            assert not already_flipped
            var = flip_variant(var, variants)
            already_flipped = True
        else:
            var.REF = new_ref
            var.ALT = [new_alt]
    return variants


def precision_edit_lifted_ref(var, ref):
    ''' adjusting ref alleles to ensure they are at the same start point as var_data '''
    ''' to be clear: origin ref is the variant's reference at the original assembly,
                 and target ref is the variant's reference at the target assembly    '''
    "Does not work for overlapping insertion/deletions"
    # ('chr18', 465979): var: GGA/G, ref:G/GGAGA should lead to G/GGA, instead leads to G/G because alt cannot change
    # G[:1] | GA[2:] | / G[:1] |(GAGA[2:4])
    # in terms of the below:
    # first condition: no prefix to be added to ref (var.REF[:0]) plus ref itself (G) plus nothing after the variant deletion in the ref [5:]
    # ALT: alt itself [G]+suffix from what is leftover in the ref.ALT[0][(ref.start-var.start+1)+len(var.ALT[0])]
    if var.start <= ref.start and var.end > ref.start:
        new_ref = var.REF[:(ref.start-var.start)] + ref.REF + var.REF[(ref.start-var.start+1)+(len(ref.ALT[0])-1):]
        new_alt = var.ALT[0][0] + ref.ALT[0][len(var.REF)-(len(new_ref)-1):]
    else:
        new_ref = var.REF[:(ref.start-var.start)] + var.REF[(ref.start-var.start+1)+len(ref.ALT[0]):]
        new_alt = var.ALT[0]
    while len(new_ref) > 1 and len(new_alt) > 1 and new_ref[-1] == new_alt[-1]:
        new_ref = new_ref[:-1]
        new_alt = new_alt[:-1]
    return new_ref, new_alt



def get_adjusted_ref_alt(var, ref):
    target_ref = ref.REF
    origin_ref = ref.ALT[0]
    var_ref = var.REF
    var_alt = var.ALT[0]
    var_base_char = var_ref[0]
    ref_base_char = target_ref[0]
    if target_ref[0] != origin_ref[0]: # if the ref change is a snp, then the initial base in the ref and var alleles change
        ref_base_char = target_ref[0]
        var_base_char = target_ref[0]
        if var_ref[0] != var_alt[0]: # if both are snps, then the ref snp just gets swapped out
            return (target_ref, var.ALT[0])
    elif var_ref[0] != var_alt[0]:
        # if the variant is a snp, but the reference change is an indel, then the reference change does not affect the base snp
        # so we can simply return the variant unchanged.
        return (var_ref, var_alt)
    var_ref = var_ref[1:]
    var_alt = var_alt[1:]
    target_ref = target_ref[1:]
    origin_ref = origin_ref[1:]
    new_ref = ref_base_char + target_ref + var_ref
    new_alt = var_base_char + origin_ref + var_alt
    while len(new_ref) > 1 and len(new_alt) > 1 and new_ref[-1] == new_alt[-1]:
        new_ref = new_ref[:-1]
        new_alt = new_alt[:-1]
    return (new_ref, new_alt)


def convert_VCF_to_intervaltree(vcf_file):
    tree_dict = dict()
    vcf = VCF(vcf_file, gts012=True, threads=2)
    for chrom in vcf.seqnames:
        tree_dict[chrom] = list()
    for var in tqdm(vcf, total=get_num_vars(vcf_file)):
        tree_dict[var.CHROM].append(var)
    sys.stderr.write('Organizing vcf containing variation between builds...\n')
    for chrom in tqdm(tree_dict.keys()):
        tree_dict[chrom] = it.IntervalTree.from_tuples([(var.start, var.end, var) for var in tree_dict[chrom]])
    return tree_dict


def iterate_over_positions(vcf):
    var = next(vcf)
    current_position = var.POS-1
    current_chrom = var.CHROM
    poslist = [var]
    for var in vcf:
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

unlifted_vcf = sys.argv[1]
ref_alignment = sys.argv[2]
out_vcf = sys.argv[3]
chain = sys.argv[4]
target_fasta = sys.argv[5]

if unlifted_vcf is None:
    raise AttributeError('need a vcf file to lift')
if ref_alignment is None:
    raise AttributeError("need a reference alignment")
if out_vcf is None:
    out_vcf = '/mnt/ssd/lalli/phasing_T2T/liftover_071224/lifted.bcf'
if chain is None:
    chain = '/dev/shm/chm13v2-hg38.over.chain'
if target_fasta is None:
    target_fasta='/dev/shm/GRCh38_full_analysis_set_plus_decoy_hla.fasta'


if out_vcf == '/dev/stdout':
    out_vcf = '-'
    mode='wbu'
else:
    mode=None # infer


def get_num_unique_positions(vcf_file: str) -> int:
    vcf = VCF(vcf_file, threads=2)
    positions = {x.POS for x in vcf}
    return len(positions)

degenerate_nucleotides='UWSMKRYBDHV*'

sys.stderr.write ('Loading target reference genome...\n')
target_genome_seq = {chrom: seq.seq.__str__().upper()
                     for chrom, seq in SeqIO.to_dict(SeqIO.parse(target_fasta, "fasta")).items()}

for chrom in target_genome_seq.keys():
    for nucleotide in degenerate_nucleotides:
        if nucleotide in target_genome_seq[chrom]:
            target_genome_seq[chrom] = target_genome_seq[chrom].replace(nucleotide, 'N')

sys.stderr.write ('Loading chainfile...\n')
ls = LiftOver(chain)

sys.stderr.write('Loading vcf containing variation between builds...\n')
ref_diffs = convert_VCF_to_intervaltree(ref_alignment)

invcf = VCF(unlifted_vcf, gts012=True, threads=2)

invcf.add_info_to_header({'ID': 'Original_Contig', 
                          'Description': 'Original contig of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Original_POS', 
                          'Description': 'Original position of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Original_REF',
                          'Description': 'Original reference sequence of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Original_ALT',
                          'Description': 'Original alt sequence of variant before liftover', 
                          'Type':'Character', 'Number': '1'})
invcf.add_info_to_header({'ID': 'Flipped_during_liftover',
                          'Description': 'REF/ALT were flipped during liftover. GTs were altered accordingly.', 
                          'Type':'Character', 'Number': '1'})

if os.path.exists(out_vcf):
    os.remove(out_vcf)

out = vcfOpen(out_vcf, invcf, mode=mode)


sys.stderr.write ('Lifting...\n')
lifted_variants_no_overlap = dict()
lifted_variants_with_overlap = dict()
unliftable=dict()
multiple_overlaps = dict()
ref_seq_problem_liftovers = dict()
ref_seq_problem_liftovers_no_overlap = dict()
reported = True
for pos, variants in tqdm(iterate_over_positions(invcf), total=get_num_unique_positions(unlifted_vcf)):
    positions = defaultdict(list)
    overlap=list()
    lifted_start=tuple()
    lifted_end=tuple()
    lifted_variants = list()
    for var in variants:
        try:
            lifted_start, lifted_end = perform_clean_liftover(var, ls, return_coordinates=True)
        except Exception:
            unliftable[pos] = variants
            break
        if len(lifted_start) == 0 or len(lifted_end) == 0:
            break
        else:
            positions[lifted_start[0]].extend((lifted_start[1], lifted_end[1]))
    if len(lifted_start) == 0 or len(lifted_end) == 0 or pos in unliftable.keys():
        unliftable[pos] = variants
        continue

    for chrom, lifted_poses in positions.items():
        overlap.extend(ref_diffs[chrom][min(lifted_poses):max(lifted_poses)])
    
    if len(overlap) > 1:
        variants = [revert_variant(var) for var in variants]
        multiple_overlaps[pos] = variants
        continue
    elif any(['N' in x.data.REF for x in overlap]):
        variants = [revert_variant(var) for var in variants]
        unliftable[pos] = variants
        continue
    elif len(overlap) == 0:
        try:
            lifted_variants.extend([check_var_ref(var) for var in variants])
        except ValueError:
            variants = [revert_variant(var) for var in variants]
            ref_seq_problem_liftovers[pos] = variants
            continue
    else:
        try:
            lifted_variants.extend([check_var_ref(x) for x in liftover_variants_adjusting_for_overlapping_ref_diff(variants, overlap[0].data)])
        except ValueError:
            variants = [revert_variant(var) for var in variants]
            ref_seq_problem_liftovers[pos] = variants
            continue
        except AssertionError:
            variants = [revert_variant(var) for var in variants]
            ref_seq_problem_liftovers[pos] = variants
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
