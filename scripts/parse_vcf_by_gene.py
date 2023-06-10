from cyvcf2 import VCF, Writer
import argparse

def parse_csq_to_dict(var, csq_fields_str):
	fields_name = csq_fields_str.split('|')
	assert fields_name[1]  == 'Consequence', "'Consequence' dose not in the second fields."

	# order of severity (more severe to less severe) http://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
	severity_order = {
		'transcript_ablation':	1,
		'splice_acceptor_variant':	2,
		'splice_donor_variant':	3,
		'stop_gained':	4,
		'frameshift_variant':	5,
		'stop_lost':	6,
		'start_lost':	7,
		'transcript_amplification':	8,
		'inframe_insertion':	9,
		'inframe_deletion':	10,
		'missense_variant':	11,
		'protein_altering_variant':	12,
		'splice_region_variant':	13,
		'splice_donor_5th_base_variant':	14,
		'splice_donor_region_variant':	15,
		'splice_polypyrimidine_tract_variant':	16,
		'incomplete_terminal_codon_variant':	17,
		'start_retained_variant':	18,
		'stop_retained_variant':	19,
		'synonymous_variant':	20,
		'coding_sequence_variant':	21,
		'mature_miRNA_variant':	22,
		'5_prime_UTR_variant':	23,
		'3_prime_UTR_variant':	24,
		'non_coding_transcript_exon_variant':	25,
		'intron_variant':	26,
		'NMD_transcript_variant':	27,
		'non_coding_transcript_variant':	28,
		'upstream_gene_variant':	29,
		'downstream_gene_variant':	30,
		'TFBS_ablation':	31,
		'TFBS_amplification':	32,
		'TF_binding_site_variant':	33,
		'regulatory_region_ablation':	34,
		'regulatory_region_amplification':	35,
		'feature_elongation':	36,
		'regulatory_region_variant':	37,
		'feature_truncation':	38,
		'intergenic_variant':	39
	}

	csq_list = var.INFO['CSQ'].split(',')

	most_severity = [x if x != "" else None for x in csq_list[0].split('|')]
	severity_level = severity_order[most_severity[1].split('&')[0]]
	for csq in csq_list:
		value = csq.split('|')
		for consequence in value[1].split('&'):
			if severity_order[consequence] < severity_level:
				most_severity = value
				most_severity[1] = consequence
				severity_level = severity_order[most_severity[1]]

	assert len(fields_name) == len(most_severity), f"Numbers of CSQ fields({len(fields_name)}) defined in header dose not equal to CSQ value({len(most_severity)}) of variant: {var.CHROM}:{var.POS}:{var.REF}:{var.ALT}"

	most_severity_dict = dict(zip(fields_name, most_severity))


	return most_severity_dict
	
parser = argparse.ArgumentParser("Parse variant from VCF by gene name list.")
parser.add_argument("--input", '-i', help='Input VCF file.')
parser.add_argument("--genes", '-g', help='Input gene list.')
parser.add_argument("--output", '-o', help='Output VCF file')

args = parser.parse_args()

invcf=VCF(args.input, strict_gt=True)
outvcf = Writer(args.output, invcf)
info_csq_str = invcf.get_header_type('CSQ')['Description']

gene_list = []
gene_dict = {}
n_genes = 0
with open(args.genes, 'r') as gene_list_fh:
	for line in gene_list_fh:
		genes = line.replace("\s+", "").strip().split(",")
		if len(genes) == 1:
			gene_dict[genes[0]] = genes[0]
		else:
			for g in genes:
				gene_dict[g] = genes[0]
		n_genes += 1
gene_list = gene_dict.keys()

print(f"Parse variants from {args.input} by {n_genes} genes of {args.genes}")
found_by_symbol = 0
found_by_gene = 0
found_by_geneinfo = 0
not_in_gene_list = 0

found_genes = {}
notfound_genes = {}
with open(args.output, 'w') as p_out_fh:
	for v in invcf:
		csq = parse_csq_to_dict(v, info_csq_str)
		geneinfo = v.INFO.get("GENEINFO")
		if csq['SYMBOL'] in gene_list:
			outvcf.write_record(v)
			found_by_symbol += 1
			found_genes[csq['SYMBOL']] = 1
		elif csq['Gene'] in gene_list:
			outvcf.write_record(v)
			found_by_gene += 1
			found_genes[csq['Gene']] = 1
		elif geneinfo and geneinfo.split(":")[0] in gene_list:
			outvcf.write_record(v)
			found_by_geneinfo += 1
			found_genes[geneinfo.split(":")[0]] = 1
		else:
			notfound_genes[csq['SYMBOL']] = 1
			not_in_gene_list += 1

print(f"#Genes: {n_genes}")
print(f"#Vars found by SYMBOL: {found_by_symbol}")
print(f"#Vars found by Gene: {found_by_gene}")
print(f"#Vars found by GENEINFO: {found_by_geneinfo}")
print(f"#Found genes: {found_genes.keys()}")
print(f"#Not Found genes: {notfound_genes.keys()}")
