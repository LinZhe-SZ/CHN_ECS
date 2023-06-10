from cyvcf2 import VCF, Writer
import argparse

from autopvs1 import AutoPVS1
import pandas as pd

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def calc_gcr(var_freq, gene_colname = 'SYMBOL'):
	var_freq_cp = var_freq.copy()
	var_freq_cp.loc[:, 'NormalFreq'] = 1 - var_freq_cp['AF']*2 # use allele frequency to calculate VCR
	if gene_colname not in var_freq_cp:
		gene_colname = 'Gene'

	gcr = var_freq_cp.groupby(gene_colname,as_index=False).NormalFreq.prod()
	gcr.loc[:, 'GCR'] =  1 - gcr.NormalFreq
	gcr.set_index(gene_colname)
	return gcr

def calc_ccr(gcr):
	cr = 1
	for i_gcr in gcr.NormalFreq:
		cr *= i_gcr
	ccr =  1 - cr
	return ccr

class VariantPredicate:
	def __init__(self, var, csq_fields_str, an_name = 'AN', ac_name = 'AC', af_name = 'AF', n_tools_cutoff = 9, missense_predict_tools = ['CADD', 'Eigen', 'REVEL', 'DANN', 'Polyphen2', 'SIFT', 'MetaSVM', 'MutationAssessor', 'PROVEAN'], output_csq_fields = ['rs_dbSNP', '1000Gp3_EAS_AC', '1000Gp3_EAS_AF', 'gnomAD_exomes_EAS_AF', 'gnomAD_genomes_EAS_AF', 'gnomAD_genomes_controls_and_biobanks_EAS_AF', 'gnomAD_exomes_controls_POPMAX_AF', 'gnomAD_genomes_POPMAX_AF', 'DisGeNET', 'LoF']) -> None:
		self.chr = var.CHROM
		self.pos = var.POS
		self.ref = var.REF
		self.alt = var.ALT[0]
		self.var = var
		
		self.an = var.INFO[an_name]
		self.ac = var.INFO[ac_name]
		self.af = var.INFO[af_name]

		self.af_dict = {}
		self.af_dict[af_name] = None
		if var.INFO[af_name] is not None:
			self.af_dict[af_name] = var.INFO[af_name]
		else:
			print(f"{af_name} dose not found in {var}")

		self.output_csq_fields = output_csq_fields
			
		self.pathogenicity = 'Unknown'
		self.missense_n_tools = None
		self.gene = None
		self.HGVSc = None
		self.consequence = None

		self.n_tools_cutoff = min(n_tools_cutoff, len(missense_predict_tools))
		self.failed_tools = None
		self.passed_tools = None

		self.csq = VariantPredicate.parse_csq_to_dict(var, csq_fields_str)

		if self.csq['SYMBOL']: self.gene = self.csq['SYMBOL']
		if self.csq['HGVSc']: self.HGVSc = self.csq['HGVSc']
		if self.csq['Consequence']: self.consequence = self.csq['Consequence']

		self.CLNSIG = None
		self.CLNREVSTAT = None
		self.CLNDN = None
		self.GENEINFO = None
		if self.var.INFO.get('CLNSIG') and self.var.INFO.get('CLNREVSTAT'):
			self.CLNSIG = self.var.INFO.get('CLNSIG')
			self.CLNREVSTAT = self.var.INFO.get('CLNREVSTAT')
			self.CLNDN = self.var.INFO.get('CLNDN')
			self.GENEINFO = self.var.INFO.get('GENEINFO')
		elif 'clinvar_clnsig' in self.csq and 'clinvar_review' in self.csq:
			if self.csq['clinvar_clnsig']: self.CLNSIG = self.csq['clinvar_clnsig']
			if self.csq['clinvar_review']: self.CLNREVSTAT = self.csq['clinvar_review']
		else:
			raise "No clinvar annotation found: {var.CHROM}:{var.POS}:{var.REF}:{var.ALT}"

		if not self.gene:
			if self.GENEINFO:
				self.gene = self.GENEINFO.split(":")[0]			
			elif self.csq['Gene']:
				self.gene = self.csq['Gene']
			elif self.csq['CANONICAL']:
				self.gene = self.csq['CANONICAL']
			else:
				print(f"Error: no gene name found: {self.var.CHROM}:{self.var.POS}:{self.var.REF}:{self.var.ALT} - {self.gene}|{self.csq['Gene']}|{self.GENEINFO}")

		# missense tools
		# !IMPORTANT: for float cutoff value, here assume the value lower than the cutoff is not damage, but maybe some other tools are not consisit to this assumption?
		default_missense_tools_cutoff = {
			'CADD_phred': 			'20', 	# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006481#:~:text=CADD%20predicts%20a%20continuous%20phred,by%20the%20authors%20i.e.%2020.
			'Eigen-PC-phred_coding':	'1', 	# https://www.theanalysisfactor.com/factor-analysis-how-many-factors/#:~:text=Eigenvalue%20%3E%201,an%20eigenvalue%20of%20%E2%89%A51.
			'REVEL_rankscore':		'0.75', 	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5065685/#:~:text=Selecting%20a%20more%20stringent%20REVEL,ESVs%20being%20classified%20as%20pathogenic.
			'DANN_score':			'0.5', 	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8682775/
			'Polyphen2_HDIV_pred':	['D'], 
			'SIFT_pred':			['D'],
			'MetaSVM_pred':			['D'],
			'MutationAssessor_pred':	['H', 'M'],
			'PROVEAN_pred':			['D'],
		}

		short_name_to_vep_field = {
			'CADD': 	'CADD_phred', 	
			'Eigen':	'Eigen-PC-phred_coding', 	
			'REVEL':	'REVEL_rankscore', 	
			'DANN':		'DANN_score', 	
			'Polyphen2':	'Polyphen2_HDIV_pred', 
			'SIFT':		'SIFT_pred',
			'MetaSVM':	'MetaSVM_pred',
			'MutationAssessor':	'MutationAssessor_pred',
			'PROVEAN':	'PROVEAN_pred',
		}

		self.missense_tools_name = []
		self.missense_tools_cutoff = {}


		for item in missense_predict_tools:
			if ':' not in item:
				tools = short_name_to_vep_field[item]
				cutoff = default_missense_tools_cutoff[tools]
			else:
				(tools, cutoff) = item.split(':')
				assert tools in self.csq, f"Missense predicte tools name '{tools}' dose not exists in CSQ"
			
			self.missense_tools_name.append(tools)
			self.missense_tools_cutoff[tools] = cutoff
		
		self.predicated_pathogenicity()
		if self.var.FILTER is not None:
			self.var.FILTER = self.var.FILTER + ',' + self.pathogenicity
		else:
			self.var.FILTER = self.pathogenicity

		if self.missense_n_tools:
			self.var.INFO['N_TOOLS'] = self.missense_n_tools

	@staticmethod
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


	def clinvar_judge(self):
		clnsig = 'UNKNOWN'
		if self.CLNREVSTAT and 'no_assertion_criteria_provided' not in self.CLNREVSTAT.split(','):
			if self.CLNSIG and self.CLNSIG.startswith('Benign') or self.CLNSIG == 'Likely_benign':
				clnsig = 'BLB'
			elif self.CLNSIG and self.CLNSIG.startswith('Pathogenic') or self.CLNSIG == 'Likely_pathogenic':
				clnsig = 'PLP'
		
		return clnsig

	def is_missense(self):
		if self.csq['Consequence'].startswith('missense_variant'):
			return True
		else:
			return False

	def missense_predict(self):
		assert self.csq['Consequence'].startswith('missense_variant'), f"Trying to predict a non-missense variant: {var}"
		
		if self.n_tools_cutoff is None:
			self.n_tools_cutoff = len(self.missense_tools_name) 

		assert self.n_tools_cutoff <= len(self.missense_tools_name), f"The n tools{self.n_tools_cutoff} cutoff is greater than number of tools{len(self.missense_tools_name)}"

		i_damage = 0
		for tool_name in self.missense_tools_name:
			cutoff = self.missense_tools_cutoff[tool_name]
			assert tool_name in self.csq, f"{tool_name} is not in CSQ annotation."
			if not self.csq[tool_name]: continue
			for val in self.csq[tool_name].split('&'):
				if not isinstance(cutoff,list):
					if float(val) >= float(cutoff):
						i_damage += 1
						self.passed_tools = f"{self.passed_tools},{tool_name}:{val}" if self.passed_tools else f"{tool_name}:{val}"
						break
					else:
						self.failed_tools = f"{self.failed_tools},{tool_name}:{val}" if self.failed_tools else f"{tool_name}:{val}"
				else:
					if val in cutoff:
						i_damage += 1
						self.passed_tools = f"{self.passed_tools},{tool_name}:{val}" if self.passed_tools else f"{tool_name}:{val}"
						break
					else:
						self.failed_tools = f"{self.failed_tools},{tool_name}:{val}" if self.failed_tools else f"{tool_name}:{val}"
		if self.passed_tools: self.passed_tools.rstrip(',')
		if self.failed_tools: self.failed_tools.rstrip(',')

		self.missense_n_tools = i_damage
		if i_damage >= self.n_tools_cutoff:
			return 'DAMAGE'
		else:
			return 'NONDAMAGE'

	def is_nonsense(self):
		for c in self.csq['Consequence'].split('&'):
			if c in ['frameshift_variant', 'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant', 'start_lost']:
				return True
		return False

	def lof_predict(self):
		if self.csq['LoF']:
			return self.csq['LoF']
		else:
			return "UNKNOWN"

	def predicated_pathogenicity(self):
		predicated = self.clinvar_judge()

		if predicated in ['PLP', 'BLB']:
			self.pathogenicity = f"CLINVAR_{predicated}"
		elif self.is_missense():
			predicated = self.missense_predict()
			self.pathogenicity = f"MISSENSE_{predicated}"
		elif self.is_nonsense():
			predicated = self.lof_predict()
			self.pathogenicity = f"NONSENSE_{predicated}"
		else:
			self.pathogenicity = "NONPATHOGENETIC"
	
	def __repr__(self):
		#af_str = "\t".join([str(self.af_dict[x]) for x in self.af_name])
		csq_str = "\t".join([str(self.csq[x]) if self.csq[x] != "" else "None" for x in self.output_csq_fields])
		format_output = f"{self.chr}\t{self.pos}\t{self.ref}\t{self.alt}\t{self.an}\t{self.ac}\t{self.af}\t{self.gene}\t{self.HGVSc}\t{self.consequence}\t{self.pathogenicity}\t{self.CLNSIG}\t{self.CLNREVSTAT}\t{self.CLNDN}\t{self.GENEINFO}\t{self.missense_n_tools}\t{self.passed_tools}\t{csq_str}"
		return format_output

parser = argparse.ArgumentParser(
	"""
	Given a VCF annotated by VEP, this script will annotation pathogenetic of each variants by the method described in this paper:
	""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("--invcf", "-i", help="Input annotation file in VCF format", required=True)
parser.add_argument("--outvcf", "-o", help="Output annotation file in VCF format")
parser.add_argument("--output-table", "-p", help="Output table file include infomation of each variant. format: CHROM\tPOS\tREF\tALT\tSYMBOL\tTYPE\tAF", required=True)
parser.add_argument("--output-gcr", "-g", help="Output table file include GCR of each gene, format: SYMBOL\tGCR", required=True)
parser.add_argument("--missense-predicate-tools", "-t", help="Tools names used for missense predicate", nargs = '*', default=['CADD', 'Eigen', 'REVEL', 'DANN', 'Polyphen2', 'SIFT', 'MetaSVM', 'MutationAssessor', 'PROVEAN'])
parser.add_argument("--n-tools-cutoff", "-n", help="At least this many tools predicate as damage will be accept as damage", default=9)
parser.add_argument("--an-name", '-an', help="A list of tag name of Alelle number which will be output in table", default='AN')
parser.add_argument("--ac-name", '-ac', help="A list of tag name of Alt Alelle count which will be output in table", default='AC')
parser.add_argument("--af-name", '-af', help="A list of tag name of Alelle freuquency which will be used to calculate allele carrying rate", default='AF')
parser.add_argument("--csq-fields", '-c', help="A list of tag names in CSQ, which will parse and output to the table file.", nargs = '+', default = ['rs_dbSNP', '1000Gp3_EAS_AC', '1000Gp3_EAS_AF', 'gnomAD_exomes_EAS_AF', 'gnomAD_genomes_EAS_AF', 'gnomAD_genomes_controls_and_biobanks_EAS_AF', 'gnomAD_exomes_controls_POPMAX_AF', 'gnomAD_genomes_POPMAX_AF', 'DisGeNET', 'LoF'])
parser.add_argument("--filter-by-clinvar-top1", '-f', action='store_true', help="Filter out variants of a gene, which AF is greater than max AF of CLINVAR_PLP of this gene", default=False)
parser.add_argument("--filter-by-autopvs1", '-v', action='store_true', help="Filter out nonsense variant by AutoPVS1", default=False)

args = parser.parse_args()

invcf=VCF(args.invcf, strict_gt=True)
invcf.add_info_to_header({'ID': 'N_TOOLS', 'Description': 'How many tools predicated a missense as damage', 'Type':'Integer', 'Number': '1'})
invcf.add_info_to_header({'ID': 'AF', 'Description': 'Allele frequency', 'Type':'Float', 'Number': 'A'})
invcf.add_filter_to_header({'ID': 'NONPATHOGENETIC', 'Description': 'Variant do not meet any condition of the algorithm.'})
invcf.add_filter_to_header({'ID': 'CLINVAR_PLP', 'Description': 'Variant can be determined as pathogenic/likey pathogenic by Clinvar.'})
invcf.add_filter_to_header({'ID': 'CLINVAR_BLB', 'Description': 'Variant can be determined as benign/likey benign by Clinvar.'})
invcf.add_filter_to_header({'ID': 'MISSENSE_DAMAGE', 'Description': 'Missense variant predicated as damage'})
invcf.add_filter_to_header({'ID': 'MISSENSE_NONDAMAGE', 'Description': 'Missense variant predicated as non-damage'})
invcf.add_filter_to_header({'ID': 'NONSENSE_HC', 'Description': 'Nonsense variant predicated as high confident LoF'})
invcf.add_filter_to_header({'ID': 'NONSENSE_LC', 'Description': 'Nonsense variant predicated as low confident LoF'})
invcf.add_filter_to_header({'ID': 'NONSENSE_UNKNOWN', 'Description': 'Failed to predicated nonsense variant'})

if args.outvcf:
	outvcf = Writer(args.outvcf, invcf)


info_csq_str = invcf.get_header_type('CSQ')['Description']
csq_fields_str = info_csq_str.split(':')[1].lstrip().rstrip('"')

max_clnaf = {}

vars_list = []

with open(args.output_table, 'w') as p_out_fh:
	for v in invcf:
		if not v.INFO.get(args.af_name):
			for k in v.INFO:
				print(f"{k}")
			continue

		var = VariantPredicate(v, csq_fields_str, args.an_name, args.ac_name, args.af_name, args.n_tools_cutoff, args.missense_predicate_tools)

		if args.outvcf:
			outvcf.write_record(var.var)

		if var.pathogenicity == 'CLINVAR_PLP':
			if var.gene not in max_clnaf: max_clnaf[var.gene] = var.af
			elif var.af > max_clnaf[var.gene]: max_clnaf[var.gene] = var.af

		if var.gene:
			print(var, file=p_out_fh)
p_out_fh.close()

table_df = pd.read_csv(args.output_table, header=None, sep="\t")
columns_name = ["Chr", "Pos", "Ref", "Alt", "AN", "AC", "AF", "Gene", "HGVSc", "Consequence", "Pathogenicity", "CLNSIG", "CLNREVSTAT", "CLNDN", "GENEINFO", "missense_n_tools", "pass_tools"]
columns_name.extend(args.csq_fields)
table_df.columns = columns_name
table_df['MoreFilter'] = 'PASS'

final_var_af = []
for idx, var in table_df.iterrows():
	morefilter = '.'
	if args.filter_by_clinvar_top1 and var.Gene in max_clnaf and var.AF > max_clnaf[var.Gene]:
		morefilter = 'CLINVARTOP1'
		if var.CLNSIG.startswith("Conflict"):
			morefilter = 'CLINVARTOP1_CONFLICT'

	if var.Gene not in max_clnaf and var.AF > 0.005:
		morefilter += ',COMMONAF'

	if args.filter_by_autopvs1 and var.Pathogenicity == 'NONSENSE_HC':
		var_input = '%s-%d-%s-%s' % (var.Chr.replace('chr', ''), var.Pos, var.Ref, var.Alt)
		apvs1 = AutoPVS1(var_input, 'GRCh38')
		if apvs1.islof:
			if str(apvs1.pvs1.strength) not in ["Strength.Strong", "Strength.VeryStrong"]:
				morefilter += f',AutoPVS1_FILTER,{str(apvs1.pvs1.strength)}'
			else:
				morefilter += f',{str(apvs1.pvs1.strength)}'
			print(f"{var_input}:{morefilter}")
		else:
			morefilter += ',AutoPVS1_NONLOF'

	if morefilter.startswith('.,'): morefilter = morefilter.lstrip(".,")
	table_df.iloc[idx, -1] = morefilter

	if var.Pathogenicity in ['CLINVAR_PLP', 'MISSENSE_DAMAGE', 'NONSENSE_HC'] and morefilter == '.':
		final_var_af.append([var.Gene, var.AF])
table_df.to_csv(args.output_table, sep="\t", index=None)

if len(final_var_af) != 0:
	with open(args.output_gcr, 'w') as g_out_fh:
		final_var_af_df = pd.DataFrame(final_var_af, columns=['SYMBOL', 'AF'] )
		gcr = calc_gcr(final_var_af_df)
		for idx, row in gcr.iterrows():
			print(f"{row['SYMBOL']}\t{row['GCR']}", file=g_out_fh)
else:
	with open(args.output_gcr, 'w') as g_out_fh:
		print(f"", file=g_out_fh)


if args.outvcf:
	outvcf.close()

invcf.close()
