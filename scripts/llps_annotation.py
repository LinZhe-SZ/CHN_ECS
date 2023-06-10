import pandas as pd
import argparse

parser = argparse.ArgumentParser("Annotate pathogenic genes with LLPS data.")
parser.add_argument("--table", '-i', help='Input pathogenic variant list with gene and protein annotation. eg."CHROM\tPOS\tREF\tALT\tSYMBOL\tConsequence\tAF\tUniprot_acc\tHGVSp_VEP\tVEP_canonical\tCLNSIG\tCLNSIGCONF')
parser.add_argument("--llps", '-l', help='Input LLPS database')
parser.add_argument("--output", '-o', help='Output VCF file')

args = parser.parse_args()

aa_mut = []
with open(args.table, 'r') as fh:
    for line in fh:
        fields = line.strip().split("\t")
        symbol = fields[4].split(",")
        cons = fields[5].split(',')
        uniprot_id = fields[7].split(',')
        hgvs_p = fields[8].split(',')
        canonical_flag = fields[9].split(',')

        for i, c in enumerate(canonical_flag):
            f = c.split("&")
            con = cons[i]

            for j, v in enumerate(f):
                if v == "YES":
                    id = uniprot_id[i].split("&")[j]
                    hgvs = hgvs_p[i].split("&")[j].split('.')[1]
                    s = symbol[i]
                    c = con.split("&")
                    
                    aa_mut.append([fields[0], fields[1], fields[2], fields[3],  fields[6], s, con, id, hgvs, fields[10], fields[11]])

aa_mut_df = pd.DataFrame(aa_mut)
aa_mut_df.columns = ["CHROM", "POS", "REF", "ALT", "AF", "SYMBOL", "Consequence", 'UniprotID', 'AAChange', 'CLNSIG', 'CLNSIGCONF']
aa_mut_df = aa_mut_df.drop_duplicates(subset=['UniprotID', 'AAChange'])

table_a_df = pd.read_excel(args.llps, sheet_name='A', header=0)
table_a_df.drop(columns=['Sequence'], inplace=True)

table_h_df = pd.read_excel(args.llps, sheet_name='H', header=0)
#table_h_df.columns = ['Protein', 'MIDs', 'LCSs', 'Mutations', 'Mendelian diseases', 'Cancers']

table_j_df = pd.read_excel(args.llps, sheet_name='J', header=0)

table_k_df = pd.read_excel(args.llps, sheet_name='K', header=0)

uid_list = table_a_df['Uniprot ID'].tolist()
pid_list = table_h_df['Protein'].tolist()

with open(args.output, 'w') as outfh:
    outfh.write("\t".join(['CHROM', 'POS', 'REF', 'ALT', 'AF', 'SYMBOL', 'Consequence', 'UniprotID', 'AAChange', 'Disease', 'Cancer', 'FoundFlags', 'CLNSIG', 'CLNSIGCONF']))
    for idx, row in aa_mut_df.iterrows():
        found = False
        if row['UniprotID'] in uid_list:
            protein_id = table_a_df.loc[table_a_df['Uniprot ID']==row['UniprotID'], 'Protein'].values[0]
            if protein_id in pid_list:
                query_result_df = table_h_df[protein_id == table_h_df['Protein']]
                
                if pd.notna(query_result_df['Mutations'].values[0]):
                    mut_list = [x.strip() for x in query_result_df['Mutations'].iloc[0].split(",")]
                    if row['AAChange'] in mut_list:
                        found = True
                        outfh.write(f"{row['CHROM']}\t{row['POS']}\t{row['REF']}\t{row['ALT']}\t{row['AF']}\t{row['SYMBOL']}\t{row['Consequence']}\t{row['UniprotID']}\t{row['AAChange']}\t{query_result_df.iloc[0]['Mendelian diseases']}\t{query_result_df.iloc[0]['Cancers']}\t{found}\t{row['CLNSIG']}\t{row['CLNSIGCONF']}\n")
        if not found:
            outfh.write(f"{row['CHROM']}\t{row['POS']}\t{row['REF']}\t{row['ALT']}\t{row['AF']}\t{row['SYMBOL']}\t{row['Consequence']}\t{row['UniprotID']}\t{row['AAChange']}\tnan\tnan\t{found}\t{row['CLNSIG']}\t{row['CLNSIGCONF']}\n")
