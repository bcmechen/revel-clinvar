import clinvar
import revel
import dosage_sensitivity
import acmg59_genes


# Load and process ClinVar data.
clinvar_data = clinvar.Clinvar('../data/raw/clinvar_20210710.vcf.gz')
clinvar_data.vcf_to_dataframe()
clinvar_data.clean_vcf()

# Load and process REVEL data.
revel_data = revel.Revel('../data/raw/revel_grch38_all_chromosomes.csv.zip')
revel_data.read_into_dataframe()

# Load and process ClinGen Dosage Sensitivity data.
dosage_data = dosage_sensitivity.DosageSensitivityMap(
    '../data/raw/ClinGen_gene_curation_list_GRCh37.tsv')
dosage_data.read_into_dataframe()
dosage_data.drop_columns([
    'cytoBand', 'Genomic Location', 'Haploinsufficiency PMID1',
    'Haploinsufficiency PMID2', 'Haploinsufficiency PMID3',
    'Haploinsufficiency PMID4', 'Haploinsufficiency PMID5',
    'Haploinsufficiency PMID6', 'Triplosensitivity Score',
    'Triplosensitivity Description', 'Triplosensitivity PMID1',
    'Triplosensitivity PMID2', 'Triplosensitivity PMID3',
    'Triplosensitivity PMID4', 'Triplosensitivity PMID5',
    'Triplosensitivity PMID6', 'Date Last Evaluated', 'Loss phenotype OMIM ID',
    'Triplosensitive phenotype OMIM ID'])

# Create a master table by joining ClinVar, REVEL, and Dosage data.
master_df = revel_data.dataframe.merge(
    clinvar_data.dataframe, how='inner',
    left_on=['chr', 'hg19_pos', 'ref', 'alt'],
    right_on=['CHROM', 'POS', 'REF', 'ALT'])

master_df = master_df.merge(dosage_data.dataframe, how='left', left_on='GENE',
                            right_on='#Gene Symbol')

# Filter ACMG59 genes.
master_df = master_df[master_df['GENE'].isin(acmg59_genes.acmg59_genes)]

# Drop columns.
master_df.drop(['CHROM', 'POS', 'REF', 'ALT', '#Gene Symbol', 'grch38_pos'],
               axis=1, inplace=True)

# Export processed data.
master_df.to_csv('../data/processed/master_table.csv.zip', index=False)
