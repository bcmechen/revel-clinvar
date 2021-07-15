import gzip
import re
import pandas as pd
import os


class Clinvar:
    """To model ClinVar data.

    Attributes:
        filepath (str): A path to a file.
        dataframe (obj): A Pandas DataFrame object.
    """

    # Class variable.
    # Fields that are commonly useful for analyzing ClinVar data.
    commonly_used_columns = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT', 'CLNHGVS', 'CLNREVSTAT',
        'CLNSIG', 'CLNVC', 'ORIGIN', 'CLNSIGCONF', 'CLASS', 'GENE',
        'MOLECULAR CONSEQUENCE', 'IS_CONFLICTING'
    ]

    def __init__(self, filepath, dataframe=None):
        """Initialize attributes."""
        self.filepath = filepath
        self.dataframe = dataframe

    @property
    def filename(self):
        """Get filename without the extension.

        Returns:
            (str): Filename without extension.
        """
        return os.path.basename(self.filepath).split('.')[0]

    @staticmethod
    def parse_vcf_meta_info(file_obj):
        """Parse ID and Description from VCF meta-information lines and
        store them in a dictionary format.

        VCF file meta-information is included after the ## string and must
        be key=value pairs. For example:
        ##fileformat=VCFv4.1
        ##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar
        Allele ID">
        ##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical
        significance for this single variant">

        Args:
            file_obj (obj): A file object from opening a VCF file.

        Returns:
            meta_info (dict): A dictionary containing ID and Description
            from ##INFO of meta-information. For example:
                {'ALLELEID': 'the ClinVar Allele ID',
                 'CLNSIG': 'Clinical significance for this single variant'}
            total_meta_rows (int): Total number of rows in meta-information.
        """
        total_meta_rows = 0
        meta_info = {}
        for line in file_obj:
            if line.startswith('##'):
                total_meta_rows = total_meta_rows + 1
            if line.startswith('##INFO='):
                # Find ID and store it as dict key.
                # Find Description and store it as dict value.
                meta_info[re.search('ID=(\w+)', line).group(1)] = re.search(
                    'Description="(.*)"', line).group(1)

        return meta_info, total_meta_rows

    def get_vcf_meta_info(self):
        """Get ID and Description from VCF meta-information in a dictionary
        format."""
        # Check if file is in gzip format.
        if 'gz' in self.filepath:
            with gzip.open(self.filepath, 'rt') as f:
                return self.parse_vcf_meta_info(f)
        else:
            with open(self.filepath) as f:
                return self.parse_vcf_meta_info(f)

    def vcf_to_dataframe(self, num_of_rows=None):
        """Read ClinVar VCF data into a Pandas DataFrame.

        VCF is a text file format. It contains meta-information lines,
        a header line, and then data lines each containing information about
        a position in the genome. File meta-information is included after
        the ## string and must be key=value pairs. The header line names the 8
        fixed, mandatory columns: #CHROM, POS, ID, REF, ALT, QUAL, FILTER,
        INFO.

        Args:
            num_of_rows (int): Number of rows of file to read.
        """
        # Include CHROM, POS, ID, REF, ALT, and INFO fields and exclude QUAL
        # and FILTER fields.
        # Set data type of the CHROM column to object.
        # Note: Can't use the comment parameter to skip the comment rows
        # because certain INFO also contains the symbol "#".
        self.dataframe = pd.read_csv(self.filepath, sep='\t',
            nrows=num_of_rows, header=self.get_vcf_meta_info()[1],
            names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO'],
            usecols=[0, 1, 2, 3, 4, 7], dtype={0: 'object'})

    @staticmethod
    def convert_info_string_to_dict(info_string):
        """A function to convert a string from the INFO field to a dictionary.

        INFO fields are encoded as a semicolon-separated series of short
        keys with optional values in the format: <key>=<data>[,data]. For
        example:

        ALLELEID=626482;CLNDISDB=MedGen:C3808739,OMIM:615120;CLNSIG=Uncertain_significance

        Args:
            info_string (str): A string from the INFO field.

        Returns:
            A dict. For example:

            {'ALLELEID': 626482,
             'CLNSIG': 'Uncertain_significance'}
        """
        info_dict = {}
        # Iterate through all key=value pairs in the string.
        for pair in info_string.split(';'):
            # Convert key=value into the dictionary format.
            info_dict[pair.split('=')[0]] = pair.split('=')[1]
        return info_dict

    def split_info_columns(self):
        """Split INFO column into individual columns using keys as column
        names.

        First convert the INFO field format from a string to a dictionary.
        Then split the INFO field into multiple columns using keys as new
        column names. For example:

        {'ALLELEID': 626482, 'CLNSIG': 'Uncertain_significance'} would
        create a column named ALLELEID and a column named CLNSIG.
        """
        # Convert the INFO field format from a string to a dictionary.
        self.dataframe['INFO'] = self.dataframe['INFO'].apply(
            self.convert_info_string_to_dict)
        # Split the INFO field into columns.
        self.dataframe = pd.concat(
            [self.dataframe, pd.DataFrame(self.dataframe['INFO'].tolist())],
            axis=1)

    def create_gene_column(self):
        """Create a new column named GENE.

        GENEINFO column split from INFO field has a gene symbol and a gene
        id separated by a colon. For example: AGRN:375790.
        This function creates a new GENE column that has only gene symbol (
        leaves out gene id).
        """
        # Check if GENEINFO column exists.
        if 'GENEINFO' in self.dataframe.columns:
            self.dataframe['GENE'] = self.dataframe['GENEINFO'].astype(
                'str').apply(lambda x: x.split(':')[0])

    def create_mol_conseq_column(self):
        """Create a new column named MOLECULAR CONSEQUENCE.

        Original MC column is a comma separated list of molecular
        consequence in the form of Sequence Ontology ID|molecular_consequence.
        This function creates a new MOLECULAR CONSEQUENCE column that has
        only molecular consequence (leaves out SO and ID). If there is more
        than one molecular consequence, the function returns 'uncertain'.
        For example:

        SO:0001819|synonymous_variant,SO:0001623|5_prime_UTR_variant returns 'uncertain'.
        SO:0001583|missense_variant returns 'missense_variant'.
        """
        # Check if MC column exists.
        if 'MC' in self.dataframe.columns:
            # Check if there is more than one molecular consequence.
            self.dataframe['MOLECULAR CONSEQUENCE'] = self.dataframe[
                'MC'].astype('str').apply(lambda x: x.split('|')[-1] if len(
                x.split('|')) == 2 else 'uncertain')

    @staticmethod
    def simplify_class_names(clnsig_name):
        """Simplify classification names used by ClinVar to Benign, VUS,
        Pathogenic, and Other.

        Args:
            clnsig_name (str): A string from the CLNSIG field.

        Returns:
            simplified_class_name (str): Benign, VUS, Pathogenic, or Other.
        """
        # Group ClinVar classification names into 3 classes.
        benign_class_names = ['Likely_benign', 'Benign',
                              'Benign/Likely_benign']
        vus_class_names = ['Uncertain_significance']
        pathogenic_class_names = ['Pathogenic', 'Likely_pathogenic',
                                  'Pathogenic/Likely_pathogenic',
                                  'risk_factor', 'Pathogenic,_risk_factor',
                                  'Likely_pathogenic,_risk_factor',
                                  'Pathogenic/Likely_pathogenic,_risk_factor']

        if clnsig_name in benign_class_names:
            simplified_class_name = 'Benign'
        elif clnsig_name in vus_class_names:
            simplified_class_name = 'VUS'
        elif clnsig_name in pathogenic_class_names:
            simplified_class_name = 'Pathogenic'
        else:
            simplified_class_name = 'Other'

        return simplified_class_name

    def create_class_column(self):
        """Create a new column named CLASS for simplified classification
        names."""
        # Check if CLNSIG column exists.
        if 'CLNSIG' in self.dataframe.columns:
            self.dataframe['CLASS'] = self.dataframe['CLNSIG'].apply(
                self.simplify_class_names)

    def create_is_conflicting_column(self):
        """Create a new column named IS_CONFLICTING that has boolean value.

        There is a key=value pair
        CLNSIG=Conflicting_interpretations_of_pathogenicity in the original
        INFO field when there are conflicts in classification for a variant.
        This function creates a new IS_CONFLICTING column that has a boolean
        value to indicate whether the variant has conflicting
        classifications (True) or not (False).
        """
        # Check if CLNSIGCONF column exists.
        if 'CLNSIGCONF' in self.dataframe.columns:
            self.dataframe['IS_CONFLICTING'] = self.dataframe[
                'CLNSIGCONF'].notna()
        else:
            self.dataframe['IS_CONFLICTING'] = False

    def drop_columns(self, column_names=None):
        """Drop columns from dataframe.

        If a list of column names are provided, drop those columns.
        If not, drop all columns that are not in Clinvar.commonly_used_columns.

        Args:
            column_names (list): A list of names of columns.
        """
        if column_names is not None:
            self.dataframe.drop(column_names, axis=1, inplace=True)
        else:
            # Check all columns in the dataframe and create a list to store
            # all column names that are not in the commonly used list.
            columns_to_drop = [x for x in
                               self.dataframe.columns if
                               x not in Clinvar.commonly_used_columns]
            self.dataframe.drop(columns_to_drop, axis=1, inplace=True)

    def export_dataframe_to_csv(self, csv_filename=None):
        """Export dataframe to a csv file.

        If a filename is provided, use the filename. If not, the csv will
        use the same filename as the VCF with _processed suffixed.

        Args:
            csv_filename (str): Filename for the csv file.
        """
        if csv_filename is not None:
            self.dataframe.to_csv('{}.csv'.format(csv_filename), index=False)
        else:
            self.dataframe.to_csv('{}_processed.csv'.format(self.filename),
                                  index=False)

    def clean_vcf(self):
        """Process and clean VCF data.

        This function processes the INFO field, creates new columns,
        and drop unwanted columns.
        The output dataframe has this list of columns:

        ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'CLNHGVS', 'CLNREVSTAT', 'CLNSIG',
       'CLNVC', 'ORIGIN', 'CLNSIGCONF', 'GENE', 'MOLECULAR CONSEQUENCE',
       'CLASS', 'IS_CONFLICTING']
        """
        # Split INFO field into columns.
        self.split_info_columns()
        # Create GENE column.
        self.create_gene_column()
        # Create MOLECULAR CONSEQUENCE column.
        self.create_mol_conseq_column()
        # Create CLASS column.
        self.create_class_column()
        # Create IS_CONFLICTING column.
        self.create_is_conflicting_column()
        # Drop all columns that are not in the commonly used list.
        self.drop_columns()

    def left_join(self, df, df_columns, clinvar_columns):
        """Join DataFrames with ClinVar data as the left table.

        Args:
            df (obj): A Pandas DataFrame object to join to ClinVar data.
            df_columns (list): Column names to join on in the other DataFrame.
            clinvar_columns (list): Column names to join on in the ClinVar
            DataFrame.

        Returns:
            (obj): A Pandas DataFrame object from joined DataFrames.
        """
        return self.dataframe.merge(df, how='left', left_on=clinvar_columns,
                                    right_on=df_columns)


def main():
    # Execute module as script.
    import sys
    # Create a Clinvar instance.
    clinvar_data = Clinvar(sys.argv[1])
    try:
        # If number of rows of file to read is provided.
        clinvar_data.vcf_to_dataframe(int(sys.argv[2]))
    except IndexError:
        # If number of rows of file to read is not provided.
        clinvar_data.vcf_to_dataframe()
    # Process and clean ClinVar data.
    clinvar_data.clean_vcf()
    # Export dataframe to a csv file.
    clinvar_data.export_dataframe_to_csv()


if __name__ == '__main__':
    main()
