import os
import pandas as pd


class Revel:
    """To model REVEL data.

    REVEL is an ensemble method for predicting the pathogenicity of missense
    variants based on a combination of scores from 13 individual tools:
    MutPred, FATHMM v2.3, VEST 3.0, PolyPhen-2, SIFT, PROVEAN,
    MutationAssessor, MutationTaster, LRT, GERP++, SiPhy, phyloP,
    and phastCons (PMID: 27666373).

    REVEL file is in comma-separated value (CSV) format and contain 7 fields:
        chr: Chromosome
        hg19_pos: Sequence position on the hg19 (GRCh37) human genome build
        (1-based coordinates)
        grch38_pos: Sequence position on the GRCh38 human genome build (
        1-based coordinates)
        ref: Reference nucleotide
        alt: Alternate nucleotide
        aaref: Reference amino acid
        aaalt: Alternate amino acid
        REVEL: REVEL score

    Attributes:
        filepath (str): A path to a file.
        dataframe (obj): A Pandas DataFrame object.
    """
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

    def read_into_dataframe(self, num_of_rows=None):
        """Read REVEL data into a pandas DataFrame.

        Args:
            num_of_rows (int): Number of rows of file to read.
        """
        # Set data type of the chr and grch38_pos column to object.
        self.dataframe = pd.read_csv(self.filepath, nrows=num_of_rows,
                                     dtype={0: 'object', 2: 'object'})

    def drop_columns(self, column_names=None):
        """Drop columns from dataframe.

        Args:
            column_names (list): A list of names of columns.
        """
        self.dataframe.drop(column_names, axis=1, inplace=True)

    def export_dataframe_to_csv(self, csv_filename=None):
        """Export dataframe to a csv file.

        If a filename is provided, use the filename. If not, the csv will
        use the same filename as the source with _processed suffixed.

        Args:
            csv_filename (str): Filename for the csv file.
        """
        if csv_filename is not None:
            self.dataframe.to_csv('{}.csv'.format(csv_filename), index=False)
        else:
            self.dataframe.to_csv('{}_processed.csv'.format(self.filename),
                                  index=False)


def main():
    # Execute module as script.
    import sys
    # Create a Revel instance.
    revel_data = Revel(sys.argv[1])
    # Read REVEL data into a pandas DataFrame.
    try:
        # If number of rows of file to read is provided.
        revel_data.read_into_dataframe(int(sys.argv[2]))
    except IndexError:
        # If number of rows of file to read is not provided.
        revel_data.read_into_dataframe()
    # Export dataframe to a csv file.
    revel_data.export_dataframe_to_csv()


if __name__ == '__main__':
    main()
