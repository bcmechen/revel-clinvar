import os
import pandas as pd


class DosageSensitivityMap:
    """To model ClinGen Dosage Sensitivity data.

    ClinGen_gene_curation_list.tsv is a tab separated file for gene curation.
    The tsv files have a header and contain all of the curation information
    found on the ClinGen Dosage Sensitivity Map web pages.

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

    @property
    def header_in_source(self):
        """Get column names.

        ClinGen_gene_curation_list.tsv contains a few information lines
        before the header line. This function finds and returns the header
        line in a list.

        Returns:
            (list): A list of column names.
        """
        with open(self.filepath) as f:
            for line in f:
                # The length of information line after split is 1.
                # The first line with length greater than 1 is the header.
                if len(line.split('\t')) > 1:
                    return [x.strip() for x in line.split('\t')]

    def read_into_dataframe(self, num_of_rows=None):
        """Read Dosage Sensitivity data into a pandas DataFrame.

        Args:
            num_of_rows (int): Number of rows of file to read.
        """
        self.dataframe = pd.read_csv(self.filepath, nrows=num_of_rows,
                                     sep='\t', header=None,
                                     names=self.header_in_source, comment='#')

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
    # Create a DosageSensitivityMap instance.
    dsm_data = DosageSensitivityMap(sys.argv[1])
    # Read Dosage Sensitivity data into a pandas DataFrame.
    try:
        # If number of rows of file to read is provided.
        dsm_data.read_into_dataframe(int(sys.argv[2]))
    except IndexError:
        # If number of rows of file to read is not provided.
        dsm_data.read_into_dataframe()
    # Export dataframe to a csv file.
    dsm_data.export_dataframe_to_csv()


if __name__ == '__main__':
    main()
