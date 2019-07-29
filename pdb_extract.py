"""
Simple command-line tool for fetching data from the RCSB Protein Data Bank.

Uses the Protein Data Bank's REST API to retrieve data on PDB structures
obtained with one of the experimental methods in EXP_METHODS.
Each of REST API keys below corresponds to a column header in the output CSV.
This is intended to be run automatically with cron, but can also be run
manually (see help menu for usage).
"""

import argparse
from collections import OrderedDict
import csv
import os, sys
import pandas as pd
import logging
import urllib.request, urllib.error, urllib.parse
import time
import datetime

DEFAULT_OUTPUT = 'PDB_extract-out.csv'
BASR_URL = 'http://www.rcsb.org/pdb/rest/'
SEARCH_URL = BASR_URL + '/search'
FETCH_URL_TEMPLATE = BASR_URL + (
    'customReport.csv?pdbids=%s&customReportColumns=%s'
    '&service=wsfile&format=csv')

logger = logging.getLogger('pdb_extract')
logger.setLevel('INFO')

# Keys corresponding to REST API fields we want data from
PDB_ID_KEY = 'structureId'
STRUCT_TITLE_KEY = 'structureTitle'
RES_KEY = 'resolution'
LIG_NAME_KEY = 'ligandName'
CLASS_KEY = 'classification'
MACRO_TYPE_KEY = 'macromoleculeType'
EMDB_KEY = 'emdbId'
PUBMED_KEY = 'pubmedId'
REL_DATE_KEY = 'releaseDate'
EXP_TECHNIQUE_KEY = 'experimentalTechnique'
UNIPROT_ID_KEY = 'uniprotAcc'
SOURCE_KEY = 'source'
LIG_SMILES_KEY = 'ligandSmiles'

# Dict for mapping REST API fields to output CSV headers
FIELDS_DICT = OrderedDict(
    [(PDB_ID_KEY, 'PDB ID'), (STRUCT_TITLE_KEY,
                              'Structure Title'), (RES_KEY, 'Resolution'),
     (LIG_NAME_KEY, 'Ligand Name(s)'), (CLASS_KEY, 'Classification'),
     (MACRO_TYPE_KEY, 'Macromol. Type'), (EMDB_KEY, 'EMDB ID'), (PUBMED_KEY,
                                                                 'Pubmed ID'),
     (REL_DATE_KEY, 'Rel. Date'), (EXP_TECHNIQUE_KEY, 'Exp. Technique'),
     (UNIPROT_ID_KEY, 'Uniprot ID'), (SOURCE_KEY,
                                      'Source'), (LIG_SMILES_KEY,
                                                  'Ligand SMILES')])

FIELDS_LIST = ','.join(FIELDS_DICT)

# Valid experimental methods to search for PDB IDs
EXP_METHODS = ('X-RAY', 'SOLUTION NMR', 'SOLID-STAE NMR',
               'ELECTRON MICROSCOPY', 'ELECTRON CRYSTALLOGRAPHY',
               'FIBER DIFFRACTION', 'NEUTRON DIFFRACTION',
               'SOLUTION SCATTERING', 'OTHER', 'HYBRID')

DEFAULT_METHOD = EXP_METHODS[3]

# Template for XML query to get PDB IDs list by experimental method
XML_TEMPLATE = (
    '<?xml version="1.0" encoding="UTF-8"?>'
    '<orgPdbQuery>'
    '<version>B0907</version>'
    '<queryType>org.pdb.query.simple.ExpTypeQuery</queryType>'
    '<description>Experimental Method Search: Experimental Method=%s</description>'
    '<mvStructure.expMethod.value>%s</mvStructure.expMethod.value>'
    '</orgPdbQuery>')


def parse_args():
    """
    Parse the command line options.

    @return:  All script arguments and options.
    @rtype:  class:`argparse.Namespace
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.HelpFormatter)

    parser.add_argument(
        '-outfile',
        help='The name of the output CSV file',
        default=DEFAULT_OUTPUT)

    parser.add_argument(
        '-method',
        help='The experimental technique used to get the crystal structure',
        choices=EXP_METHODS,
        default=DEFAULT_METHOD)

    parser.add_argument(
        '-min_res',
        help='The minimum resolution for the crystal structure (Angstroms)',
        type=float)

    parser.add_argument(
        '-condensed',
        help='Write a condensed output CSV file in addition to the default '
        'full format',
        action='store_true')

    args = parser.parse_args()
    outfile_ext = os.path.splitext(args.outfile)[1]
    if not outfile_ext == '.csv':
        parser.error("Output file must have a *.csv extension (CSV format).")

    return args


def query_rcsb_db(url, data=None):
    """
    Query RCSB database with the URL and return the response.
    """

    request = urllib.request.Request(url, data=data)

    try:
        response = urllib.request.urlopen(request)
    except urllib.error.URLError as url_err:
        logger.error(url_err)
        sys.exit(1)

    return response


def get_pdbs_from_method(exp_method):
    """
    Get the list of PDB IDs for the given experimental method.

    @return PDB ID list
    @rtype: list
    """

    xml_str = XML_TEMPLATE % (exp_method, exp_method)
    # logger.info("Getting PDB IDs for experimental method %s..." % exp_method)
    pdb_ids_response = query_rcsb_db(SEARCH_URL, data=xml_str.encode('UTF-8'))
    pdb_ids_str = pdb_ids_response.read()
    # pdbs_list = pdb_ids_str.splitlines()
    pdbs_list = [pdb.decode('UTF-8') for pdb in pdb_ids_str.splitlines()]
    logger.info("Found %i PDB IDs for experimental method %s" %
                (len(pdbs_list), exp_method))

    return pdbs_list


def get_dataframe_from_pdbs(pdb_id_list, field_names):
    """
    Create a pandas dataframe from with data returned by querying the RCSB with
    the PDB IDs list and fieldnames list.

    @return: The pandas dataframe
    @rtype: <pandas.core.frame.DataFrame>
    """

    # Build data_str with PDB data
    data_str = ''
    logger.info("Fetching PDB data...")
    for pdb_num in pdb_id_list:
        fetch_url = FETCH_URL_TEMPLATE % (pdb_num, field_names)
        pdb_data = query_rcsb_db(fetch_url).read()
        pdb_data = pdb_data.decode('UTF-8')
        data_str += pdb_data

    # Build rows of PDB data for creating the DataFrame
    data_rows = list(csv.reader(data_str.splitlines()))

    # Create the DataFrame object
    df = pd.DataFrame(data_rows)

    if df.empty:
        raise Exception("The DataFrame object is empty.")

    # Do some cleanup
    df = df.drop_duplicates()
    # Set column names to the REST API keys
    df.columns = df.iloc[0]
    # Slice off the first and last row (the first it's a duplicate of the column
    # names and the last is blank)
    df = df[1:-1]
    df.rename(columns=FIELDS_DICT, inplace=True)

    return df


def get_condensed_df(col_names, d_frame):
    """
    Get a 'condensed' dataframe of the PDB results.  Concantenates differences
    in column data within each group of identical PDBs, then removes duplicate
    rows from the dataframe.

    @return: The condensed DataFrame
    @rtype: <pandas.core.frame.DataFrame>
    """

    concat_df = concat_column_data(col_names, d_frame)
    concat_df.drop('chainId', axis=1, inplace=True)
    concat_df.reset_index(inplace=True, level=0, drop=True)
    condensed_df = concat_df.drop_duplicates()

    return condensed_df


def concat_column_data(col_names, d_frame):
    """
    Concatenate unique column data in the dataframe for each identical PDB ID.

    @param col_names: The names of column headers containing data to
    concatenate.

    @type col_names: List

    @return:

    PDB ID | Column Name |             PDB ID | Column Name |
    ----------------------      -->    --------------------------
    1F6H   | Value_A                   1F6H   | Value_A Value_B
    1F6H   | Value_B                   1F6H   | Value_A Value_B

    @rtype: <pandas.core.frame.DataFrame>
    """

    # Group by PDB first
    grouped = d_frame.groupby(FIELDS_DICT[PDB_ID_KEY])
    results = []

    # FIXME: With the latest version of pandas, we could avoid having to
    # iterate through each group.  We are currently using 0.14.0.
    for name, grp in grouped:
        for col_nm in col_names:
            grp[col_nm] = grp[col_nm].astype(str)
            unique_entries = list(grp[col_nm].unique())
            join_char = ' '
            if col_nm == FIELDS_DICT[LIG_NAME_KEY]:
                join_char = ' | '
            # Don't include blank entries in the replaced column
            grp[col_nm] = join_char.join([_f for _f in unique_entries if _f])
            results.append(grp)

    concat_df = pd.concat(results)

    return concat_df


def filter_by_res(d_frame, res_val, res_col_header):
    """
    Filter the dataframe by resolution.

    @param res_val: The minimum resolution required to retain a row of data
    @type: res_val: float

    @param res_col_header: The name of the column header containing the resolution
    @type res_col_header: str

    @return: The filtered DataFrame object
    @rtype: <pandas.core.frame.DataFrame>
    """

    d_frame[res_col_header] = pd.to_numeric(d_frame[res_col_header])
    logger.info(
        "Filtering out PDB structures with resolution > %0.1f Angstroms" %
        res_val)
    res_mask = d_frame[res_col_header] <= res_val
    d_frame = d_frame[res_mask]

    return d_frame

def main():
    cmd_args = parse_args()
    # Get the list of PDB IDs
    pdb_ids = get_pdbs_from_method(cmd_args.method)
    fieldnames = FIELDS_LIST
    # Get the pandas dataframe
    df = get_dataframe_from_pdbs(pdb_ids, fieldnames)

    if cmd_args.min_res:
        df = filter_by_res(df, cmd_args.min_res, FIELDS_DICT[RES_KEY])

    # Write full output to CSV
    tmp_outfile = os.path.splitext(cmd_args.outfile)[0] + '_tmp.csv'
    df.to_csv(tmp_outfile, index=False)

    # Check that output was written and contains data before overwriting the
    # previous CSV file
    if not os.path.exists(tmp_outfile):
        logger.error("Output CSV %s file was not written." % cmd_args.outfile)
        sys.exit(1)

    filesize = os.path.getsize(tmp_outfile)
    if filesize == 0:
        logger.error("Output CSV %s is empty." % cmd_args.outfile)
        sys.exit(1)

    most_recent_date = pd.to_datetime(df[FIELDS_DICT[REL_DATE_KEY]]).max()
    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime(
        '%Y-%m-%d %H:%M:%S')

    os.rename(tmp_outfile, cmd_args.outfile)
    logger.info("Output CSV written to: %s" % cmd_args.outfile)
    logger.info("Date and time run: %s\n"
                "Date of most recently entered PDB data: %s\n" %
                (timestamp, str(most_recent_date)))

    # Write condensed CSV output file
    if cmd_args.condensed:
        pubmed_col = FIELDS_DICT[PUBMED_KEY]
        lig_col = FIELDS_DICT[LIG_NAME_KEY]
        source_col = FIELDS_DICT[SOURCE_KEY]
        mac_type_col = FIELDS_DICT[MACRO_TYPE_KEY]
        uniprot_col = FIELDS_DICT[UNIPROT_ID_KEY]

        cols_to_concat = [
            pubmed_col, lig_col, source_col, mac_type_col, uniprot_col
        ]
        condensed_df = get_condensed_df(cols_to_concat, df)
        out_name = 'CONDENSED_' + cmd_args.outfile
        logger.info("Writing condensed CSV output to %s" % out_name)

        # Move QED columns next to ligand columnds for final output.
        export_cols = condensed_df.columns.tolist()[:-2]
        qed_items = condensed_df.columns.tolist()[-2:]
        export_cols[5:5] = qed_items
        condensed_df.to_csv(out_name, index=False, columns=export_cols)


if __name__ == '__main__':
    main()
