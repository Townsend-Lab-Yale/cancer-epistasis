import os
import re
import pandas as pd

from locations import location_data
from locations import gene_coordinates_file
from locations import default_db_file_name

## Some common genes
## Positions in GRCh38 coordinates. Source: http://grch37.ensembl.org/index.html
genes = pd.read_excel(gene_coordinates_file,
                      index_col=0, skiprows=1)


def db_full_name(db_file):
    """Return full path file name of `db_file'.

    :typed db_file: str
    :param db_file: Name of the MAF file containing the data. It
        should be located in the :const:`location_data`. Data files
        names available thus far are:

            - 'cca_ihc.txt' alias 'ihc'
            - 'melanoma.maf.txt' alias 'melanoma'
            - 'rectal_adenocarcinoma.maf.txt' alias 'rectal'
            - 'prostate_adenocarcinoma.maf.txt' alias 'prostate'
            - 'ucec.maf.txt' alias 'ucec'
            - 'tcga.luad.maf.txt' alias 'luad'
            - 'tcga.lusc.maf.txt' alias 'lusc'

         This repository only contains LUAD data, but data for other
         cancer types can be downloaded from TCGA and placed in the
         data directory.

    """
    db_file_name = db_file.lower()
    if db_file_name in "ihc":
        db_file_name = "cca_ihc.txt"
    if db_file_name == "melanoma":
        db_file_name = 'melanoma.maf.txt'
    if db_file_name == "rectal":
        db_file_name = 'rectal_adenocarcinoma.maf.txt'
    if db_file_name == "prostate":
        db_file_name = 'prostate_adenocarcinoma.maf.txt'
    if db_file_name == "ucec":
        db_file_name = 'ucec.maf.txt'
    if db_file_name == 'luad':
        db_file_name = 'tcga.luad.maf.txt'
    if db_file_name == 'lusc':
        db_file_name = 'tcga.lusc.maf.txt'

    if "/" not in db_file_name:
        db_file_name = os.path.join(location_data, db_file_name)

    return db_file_name


def filter_db_by_mutation(db=default_db_file_name,
                          tumor_col_name=None,
                          patient_id_col_name=None,
                          restrict_to_snvs=True,
                          clear_silent=True):
    """Build a data frame with patient IDs and mutations.

    :type db: str or list
    :param db: Name of the MAF file (or files if a db is a list)
        containing the data. It should be located in the
        :const:`location_data`. Data files names available thus far
        are:

            - 'cca_ihc.txt' alias 'ihc'
            - 'melanoma.maf.txt' alias 'melanoma'
            - 'rectal_adenocarcinoma.maf.txt' alias 'rectal'
            - 'prostate_adenocarcinoma.maf.txt' alias 'prostate'
            - 'ucec.maf.txt' alias 'ucec'
            - 'tcga.luad.maf.txt' alias 'luad'
            - 'tcga.lusc.maf.txt' alias 'lusc'

         This repository only contains LUAD data, but data for other
         cancer types can be downloaded from TCGA and placed in the
         data directory.

    :type tumor_col_name: str
    :param tumor_col_name: Name of the column in the data that
        contains the tumor allele. If None, try to infer it from the
        data.

    :type patient_id_col_name: str
    :param patient_id_col_name: Name of the column in the data that
        contains the patient identification. If None, try to infer it
        from the data.

    :type restrict_to_snvs: bool
    :param restrict_to_snvs: If True (default), restrict to only
        single nucleotide variants. It is important because CES only
        can calculate mutation rates for SNVs.

    :type clear_silent: bool
    :param clear_silent: Remove silent mutations from the data. Leave
        as True (default), to later provide ranges for a gene.

    :rtype: pandas.core.frame.DataFrame
    :return: Pandas DataFrame with patient IDs and mutations.

    """

    if isinstance(db, list):
        db_full = [db_full_name(x) for x in db]
        dfs = [pd.read_csv(x, sep="\t", comment="#", low_memory=False)
               for x in db_full]
        dfs = [df.assign(Source=db_name)
               for df, db_name in zip(dfs, db)]
        data = pd.concat(dfs, ignore_index=True)
    else:
        db_full = db_full_name(db)
        data = pd.read_csv(db_full, sep="\t", comment="#", low_memory=False)
        data = data.assign(Source=db)

    if restrict_to_snvs:
        ## TODO: When using the maf file out of CES change to 'snv'
        data = data[data['Variant_Type'] == 'SNP']

    if clear_silent:
        data = data[data['Variant_Classification'] != 'Silent']

    if tumor_col_name is None:
        if 'Tumor_Seq_Allele2' in data.columns:
            tumor_col_name = 'Tumor_Seq_Allele2'
        elif 'Tumor_Allele' in data.columns:
            tumor_col_name = 'Tumor_Allele'
        else:
            raise Exception("Unknown tumor allele. "
                            "Provide variable 'tumor_col_name'")

    if patient_id_col_name is None:
        if 'Tumor_Sample_Barcode' in data.columns:
            patient_id_col_name = 'Tumor_Sample_Barcode'
        elif 'Unique_Patient_Identifier' in data.columns:
            patient_id_col_name = 'Unique_Patient_Identifier'
        else:
            raise Exception("Unknown patient identifier. "
                            "Provide variable 'patient_id_col_name'")

    # removes 'chr' when there:
    data['Chromosome'] = data.apply(
        lambda x: x['Chromosome'].split("chr")[-1], axis=1)
    # includes new column with unique description of mutation
    data['Mutation'] = data.apply(
        lambda x: "{}:{} {}>{}".format(
            x['Chromosome'],
            x['Start_Position'],
            x['Reference_Allele'],
            x[tumor_col_name]),
        axis=1)

    cols_except_id =  (['Chromosome', 'Start_Position', 'Mutation', 'Source']
                       + (['Variant_Classification'] if not clear_silent else []))
    data = data[[patient_id_col_name] + cols_except_id]
    data.columns = ['Patient ID'] + cols_except_id

    return data
