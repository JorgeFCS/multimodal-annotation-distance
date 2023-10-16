"""Tool for determinating distances between multimodal annotations.
"""

import argparse
import configparser
from os import listdir
from os.path import isfile, join, isdir
from typing import List, Literal

import pandas as pd
import pympi
from tqdm import tqdm

Mode = Literal['span', 'point_comparison']

__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Development"


def parse_settings(config: argparse.Namespace) -> dict:
    """
    Parses the configuration file and returns a dictionary with the proper values. Realizes pertinent
    validations on the configuration settings.

    :param config: the namespace containing the configuration file path.
    :type config: argparse.Namespace

    :return: dictionary containing the parsed settings.
    :rtype: dictionary
    """
    # Initializing settings dictionary.
    settings_dict = {}
    # Parsing configuration file.
    settings = configparser.ConfigParser()
    settings.read(config.config_file_path)

    # Loading settings.
    settings_dict['dir_path'] = settings['MAIN_SETTINGS'].get('DIR_PATH')
    # Validating that the directory exists.
    if not isdir(settings_dict['dir_path']):
         raise ValueError('The specified directory does not exist!')
    settings_dict['save_path'] = settings['MAIN_SETTINGS'].get('SAVE_PATH')
    if not isdir(settings_dict['save_path']):
         raise ValueError('The specified directory does not exist!')
    settings_dict['results_file_name'] = settings['MAIN_SETTINGS'].get('RESULTS_FILE_NAME')
    settings_dict['ref_tier'] = settings['MAIN_SETTINGS'].get('REFERENCE_TIER')
    settings_dict['search_value'] = settings['MAIN_SETTINGS'].get('SEARCH_VALUE')
    settings_dict['comparison_tiers'] = settings['MAIN_SETTINGS'].get('COMPARISON_TIERS').split(',')
    settings_dict['mode'] = settings['MAIN_SETTINGS'].get('MODE')
    settings_dict['buffer'] = settings['MAIN_SETTINGS'].getint('BUFFER')

    # Validating that the reference tier is not part of the comparison tiers.
    if settings_dict['ref_tier'] in settings_dict['comparison_tiers']:
        raise ValueError('The reference tier cannot be part of the comparison tiers!')

    return settings_dict


def get_file_paths(dir_path: str) -> List[str]:
    """
    Gets the paths for all EAF files in a given directory.

    :param dir_path: Path of the directory in which to search.
    :type dir_path: str

    :return: List of paths for the found EAF files.
    :rtype: List[string]
    """
    file_list = [join(dir_path, file) for file in listdir(dir_path) if isfile(join(dir_path, file)) and '.eaf' in file]
    
    return file_list


def construct_columns(ref_tier: str,
                      search_value: str,
                      mode: Mode) -> List[str]:
    """
    Constructs the column names for the results dataframe.

    :param ref_tier: Reference tier against which the `comparison_tiers`
        will be compared.
    :type ref_tier: string
    :param search_value: Particular value of `ref_tier` for which to search
        for overlaps.
    :type search_value: string
    :param mode: Either `span` or `point_comparison` for both possible operation
        modalities.
    :type mode: string

    :return: List containing the column names for the results dataframe.
    :rtype: List[string]
    """
    # Base column names.
    columns = [ref_tier + "_" + search_value + "_start_ms",
               ref_tier + "_" + search_value + "_end_ms",
               ref_tier + "_" + search_value + "_duration",
               "Buffer_ms"]
    if mode == "span":
        # Column names for the different `comparison tiers`.
        columns.extend(["Tier", "Value", "Begin_ms", "End_ms",
                        "Duration", "Overlap_time", "Overlap_ratio",
                        "Diff_start", "Diff_end"])
    else:
        # Column names for the different `comparison tiers`.
        columns.extend(["Tier", "Value", "Begin_ms", "End_ms",
                        "Duration"])

    return columns


def get_sym_tier_parent(eaf_file: pympi.Elan.Eaf,
                        ref_tier: str) -> str:
    """
    Function that obtains the ID of the first parent of a symbolic tier that is
    not a symbolic tier itsef.

    :param eaf_file: the Eaf object that contains the file information.
    :type eaf_file: pympi.Elan.Eaf
    :param ref_tier: reference tear to check.
    :type ref_tier: string

    :return: A string representing the name of the found parent tier.
    :rtype: string
    """
    if eaf_file.tiers[ref_tier][1]:
        parent_tier = eaf_file.tiers[ref_tier][2]['PARENT_REF']
        return get_sym_tier_parent(eaf_file, parent_tier)
    else:
        return ref_tier


def get_annotation_data_for_sym_tier(eaf_file: pympi.Elan.Eaf,
                                     ref_tier: str,
                                     search_value: str) -> List:
    """
    Function that obtains the annotation data for a symbolic tier that matches a given
    search value.

    :param eaf_file: the Eaf object that contains the file information.
    :type eaf_file: pympi.Elan.Eaf
    :param ref_tier: reference tear to check.
    :type ref_tier: string
    :param search_value: the value within the tier that we want to find.
    :type search_value: string

    :return: The list containing the found matches.
    :rtype: List
    """
    # Obtain the name of the first parent that is not a symbolic one.
    parent_tier_name = get_sym_tier_parent(eaf_file,
                                           ref_tier)
    # Getting the time intervals for the parent's data.
    ref_search_intervals = []
    # We only add it if the reference tier contains our specific search value.
    for annotation in eaf_file.get_annotation_data_for_tier(parent_tier_name):
        # Taking the start and end times of the parent tier to search for the
        # symbolic reference tier.
        for match in eaf_file.get_annotation_data_between_times(ref_tier,
                                                                annotation[0],
                                                                annotation[1]):
            # We only add it if the reference tier contains our specific search value.
            if search_value in match[2]:
                ref_search_intervals.append(match)

    return ref_search_intervals


def get_point_comparison_matches(file_paths: List[str],
                                 ref_tier: str,
                                 search_value: str,
                                 comparison_tiers: List[str], 
                                 buffer: int = 0) -> pd.DataFrame:
    """
    Function that searches for matches between a given `ref_value` in 
    a reference tier `ref_tier` and one or more `comparison_tiers` within a time frame given
    by the ending time of the found `ref_values`.
    Returns the results in tabular format.

    :param file_paths: List of EAF files for analysis.
    :type file_paths: List[string]
    :param ref_tier: Reference tier against which the `comparison_tiers`
        will be compared.
    :type ref_tier: string
    :param search_value: Particular value of `ref_tier` for which to search
        for overlaps.
    :type search_value: string
    :param comparison_tiers: List of tiers to compare against `ref_tier`.
    :type comparison_tiers: List[string]
    :param buffer: Time buffer with which to extend the time frame start and end times.
        Deafults to 0.
    :type buffer: integer, optional

    :return: Dataframe containing the computed results in tabular format.
    :rtype: pandas Dataframe
    """
    # Constructing columns
    columns = construct_columns(ref_tier, search_value, 'point_comparison')
    # Initializing results dataframe.
    results_df = pd.DataFrame([], columns=columns)

    for file_path in file_paths:
        print("Analyzing file ", file_path, "...")
        # Reading EALN file.
        eaf = pympi.Elan.Eaf(file_path)

        # Get the elements in ref_tier that match the search_value.
        if eaf.tiers[ref_tier][1]:
            # In case of a symbolic reference tier.
            ref_search_intervals = get_annotation_data_for_sym_tier(eaf, ref_tier, search_value)
        else:
            ref_search_intervals = []
            for annotation in eaf.get_annotation_data_for_tier(ref_tier):
                # We only add it if the reference tier contains our specific search value.
                if search_value in annotation[2]:
                    ref_search_intervals.append(annotation)

        for annotation in tqdm(ref_search_intervals):
            # Constructing first part of the new entry with start time, end time and 
            # duration.
            row_base = [annotation[0], annotation[1], (annotation[1] + buffer) - (annotation[0] - buffer), buffer]
            for tier in comparison_tiers:
                # Finding matches taking as reference only the ending time of the found
                # reference annotations.
                matches = eaf.get_annotation_data_between_times(tier,
                                                                annotation[1] - buffer,
                                                                annotation[1] + buffer)
                for match in matches:
                    # Saving results.
                    overlap_time = min(annotation[1] + buffer, match[1]) - max(annotation[0] - buffer, match[0])
                    if overlap_time > 0:
                        duration = match[1] - match[0]
                        row_part = [tier, match[2], match[0], match[1], duration]
                        row = row_base + row_part
                        results_df.loc[len(results_df.index)] = row
        print("... done.")

    return results_df 


def get_span_overlaps(file_paths: List[str],
                      ref_tier: str,
                      search_value: str,
                      comparison_tiers: List[str],
                      buffer: int = 0) -> pd.DataFrame:
    """
    Function that searches for overlaps between a given `ref_value` in 
    a reference tier `ref_tier` and one or more `comparison_tiers`. Computes additional
    data such as overlap ratio, overlap time, etc., and returns the results in tabular format.

    :param file_paths: List of EAF files for analysis.
    :type file_paths: List[string]
    :param ref_tier: Reference tier against which the `comparison_tiers`
        will be compared.
    :type ref_tier: string
    :param search_value: Particular value of `ref_tier` for which to search
        for overlaps.
    :type search_value: string
    :param comparison_tiers: List of tiers to compare against `ref_tier`.
    :type comparison_tiers: List[string]
    :param buffer: Time buffer with which to extend the time frame start and end times.
        Deafults to 0.
    :type buffer: integer, optional

    :return: Dataframe containing the computed results in tabular format.
    :rtype: pandas Dataframe
    """
    # Constructing columns
    columns = construct_columns(ref_tier, search_value, 'span')
    # Initializing results dataframe.
    results_df = pd.DataFrame([], columns=columns)

    for file_path in file_paths:
        print("Analyzing file ", file_path, "...")
        # Reading EALN file.
        eaf = pympi.Elan.Eaf(file_path)
        # Get the elements in ref_tier that match the search_value.
        if eaf.tiers[ref_tier][1]:
            # In case of a symbolic reference tier.
            ref_search_intervals = get_annotation_data_for_sym_tier(eaf, ref_tier, search_value)
        else:
            ref_search_intervals = []
            for annotation in eaf.get_annotation_data_for_tier(ref_tier):
                # We only add it if the reference tier contains our specific search value.
                if search_value in annotation[2]:
                    ref_search_intervals.append(annotation)

        for annotation in tqdm(ref_search_intervals):
            # Constructing first part of the new entry with start time, end time and 
            # duration.
            row_base = [annotation[0], annotation[1], (annotation[1] + buffer) - (annotation[0] - buffer), buffer]
            for tier in comparison_tiers:
                potential_matches = eaf.get_annotation_data_between_times(tier,
                                                                          annotation[0] - buffer,
                                                                          annotation[1] + buffer)
                # Validate that it is not the immediate previous - next entry.
                for match in potential_matches:
                     if not (annotation[0] - buffer) == match[1] and not (annotation[1] + buffer) == match[0]:
                        # Computing overlap data.
                        overlap_time = min(annotation[1] + buffer, match[1]) - max(annotation[0] - buffer, match[0])
                        if overlap_time > 0:
                            duration = match[1] - match[0]
                            overlap_ratio = round(overlap_time / ((annotation[1] + buffer) - (annotation[0] - buffer)), 3)
                            start_diff = (annotation[0] - buffer) - match[0]
                            end_diff = (annotation[1] + buffer) - match[1]
                            row_part = [tier, match[2], match[0], match[1],
                                        duration, overlap_time, overlap_ratio,
                                        start_diff, end_diff]
                            row = row_base + row_part
                            results_df.loc[len(results_df.index)] = row
        print("... done.")

    return results_df


def main(config: argparse.Namespace) -> None:
    # Parsing settings.
    settings_dict = parse_settings(config)

    # Getting all EAF file paths in the specified directory.
    file_paths = get_file_paths(settings_dict['dir_path'])

    #file_paths = ["./ELAN-data/bgest_001.eaf", "./ELAN-data/bgest_002.eaf"]
    # Getting the results according to the specified modality.
    if settings_dict['mode'] == "span":
        results_df = get_span_overlaps(file_paths,
                                       settings_dict['ref_tier'],
                                       settings_dict['search_value'],
                                       settings_dict['comparison_tiers'],
                                       settings_dict['buffer'])
    elif settings_dict['mode'] == "point_comparison":
        results_df = get_point_comparison_matches(file_paths,
                                                  settings_dict['ref_tier'],
                                                  settings_dict['search_value'],
                                                  settings_dict['comparison_tiers'],
                                                  settings_dict['buffer'])
    else:
        raise ValueError('Invalid mode.')
    # Saving CSV file.
    results_df.to_csv(settings_dict['save_path'] + settings_dict['results_file_name'],
                      index=False,
                      encoding="utf-8-sig")


if __name__ == '__main__':
    # Initializing the argument parser.
    parser = argparse.ArgumentParser(description='Multimodal annotation distance tool.')
    parser.add_argument('--config-file-path',
                        help='Path to the configuration .ini file.',
                        action='store', required=True, type=str)
    # Parsing the arguments.
    config = parser.parse_args()
    # Calling the main function.
    main(config)
