"""Tool for determinating distances between multimodal annotations.

Based on the R code originally proposed by Camila Barros and 

"""

from os import listdir
from os.path import isfile, join, isdir
from typing import List

import pandas as pd
import pympi
from tqdm import tqdm

__author__ = "Jorge Ciprian"
__credits__ = ["Jorge Ciprian", "Camila Barros"]
__license__ = ""  # TBD
__version__ = "0.1.0"
__status__ = "Development"


def main():
    search_value = "stroke"
    ref_tier = "GE-Phase"
    comparison_tiers = ["GE-Phrase", "NTB"]
    dir_path = "./ELAN-data"
    file_paths = get_file_paths(dir_path)
    #file_paths = ["./ELAN-data/bgest_001.eaf", "./ELAN-data/bgest_002.eaf"]
    get_overlaps(file_paths, ref_tier, search_value, comparison_tiers)


def get_file_paths(dir_path: str) -> List[str]:
    """
    Gets the paths for all EAF files in a given directory.

    :param dir_path: Path of the directory in which to search.
    :type dir_path: str

    :return: List of paths for the found EAF files.
    :rtype: List[string]
    """
    # Validating that the directory exists.
    assert isdir(dir_path), "The specified directory does not exist!"

    file_list = [join(dir_path, file) for file in listdir(dir_path) if isfile(join(dir_path, file)) and '.eaf' in file]
    
    return file_list



def construct_columns(ref_tier: str,
                      search_value: str) -> List[str]:
    """
    Constructs the column names for the results dataframe.

    :param ref_tier: Reference tier against which the `comparison_tiers`
        will be compared.
    :type ref_tier: string
    :param search_value: Particular value of `ref_tier` for which to search
        for overlaps.
    :type search_value: string

    :return: List containing the column names for the results dataframe.
    :rtype: List[string]
    """
    # Base column names.
    columns = [ref_tier + "_" + search_value + "_start_ms",
               ref_tier + "_" + search_value + "_end_ms",
               ref_tier + "_" + search_value + "_duration"]
    # Column names for the different `comparison tiers`.
    columns.extend(["Tier", "Value", "Begin_ms", "End_ms",
                    "Duration", "Overlap_time", "Overlap_ratio",
                    "Diff_start", "Diff_end"])
    return columns


def get_overlaps(file_paths: List[str],
                 ref_tier: str,
                 search_value: str,
                 comparison_tiers: List[str]) -> pd.DataFrame:
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

    :return: Dataframe containing the computed results in tabular format.
    :rtype: pandas Dataframe
    """
    # Constructing columns
    columns = construct_columns(ref_tier, search_value)
    # Initializing results dataframe.
    results_df = pd.DataFrame([], columns=columns)

    for file_path in file_paths:
        print("Analyzing file ", file_path, "...")
        # Reading EALN file.
        eaf = pympi.Elan.Eaf(file_path)
        # Get the elements in ref_tier that match the search_value.
        ref_search_intervals = []
        for annotation in eaf.get_annotation_data_for_tier(ref_tier):
            if search_value in annotation:
                ref_search_intervals.append(annotation)

        for annotation in tqdm(ref_search_intervals):
            row_base = [annotation[0], annotation[1], annotation[1] - annotation[0]]
            for tier in comparison_tiers:
                potential_matches = eaf.get_annotation_data_between_times(tier,
                                                                        annotation[0],
                                                                        annotation[1])
                # Validate that it is not the immediate previous - next entry.
                for match in potential_matches:
                    if not annotation[0] == match[1] and not annotation[1] == match[0]:
                        # Computing overlap data.
                        duration = match[1] - match[0]
                        overlap_time = min(annotation[1], match[1]) - max(annotation[0], match[0])
                        overlap_ratio = round(overlap_time / (annotation[1] - annotation[0]), 3)
                        row_part = [tier, match[2], match[0], match[1],
                                    duration, overlap_time, overlap_ratio,
                                    annotation[0] - match[0],
                                    annotation[1] - match[1]]
                        row = row_base + row_part
                        results_df.loc[len(results_df.index)] = row
        print("... done.")
    results_df.to_csv("./results.csv", index=False, encoding="utf-8-sig") 


def test_experiments():
    file_path = "./ELAN-data/bgest_001.eaf"
    # Initializing the ELAN file.
    eaf = pympi.Elan.Eaf(file_path)
    print(dir(eaf))
    print("*"*60)
    #print(eaf.__dict__)
    #print("*"*60)
    #print("TIER keys:")
    #print(eaf.tiers.keys())
    #print("*"*60)
    print("Tier names:")
    print(eaf.get_tier_names())
    print("*"*60)
    print("GE-Phrase tier keys")
    #print(eaf.tiers['GE-Phrase'])
    for i in range(4):
        print(eaf.tiers['GE-Phrase'][i])
        print("-"*60)
    print("Length of first dict: ", len(eaf.tiers['GE-Phrase'][0]))
    print("*"*60)
    print("GE-Units tier keys")
    #print(eaf.tiers['GE-Units'])
    for i in range(4):
        print(eaf.tiers['GE-Units'][i])
        print("-"*60)
    print("Length of first dict: ", len(eaf.tiers['GE-Units'][0]))
    print("*"*60)
    print("Timeslots:")
    print(eaf.timeslots)
    print("*"*60)
    print("REFERENCE TEXT FILES: ")
    # To explore stuff: load the TXT files as dataframes to check them out.
    #print("Transcription file")
    #pd_transcription = pd.read_csv('./text-data-samples/transc/bgest_001_transc.txt', sep="\t")
    #print(pd_transcription)
    #print(pd_transcription.columns)
    #print("*"*60)
    print("Gesture file")  # The last column in the gesture file is filled with NaN values; can be safely removed.
    pd_gesture = pd.read_csv('./text-data-samples/gesture/bgest_001_gesture.txt', sep='\t')
    print(pd_gesture)
    #print(pd_gesture.columns)
    #selected_rows = pd_gesture[pd_gesture['Unnamed: 13'].isnull()]
    #print(selected_rows)

if __name__ == '__main__':
    main()