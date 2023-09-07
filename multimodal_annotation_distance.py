"""Tool for determinating distances between multimodal annotations.

Based on the R code originally proposed by Camila Barros and 

"""

import pandas as pd
import pympi

__author__ = "Jorge Ciprian"
__credits__ = ["Jorge Ciprian", "Camila Barros"]
__license__ = ""  # TBD
__version__ = "0.1.0"
__status__ = "Development"

def main():
    get_overlap()

def get_overlap():
    file_path = "./ELAN-data/bgest_001.eaf"
    search_value = "stroke"
    ref_tier = "GE-Phase"
    comparison_tiers = ["GE-Phrase", "NTB"]

    # Reading EALN file.
    eaf = pympi.Elan.Eaf(file_path)
    results_df = pd.DataFrame([], columns=['Begin_Time_ms',
                                           'End_Time_ms',
                                           'GE-Phrase',
                                           'GE-Phase'])
    # Get the elements in ref_tier that match the search_value.
    ref_search_intervals = []
    for annotation in eaf.get_annotation_data_for_tier(ref_tier):
        if search_value in annotation:
            ref_search_intervals.append(annotation)
    
    for annotation in ref_search_intervals:
        print("Reference interval times for ", annotation[2] ,": ", annotation[0], " - ", annotation[1])
        #print("Found matches:")
        for tier in comparison_tiers:
            potential_matches = eaf.get_annotation_data_between_times(tier,
                                                                    annotation[0],
                                                                    annotation[1])
            # Validate that it is not the immediate previous - next entry.
            for potential_match in potential_matches:
                if not annotation[0] == potential_match[1] and not annotation[1] == potential_match[0]:
                    print("Found matches:")
                    print(potential_match)
            


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