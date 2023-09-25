import pympi



"""

parentTier = self.tiers[self.annotations[ref_id]]
while 'PARENT_REF' in parentTier[2] and parentTier[2]['PARENT_REF'] and len(parentTier[2]) > 0:
    ref_id = parentTier[1][ref_id][0]
    parentTier = self.tiers[self.annotations[ref_id]]

return parentTier[0][ref_id]
"""


def get_sym_tier_parent(eaf_file: pympi.Elan.Eaf,
                        ref_tier: str):
    if eaf_file.tiers[ref_tier][1]:
        parent_tier = eaf_file.tiers[ref_tier][2]['PARENT_REF']
        return get_sym_tier_parent(eaf_file, parent_tier)
    else:
        return ref_tier

def get_parent_aligned_annotation(eaf, ref_id):
    parentTier = eaf.tiers[eaf.annotations[ref_id]]
    while 'PARENT_REF' in parentTier[2] and parentTier[2]['PARENT_REF'] and len(parentTier[2]) > 0:
        print("Parent tier [1]: ", parentTier[1])
        print("Parent tier [2]: ", parentTier[2])
        ref_id = parentTier[1][ref_id][0]
        parentTier = eaf.tiers[eaf.annotations[ref_id]]
    
    return parentTier[0][ref_id]

def get_annotation_data_for_sym_tier(eaf_file: pympi.Elan.Eaf,
                                     ref_tier: str,
                                     search_value: str):
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
            if search_value in match:
                ref_search_intervals.append(match)

    return ref_search_intervals

eaf = pympi.Elan.Eaf('./ELAN-data/bgest_001.eaf')

#tier = 'NTB'
tier = 'InfoStructure'
#tier = 'left - GE-Movement'
#tier = 'GE-Phase'
search_value = 'PAR'

parent_tier = get_sym_tier_parent(eaf, tier)
print("PARENT TIER: ", parent_tier)
annotation_data = get_annotation_data_for_sym_tier(eaf,
                                                   tier,
                                                   search_value)
print("ANNOTATION DATA: ", annotation_data)

#annotations = eaf.get_annotation_data_for_tier(tier)
#print(annotations)
#print(type(annotations[0]))

#print(eaf.tiers[tier][1])
#print("*"*70)
#print(eaf.tiers[tier][0])

#print(type(eaf.tiers))

#print("-"*70)
#print(eaf.tiers[tier].items())
#parentTier = eaf.tiers[eaf.annotations['a975']]


#annotations = eaf.get_annotation_data_between_times(tier,
#                                                    9938,
#                                                    11306)

#print(annotations)

#for aid, (ref, value, _, _) in eaf.tiers[tier][1].items():
#    print("REF: ", ref)
#    print("VALUE: ", value)
#

#print("MAIN REFERENCE TIER: ")
#print("TIER [0]: ", eaf.tiers[tier][0])
#print("TIER [1]: ", eaf.tiers[tier][1])
#print("TIER [2]: ", eaf.tiers[tier][2])
"""
if eaf.tiers[tier][1]:
    bucket = []
    for aid, (ref, value, prev, _) in eaf.tiers[tier][1].items():
        parentTier = get_parent_aligned_annotation(eaf, ref)
        print("Parent tier: ", parentTier)
        #print("ITEM: ", ref)
        #print(eaf.tiers[eaf.annotations[ref]])
        #break


#print(eaf.annotations['a975'])
#print(eaf.tiers[eaf.annotations['a975']])
print("#"*70)



"""