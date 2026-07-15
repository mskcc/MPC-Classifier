import pandas as pd
import logging
import sys
import yaml
from datetime import datetime
import os 
import numpy as np

#Add logging configuration
logging.basicConfig(
format='%(asctime)s - %(levelname)s - %(message)s',
level=logging.INFO
)

# ─── Load yaml configuration ─────────────────────────────────────────────────────────
with open('mpc_config.yaml', 'r') as f:
    config = yaml.safe_load(f)


SYSTEMIC_CANCERS = config.get('systemic_cancers', []) # e.g. ['Leukemia', 'Lymphoma']
TOPOLOGY_GROUPS = config.get('topology_groups', {}) # e.g. {'GI': ['C181']}
RAW_MORPHOLOGY_GROUPS = config.get('morphology_groups', {}) # e.g. {'Adeno': ['8140/3', '8141/3']}
PAIRED_ORGANS = config.get('paired_organs', {}) # e.g. {'Lung': ['C341', 'C342']}
EXCEPTIONS = config.get('exceptions', [])# e.g. ['C18', 'C44']
COLLAPSED_CANCER_TYPES = config.get("collapsed_cancer_types", {})
SYNDROMIC_PATTERNS = config.get("syndromic_patterns", {})
DATE_DIFFERENCES = config.get("date_differences", {})
ADJACENCY_SCORES = config.get("adjacency_scores", {})

# ─── Helpers ────────────────────────────────────────────────────────────────────
def expand_code_range(code_range):
    """Expand a single code or range and append M"""
    if isinstance(code_range, str) and '-' in code_range:
        start, end = map(int, code_range.split('-'))
        return [f"M{i}" for i in range(start, end + 1)]
    elif isinstance(code_range, (int, str)):
        return [f"M{int(code_range)}"]
    else:
        raise ValueError(f"Unsupported format: {code_range}")


def expand_morphology_groups(raw_groups):
    """Expand morphology code ranges and individual codes into
       their corresponding list of formatted morphology codes."""
    try:
        expanded = {}
        for group, code_list in raw_groups.items():
            expanded[group] = []
            for item in code_list:
                expanded[group].extend(expand_code_range(item))
        return expanded
    except Exception:
        logger.exception(f"Failed to expand morphology group '%s'.", group)
        raise

# Expand raw morphology group config into full list of formatted codes
MORPHOLOGY_GROUPS = expand_morphology_groups(RAW_MORPHOLOGY_GROUPS)


def get_morpho_group(code):
    """Get morphology group for a given code"""
    code = code.split("/")[0]
    for grp, codes in MORPHOLOGY_GROUPS.items():
        if code in codes:
            return grp
    if code != 'UNKNOWN':
        logging.warning(f"Morphology code: {code}: not found")
        return None


def get_topo_group(code):
    """Get topology group for a given code"""
    # Remove dot and merge parts, then uppercase
    merged_code = code.replace('.', '').upper()
    for grp, codes in TOPOLOGY_GROUPS.items():
        if merged_code in codes:
            return grp
    
    logging.warning(f"Topology code: {code}: not found")
    return None


def get_collapsed_cancer_types(codes, morphologies_desc):
    """Get cancer for a given code"""
    # Remove dot and merge parts, then uppercase
    cancer_types = []
    for i, code in enumerate(codes):
        merged_code = code.replace('.', '').upper()

        #string trim for all morphology descriptions and lower case
        morpho_desc = morphologies_desc[i].strip().lower()

        #if the morphology is of lymphoma, leukmia, myeloid call the cancer types that
        if morpho_desc in ["myeloid", "b_cell_neoplasms", "t_cell_and_nk_cell_neoplasms", "hodgkin_lymphoma", "mast_cell_tumours", "erdheim-chester", "histiocytes_and_accessory_lymphoid_cells", "unspecified_types"]:
             cancer_types.append("HEMATOPOIETIC")
        else:
            matched = False
            for grp, grp_codes in COLLAPSED_CANCER_TYPES.items():
                if merged_code in grp_codes:
                    cancer_types.append(grp)
                    matched = True
                    break
            if not matched:
                print(f"Topology code: {code}: not found in collapsed cancer types")
                cancer_types.append("OTHER")

    return cancer_types
    

def is_systemic_cancer(code):
    """Check if a code is a systemic cancer"""
    grp = get_morpho_group(code)
    return grp in SYSTEMIC_CANCERS

def find_paired_organ(code):
    """Find paired organ for a given code"""
    for organ, codes in PAIRED_ORGANS.items():
        if code in codes:
            return organ
    return None


def expand_syndromic_morphologies(morphology_ranges):
    """Expand a list of morphology ranges to list of strings with M prefix."""
    expanded = []
    for code_range in morphology_ranges:
        expanded.extend(expand_code_range(code_range))
    return expanded


def is_synchronous(pid, start_dates):

    ''' Determine if cancer pairs are synchronous or not metachronous'''
    if len(start_dates) < 2:
        logging.warning(f"Start dates are not of length 2 for {pid}")
        return None

    #look through all pairs of start dates
    start_dates_sorted = sorted(start_dates)

    for i in range(len(start_dates_sorted) - 1):   
        if abs(start_dates_sorted[i + 1] - start_dates_sorted[i])< 183:
            return 'synchronous'

    return 'metachronous'



def has_missing_critical_info(r):
    for f in ('Patient_ID','TOPOLOGY','MORPHOLOGY'):
        if pd.isna(r[f]) or r[f] == '' or r[f] == 'UNKNOWN' or r[f] == None:
            return True
    return False



# #Extract info from the format of the tumor registry
def extract_info(df):
    df['MORPHOLOGY'] = df['DX_DESCRIPTION'].str.extract(r'\(\s*(M(?:\d+|[A-Z]{2})/\d+)\s*\|')
    df['MORPHOLOGY'] = df['MORPHOLOGY'].fillna("UNKNOWN")
    df['TOPOLOGY'] = df['DX_DESCRIPTION'].str.extract(r'\(.*\|\s*(C\d+)\s*\)')
    df['TOPOLOGY'] = df['TOPOLOGY'].fillna("UNKNOWN")
    df['Patient_ID'] = df['PATIENT_ID']
    df['START_DATE'] = df['START_DATE'] #numeric offset of dates. Can be offset from anything (diagnosis, birth, etc.) as long as it's consistent within patients and allows us to determine order and time differences between tumors.
    df['SUMMARY'] = df['SUMMARY'] if 'SUMMARY' in df.columns else pd.NA
    df['STAGE_CDM_DERIVED_GRANULAR'] = df['STAGE_CDM_DERIVED_GRANULAR']
    df['LATERALITY'] = df['LATERALITY'] if 'LATERALITY' in df.columns else pd.NA
    df['STAGE_CDM_DERIVED_GRANULAR'] = pd.to_numeric(
        df['STAGE_CDM_DERIVED_GRANULAR'], errors='coerce'
    ).fillna(-1)
    return df

#check for colon-rectum or prostate-bladder exceptions
def check_cr_pb_exceptions(primaries, topos, morphos, start_dates, morphologies_descrptions, stage_data, summary_data, sides):
       
        dedup_reasons = []
        # Check if the lengths of all lists are equal
        if len(primaries) != len(topos) or len(primaries) != len(morphos) or len(primaries) != len(start_dates) or len(primaries) != len(morphologies_descrptions):
            logging.warning(f"cr-br fields not the same length")
            return primaries, topos, morphos, start_dates, morphologies_descrptions, dedup_reasons, stage_data, summary_data, sides


        # First collect all PROSTATE primaries with adenocarcinoma morphology and their start date:
        prostate_index = [
            i for i, p in enumerate(primaries)
            if len(p) >= 2 and p[0] == 'PROSTATE' and p[1] == 'adenocarcinomas'
            ]

        # If there are any prostate adenocarcinomas, we need to check for bladder adenocarcinomas
        # that have the same start date as the prostate adenocarcinoma.
        # If so, we will remove the bladder adenocarcinoma.
        if prostate_index:
            prostate_start_time = [start_dates[i] for i in prostate_index][0] #can't have more than 1 prostate adenocarcinoma
            primaries_to_remove = ("BLADDER,_BLADDER_NECK,_URACHUS_,_URETERIC_ORIFIC", 'adenocarcinomas')
            bladder_adenocarcinomas_indices = [
                i for i, p in enumerate(primaries)
                if len(p) >= 2 and p[0] == primaries_to_remove[0] and p[1] == primaries_to_remove[1]
            ]
            # Now we need to check if any of the bladder adenocarcinomas have the same start date as the prostate adenocarcinoma
            indices_to_remove = []
            for index in bladder_adenocarcinomas_indices:
                if start_dates[index] == prostate_start_time:
                    indices_to_remove.append(index)
       
            # Remove the bladder adenocarcinomas that have the same start date as the prostate adenocarcinoma
            for index in sorted(set(indices_to_remove), reverse=True):
                del primaries[index]
                del topos[index]
                del morphos[index]
                del start_dates[index]
                del stage_data[index]
                del summary_data[index]
                del morphologies_descrptions[index]
                del sides[index]

            # If we removed any indices, we add a reason for deduplication
            if indices_to_remove:
                dedup_reasons.append("PROSTATE_BLADDER")


        # First collect all RECTUM primaries with their morphology and start date:
        rectum_entries = [
        (i, p[1], start_dates[i])  # index, morphology, start_date
        for i, p in enumerate(primaries)
        if len(p) >= 2 and p[0] == 'RECTUM'
        ]

        # Collect all COLON primaries with morphology and start date: 
        colon_entries = [
            (i, p[1], start_dates[i])  # index, morphology, start_date
            for i, p in enumerate(primaries)
            if len(p) >= 2 and p[0] in TOPOLOGY_GROUPS['ASCENDING_COLON_,_CECUM_,_COLON,_DESCENDING_COLON_,_RECTOSIGMOID_JUNCTION_,_SIGMOID_COLON_,_SPLENIC_FLEXURE_OF_COLON_,_TRANSVERSE_COLON']
        ]
 
        # If there are both rectum and colon entries, we need to check for duplicates.
        # If the morphology and start date are the same, we will remove the one with the lower stage.
        indices_to_remove = []
        if rectum_entries and colon_entries:
            for r_index, r_morph, r_date in rectum_entries:
                for c_index, c_morph, c_date in colon_entries:
                    if r_morph == c_morph and r_date == c_date:
                        c_stage = stage_data[c_index]  # stage for colon entry
                        r_stage = stage_data[r_index]  # stage for rectum entry
                        if c_stage < r_stage:
                            indices_to_remove.append(c_index)
                        elif r_stage < c_stage:
                            indices_to_remove.append(r_index)
                        else:
                            # If stages are equal, you can decide what to do — keep both
                            pass

            # Remove the duplicates from the lists
            for index in sorted(indices_to_remove, reverse=True):
                del primaries[index]
                del topos[index]
                del morphos[index]
                del start_dates[index]
                del morphologies_descrptions[index]
                del stage_data[index]
                del summary_data[index]
                del sides[index]

            # If we removed any indices, we add a reason for deduplication
            if indices_to_remove:
                dedup_reasons.append("COLON_RECTUM")

        return primaries, topos, morphos, start_dates, morphologies_descrptions, dedup_reasons, stage_data, summary_data, sides


# Normalize side values to 0 if missing, otherwise keep original value.  
def normalize_side(side_val):
    if pd.isna(side_val):
        return 0
    return side_val


    
# ─── Core classifier ────────────────────────────────────────────────────────────
def mpc_classifier(df, output_path_mpc, output_path_sp, no_pc_path, output_dir):
    if not MORPHOLOGY_GROUPS or not TOPOLOGY_GROUPS:
        raise ValueError("Topology or Morphology groups not defined in config.")

    df = extract_info(df)

    #Get all patient IDs for later use to see which patients were excluded
    all_patient_ids = set(df['Patient_ID'].unique())

    #only keep primary sites
    cond1 = df['MORPHOLOGY'].str.contains('/3', na=False)

    # allow in-situ DCIS breast cancers
    cond2 = df['MORPHOLOGY'].str.contains('/2', na=False) & df['DX_DESCRIPTION'].str.contains('DUCTAL CARCINOMA IN SITU|INTRADUCTAL CA', na=False, case=False) & df['TOPOLOGY'].str.startswith('C50', na=False)

    # Combine conditions
    df = df[cond1 | cond2 ].copy()

    results = []
    time_exception_keys = [] #store patient IDs with time exceptions

    # Group by Patient_ID and process each group
    for pid, group in df.groupby('Patient_ID'):

        # Initialize variables to store for each patient
        primaries = []  
        valid_topos = [] 
        valid_morphos = []
        morphologies_descrptions = []
        start_dates = []
        stage_data = []
        summary_data = []
        desc = []
        sides = []
        time_exception_added = False

        
        # Sort the group by START_DATE to ensure chronological order
        group_sorted = group.sort_values(by='START_DATE', ascending=True)
        
        # Iterate through each row in the sorted group
        for _, r in group_sorted.iterrows():
            pid = r['Patient_ID']
            topo = r['TOPOLOGY']
            topo = topo.replace('.', '').upper() #standardize code by removing dots and uppercasing
            morph = r['MORPHOLOGY']
            morph = morph.replace('.', '').upper() #standardize code by removing dots and uppercasing
            start_date = r['START_DATE']
            stage = r['STAGE_CDM_DERIVED_GRANULAR']
            summary = r['SUMMARY'] if 'SUMMARY' in r else ''
            side = normalize_side(r['LATERALITY'] if 'LATERALITY' in r else None)


            morph_grp = get_morpho_group(morph) #get morpho group
            topo_grp = get_topo_group(topo) #get topo group
            paired_organ = find_paired_organ(topo) #check if paired organ

              # Skip this tumor if it doesn’t have a valid morpho/topo group
            if morph_grp is None or topo_grp is None:
                continue
                
            if is_systemic_cancer(morph):
                key = (morph_grp,) #systemic cancers only rely on morphology.

            # Skin and colon rely on full code, not group
            elif topo_grp in EXCEPTIONS:
                key = (topo, morph_grp)

            #handle paired organs
            elif paired_organ:
                key = (topo_grp, morph_grp, side)
                
        
            #morphology and topology are the key
            else:
                key = (topo_grp, morph_grp)

            
            # Check if the key already exists in primaries and if not, add it
            if key not in primaries:
                primaries.append(key)
                valid_topos.append(topo) 
                valid_morphos.append(morph)
                start_dates.append(start_date)
                morphologies_descrptions.append(morph_grp)
                stage_data.append(stage)
                summary_data.append(summary)
                sides.append(side)

            # If the key already exists, check for time exceptions 
            else:
                # Check if the existing key is a duplicate based on time exceptions
                if  (abs(start_dates[primaries.index(key)] - start_date) > 5475) and not time_exception_added:
                    new_key = key + ('DUP',)
                    primaries.append(new_key)
                    valid_topos.append(topo) 
                    valid_morphos.append(morph)
                    start_dates.append(start_date)
                    morphologies_descrptions.append(morph_grp)
                    stage_data.append(stage)
                    summary_data.append(summary)
                    sides.append(side)
                    time_exception_added = True
                    time_exception_keys.append(pid)

                # Not in exceptions, so we need to add the DEDUP DESC
                else:
                    description = ''
                    # add to desc based on type of key
                    if is_systemic_cancer(morph):
                        description = "SYSTEMIC"
                    elif paired_organ:
                        description = "PAIRED_ORGAN_BILATERALITY"
                    elif len(key) == 2:
                        description = "MORPHO_TOPO"
                    elif len(key) == 3:
                        description = "MORPHO_TOPO_SIDE"
                    else:
                        description = "DEDUPLICATED_UNKNOWN"

                    desc.append(description)


          
        if len(primaries) == 0:
            logging.warning(f"Excluding {pid}: no valid tumors after filtering")
            continue

        #Check for exceptions like colon-rectum or prostate-bladder where we need to remove one of the primaries based on start date and stage / morphology
        primaries, valid_topos, valid_morphos, start_dates, morphologies_descrptions, dedup_reasons, stage_data, summary_data, sides = check_cr_pb_exceptions(primaries, valid_topos, valid_morphos, start_dates, morphologies_descrptions, stage_data, summary_data, sides)


        desc.extend(dedup_reasons)
        mpc_count = len(primaries) #count unique combinations
        cancer_types = get_collapsed_cancer_types(valid_topos, morphologies_descrptions) 


        
        result = {
        'Patient_ID': pid,
        'START_DATES': start_dates,
        'MPC_COUNT': mpc_count,
        'TOPOLOGIES': list(valid_topos),
        'MORPHOLOGIES': list(valid_morphos),
        'CANCER_TYPES': list(cancer_types),
        'STAGE': list(stage_data),
        'SUMMARY': list(summary_data),
        'MORPHOLOGIES_DESC': list(morphologies_descrptions),
        'SIDES': list(sides),
        'DESC_DEDUP': desc,
        }

        results.append(result)
    
    if not results:
        print("No patients passed filtering.")
        return    

    results_df = pd.DataFrame(results)
    


    included_patient_ids = set(results_df['Patient_ID'].unique())
    excluded_patient_ids = all_patient_ids - included_patient_ids

    df_sp = results_df[results_df['MPC_COUNT'] == 1]
    df_mpc_original = results_df[results_df['MPC_COUNT'] > 1]
    df_mpc = df_mpc_original.copy()

    # Add SYNCHRONOUS column
    for _, row in df_mpc.iterrows():
        pid = row['Patient_ID']
        start_dates = row['START_DATES']
        synchronous = is_synchronous(pid, list(start_dates))

        # Add SYNCHRONOUS column
        if synchronous is not None:
            df_mpc.loc[df_mpc['Patient_ID'] == pid, 'SYNCHRONOUS'] = synchronous

    df_sp.to_csv(output_path_sp, sep='\t', index=False)
    df_mpc.to_csv(output_path_mpc, sep='\t', index=False)


    with open(no_pc_path, "w") as file:
            for pid in excluded_patient_ids:
                file.write(f"{pid}\n")


    with open(f'{output_dir}/start_date_exceptions.txt', "w") as file:
            for pid in time_exception_keys:
                file.write(f"{pid}\n")

    print(f'{len(df_mpc)} Patients classified with MPC')
    print(f'{len(df_sp)} Patients classified with SPC')
    print(f'{len(excluded_patient_ids)} Patients excluded')





if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python mpc.py <input.tsv>" " <output_name> <output_dir (optional)>")
        sys.exit(1)

    df = pd.read_csv(sys.argv[1], sep='\t', dtype={'PATIENT_ID': str, 'START_DATE': int, 'STAGE_CDM_DERIVED': str}, low_memory=False)



    date = datetime.now().date()

    if len(sys.argv) < 4:
        print("No output directory provided, using default 'outputs'")
        output_dir = 'outputs'
    else:
        output_dir = sys.argv[3]
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    run_name = sys.argv[2]
    output_path_mpc = f'{output_dir}/MPC_{run_name}_{date}.tsv'
    output_path_sp = f'{output_dir}/SP_{run_name}_{date}.tsv'
    no_pc_path = f'{output_dir}/NO_PC_{run_name}_{date}.txt'
    mpc_classifier(df, output_path_mpc,output_path_sp, no_pc_path, output_dir)
