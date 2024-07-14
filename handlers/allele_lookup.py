from typing import Dict, Union, List

import tidytcells as tt

from functions.text import slugify


def find_allele_match(raw_input:str, alleles):
    allele_slug = None
    allele_group = None
    allele_number = None
    allele_data = None
    locus = None

    suggestion = False
    match = False
    clean_input = tt.mh.standardize(raw_input, precision='protein')

    # it tidytcells can't clean the input, it will be set to None, it may be that someone has typed an extra character
    if clean_input is None:
        clean_input = tt.mh.standardize(raw_input[:-1], precision='protein')
        if clean_input is not None:
            allele_number = clean_input
            suggestion = True
    else:
        allele_number = clean_input
        match = True

    if not ':' in clean_input:
        allele_number = clean_input + ':01'  
        suggestion = True
        match = False

    if clean_input is None:
        return None
    else:
        allele_slug = slugify(allele_number)
        allele_group = clean_input.split(':')[0]
        locus = clean_input.split('*')[0]
        locus_slug = slugify(locus)
    
        if allele_slug in alleles[locus_slug]:
            allele_data = alleles[locus_slug][allele_slug]
        else:
            allele_data = None  

    return {
        'allele_slug': allele_slug,
        'allele_group': allele_group,
        'allele_group_slug': slugify(allele_group),
        'allele_number': allele_number,
        'locus': locus, 
        'locus_slug': locus_slug,
    }, suggestion, match



def allele_lookup(request_data:Dict, app_data:Dict) -> Dict:
    """
    This function takes a request object and returns a response object.
    """
    # we'll initialise the locus and allele_data variables to None


    raw_input = request_data['allele_number_query']

    if raw_input:
        allele_info, suggestion, match = find_allele_match(raw_input, app_data.copy()['protein_alleles'])
    else:
        allele_info = None
        suggestion = False
        match = False
        

    return {
        'request_data': request_data,
        'match_type': {
            'suggestion': suggestion,
            'match': match
        },
        'allele_info': allele_info
    }