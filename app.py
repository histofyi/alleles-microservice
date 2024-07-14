from typing import Dict, List, Tuple, Union
from flask import Flask, request, url_for, redirect

import os
import json
import toml

import requests
import py3Dmol

import tidytcells as tt

from io import StringIO
#from Bio.PDB import PDBParser, PDBIO, Select

from functions.decorators import templated
from functions.templating import render
from functions.text import slugify
from functions.forms import get_request_data

import handlers

import sys


def build_pmbec_matrix():
    with open ('data/pmbec_covariance_matrix.mat', 'r') as f:
        raw_matrix = f.read()

    aa_list = []
    new_matrix = {}

    max_val = 0.0
    min_val = 10.0

    for i, line in enumerate(raw_matrix.split('\n')):   
        if i == 0:
            aa_list = line.split()
        else:
            elements = [element for element in line.split()]
            if len(elements) > 1:
                new_matrix[elements[0]] = {}
                for j, element in enumerate(elements[1:]):
                    new_matrix[elements[0]][aa_list[j]] = element
                    if float(element) > max_val:
                        max_val = float(element)
                    elif float(element) < min_val:
                        min_val = float(element)


    # now perform min max normalisation of the matrix
    normalised_matrix = {}

    for aa in new_matrix:
        normalised_matrix[aa] = {}
        for sub in new_matrix[aa]:
            normalised_matrix[aa][sub] = round((float(new_matrix[aa][sub]) - min_val) / (max_val - min_val) *10, 1)

    return normalised_matrix

matrix = build_pmbec_matrix()


page_size = 25

pockets = {
        "a": ["5","59","63","66","159","163","167","171"],
        "b": ["7","9","24","25","33","34","45","60","67","70"],
        "c": ["73","74"],
        "d": ["99","114","155","156"],
        "e": ["97","114","147","152"],
        "f": ["77","80","81","84","95","116","123","143","146","147"]
}


def map_pocket(position:int) -> str:
    for pocket in pockets:
        if str(position) in pockets[pocket]:
            return pocket
    return 'o'

netmhcpan_pocket_residues = [7,9,24,45,59,62,63,66,67,69,70,73,74,76,77,80,81,84,95,97,99,114,116,118,143,147,150,152,156,158,159,163,167,171]

netmhc_pocket_labels = [map_pocket(position) for position in netmhcpan_pocket_residues]



def zero_pad(number:int) -> str:
    if number < 10:
        return f"0{number}"
    else:
        return str(number)


def pagination(records:List, page_size:int, page:int) -> Tuple[List, int]:
    start = (page - 1) * page_size
    end = page * page_size
    page_count = (len(records) // page_size) + 1
    return records[start:end], page_count


def load_json_data(dataset_name:str) -> Dict:
    """
    This is the function which loads the generated datasets which are used by the site.

    By loading them in here, we can reduce S3 calls and speed the app up significantly.
    """
    filename = f"data/{dataset_name}.json"
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            return json.load(f)
    else:
        return {}


def load_json_data_folder(dataset_name:str) -> Dict:
    """
    This is the function which loads the generated datasets which are used by the site.

    By loading them in here, we can reduce S3 calls and speed the app up significantly.
    """
    folder_name = f"data/{dataset_name}"
    if os.path.exists(folder_name):
        data = {}
        for file in os.listdir(folder_name):
            with open(f"{folder_name}/{file}", 'r') as f:
                data[file.replace('.json', '')] = json.load(f)
        return data
    else:
        return {}


def load_pandas_data(dataset_name:str) -> Dict:
    pass


def process_species_count(data:Dict) -> Dict:
    completed_species = len(data['species'])
    ipd_species = []
    for locus in data['protein_alleles']:
        species_stem = locus.split('_')[0]
        if species_stem not in ipd_species:
            ipd_species.append(species_stem)
    return {
        'completed_species': completed_species,
        'uncompleted_species': len(ipd_species) - completed_species,
        'total_species': len(ipd_species),
    }

def process_locus_count(data:Dict) -> Dict:
    completed_loci = [len(data['species'][species]['loci']) for species in data['species']]
    ipd_loci = len(data['protein_alleles'].keys())
    return {
        'completed_loci': sum(completed_loci),
        'uncompleted_loci': ipd_loci - sum(completed_loci),
        'total_loci': ipd_loci
    }


def process_allele_group_count(data:Dict, locus:str) -> Dict:
    allele_groups = {}
    for allele in data['protein_alleles'][locus]:
        allele_group = '_'.join(allele.split('_')[:3])
        if allele_group not in allele_groups:
            allele_groups[allele_group] = {'allele_count':0, 'alleles':[], 'variability':'', 'motifs':[], 'motif_count':0, 'structure_count':0, 'structures':[]}
        allele_groups[allele_group]['allele_count'] += 1
        allele_groups[allele_group]['alleles'].append(allele)
    return {
        'allele_group_count': len(allele_groups),
        'allele_count': len(data['protein_alleles'][locus]),
        'allele_groups': allele_groups
    }




def create_app():
    """
    Creates an instance of the Flask app, and associated configuration and blueprints registration for specific routes. 

    Configuration includes

    - Relevant secrets stored in the config.toml file
    - Storing in configuration a set of credentials for AWS (decided upon by the environment of the application e.g. development, live)
    
    Returns:
            A configured instance of the Flask app

    """
    app = Flask(__name__)

    app.config.from_file('config.toml', toml.load)
    # removing whitespace from templated returns    
    app.jinja_env.trim_blocks = True
    app.jinja_env.lstrip_blocks = True

    json_datasets = ['species','sets', 'peptide_length_distributions', 'simplified_motifs', 'sorted_amino_acid_distributions', 'polymorphisms_and_motifs','hla_class_i_variability']
    
    json_dataset_folders = ['protein_alleles','pocket_pseudosequences', 'gdomain_sequences', 'allele_groups', 'reference_alleles']

    pandas_datasets = []

    app.data = {'stats':{}}

    with open('forms.json', 'r') as f:
        app.data['forms'] = json.load(f)

    for dataset in json_datasets:
        app.data[dataset] = load_json_data(dataset)

    for dataset in json_dataset_folders:
        app.data[dataset] = load_json_data_folder(dataset)

    for dataset in pandas_datasets:
        app.data[dataset] = load_pandas_data(dataset)

    app.data['stats']['species'] = process_species_count(app.data)
    app.data['stats']['loci'] = process_locus_count(app.data)

    app.data['stats']['allele_groups'] = {}

    for locus in app.data['species']['homo_sapiens']['loci']:
        app.data['stats']['allele_groups'][locus] = process_allele_group_count(app.data, locus)

    app.data['stats']['motifs'] = len(app.data['sorted_amino_acid_distributions'].keys())

    dataset_size = 0

    for item in app.data:
        item_size = sys.getsizeof(app.data[item])
        dataset_size += item_size

    print (f"Data held in memory for app.data is {round(dataset_size / 1024, 1)}MB")
        
    return app


app = create_app()


@app.template_filter()
def polymorphism_information(polymorphism:str) -> str:
    locus_slug = polymorphism.split('|')[0]
    locus = locus_slug.replace('_','-').upper()
    polymorphism_details = polymorphism.split('|')[1].split('_')
    amino_acid = polymorphism_details[0]
    position = polymorphism_details[1]
    position_information = app.data['hla_class_i_variability'][locus_slug]['variability'][position]
    index = position_information['labels'].index(amino_acid)
    rarity = position_information['rarities'][index]
    if rarity == 'majority':
        info_string = f"{amino_acid} is found at position {position} in the majority of {locus} alleles."
    elif rarity == 'unique':    
        info_string = f"{amino_acid} is uniquely found at position {position} in this particular {locus} allele."
    elif 'only' in rarity:    
        info_string = f"{amino_acid} is found at position {position} in {rarity.replace('_',' ')} {locus} alleles."
    else:
        info_string = f"{amino_acid} is {rarity.replace('_',' ')}ly found at position {position} of {locus}."
    return info_string
    


@app.template_filter()
def substitution_effect(substitution:str) -> str:
    from_aa = substitution[0]
    to_aa = substitution[1]
    val = matrix[from_aa][to_aa]
    # now write a categorical label for the value
    if val < 1:
        val_name = 'strongly non-conservative'
    elif val <= 3:
        val_name = 'non-conservative'
    elif val <= 4:
        val_name = 'neutral'
    elif val <= 6:
        val_name = 'conservative'
    else:
        val_name = 'highly conservative'

    info_string = f"The change from {from_aa} to {to_aa} is a {val_name} one [{val}]."
    return info_string


@app.template_filter()
def structure_count_display(structure_count:int) -> str:
    if structure_count == 0:
        return 'No structures'
    elif structure_count == 1:
        return '1 structure'
    else:
        return f"{structure_count} structures"


@app.template_filter()
def deslugify_locus(text:str) -> str:
    return text.replace('_', '-').upper()


@app.template_filter()
def deslugify_allele_group(text:str) -> str:
    elements = text.split('_')
    return f"{elements[0]}-{elements[1]}*{elements[2]}".upper()


@app.template_filter()
def deslugify_allele(text:str) -> str:
    elements = text.split('_')
    return f"{elements[0]}-{elements[1]}*{elements[2]}:{elements[3]}".upper()


@app.template_filter()
def mhc_flurry_url(allele_number:str) -> str:
    return f"https://openvax.github.io/mhcflurry-motifs/{allele_number.upper().replace('_','-')}.html"


@app.template_filter()
def display_simple_motif(motif:Dict) -> str:
    spacer = '<span>.</span>'
    motif_string = "<table width='90%'>"
    rows = [0,1]
    motif_string += "<tr>"
    for position in motif:
        if motif[position] != []:
            motif_string += "<td width='10%' class='motif-bar'>|</td>"
        else:
            motif_string += "<td width='10%'></td>"
    motif_string += "</tr>"
    for row in rows:
        motif_string += "<tr>"
        for position in motif:
            motif_string += "<td width='10%'>"
            if motif[position] == [] and row == 0:
                motif_string += "<span class='motif-spacer'>.</span>"
            elif len(motif[position]) > row:
                motif_string += f"<span class='motif-amino-acid simplified-{motif[position][row]['grade']}-frequency'>{motif[position][row]['amino_acid']}</span>"
            motif_string += "</td>"
        motif_string += "</tr>"
    motif_string += "</table>"
    return motif_string


@app.template_filter()
def slugify_this(text:str) -> str:
    return slugify(text)


def add_prototype_message(message_type:str, text:str) -> str:
    return f"<div><{message_type}><strong>{message_type}:</strong> {text}</{message_type}></div><br />"

@app.template_filter()
def aside(text:str) -> str:
    return add_prototype_message('aside', text)

@app.template_filter()
def todo(text:str) -> str:
    return add_prototype_message('todo', text)

@app.template_filter()
def question(text:str) -> str:
    return add_prototype_message('question', text)

@app.template_filter()
def placeholder(text:str) -> str:
    return add_prototype_message('placeholder', text)




@app.route('/alleles/')
@app.route('/alleles')
@templated('alleles_home')
def alleles_home(api=False):
    """
    This is the handler for the alleles homepage. 
    """
    data = app.data.copy()
    return {
        'species':data['species'],
        'stats':data['stats'], 
        'form': app.data['forms']['allele_number_lookup']
    }


@app.route('/alleles/lookup/', methods=['GET', 'POST'])
@app.route('/alleles/lookup', methods=['GET', 'POST'])
def alleles_lookup(api=False):
    suggestions = []
    error = None
    request_data = get_request_data(request, app.data['forms']['allele_number_lookup'])
    response_dict = handlers.allele_lookup(request_data, app.data)
    if response_dict['match_type']['match']:
        return redirect(url_for('allele_page', allele=response_dict['allele_info']['allele_slug']))
    else:
        raw_input = response_dict['request_data']['allele_number_query']
        if response_dict['match_type']['suggestion']:
            error = f"Nothing matched your query <strong>\"{raw_input}\"</strong> exactly."
            if response_dict['allele_info']['allele_group'] is not None:
                suggestions.append({'type':'allele_group', 'id':response_dict['allele_info']['allele_group'], 'slug':response_dict['allele_info']['allele_group_slug']})
            if response_dict['allele_info']['allele_number'] is not None:
                suggestions.append({'type':'allele_number', 'id':response_dict['allele_info']['allele_number'], 'slug':response_dict['allele_info']['allele_slug']})

        elif response_dict['request_data']['allele_number_query'] is None:
            error = f"You didn't enter any text in the search box. Please try again."
        response_dict['raw_input'] = raw_input
        response_dict['form'] = app.data['forms']['allele_number_lookup']
        response_dict['error'] = error
        response_dict['suggestions'] = suggestions
        response_dict['static_route'] = app.config['STATIC_ROUTE']
        return render("lookup", response_dict)



@app.route('/alleles/search/', methods=['GET', 'POST'])
@app.route('/alleles/search', methods=['GET', 'POST'])
@templated('advanced_search')
def alleles_search(api=False):
    """
    This is the handler that performs advanced search for alleles

    Args:
        None
    The arguments are provided either as querystring or post variables
    """
    return {'search_term':None}    




@app.route('/alleles/species/<string:species_stem>/')
@app.route('/alleles/species<string:species_stem>')
@templated('alleles_species')
def species_page(species_stem, api=False):
    """
    This is the handler for the species page, it provides a list of loci

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
    """
    data = app.data.copy()

    if species_stem not in data['species']:
        return {
            'error': f"Species {species_stem} not found",
            'code': 404
        }
    else:
        loci = data['species'][species_stem]['loci']

        raw_allele_groups = {}

        locus_stats = {}

        for locus in loci:
            locus_group = data['stats']['allele_groups'][locus]['allele_groups']
            
            allele_count = sum([locus_group[allele_group]['allele_count'] for allele_group in locus_group])
            locus_stats[locus] = {
                'allele_group_count': len(locus_group),
                'allele_count': allele_count
            }
        

        return {
            'species': species_stem,
            'loci': locus_stats
        }


@app.route('/alleles/locus/<string:locus>/')
@app.route('/alleles/locus/<string:locus>')
@templated('alleles_locus')
def locus_page(locus, api=False):
    """
    This is the handler for the locus page, it provides a list of allele groups

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
        locus (string): the slugified locus e.g. hla_a
    """

    data = app.data.copy()

    raw_allele_groups = data['allele_groups'][locus]



    raw_structure_sets = data['sets']['allele_groups']
    raw_motifs = data['simplified_motifs']

    allele_group_summary = {}
    allele_group_count = len(raw_allele_groups)
    allele_count = 0
    structure_count = 0
    motif_count = 0

    for allele_group in raw_allele_groups:
        allele_group_summary[allele_group] = {
            'allele_count': len(raw_allele_groups[allele_group]),
            'structure_count': 0,
            'motif_count': 0
        }
        allele_count += len(raw_allele_groups[allele_group])

    for allele_group in raw_structure_sets:
        if allele_group in allele_group_summary:
            allele_group_summary[allele_group]['structure_count'] = raw_structure_sets[allele_group]['count']
            structure_count += raw_structure_sets[allele_group]['count']

    for motif in raw_motifs:
        allele_group = '_'.join(motif.split('_')[0:3])
        if allele_group in allele_group_summary:
            allele_group_summary[allele_group]['motif_count'] += 1
            motif_count += 1
    
    return {
        'locus': locus,
        'allele_groups': allele_group_summary,
        'allele_group_count': allele_group_count,
        'allele_count': allele_count, 
        'structure_count': structure_count,
        'motif_count': motif_count,
        'alt_text': '',
        'page_size': 25,
        'page_url': url_for('locus_page', locus=locus)
    }


@app.route('/alleles/allele_group/<string:allele_group>/')
@app.route('/alleles/allele_group/<string:allele_group>')
@templated('alleles_allele_group')
def allele_group_page(allele_group, api=False):
    """
    This is the handler for the allele_group page, it provides a list of alleles and information on polymorphisms/features

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
        locus (string): the slugified locus e.g. hla_a
        allele_group (string): the slugified allele group e.g. hla_a_01
    """
    locus = '_'.join(allele_group.split('_')[0:2])

    data = app.data.copy()
    
    allele_group_data = data['allele_groups'][locus][allele_group]

    reference_allele = data['reference_alleles'][locus]['allele_groups'][allele_group]

    sort_order = sorted([int(allele.split('_')[3]) for allele in allele_group_data])

    polymorphisms_and_motifs = data['polymorphisms_and_motifs']

    allele_count = len(allele_group_data)

    if allele_group in data['sets']['allele_groups']:
        structure_count = data['sets']['allele_groups'][allele_group]['count']
    else:
        structure_count = 0


    
    reference_allele_info = polymorphisms_and_motifs[locus][reference_allele]
    reference_allele_info['allele'] = reference_allele

    if reference_allele in data['sets']['alleles']:
        reference_allele_info['structure_count'] = data['sets']['alleles'][reference_allele]['count']


    raw_alleles = [f"{allele_group}_{zero_pad(number)}" for number in sort_order]

    experimental_motif_count = 0
    inferred_motif_count = 0

    for allele in raw_alleles:
        if allele in polymorphisms_and_motifs[locus]:
            if 'motif_type' in polymorphisms_and_motifs[locus][allele]:
                if polymorphisms_and_motifs[locus][allele]['motif_type'] == 'experimental':
                    experimental_motif_count += 1
                elif polymorphisms_and_motifs[locus][allele]['motif_type'] == 'infered':
                    inferred_motif_count += 1

            

    alleles = []

    raw_alleles = [allele for allele in raw_alleles if allele != reference_allele]

    current_page = 1

    if 'page_number' in request.args:
        current_page = int(request.args['page_number'])

    paged_alleles, page_count = pagination(raw_alleles, page_size, current_page)

    pages = [i for i in range(1, page_count + 1)]

    for allele in paged_alleles:
        allele_info = polymorphisms_and_motifs[locus][allele]
        allele_info['allele'] = allele
        if allele in data['sets']['alleles']:
            allele_info['structure_count'] = data['sets']['alleles'][allele]['count']
        alleles.append(allele_info)
    

    return {
        'locus': locus,
        'allele_group': allele_group,
        'reference_allele_info': reference_allele_info,
        'alleles': alleles,
        'known_motif_count': experimental_motif_count,
        'inferred_motif_count': inferred_motif_count,
        'structure_count': structure_count,
        'allele_count': allele_count,
        'page_size': page_size,
        'page_count': page_count,
        'pages': pages,
        'current_page': current_page,
        'page_url': url_for('allele_group_page', allele_group=allele_group)
    }


@app.route('/alleles/allele/<string:allele>/')
@app.route('/alleles/allele/<string:allele>')
@templated('alleles_allele')
def allele_page(allele, api=False):
    """
   This is the handler for the allele page, it provides information on that allele, including an ESMFold prediction

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
        locus (string): the slugified locus e.g. hla_a
        allele (string): the slugified allele number e.g. hla_a_01_01
    """
    locus = '_'.join(allele.split('_')[0:2])
    allele_group = '_'.join(allele.split('_')[0:3])

    data = app.data.copy()

    allele_data = data['protein_alleles'][locus][allele]
    gdomain_matches = data['gdomain_sequences'][locus][allele_data['gdomain_sequence']]['alleles']


    pocket_pseudosequence_matches = data['pocket_pseudosequences'][locus][allele_data['pocket_pseudosequence']]['alleles']
    pocket_pseudosequence_match_alleles = []
    cleaned_pocket_pseudosequence_matches = {}

    for pocket_pseudosequence_match in pocket_pseudosequence_matches:

        match_allele_slug = slugify(pocket_pseudosequence_match['protein_allele_name'])
        if match_allele_slug not in pocket_pseudosequence_match_alleles and match_allele_slug != allele:
            pocket_pseudosequence_match_alleles.append(match_allele_slug)
            cleaned_pocket_pseudosequence_matches[match_allele_slug] = pocket_pseudosequence_match


    if allele in data['sets']['alleles']:
        structures = data['sets']['alleles'][allele]['members']
    else:
        structures = []

    # Put something in here for inferred motif
    polymorphisms_and_motifs = data['polymorphisms_and_motifs']

    motif_allele = None
    motif_type = None

    if allele in polymorphisms_and_motifs[locus]:
        print (polymorphisms_and_motifs[locus][allele])
        if 'motif_type' in polymorphisms_and_motifs[locus][allele]:
            if polymorphisms_and_motifs[locus][allele]['motif_type'] == 'experimental':
                motif_allele = allele
                motif_type = 'experimental'
            elif polymorphisms_and_motifs[locus][allele]['motif_type'] == 'infered':
                motif_allele = polymorphisms_and_motifs[locus][allele]['motif_allele']
                motif_type = 'infered'
    if motif_allele:
        if motif_allele in data['sorted_amino_acid_distributions']:
            raw_motif = data['sorted_amino_acid_distributions'][motif_allele]
            processed_motif = data['sorted_amino_acid_distributions'][motif_allele]['9']
        else:
            processed_motif = None
    else:
        processed_motif = None

    reference_allele = data['reference_alleles'][locus]['allele_groups'][allele_group]
    if allele != reference_allele:
        polymorphisms = data['polymorphisms_and_motifs'][locus][allele]['polymorphisms']

        netmhcpan_polymorphisms = [polymorphism['position'] for polymorphism in polymorphisms['binding_pocket']]

        abd_polymorphisms = []

        for polymorphism in polymorphisms['abd']:
            if polymorphism['position'] not in netmhcpan_polymorphisms:
                print (polymorphism)
                abd_polymorphisms.append(polymorphism)

        polymorphisms['abd'] = abd_polymorphisms
    else:
        polymorphisms = None
        netmhcpan_polymorphisms = None
    reference_allele = data['reference_alleles'][locus]['allele_groups'][allele_group]

    polymorphism_view = polymorphism_structure_viewer(locus, allele, reference_allele, netmhcpan_polymorphisms)

    return {
        'locus': locus,
        'allele_group': allele_group,
        'allele': allele,
        'allele_data': allele_data,
        'motif': '',
        'pocket_pseudosequence_matches': cleaned_pocket_pseudosequence_matches,
        'pocket_pseudosequence_match_count': len(cleaned_pocket_pseudosequence_matches),
        'pocket_pseudosequence_positions': netmhcpan_pocket_residues,
        'netmhc_pocket_labels': netmhc_pocket_labels,
        'netmhcpan_polymorphisms': netmhcpan_polymorphisms,
        'polymorphisms': polymorphisms,
        'reference_allele': reference_allele,
        'structures': structures,
        'structure_count': len(structures),
        'processed_motif': processed_motif,
        'motif_allele': motif_allele,
        'motif_type': motif_type,
        'polymorphism_view': polymorphism_view,
        'page_size': 25,
        'page_url': url_for('allele_page', allele=allele)
    }




def polymorphism_structure_viewer(locus:str, allele_slug:str, reference_allele_slug:str, polymorphisms:List) -> str:
    
    reference_url = f"https://coordinates.histo.fyi/predictions/view/class_i/{locus}/{reference_allele_slug}_canonical.pdb"
    print (reference_url)
    
    allele_url = f"https://coordinates.histo.fyi/predictions/view/class_i/{locus}/{allele_slug}_canonical.pdb"
    print (allele_url)

    view = py3Dmol.view(width=800, height=500)

    reference_request = requests.get(reference_url)
    reference_structure = reference_request.text

    #print (len(reference_structure))

    allele_request = requests.get(allele_url)
    allele_structure = allele_request.text   

    #print (len(allele_structure)) 

    #reference_polymporphisms = extract_polymorphic_residues(reference_structure, polymorphisms)
    #allele_polymorphisms = extract_polymorphic_residues(allele_structure, polymorphisms)

    #print (len(reference_polymporphisms))
    #print (len(allele_polymorphisms))

    #print (allele_polymorphisms)

    view.addModelsAsFrames(reference_structure)
    view.addModelsAsFrames(reference_structure)
    view.addModelsAsFrames(allele_structure)


    view.setStyle({'model': 0}, {"cartoon": {'colorscheme': 'grey'}})
    view.setStyle({'model': 1}, {"stick": {'colorscheme': 'greyCarbon'}})
    view.setStyle({'model': 2}, {"stick": {'colorscheme': 'yellowCarbon'}})
    view.zoomTo()

    return view.write_html()


@app.route('/alleles/identifier/<string:datasource>/<string:identifier>/')
@app.route('/alleles/identifier/<string:datasource>/<string:identifier>')
def allele_identifier_page(datasource, identifier, api=False):
    """
   This is the handler for the allele page, it provides information on that allele, including an ESMFold prediction. This version takes a datasource and identifier combination

    Args:
        datasource (string): the slugified datasource  e.g. ipd_imgt
        identifier (string): the slugified identifier e.g. hla00001
    """
    return {}

