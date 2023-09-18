from typing import Dict, List, Tuple, Union
from flask import Flask, request, url_for

import os
import json
import toml



from functions.decorators import templated
from functions.text import slugify

import sys

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

    json_datasets = ['species','sets', 'peptide_length_distributions', 'simplified_motifs', 'sorted_amino_acid_distributions', 'polymorphisms_and_motifs']
    
    json_dataset_folders = ['protein_alleles','pocket_pseudosequences', 'gdomain_sequences', 'allele_groups', 'reference_alleles']

    pandas_datasets = []

    app.data = {'stats':{}}

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
        'stats':data['stats']
    }


@app.route('/alleles/lookup/')
@app.route('/alleles/lookup')
def alleles_lookup(api=False):
    """
    This is the handler that performs lookups for alleles

    Args:

    The arguments are provided either as querystring or post variables
    """
    return "lookup"


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

    for allele in raw_alleles[0:25]:
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
        'page_size': 25,
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


    if allele in data['sorted_amino_acid_distributions']:
        raw_motif = data['sorted_amino_acid_distributions'][allele]
        processed_motif = data['sorted_amino_acid_distributions'][allele]['9']
    else:
        processed_motif = None




        
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
        'structures': structures,
        'structure_count': len(structures),
        'processed_motif': processed_motif,
        'page_size': 25,
        'page_url': url_for('allele_page', allele=allele)
    }


@app.route('/alleles/identifier/<string:datasource>/<string:identifier>/')
@app.route('/alleles/identifier/<string:datasource>/<string:identifier>')
def allele_identifier_page(datasource, identifier, api=False):
    """
   This is the handler for the allele page, it provides information on that allele, including an ESMFold prediction. This version takes a datasource and identifier combination

    Args:
        datasource (string): the slugified datasource  e.g. ipd_imgt
        identifier (string): the slugified identifier e.g. hla00001
    """
    return f'{datasource=}:{identifier=}'

