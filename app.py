from typing import Dict, List, Tuple, Union
from flask import Flask, request, url_for

import os
import json
import toml

import base64
from io import BytesIO

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np

from functions.decorators import templated
from functions.text import slugify



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


def top_n(dataset:Dict, n:int=10):
    # Sort the dictionary by percentage in descending order
    sorted_data = sorted(dataset.items(), key=lambda x: x[1]['percent'], reverse=True)
    
    # Extract the top 10 labels and values
    labels = [item[0] for item in sorted_data[:n]]
    values = [item[1]['percent'] for item in sorted_data[:n]]
    
    # Extract the rest of the labels and values for others category
    others = [item[0] for item in sorted_data[n:]]
    others_values = [item[1]['percent'] for item in sorted_data[n:]]
    
    return labels, values, others, others_values



def generate_allele_group_pie_chart(allele_groups:Dict, allele_count:int, locus:str) -> Tuple[str, str, str]:
    labels = []
    values = []
    others = []

    for allele_group in allele_groups:
        allele_groups[allele_group]['percent'] = allele_groups[allele_group]['allele_count']*100/allele_count

    labels, values, others, others_values = top_n(allele_groups, 9)

    labels = [f"{deslugify_allele_group(allele_group)} [{round(allele_groups[allele_group]['percent'], 1)}%]" for allele_group in labels]

    others_percent = sum(others_values)

    if len(others) > 0:
        labels.append('Others')
        values.append(others_percent)
    alt_text = f"Donut chart of the allele group distribution for the {locus} locus. Alleles shown indvidually are {labels}. Alleles shown in the 'Others' category are {others}."
    figsize = 15

    fig = Figure()
    fig.set_figwidth(figsize+4)
    fig.set_figheight(figsize-3)
    ax = fig.subplots()
    bbox_props = dict(boxstyle="square,pad=0.2", fc="w", ec="k", lw=0)
    wedges, texts = ax.pie(values, textprops={'fontsize': 30}, wedgeprops=dict(width=0.3), startangle=-260, counterclock=False)

    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = f"angle,angleA=0,angleB={ang}"
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),horizontalalignment=horizontalalignment, **kw, size=32)

    # Save it to a temporary buffer.
    png = BytesIO()
    svg = BytesIO()
    fig.savefig(png, format="png")


    fig.savefig(svg, format="svg")


    png_data = base64.b64encode(png.getbuffer()).decode("ascii")    
    svg_data = base64.b64encode(svg.getvalue()).decode("ascii")

    return png_data, svg_data, alt_text



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

    json_datasets = ['species','core','sets', 'peptide_length_distributions', 'simplified_motifs', 'sorted_amino_acid_distributions']
    
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

    app.files = {}

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
    
    pie_chart, svg_pie_chart, alt_text = generate_allele_group_pie_chart(allele_group_summary, allele_count, locus)


    return {
        'locus': locus,
        'allele_groups': allele_group_summary,
        'allele_group_count': allele_group_count,
        'allele_count': allele_count, 
        'structure_count': structure_count,
        'motif_count': motif_count,
        'svg_pie_chart': svg_pie_chart,
        'alt_text': alt_text,
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
    
    known_motifs = {}

    all_allele_group_alleles = data['stats']['allele_groups'][locus]['allele_groups'][allele_group]['alleles']
    paged_alleles = all_allele_group_alleles[1:26]
    reference_allele = data['stats']['allele_groups'][locus]['allele_groups'][allele_group]['alleles'][0]
    i = 0
    allele_group_pocket_sequence = data['protein_alleles'][locus][reference_allele]['pocket_pseudosequence']
    
    allele_dict = {}
    
    allele_dict[reference_allele] = data['protein_alleles'][locus][reference_allele]
    allele_dict[reference_allele]['pocket_polymorphisms'] = 'Reference'

    for allele in paged_alleles:
        allele_dict[allele] = data['protein_alleles'][locus][allele]

        current_allele_pocket_sequence = data['protein_alleles'][locus][allele]['pocket_pseudosequence']
        if current_allele_pocket_sequence == allele_group_pocket_sequence:
            allele_dict[allele]['pocket_polymorphisms'] = None
        else:
            polymorphisms = {}
            for i in range(0, len(current_allele_pocket_sequence)):
                if current_allele_pocket_sequence[i] != allele_group_pocket_sequence[i]:
                    position = netmhcpan_pocket_residues[i]
                    pocket = map_pocket(position)
                    polymorphisms[position] = {'from': allele_group_pocket_sequence[i], 'to': current_allele_pocket_sequence[i], 'pocket': pocket, 'position': position}
            allele_dict[allele]['pocket_polymorphisms'] = polymorphisms

    for allele in allele_dict:     
        if allele in data['sets']['alleles']:
            allele_dict[allele]['structure_count'] = data['sets']['alleles'][allele]['count']
        else:
            allele_dict[allele]['structures'] = None
            allele_dict[allele]['structure_count'] = None

    # THIS CURRENTLY ONLY WORKS ON THE PAGED ALLELES, DO IN THE PIPELINE FOR ALL ALLELES
    motif_index = 1

    known_motif_count = 0
    inferred_motif_count = 0
    structure_count = 0

    if allele_group in data['sets']['allele_groups']:
        structure_count = data['sets']['allele_groups'][allele_group]['count']

    for allele in allele_dict:
        if allele in data['sorted_amino_acid_distributions']:
            simplified_motif = [None,['V','L'],None,None,None,None,None,None,['L','V']]
            motif_entry = {'type': 'known', 'simplified_motif': simplified_motif, 'motif_index': motif_index, 'motif_allele': allele}
            allele_dict[allele]['motif'] = motif_entry
            pocket_pseudosequence = data['protein_alleles'][locus][allele]['pocket_pseudosequence']
            if not pocket_pseudosequence in known_motifs:
                known_motifs[pocket_pseudosequence] = []
            known_motifs[pocket_pseudosequence].append(motif_entry)
            known_motif_count += 1


    for allele in allele_dict:
        if 'motif' not in allele_dict[allele]:
            pocket_pseudosequence = data['protein_alleles'][locus][allele]['pocket_pseudosequence']
            if pocket_pseudosequence in known_motifs:
                motif_entry = known_motifs[pocket_pseudosequence][0]
                motif_entry['type'] = 'inferred'
                allele_dict[allele]['motif'] = motif_entry
                inferred_motif_count += 1

            

    return {
        'locus': locus,
        'allele_group': allele_group,
        'alleles': allele_dict,
        'known_motif_count': known_motif_count,
        'inferred_motif_count': inferred_motif_count,
        'structure_count': structure_count,
        'allele_count': len(all_allele_group_alleles),
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




def build():
    pass