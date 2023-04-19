from flask import Flask
import json


def create_app():
    app = Flask(__name__)
    app.data = {}
    return app


app = create_app()


@app.before_first_request
def load_data():
    """
    This is the function which loads the generated datasets which are used by the site.

    By loading them in here, we can reduce S3 calls and speed the app up significantly.
    """


@app.route('/alleles/')
@app.route('/alleles')
def alleles_home(api=False):
    """
    This is the handler for the alleles homepage. 
    """
    return app.data


@app.route('/alleles/lookup/')
@app.route('/alleles/lookup')
def alleles_lookup(api=False):
    """
    This is the handler that performs lookups for alleles

    Args:

    The arguments are provided either as querystring or post variables
    """
    return "lookup"


@app.route('/alleles/<string:species_stem>/')
@app.route('/alleles/<string:species_stem>')
def species_page(species_stem, api=False):
    """
    This is the handler for the species page, it provides a list of loci

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
    """
    return f'{species_stem=}'


@app.route('/alleles/<string:species_stem>/<string:locus>/')
@app.route('/alleles/<string:species_stem>/<string:locus>')
def locus_page(species_stem, locus, api=False):
    """
    This is the handler for the locus page, it provides a list of allele groups

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
        locus (string): the slufiied locus e.g. hla_a
    """
    return f'{species_stem=}:{locus=}'


@app.route('/alleles/<string:species_stem>/<string:locus>/<string:allele_group>/')
@app.route('/alleles/<string:species_stem>/<string:locus>/<string:allele_group>')
def allele_group_page(species_stem, locus, allele_group, api=False):
    """
    This is the handler for the allele_group page, it provides a list of alleles and information on polymorphisms/features

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
        locus (string): the slugified locus e.g. hla_a
        allele_group (string): the slugified allele group e.g. hla_a_01
    """
    return f'{species_stem=}:{locus=}:{allele_group=}'


@app.route('/alleles/<string:species_stem>/<string:locus>/<string:allele_group>/<string:allele>/')
@app.route('/alleles/<string:species_stem>/<string:locus>/<string:allele_group>/<string:allele>')
def allele_page(species_stem, locus, allele_group, allele, api=False):
    """
   This is the handler for the allele page, it provides information on that allele, including an ESMFold prediction

    Args:
        species_stem (string): the slugified MHC species stem  e.g. hla
        locus (string): the slugified locus e.g. hla_a
        allele (string): the slugified allele number e.g. hla_a_01_01
    """
    return f'{species_stem=}:{locus=}:{allele_group=}:{allele=}'


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

