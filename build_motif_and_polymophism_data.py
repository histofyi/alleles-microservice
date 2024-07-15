from typing import Dict, List, Union

import json

def map_pocket(position:int) -> str:
    for pocket in pockets:
        if str(position) in pockets[pocket]:
            return pocket
    return 'o'


def check_position_for_polymorphism(reference, test):
    if test == '-':
        return False
    if reference == test:
        return False
    else:
        return True




# constants

pockets = {
        "a": ["5","59","63","66","159","163","167","171"],
        "b": ["7","9","24","25","33","34","45","60","67","70"],
        "c": ["73","74"],
        "d": ["99","114","155","156"],
        "e": ["97","114","147","152"],
        "f": ["77","80","81","84","95","116","123","143","146","147"]
}

netmhcpan_pocket_residues = [7,9,24,45,59,62,63,66,67,69,70,73,74,76,77,80,81,84,95,97,99,114,116,118,143,147,150,152,156,158,159,163,167,171]
netmhc_pocket_labels = [map_pocket(position) for position in netmhcpan_pocket_residues]











def build_sequence_polymorphism_data(reference_sequence, test_sequence):
    i = 0
    polymorphisms = []
    for position in reference_sequence:
        if i >= len(test_sequence):
            break
        if check_position_for_polymorphism(reference_sequence[i], test_sequence[i]):
            polymorphism = {'position': i + 1, 'from': reference_sequence[i], 'to': test_sequence[i]}
            polymorphisms.append(polymorphism)
        i += 1
    return polymorphisms
    

def build_allele_polymorphism_data(reference_pocket_pseudosequence, reference_cytoplasmic_sequence, pocket_pseudosequence, cytoplasmic_sequence):

    binding_pocket_polymorphisms = build_sequence_polymorphism_data(reference_pocket_pseudosequence, pocket_pseudosequence)

    for polymorphism in binding_pocket_polymorphisms:
        polymorphism['position'] = netmhcpan_pocket_residues[int(polymorphism['position']) - 1]
        polymorphism['pocket'] = map_pocket(polymorphism['position'])

    cytoplasmic_polymorphisms = build_sequence_polymorphism_data(reference_cytoplasmic_sequence, cytoplasmic_sequence)

    abd_polymorphisms = []
    non_abd_polymorphisms = []

    if cytoplasmic_polymorphisms != []:
        for polymorphism in cytoplasmic_polymorphisms:
            if polymorphism not in abd_polymorphisms:
                if polymorphism['position']  > 180:
                    non_abd_polymorphisms.append(polymorphism)
                else:
                    abd_polymorphisms.append(polymorphism)
    
    return {'binding_pocket': binding_pocket_polymorphisms, 'abd': abd_polymorphisms, 'non-abd': non_abd_polymorphisms}
        


def build_motif_and_polymophism_data(locus:str) -> Dict:

    reference_alleles = {}
    allele_groups = {}
    protein_alleles = {}
    pocket_polymorphisms = {}
    alleles = {}

    pocket_pseudosequences = {}
    locus_motifs = []

    with open("data/allele_groups/" + locus + ".json", "r") as allele_groups_file:
        allele_groups = json.load(allele_groups_file)

    with open("data/reference_alleles/" + locus + ".json", "r") as reference_alleles_file:
        reference_alleles = json.load(reference_alleles_file)['allele_groups']

    with open("data/protein_alleles/" + locus + ".json", "r") as protein_alleles_file:
        protein_alleles = json.load(protein_alleles_file)

    with open("data/simplified_motifs.json", "r") as motifs_file:
        motifs = json.load(motifs_file)


    # first we run through the motifs and build a dictionary of pocket pseudosequences relating to the motifs, and a list of the alleles that have motifs
    # we do this first as some alleles with higher allele numbers may have motifs vs lower motif numbers with the same pseudosequence
    for allele in motifs:
        if locus in allele:
            pocket_pseudosequence = protein_alleles[allele]['pocket_pseudosequence']
            pocket_pseudosequences[pocket_pseudosequence] = allele
            locus_motifs.append(allele)

    # now we run through the alleles and see if they match the reference allele, or any of the motifs/psuedosequences
    for allele_group in allele_groups:

        # first check if the allele group has a reference allele
        if allele_group in reference_alleles:
            # if it does, we set the reference allele to the first allele in the group
            reference_allele = reference_alleles[allele_group]

        # then we iterate through the alleles in the allele group
        for allele in allele_groups[allele_group]:
            
            # we create a dictionary for the allele
            alleles[allele] = {}
            
            # if the allele is in the motifs, we set the motif to the allele and mark it as experimental
            if allele in locus_motifs:
                alleles[allele]['motif_allele'] = allele
                alleles[allele]['motif_type'] = 'experimental'

            pocket_pseudosequence = protein_alleles[allele]['pocket_pseudosequence']
            cytoplasmic_sequence = protein_alleles[allele]['canonical_sequence']
            
            # if the allele matches the reference allele, we set the reference allele property to true
            if allele == reference_allele:
                alleles[allele]['reference'] = True
                print (f"Reference allele {allele}")

                # we set the reference pseudosequence to the first pseudosequence in the reference allele, this is the sequence that polymorphisms are compared to
                reference_pocket_pseudosequence = pocket_pseudosequence
                reference_cytoplasmic_sequence = cytoplasmic_sequence
            else:

                # if the allele does not match the reference allele, we check if the pseudosequence matches the reference pseudosequence
                if pocket_pseudosequence == reference_pocket_pseudosequence:

                    # if it does we set the matches property to the reference allele
                    alleles[allele]['matches'] = reference_allele

                    # we then check if the reference allele is in the motifs
                    if reference_allele in locus_motifs and 'motif' not in alleles[allele]:
                        # if it is, we set the motif to the reference allele and mark it as infered
                        alleles[allele]['motif_allele'] = reference_allele
                        alleles[allele]['motif_type'] = 'infered'
                    # we then check if the reference pseudosequence is in the pocket pseudosequences
                    elif reference_pocket_pseudosequence in pocket_pseudosequences:
                        if pocket_pseudosequences[reference_pocket_pseudosequence] in locus_motifs and 'motif' not in alleles[allele]:
                            alleles[allele]['motif_allele'] = pocket_pseudosequences[reference_pocket_pseudosequence]
                            alleles[allele]['motif_type'] = 'infered'

                    alleles[allele]['polymorphisms'] = build_allele_polymorphism_data(reference_pocket_pseudosequence, reference_cytoplasmic_sequence, pocket_pseudosequence, cytoplasmic_sequence)
                    print (f"Matches reference allele {reference_allele}")

                # next we check if the pseudosequence matches any that have been previously seen
                elif pocket_pseudosequence in pocket_pseudosequences:
                    match = pocket_pseudosequences[pocket_pseudosequence]
                    alleles[allele]['matches'] = match
                    
                    # we then check if the match is in the motifs
                    if match in locus_motifs:
                        alleles[allele]['motif_allele'] = match
                        if allele == match:
                            alleles[allele]['motif_type'] = 'experimental'
                        else:
                            alleles[allele]['motif_type'] = 'infered'
                        
                    alleles[allele]['polymorphisms'] = build_allele_polymorphism_data(reference_pocket_pseudosequence, reference_cytoplasmic_sequence, pocket_pseudosequence, cytoplasmic_sequence)
                    print (f"Matches pseudosequence {pocket_pseudosequences[pocket_pseudosequence]}")
                else:
                    # we add the pseudosequence to the dictionary of pseudosequences
                    pocket_pseudosequences[pocket_pseudosequence] = allele
                    alleles[allele]['matches'] = None   
                    alleles[allele]['polymorphisms'] = build_allele_polymorphism_data(reference_pocket_pseudosequence, reference_cytoplasmic_sequence, pocket_pseudosequence, cytoplasmic_sequence)
                    print (f"New pseudosequence {allele}")
            if 'motif_allele' in alleles[allele]:
                alleles[allele]['motif'] = motifs[alleles[allele]['motif_allele']]

    return alleles



def main():

    loci = ['hla_a','hla_b','hla_c', 'hla_e', 'hla_f', 'hla_g']

    allele_polymorphism_and_motif_data = {}

    for locus in loci:
        allele_polymorphism_and_motif_data[locus] = build_motif_and_polymophism_data(locus)

    with open("data/allele_polymorphism_and_motif_data.json", "w") as f:
        json.dump(allele_polymorphism_and_motif_data, f, indent=4)




main()