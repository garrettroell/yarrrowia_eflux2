import pandas as pd
import straindesign as sd


def get_eflux2_biomass_fluxes(
    model=None, 
    substrate='',
    reaction_transcript_dictionary_list=[]
    ):

    biomass_flux_list = []

    for reaction_transcript_dictionary in reaction_transcript_dictionary_list:
        biomass_flux = get_eflux2_biomass_flux(
            model=model, 
            substrate=substrate,
            reaction_transcript_dictionary=reaction_transcript_dictionary
        )
        biomass_flux_list.append(biomass_flux)

    return biomass_flux_list

def get_eflux2_biomass_flux(
        model=None, 
        substrate='', 
        reaction_transcript_dictionary={}
):
    # make a copy of the model
    model = model.copy()

    # update the uptake reaction id based on the substrate
    if substrate == 'glucose':
        uptake_reaction_id = 'EX_glc_e'
        biomass_reaction_id = 'biomass_glucose'
    elif substrate == 'glycerol':
        uptake_reaction_id = 'EX_glyc_e'
        biomass_reaction_id = 'biomass_glucose'
    else:
        uptake_reaction_id = 'EX_ocdcea_e'
        biomass_reaction_id = 'biomass_oil'

    # set the flux bounds for each reaction using the transcriptomics data    
    for r in model.reactions:
        if 'EX_' not in r.id:
            # if r.gene_reaction_rule:
            if r.lower_bound < 0.0:
                r.lower_bound = -reaction_transcript_dictionary[r.id]

            if r.upper_bound > 0.0:
                r.upper_bound = reaction_transcript_dictionary[r.id]

    # update the media to minimal medium with the specified sole carbon source
    medium = model.medium
    medium['EX_glc_e'] = 1000 if substrate == 'glucose' else 0
    medium['EX_glyc_e'] = 1000 if substrate == 'glycerol' else 0
    medium['EX_ocdcea_e'] = 1000 if substrate == 'oleic_acid' else 0
    medium['EX_h2o_e'] = 10000
    medium['EX_h_e'] = 10000
    medium['EX_nh4_e'] = 10000
    medium['EX_o2_e'] = 10000
    medium['EX_pi_e'] = 10000
    medium['EX_so4_e'] = 10000
    medium['trehalose_c_tp'] = 0
    model.medium = medium

    # find the optimal solution
    eflux2_solution = sd.fba(model, obj=biomass_reaction_id, obj_sense='maximize', pfba=1)

    uptake_flux = eflux2_solution[uptake_reaction_id]
    biomass_flux = eflux2_solution[biomass_reaction_id]

    # calculate the scale factor to normalize to 100 uptake
    scale_factor = -100 / uptake_flux

    normalized_biomass_flux = scale_factor * biomass_flux

    return normalized_biomass_flux


