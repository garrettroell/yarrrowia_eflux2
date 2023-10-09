import pandas as pd
import numpy as np
import straindesign as sd
from is_active_in_fva import is_active_in_fva

# this function returns a dataframe with rows for every GSM reaction
# each row has cols for reaction_id, reaction_name, full_reaction, flux, trans_LB, trans_UB, box_LB, box_UB, whisker_LB, and whisker_UB
def get_eflux2_gsm_df(model=None, substrate='', reaction_transcript_dictionary={}, box_biomass_flux=0, whisker_biomass_flux=0):
    # make a copy of the model
    model = model.copy()

    # update the uptake reaction id based on the substrate
    if substrate == 'glucose':
        uptake_reaction_id = 'EX_glc_e'
    elif substrate == 'glycerol':
        uptake_reaction_id = 'EX_glyc_e'
    else:
        uptake_reaction_id = 'EX_ocdcea_e'

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
    eflux2_solution = sd.fba(model, obj='biomass_C', obj_sense='maximize', pfba=1)

    uptake_flux = eflux2_solution[uptake_reaction_id]
    biomass_flux = eflux2_solution['biomass_C']

    # calculate the scale factor to normalize to 100 uptake
    scale_factor = -100 / uptake_flux

    # print data before normalization
    print(f'Pre-normalized {substrate} uptake flux: {uptake_flux}.')
    print(f'Maximum biomass flux: {biomass_flux}.')

    # print data after normalization
    print(f'The number of active reactions in pFBA: {sum([abs(scale_factor * flux) > 0.1 for flux in eflux2_solution.fluxes.values()])}.')
    print()

    # Run FVA with largweer biomass flux constraint
    print(f'Running FVA constraining the biomass to {100 * box_biomass_fraction}% of maximum using the constraint:')
    print(f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {box_biomass_fraction * biomass_flux}')
    # run E-Flux2 FVA with larger biomass flux constraint 
    fva_box_solution = sd.fva(
      model, 
      constraints=f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {box_biomass_fraction * biomass_flux}',
    )

    # normalize the fluxes to 100
    fva_box_solution = fva_box_solution * scale_factor

    # print data after normalization
    print(f'The number of active reactions in {100 * box_biomass_fraction}% biomass FVA: {sum([is_active_in_fva(row) for _, row in fva_box_solution.iterrows()])}')
    print()

    # Run FVA with smaller biomass flux constraint
    print(f'Running FVA constraining the biomass to {100 * whisker_biomass_fraction}% of maximum using the constraint:')
    print(f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {whisker_biomass_fraction * biomass_flux}')
    fva_whisker_solution = sd.fva(
      model, 
      constraints=f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {whisker_biomass_fraction * biomass_flux}',
    )

    # normalize the fluxes to 100
    fva_whisker_solution = fva_whisker_solution * scale_factor

    # print data after normalization
    print(f'The number of active reactions in {100 * whisker_biomass_fraction}% biomass FVA: {sum([is_active_in_fva(row) for _, row in fva_whisker_solution.iterrows()])}')

    # box and whisker labels
    box_label_lb = f'{int(100 * box_biomass_fraction)}_LB'
    box_label_ub = f'{int(100 * box_biomass_fraction)}_UB'
    whisker_label_lb = f'{int(100 * whisker_biomass_fraction)}_LB'
    whisker_label_ub = f'{int(100 * whisker_biomass_fraction)}_UB'

    # make a list of dictionaries with the reaction id, name, flux, and absolute flux
    reactions = []
    for reaction_id, flux in eflux2_solution.fluxes.items():
      # add the reaction info to the list of dictionaries
      reactions.append({
        'reaction_id': reaction_id,
        'reaction_name': model.reactions.get_by_id(reaction_id).name,
        'full_reaction': model.reactions.get_by_id(reaction_id).reaction,
        'flux': scale_factor * flux,
        'trans_LB': scale_factor * model.reactions.get_by_id(reaction_id).lower_bound,
        'trans_UB': scale_factor * model.reactions.get_by_id(reaction_id).upper_bound,
        box_label_lb: fva_box_solution.loc[reaction_id, 'minimum'],
        box_label_ub: fva_box_solution.loc[reaction_id, 'maximum'],
        whisker_label_lb: fva_whisker_solution.loc[reaction_id, 'minimum'],
        whisker_label_ub: fva_whisker_solution.loc[reaction_id, 'maximum'],
        'absolute_flux': abs(flux), # use for sorting, then drop
      })

    # make a dataframe from the list of dictionaries
    eflux2_df = pd.DataFrame(reactions)

    # sort the dataframe by absolute flux
    eflux2_df = eflux2_df.sort_values(by=['absolute_flux'], ascending=False)

    # drop the absolute flux column
    eflux2_df = eflux2_df.drop(columns=['absolute_flux'])

    return eflux2_df


# this function returns a dataframe with rows for every GSM reaction
# each row has cols for reaction_id, reaction_name, full_reaction, flux, trans_LB, trans_UB, box_LB, box_UB, whisker_LB, and whisker_UB
# def get_eflux2_gsm_df(model=None, substrate='', reaction_transcript_dictionary={}, box_biomass_fraction, whisker_biomass_fraction):
#     # make a copy of the model
#     model = model.copy()

#     # update the uptake reaction id based on the substrate
#     if substrate == 'glucose':
#         uptake_reaction_id = 'EX_glc_e'
#     elif substrate == 'glycerol':
#         uptake_reaction_id = 'EX_glyc_e'
#     else:
#         uptake_reaction_id = 'EX_ocdcea_e'

#     # set the flux bounds for each reaction using the transcriptomics data    
#     for r in model.reactions:
#         if 'EX_' not in r.id:
#             # if r.gene_reaction_rule:
#             if r.lower_bound < 0.0:
#                 r.lower_bound = -reaction_transcript_dictionary[r.id]

#             if r.upper_bound > 0.0:
#                 r.upper_bound = reaction_transcript_dictionary[r.id]

#     # update the media to minimal medium with the specified sole carbon source
#     medium = model.medium
#     medium['EX_glc_e'] = 1000 if substrate == 'glucose' else 0
#     medium['EX_glyc_e'] = 1000 if substrate == 'glycerol' else 0
#     medium['EX_ocdcea_e'] = 1000 if substrate == 'oleic_acid' else 0
#     medium['EX_h2o_e'] = 10000
#     medium['EX_h_e'] = 10000
#     medium['EX_nh4_e'] = 10000
#     medium['EX_o2_e'] = 10000
#     medium['EX_pi_e'] = 10000
#     medium['EX_so4_e'] = 10000
#     medium['trehalose_c_tp'] = 0
#     model.medium = medium

#     # find the optimal solution
#     eflux2_solution = sd.fba(model, obj='biomass_C', obj_sense='maximize', pfba=1)

#     uptake_flux = eflux2_solution[uptake_reaction_id]
#     biomass_flux = eflux2_solution['biomass_C']

#     # calculate the scale factor to normalize to 100 uptake
#     scale_factor = -100 / uptake_flux

#     # print data before normalization
#     print(f'Pre-normalized {substrate} uptake flux: {uptake_flux}.')
#     print(f'Maximum biomass flux: {biomass_flux}.')

#     # print data after normalization
#     print(f'The number of active reactions in pFBA: {sum([abs(scale_factor * flux) > 0.1 for flux in eflux2_solution.fluxes.values()])}.')
#     print()

#     # Run FVA with largweer biomass flux constraint
#     print(f'Running FVA constraining the biomass to {100 * box_biomass_fraction}% of maximum using the constraint:')
#     print(f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {box_biomass_fraction * biomass_flux}')
#     # run E-Flux2 FVA with larger biomass flux constraint 
#     fva_box_solution = sd.fva(
#       model, 
#       constraints=f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {box_biomass_fraction * biomass_flux}',
#     )

#     # normalize the fluxes to 100
#     fva_box_solution = fva_box_solution * scale_factor

#     # print data after normalization
#     print(f'The number of active reactions in {100 * box_biomass_fraction}% biomass FVA: {sum([is_active_in_fva(row) for _, row in fva_box_solution.iterrows()])}')
#     print()

#     # Run FVA with smaller biomass flux constraint
#     print(f'Running FVA constraining the biomass to {100 * whisker_biomass_fraction}% of maximum using the constraint:')
#     print(f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {whisker_biomass_fraction * biomass_flux}')
#     fva_whisker_solution = sd.fva(
#       model, 
#       constraints=f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {whisker_biomass_fraction * biomass_flux}',
#     )

#     # normalize the fluxes to 100
#     fva_whisker_solution = fva_whisker_solution * scale_factor

#     # print data after normalization
#     print(f'The number of active reactions in {100 * whisker_biomass_fraction}% biomass FVA: {sum([is_active_in_fva(row) for _, row in fva_whisker_solution.iterrows()])}')

#     # box and whisker labels
#     box_label_lb = f'{int(100 * box_biomass_fraction)}_LB'
#     box_label_ub = f'{int(100 * box_biomass_fraction)}_UB'
#     whisker_label_lb = f'{int(100 * whisker_biomass_fraction)}_LB'
#     whisker_label_ub = f'{int(100 * whisker_biomass_fraction)}_UB'

#     # make a list of dictionaries with the reaction id, name, flux, and absolute flux
#     reactions = []
#     for reaction_id, flux in eflux2_solution.fluxes.items():
#       # add the reaction info to the list of dictionaries
#       reactions.append({
#         'reaction_id': reaction_id,
#         'reaction_name': model.reactions.get_by_id(reaction_id).name,
#         'full_reaction': model.reactions.get_by_id(reaction_id).reaction,
#         'flux': scale_factor * flux,
#         'trans_LB': scale_factor * model.reactions.get_by_id(reaction_id).lower_bound,
#         'trans_UB': scale_factor * model.reactions.get_by_id(reaction_id).upper_bound,
#         box_label_lb: fva_box_solution.loc[reaction_id, 'minimum'],
#         box_label_ub: fva_box_solution.loc[reaction_id, 'maximum'],
#         whisker_label_lb: fva_whisker_solution.loc[reaction_id, 'minimum'],
#         whisker_label_ub: fva_whisker_solution.loc[reaction_id, 'maximum'],
#         'absolute_flux': abs(flux), # use for sorting, then drop
#       })

#     # make a dataframe from the list of dictionaries
#     eflux2_df = pd.DataFrame(reactions)

#     # sort the dataframe by absolute flux
#     eflux2_df = eflux2_df.sort_values(by=['absolute_flux'], ascending=False)

#     # drop the absolute flux column
#     eflux2_df = eflux2_df.drop(columns=['absolute_flux'])

#     return eflux2_df