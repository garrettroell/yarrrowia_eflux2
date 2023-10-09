import pandas as pd
import straindesign as sd
from is_active_in_fva import is_active_in_fva

# this function returns a dataframe with rows for every GSM reaction
# each row has cols for reaction_id, reaction_name, full_reaction, flux, trans_LB, trans_UB, box_LB, box_UB, whisker_LB, and whisker_UB
def get_pfba_gsm_df(model=None, substrate='', box_biomass_fraction=0.90, whisker_biomass_fraction=0.95):
    # make a copy of the model
    model = model.copy()

    # update the uptake reaction id based on the substrate
    if substrate == 'glucose':
        uptake_reaction_id = 'EX_glc_e'
    elif substrate == 'glycerol':
        uptake_reaction_id = 'EX_glyc_e'
    else:
        uptake_reaction_id = 'EX_ocdcea_e'

    # update the media to minimal medium with the specified sole carbon source
    medium = model.medium
    medium['EX_glc_e'] = 100 if substrate == 'glucose' else 0
    medium['EX_glyc_e'] = 100 if substrate == 'glycerol' else 0
    medium['EX_ocdcea_e'] = 10 if substrate == 'oleic_acid' else 0
    medium['EX_h2o_e'] = 10000
    medium['EX_h_e'] = 10000
    medium['EX_nh4_e'] = 10000
    medium['EX_o2_e'] = 10000
    medium['EX_pi_e'] = 10000
    medium['EX_so4_e'] = 10000
    medium['trehalose_c_tp'] = 0
    model.medium = medium

    # find the optimal solution
    uptake_flux = -10 if substrate == 'oleic_acid' else -100
    pfba_solution = sd.fba(model, constraints=f'{uptake_reaction_id} = {uptake_flux}', obj='biomass_C', obj_sense='maximize', pfba=1)
    biomass_flux = pfba_solution['biomass_C']
    box_biomass_flux = box_biomass_fraction * biomass_flux
    whisker_biomass_flux = whisker_biomass_fraction * biomass_flux

    # capitalize the substrate
    print(f'{substrate.capitalize()} uptake flux: {uptake_flux}.')
    print(f'Maximum biomass flux: {biomass_flux}.')
    print(f'The number of active reactions in pFBA: {sum([abs(flux) > 0.1 for flux in pfba_solution.fluxes.values()])}.')
    print()

    # run FVA with larger biomass flux constraint 
    print(f'Running FVA constraining the biomass to {100 * box_biomass_fraction}% of maximum using the constraint:')
    print(f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {box_biomass_flux}')
    fva_box_solution = sd.fva(
      model, 
      constraints=f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {box_biomass_flux}',
    )

    print(f'The number of active reactions in {100 * box_biomass_fraction}% biomass FVA: {sum([is_active_in_fva(row) for _, row in fva_box_solution.iterrows()])}')
    print()

    # Run FVA with smaller biomass flux constraint
    print(f'Running FVA constraining the biomass to {100 * whisker_biomass_fraction}% of maximum using the constraint:')
    print(f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {whisker_biomass_flux}')
    fva_whisker_solution = sd.fva(
      model, 
      constraints=f'{uptake_reaction_id} = {uptake_flux}, biomass_C >= {whisker_biomass_flux}',
    )

    print(f'The number of active reactions in {100 * whisker_biomass_fraction}% biomass FVA: {sum([is_active_in_fva(row) for _, row in fva_whisker_solution.iterrows()])}')
    
    # box and whisker labels
    box_label_lb = f'{int(100 * box_biomass_fraction)}_LB'
    box_label_ub = f'{int(100 * box_biomass_fraction)}_UB'
    whisker_label_lb = f'{int(100 * whisker_biomass_fraction)}_LB'
    whisker_label_ub = f'{int(100 * whisker_biomass_fraction)}_UB'

    # make a list of dictionaries with the reaction id, name, flux, and absolute flux
    reactions = []
    for reaction_id, flux in pfba_solution.fluxes.items():
      # add the reaction info to the list of dictionaries
      reactions.append({
        'reaction_id': reaction_id,
        'reaction_name': model.reactions.get_by_id(reaction_id).name,
        'full_reaction': model.reactions.get_by_id(reaction_id).reaction,
        'flux': 10 * flux if substrate == 'oleic_acid' else flux,
        box_label_lb: 10 * fva_box_solution.loc[reaction_id, 'minimum'] if substrate == 'oleic_acid' else fva_box_solution.loc[reaction_id, 'minimum'],
        box_label_ub: 10 * fva_box_solution.loc[reaction_id, 'maximum'] if substrate == 'oleic_acid' else fva_box_solution.loc[reaction_id, 'maximum'],
        whisker_label_lb: 10 * fva_whisker_solution.loc[reaction_id, 'minimum'] if substrate == 'oleic_acid' else fva_whisker_solution.loc[reaction_id, 'minimum'],
        whisker_label_ub: 10 * fva_whisker_solution.loc[reaction_id, 'maximum'] if substrate == 'oleic_acid' else fva_whisker_solution.loc[reaction_id, 'maximum'],
        'absolute_flux': abs(flux), # use for sorting, then drop
      })

    # make a dataframe from the list of dictionaries
    pfba_df = pd.DataFrame(reactions)

    # sort the dataframe by absolute flux
    pfba_df = pfba_df.sort_values(by=['absolute_flux'], ascending=False)

    # drop the absolute flux column
    pfba_df = pfba_df.drop(columns=['absolute_flux'])

    return pfba_df, box_biomass_flux, whisker_biomass_flux