import pandas as pd
import straindesign as sd

# this function returns a dataframe with rows for every GSM reaction
# each row has cols for reaction_id, reaction_name, full_reaction, flux, trans_LB, trans_UB, box_LB, box_UB, whisker_LB, and whisker_UB
def get_pfba_fva_df(
        model=None, 
        substrate='', 
        biomass_cutoff=0
):
    # make a copy of the model
    model = model.copy()

    # update the uptake reaction id based on the substrate
    if substrate == 'glucose':
        uptake_reaction_id = 'EX_glc_e'
        biomass_reaction_id = 'biomass_glucose'
        other_biomass_id = 'biomass_oil'
    elif substrate == 'glycerol':
        uptake_reaction_id = 'EX_glyc_e'
        biomass_reaction_id = 'biomass_glucose'
        other_biomass_id = 'biomass_oil'
    else:
        uptake_reaction_id = 'EX_ocdcea_e'
        biomass_reaction_id = 'biomass_oil'
        other_biomass_id = 'biomass_glucose'

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
    biomass_cutoff = 0.1 * biomass_cutoff if substrate == 'oleic_acid' else biomass_cutoff

    # define constraint string
    constraint_string = f'{uptake_reaction_id} = {uptake_flux}, {biomass_reaction_id} >= {biomass_cutoff}, {other_biomass_id} = 0, biomass_C = 0, biomass_N = 0'

    # run pFBA to get the pFBA flux values
    pfba_solution = sd.fba(model, constraints=constraint_string, obj=biomass_reaction_id, obj_sense='maximize', pfba=1)

    print(f'Running pFBA FVA with the constraints: {constraint_string}:')

    # run FVA to the get the pFBA flux ranges
    pfba_fva_solution = sd.fva(
      model, 
      constraints=constraint_string,
    )
    print('ran pfba fva')

    # maybe scale the fluxes to 100 uptake for oleic acid here
    if substrate == 'oleic_acid':
        pfba_fva_solution = pfba_fva_solution * 10

    flux_col_label = f'{substrate}_GSM_flux'
    lb_col_label = f'{substrate}_GSM_LB'
    ub_col_label = f'{substrate}_GSM_UB'

    # make a list of dictionaries with the reaction id, name, flux, and absolute flux
    reactions = []
    for reaction_id, flux in pfba_solution.fluxes.items():
      # add the reaction info to the list of dictionaries
      reactions.append({
        'reaction_id': reaction_id,
        'reaction_name': model.reactions.get_by_id(reaction_id).name,
        'full_reaction': model.reactions.get_by_id(reaction_id).reaction,
        flux_col_label: 10 * flux if substrate == 'oleic_acid' else flux,
        lb_col_label: pfba_fva_solution.loc[reaction_id, 'minimum'],
        ub_col_label: pfba_fva_solution.loc[reaction_id, 'maximum'],
      })

    # make a dataframe from the list of dictionaries
    pfba_df = pd.DataFrame(reactions)

    return pfba_df