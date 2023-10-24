import pandas as pd
import straindesign as sd


def get_eflux2_fva_fluxes(
    model=None, 
    substrate='',
    reaction_transcript_dictionary_list=[],
    biomass_cutoff=0
):
    fva_df_list = []

    # run pfba here
    # (Please help me add pfba here)
    pfba_fva_df = get_pfba_fva_df(
        model=model, 
        substrate=substrate,
        biomass_cutoff=biomass_cutoff
    )

    # Rename columns based on the substrate and the index
    pfba_fva_df = pfba_fva_df.rename(columns={
        'flux': f'{substrate}_pfba_flux',
        'pfba_LB': f'{substrate}_pfba_LB',
        'pfba_UB': f'{substrate}_pfba_UB'
    })

    fva_df_list.append(pfba_fva_df)

    for idx, reaction_transcript_dictionary in enumerate(reaction_transcript_dictionary_list):
        fva_df = get_eflux2_fva_df(
            model=model, 
            substrate=substrate,
            reaction_transcript_dictionary=reaction_transcript_dictionary,
            biomass_cutoff=biomass_cutoff
        )
        
        # Rename columns based on the substrate and the index
        fva_df = fva_df.rename(columns={
            'flux': f'{substrate}_{idx+1}_eflux2_flux',
            'trans_LB': f'{substrate}_{idx+1}_trans_LB',
            'trans_UB': f'{substrate}_{idx+1}_trans_UB',
            'eflux2_LB': f'{substrate}_{idx+1}_eflux2_LB',
            'eflux2_UB': f'{substrate}_{idx+1}_eflux2_UB'
        })

        fva_df_list.append(fva_df)

    # Combine all dataframes into a single dataframe
    combined_df = pd.concat(fva_df_list, axis=1)

    # Assuming 'reaction_id', 'reaction_name', and 'full_reaction' are common columns across all dataframes
    # and they are identical, we can drop the duplicates and keep only one set of these columns
    common_cols = ['reaction_id', 'reaction_name', 'full_reaction']
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

    # sort the rows by absolute value of {substrate}_pfba_flux
    combined_df = combined_df.reindex(combined_df[f'{substrate}_pfba_flux'].abs().sort_values(ascending=False).index)
    
    return combined_df

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
    elif substrate == 'glycerol':
        uptake_reaction_id = 'EX_glyc_e'
        biomass_reaction_id = 'biomass_glucose'
    else:
        uptake_reaction_id = 'EX_ocdcea_e'
        biomass_reaction_id = 'biomass_oil'

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

    pfba_solution = sd.fba(model, constraints=f'{uptake_reaction_id} = {uptake_flux}', obj=biomass_reaction_id, obj_sense='maximize', pfba=1)

    constraint_string = f'{uptake_reaction_id} = {uptake_flux}, {biomass_reaction_id} >= {biomass_cutoff}'

    print(f'Running pFBA FVA with the constraints: {constraint_string}:')

    # run pFBA FVA with larger biomass flux constraint 
    pfba_fva_solution = sd.fva(
      model, 
      constraints=constraint_string,
    )

    # maybe scale the fluxes to 100 uptake for oleic acid here
    if substrate == 'oleic_acid':
        pfba_fva_solution = pfba_fva_solution * 10

    # make a list of dictionaries with the reaction id, name, flux, and absolute flux
    reactions = []
    for reaction_id, flux in pfba_solution.fluxes.items():
      # add the reaction info to the list of dictionaries
      reactions.append({
        'reaction_id': reaction_id,
        'reaction_name': model.reactions.get_by_id(reaction_id).name,
        'full_reaction': model.reactions.get_by_id(reaction_id).reaction,
        'flux': 10 * flux if substrate == 'oleic_acid' else flux,
        'pfba_LB': pfba_fva_solution.loc[reaction_id, 'minimum'],
        'pfba_UB': pfba_fva_solution.loc[reaction_id, 'maximum'],
      })

    # make a dataframe from the list of dictionaries
    pfba_df = pd.DataFrame(reactions)

    return pfba_df


def get_eflux2_fva_df(
    model=None, 
    substrate='', 
    reaction_transcript_dictionary={},
    biomass_cutoff=0
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

    # calculate the scale factor to normalize cutoff to 100 uptake
    normalized_fva_biomass_cutoff = -(uptake_flux * biomass_cutoff) / 100

    # calculate the scale factor to normalize to 100 uptake
    scale_factor = -100 / uptake_flux

    constraint_string = f'{uptake_reaction_id} = {uptake_flux}, {biomass_reaction_id} >= {normalized_fva_biomass_cutoff}'

    print(f'Running {substrate} E-Flux2 FVA with the constraints: {constraint_string}:')

    # run E-Flux2 FVA with larger biomass flux constraint 
    eflux_2_fva_solution = sd.fva(
      model, 
      constraints=constraint_string,
    )

    # normalize the fluxes to 100
    eflux_2_fva_solution = eflux_2_fva_solution * scale_factor

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
        'eflux2_LB': eflux_2_fva_solution.loc[reaction_id, 'minimum'],
        'eflux2_UB': eflux_2_fva_solution.loc[reaction_id, 'maximum'],
      })

    # make a dataframe from the list of dictionaries
    eflux2_df = pd.DataFrame(reactions)

    return eflux2_df