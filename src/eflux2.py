import pandas as pd
import cobra
from get_reaction_transcript_dictionary import get_reaction_transcript_dictionary

def eflux2(
    model=None, 
    reaction_transcript_dictionary={}, 
    medium={},
    carbon_uptake_reaction='EX_glc(e)', 
    carbon_uptake_flux=100
):
    """Calculate eflux2 for a given a gene transcript dictionary"""
    with model:
        model.medium = medium

        # set the flux bounds for each reaction using the transcriptomics data    
        for r in model.reactions:
            if 'EX_' not in r.id:
                # if r.gene_reaction_rule:
                if r.lower_bound < 0.0:
                    r.lower_bound = -reaction_transcript_dictionary[r.id]

                if r.upper_bound > 0.0:
                    r.upper_bound = reaction_transcript_dictionary[r.id]

        # find the optimal solution
        eflux2_solution = cobra.flux_analysis.pfba(model)

        # get the carbon uptake flux value
        uptake_flux = -1 * eflux2_solution.fluxes[carbon_uptake_reaction]

        # scale all fluxes, so uptake flux is 100
        for reaction_id, flux in eflux2_solution.fluxes.items():
            eflux2_solution.fluxes[reaction_id] = flux * (carbon_uptake_flux / uptake_flux)

        return eflux2_solution

    #     # make a list of dictionaries with the reaction id, name, flux, and absolute flux
    #     reactions = []
    #     for reaction_id, flux in eflux2_solution.fluxes.items():
    #       reactions.append({
    #         'reaction_id': reaction_id,
    #         'reaction_name': model.reactions.get_by_id(reaction_id).name,
    #         'full_reaction': model.reactions.get_by_id(reaction_id).reaction,
    #         'flux': flux,
    #         'absolute_flux': abs(flux), # use for sorting, then drop
    #       })

    #     # make a dataframe from the list of dictionaries
    #     eflux2_solution_df = pd.DataFrame(reactions)

    #     # sort the dataframe by absolute flux
    #     eflux2_solution_df = eflux2_solution_df.sort_values(by=['absolute_flux'], ascending=False)

    #     # drop the absolute flux column
    #     eflux2_solution_df = eflux2_solution_df.drop(columns=['absolute_flux'])

    #     # get the carbon uptake flux value
    #     uptake_row = eflux2_solution_df[eflux2_solution_df['reaction_id'] == carbon_uptake_reaction]
    #     uptake_flux = -1 * uptake_row['flux'].iloc[0]

    #     # scale all fluxes, so uptake flux is 100
    #     eflux2_solution_df['flux'] = eflux2_solution_df['flux'] * (carbon_uptake_flux / uptake_flux)

    # return eflux2_solution_df