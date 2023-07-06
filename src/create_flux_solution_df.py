import pandas as pd

def create_flux_solution_df(model, solution):
    # make a list of dictionaries with the reaction id, name, flux, and absolute flux
    reactions = []
    for reaction_id, flux in solution.fluxes.items():
      reactions.append({
        'reaction_id': reaction_id,
        'reaction_name': model.reactions.get_by_id(reaction_id).name,
        'full_reaction': model.reactions.get_by_id(reaction_id).reaction,
        'flux': flux,
        'absolute_flux': abs(flux), # use for sorting, then drop
      })

    # make a dataframe from the list of dictionaries
    solution_df = pd.DataFrame(reactions)

    # sort the dataframe by absolute flux
    solution_df = solution_df.sort_values(by=['absolute_flux'], ascending=False)

    # drop the absolute flux column
    solution_df = solution_df.drop(columns=['absolute_flux'])

    return solution_df
    