import pandas as pd

def reaction_id_to_flux(reaction_id, solution_df):
    # 
    is_reverse = reaction_id.startswith('reverse_')

    # if the reaction is a reverse reaction, remove the reverse_ from the reaction id
    reaction_id = reaction_id.replace('reverse_', '')

    # get the row from the solution dataframe that corrosponds to the reaction id
    reaction_row = solution_df[solution_df['reaction_id'] == reaction_id].iloc[0]

    # get the flux value from the reaction row
    flux_value = reaction_row['flux']

    if is_reverse:
        return -1 * flux_value
    else:
        return flux_value

# This function takes in a row from a central flux dataframe and a string 
# of reaction ids, and it returns a flux value that corrosponds to the reaction ids.
def reaction_ids_to_flux_value(solution, reaction_ids):
    total_flux = 0
    for x in [x.strip(' ') for x in reaction_ids.split(' or ')]:
        and_split = [y.strip(' ') for y in x.split(' and ')]
        total_flux += min([reaction_id_to_flux(v, solution) for v in and_split])
            
    return total_flux
    
# This function takes in a 13C flux dataframe and a cobra solution. It determines the 
# flux value using a reaction ids string that maps the 13C-MFA reaction to GSM reactions/
# It uses the reaction_ids_to_flux_value function to calculate the GSM flux.
def add_column_to_13c_flux_df(central_flux_df, solution, column_name):
    updated_df = central_flux_df.copy()
    
    # create a blank list to hold values to add to the column
    column_values = []
    
    # loop over rows in the central flux dataframe
    for _, row in central_flux_df.iterrows():
        # Add the flux value for each row to the column values list
        reaction_ids = row['reaction_ids']
        flux_value = reaction_ids_to_flux_value(solution, reaction_ids)
        column_values.append(flux_value)

    # add the column to the dataframe
    updated_df[column_name] = column_values
    
    return updated_df