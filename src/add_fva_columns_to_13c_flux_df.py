import pandas as pd

def add_fva_columns_to_13c_flux_df(central_flux_df, genome_scale_df, LB_column_name, UB_column_name):
    updated_df = central_flux_df.copy()

    LB_list = []
    UB_list = []
    
    # loop over rows in the central flux dataframe
    for _, row in central_flux_df.iterrows():
        
        # check if row['reaction_ids'] is nan
        if pd.isna(row['reaction_ids']):
            # if it is, add a blank values to the bounds lists
            LB_list.append('')
            UB_list.append('')
        else:
            # Add the flux value for each row to the column values list
            reaction_ids = row['reaction_ids']
            
            # get the lower and upper bounds
            LB, UB = reaction_ids_to_bounds(genome_scale_df, reaction_ids, LB_column_name, UB_column_name)


            # add the lower and upper bounds to the lists
            LB_list.append(LB)
            UB_list.append(UB)

    # add the column to the dataframe
    updated_df[LB_column_name] = LB_list
    updated_df[UB_column_name] = UB_list
    
    return updated_df

def reaction_ids_to_bounds(genome_scale_df, reaction_ids, LB_column_name, UB_column_name):
    total_lower_bound = 0
    total_upper_bound = 0
    for x in [x.strip(' ') for x in reaction_ids.split(' or ')]:
        and_split = [y.strip(' ') for y in x.split(' and ')]
        lower_bounds = [reaction_id_to_bounds(v, genome_scale_df, LB_column_name, UB_column_name)[0] for v in and_split]
        upper_bounds = [reaction_id_to_bounds(v, genome_scale_df, LB_column_name, UB_column_name)[1] for v in and_split]
        total_lower_bound += min(lower_bounds)
        total_upper_bound += max(upper_bounds)
            
    return total_lower_bound, total_upper_bound

def reaction_id_to_bounds(reaction_id, genome_scale_df, LB_column_name, UB_column_name):
    
    # check if the reaction is a reverse reaction 
    is_reverse = reaction_id.startswith('reverse_')

    # if the reaction is a reverse reaction, remove the reverse_ from the reaction id
    reaction_id = reaction_id.replace('reverse_', '')

    # get the row from the genome-scale dataframe that corresponds to the reaction id
    reaction_row = genome_scale_df[genome_scale_df['reaction_id'] == reaction_id].iloc[0]

    # get the flux value from the reaction row
    lower_bound = reaction_row[LB_column_name]
    upper_bound = reaction_row[UB_column_name]

    if is_reverse:
        return -1 * upper_bound, -1 * lower_bound
    else:
        return lower_bound, upper_bound