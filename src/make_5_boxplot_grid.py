import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path
import matplotlib.patches
import numpy as np

plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['grid.color'] = 'black'
plt.rcParams['lines.color'] = 'black'

def make_5_boxplot_grid(df, substrate, box_percent_string, whisker_percent_string):
    # make a figure with subplots
    _, ax = plt.subplots(11, 4, figsize=(20, 50))
    ax = ax.flatten()  # Flatten the axes for easier indexing

    subplot_index = 0
    for _, row in df.iterrows():
        reaction = row['Equation']

        if row['pathway'] == 'biomass_formation':
            reaction = 'Biomass formation'

        # generate column names based on the substrate
        flux_col = f'{substrate}_flux'
        lb_col = f'{substrate}_LB'
        ub_col = f'{substrate}_UB'

        # construct mfa dictionary
        mfa_stats = {
            'median': pd.to_numeric(row[flux_col], errors='coerce'), 
            'q1': pd.to_numeric(row[lb_col], errors='coerce'), 
            'q3': pd.to_numeric(row[ub_col], errors='coerce'), 
            'lower_whisker': pd.to_numeric(row[lb_col], errors='coerce'),
            'upper_whisker': pd.to_numeric(row[ub_col], errors='coerce'),
        }
    
        # generate column names based on the substrate
        pfba_flux_col = f'{substrate}_pfba_flux'
        pfba_box_lb_col = f'{substrate}_pfba_{box_percent_string}_LB'
        pfba_box_ub_col = f'{substrate}_pfba_{box_percent_string}_UB'
        pfba_whisker_lb_col = f'{substrate}_pfba_{whisker_percent_string}_LB'
        pfba_whisker_ub_col = f'{substrate}_pfba_{whisker_percent_string}_UB'

        # construct pfba dictionary
        pfba_stats = {
            'median': pd.to_numeric(row[pfba_flux_col], errors='coerce'),
            'q1': pd.to_numeric(row[pfba_box_lb_col], errors='coerce'),
            'q3': pd.to_numeric(row[pfba_box_ub_col], errors='coerce'),
            'lower_whisker': pd.to_numeric(row[pfba_whisker_lb_col], errors='coerce'),
            'upper_whisker': pd.to_numeric(row[pfba_whisker_ub_col], errors='coerce'),
        }

        # generate column names based on the substrate
        eflux2_flux_col = f'{substrate}_1_eflux2_flux'
        eflux2_box_lb_col = f'{substrate}_1_eflux2_{box_percent_string}_LB'
        eflux2_box_ub_col = f'{substrate}_1_eflux2_{box_percent_string}_UB'
        eflux2_whisker_lb_col = f'{substrate}_1_eflux2_{whisker_percent_string}_LB'
        eflux2_whisker_ub_col = f'{substrate}_1_eflux2_{whisker_percent_string}_UB'
        
        # construct eflux2 dictionary
        eflux2_1_stats = {
            'median': pd.to_numeric(row[eflux2_flux_col], errors='coerce'),
            'q1': pd.to_numeric(row[eflux2_box_lb_col], errors='coerce'),
            'q3': pd.to_numeric(row[eflux2_box_ub_col], errors='coerce'),
            'lower_whisker': pd.to_numeric(row[eflux2_whisker_lb_col], errors='coerce'),
            'upper_whisker': pd.to_numeric(row[eflux2_whisker_ub_col], errors='coerce'),
        }

        # generate column names based on the substrate
        eflux2_flux_col = f'{substrate}_2_eflux2_flux'
        eflux2_box_lb_col = f'{substrate}_2_eflux2_{box_percent_string}_LB'
        eflux2_box_ub_col = f'{substrate}_2_eflux2_{box_percent_string}_UB'
        eflux2_whisker_lb_col = f'{substrate}_2_eflux2_{whisker_percent_string}_LB'
        eflux2_whisker_ub_col = f'{substrate}_2_eflux2_{whisker_percent_string}_UB'

        # construct eflux2 dictionary
        eflux2_2_stats = {
            'median': pd.to_numeric(row[eflux2_flux_col], errors='coerce'),
            'q1': pd.to_numeric(row[eflux2_box_lb_col], errors='coerce'),
            'q3': pd.to_numeric(row[eflux2_box_ub_col], errors='coerce'),
            'lower_whisker': pd.to_numeric(row[eflux2_whisker_lb_col], errors='coerce'),
            'upper_whisker': pd.to_numeric(row[eflux2_whisker_ub_col], errors='coerce'),
        }

        # generate column names based on the substrate
        eflux2_flux_col = f'{substrate}_3_eflux2_flux'
        eflux2_box_lb_col = f'{substrate}_3_eflux2_{box_percent_string}_LB'
        eflux2_box_ub_col = f'{substrate}_3_eflux2_{box_percent_string}_UB'
        eflux2_whisker_lb_col = f'{substrate}_3_eflux2_{whisker_percent_string}_LB'
        eflux2_whisker_ub_col = f'{substrate}_3_eflux2_{whisker_percent_string}_UB'

        # construct eflux2 dictionary
        eflux2_3_stats = {
            'median': pd.to_numeric(row[eflux2_flux_col], errors='coerce'),
            'q1': pd.to_numeric(row[eflux2_box_lb_col], errors='coerce'),
            'q3': pd.to_numeric(row[eflux2_box_ub_col], errors='coerce'),
            'lower_whisker': pd.to_numeric(row[eflux2_whisker_lb_col], errors='coerce'),
            'upper_whisker': pd.to_numeric(row[eflux2_whisker_ub_col], errors='coerce'),
        }

        # Check if any values are NaN
        mfa_bounds_have_nan = any(pd.isna(list(mfa_stats.values())))
        pfba_bounds_have_nan = any(pd.isna(list(pfba_stats.values())))
        eflux2_bounds_have_nan = any(pd.isna(list(eflux2_1_stats.values())))
        
        if mfa_bounds_have_nan or pfba_bounds_have_nan or eflux2_bounds_have_nan:
            continue
        else:
            show_y_label = subplot_index % 4 == 0
            make_boxplot(ax[subplot_index], mfa_stats, pfba_stats, eflux2_1_stats, eflux2_2_stats, eflux2_3_stats, show_y_label, title=reaction,)
            subplot_index += 1

    # Show the plot
    plt.show()

def make_boxplot(ax, mfa_stats, pfba_stats, eflux2_1_stats, eflux2_2_stats, eflux2_3_stats, show_y_label, title=''):
    # Define boxprops for 13C-MFA (black lines, transparent fill)
    boxprops = dict(facecolor=(0, 0, 0, 0), linewidth=2, edgecolor='black') 
    medianprops = dict(linewidth=2, color='red')
    whiskerprops = dict(color='black') 
    capprops = dict(color='black') 

    for index, box_data in enumerate([mfa_stats, pfba_stats, eflux2_1_stats, eflux2_2_stats, eflux2_3_stats]):
        ax.bxp(
            [{
            'med': box_data['median'],
            'q1': box_data['q1'],
            'q3': box_data['q3'],
            'cilo': box_data['q1'],  # not actually used, but necessary for bxp to work
            'cihi': box_data['q3'],  # not actually used, but necessary for bxp to work
            'whislo': box_data['lower_whisker'],  # set to q1 because there are no whiskers
            'whishi': box_data['upper_whisker'],  # set to q3 because there are no whiskers
            'fliers': []  # no outliers
            }], 
            positions=[index + 1], 
            patch_artist=True, 
            widths=0.5, 
            boxprops=boxprops, 
            medianprops=medianprops, 
            whiskerprops=whiskerprops,
            capprops=capprops
        )

    # conditionally set y axis
    if show_y_label:
        ax.set_ylabel('Flux (normalized to 100 uptake flux)')

    # Set the x-axis labels and title
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xticklabels(['13C-MFA', 'pFBA', 'EF2 1', 'EF2 2', 'EF2 3'])
    ax.set_title(title)