import pandas as pd
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox, TextArea)
from IPython.display import display


# This function takes in a dataframe of flux values, and makes a scatter plot of
#  predicted values vs. measured values
def flux_prediction_scatterplot(
        fluxes_df=None,
        substrate='',
        method='',
        range=None,
):
    # define column names
    prediction_column_name = f'{method} Flux'
    prediction_std_column_name = f'{method} Flux Std'
    
    # define plot area
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # create a list of pathways from the dataframe, this is sorted so the colors match on all plots
    pathway_list = sorted(list(set(fluxes_df['pathway'])))
    
    # loop over each pathway
    for pathway in pathway_list:
        # make a pathway specific filtered dataframe
        pathway_df = fluxes_df[fluxes_df['pathway'] == pathway]
        
        # make lists of measured and predicted fluxes
        measured_flux_list = list(pathway_df['13C_mfa_flux'])
        predicted_flux_list = list(pathway_df[prediction_column_name])
            
        # add the points for each pathway to the plot
        sc = ax.scatter(measured_flux_list, predicted_flux_list, label=pathway)

        # add labels to each point
        # for i, label in enumerate(list(pathway_df['reaction_ids'])):
        #     ax.annotate(label, (measured_flux_list[i], predicted_flux_list[i]))
        
        # this is incorrect (update later with custom upper and lower bound errors)
        # calculate horizontal error bars
        measured_upper_std_list = list(pathway_df['13C_mfa_flux_95%_upper_bound']-pathway_df['13C_mfa_flux'])
        measured_lower_std_list = list(pathway_df['13C_mfa_flux']-pathway_df['13C_mfa_flux_95%_lower_bound'])
        
        # if prediction standard deviations are included, then add vertical and horizontal error bars
        if prediction_std_column_name in pathway_df.columns:
            # make a list of predicted flux standard deviations
            predicted_std_list = list(pathway_df[prediction_std_column_name])
            
            ax.errorbar(
                measured_flux_list, 
                predicted_flux_list, 
                xerr=(measured_lower_std_list, measured_upper_std_list),
                yerr=[1.9*std for std in predicted_std_list], 
                ecolor="gray", 
                ls='none', 
                alpha=0.8
            )
            
        # if prediction standard deviations are NOT included, then only add horizontal error bars
        else: 
            ax.errorbar(
                measured_flux_list, 
                predicted_flux_list, 
                xerr=(measured_lower_std_list, measured_upper_std_list),
                ecolor="gray", 
                ls='none', 
                alpha=0.8
            )

    # add 45 degree dashed line
    x = np.linspace(*ax.get_xlim())
    ax.plot(x, x, ls="--", c=".3")

    # calculate r-squared values using all the values
    _, _, r_squared, _, _ = linregress(fluxes_df.loc[:, '13C_mfa_flux'], fluxes_df.loc[:, prediction_column_name])

    # calculate r-squared value without ATPM
    fluxes_df_no_atpm = fluxes_df.copy()
    fluxes_df_no_atpm = fluxes_df_no_atpm[fluxes_df_no_atpm['reaction_ids'] != 'ATPM']
    _, _, r_squared_no_atpm, _, _ = linregress(fluxes_df_no_atpm.loc[:, '13C_mfa_flux'], fluxes_df_no_atpm.loc[:, prediction_column_name])


    # add labels to the plot
    plt.title(r''+ r"$\bf{" + str(substrate) + r"\;" + str(method) + "}$"  + ': ' + f"$R^2$={r_squared_no_atpm:.2F}*", fontsize=24)
    # plt.title(r''+ r"$\bf{" + str(substrate) + r"\;" + str(method) + "}$"  + ': ' + f"$R^2$={r_squared_no_atpm:.2F} ({r_squared:.2F})", fontsize=24)
    plt.ylabel(f'{method} Flux', fontsize=22)
    plt.xlabel('13C-MFA Flux', fontsize=22)
    plt.legend(fontsize=10)
    
    # add styles to the plot
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # only show the range of fluxes specified
    if range:
        plt.xlim(range)
        plt.ylim(range)
    
    # save and show the plot
    plt.savefig(f'../figures/{substrate}_{method}_scatterplot.png', bbox_inches='tight')
    plt.show()