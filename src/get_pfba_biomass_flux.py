import pandas as pd
import straindesign as sd

# this function returns a dataframe with rows for every GSM reaction
# each row has cols for reaction_id, reaction_name, full_reaction, flux, trans_LB, trans_UB, box_LB, box_UB, whisker_LB, and whisker_UB
def get_pfba_biomass_flux(model=None, substrate=''):
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
        biomass_reaction_id = 'biomass_oleic_acid'

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
    pfba_solution = sd.fba(model, constraints=f'{uptake_reaction_id} = {uptake_flux}', obj=biomass_reaction_id, obj_sense='maximize', pfba=1)
    biomass_flux = pfba_solution[biomass_reaction_id]

    return biomass_flux if substrate != 'oleic_acid' else 10 * biomass_flux