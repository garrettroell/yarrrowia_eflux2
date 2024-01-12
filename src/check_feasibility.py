# Define a function to check feasibility
def check_feasibility(gsm_bounds, mfa_bounds):
    
    gsm_lb = float(gsm_bounds[0])
    gsm_ub = float(gsm_bounds[1])
    mfa_lb = float(mfa_bounds[0])
    mfa_ub = float(mfa_bounds[1])

    if mfa_lb >= gsm_lb and mfa_ub <= gsm_ub:
        return 'fully feasible'
    elif mfa_ub < gsm_lb or mfa_lb > gsm_ub:
        return 'not feasible'
    else:
        return 'partially feasible'


def add_mfa_bound_feasibility_column(full_central_rxn_df, substrate):
    mfa_bound_feasibilities = []
    gsm_lb, gsm_ub = f'{substrate}_GSM_LB', f'{substrate}_GSM_UB'
    mfa_lb, mfa_ub = f'{substrate}_LB', f'{substrate}_UB'
    feasibility_column = f'{substrate}_mfa_bound_feasibility'

    # loop over flux rows
    for index, row in full_central_rxn_df.iterrows():
        # get the MFA and GSM flux bounds
        gsm_flux_bounds = [row[gsm_lb], row[gsm_ub]]
        mfa_flux_bounds = [row[mfa_lb], row[mfa_ub]]

        gsm_flux_bounds = pd.to_numeric(gsm_flux_bounds, errors='coerce')
        mfa_flux_bounds = pd.to_numeric(mfa_flux_bounds, errors='coerce')

        if any(pd.isna(gsm_flux_bounds)) or any(pd.isna(mfa_flux_bounds)):
            mfa_bound_feasibilities.append('')
            continue

        mfa_bound_feasibility = check_feasibility(gsm_flux_bounds, mfa_flux_bounds)
        mfa_bound_feasibilities.append(mfa_bound_feasibility)

    full_central_rxn_df[feasibility_column] = mfa_bound_feasibilities

    return full_central_rxn_df