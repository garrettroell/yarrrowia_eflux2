# import numpy as np
# from optlang.symbolics import add
# import pandas as pd
# # import cplex

# #################################################################
# ########Functions needed for all prediction methods##############
# #################################################################
# """
#     Creates dictionary of isozymes by parsing GPR:
#     Parse GPR into a dict containing isozymes (separated by 'or'). Each isozyme has a set of subunits (separated by 'and') 'and' and 'or' can occur at the same time, or can occur by itself.
    
#         Parameters
#         ----------
#         model : cobrapy model.
        
        
#         Returns
#         -------
#         gpr_dict: dictionary with isozymes.
        
# """
# #Code only works for GPRs written in disjunctive normal form (DNF). Majority of models have them in DNF but there are some exceptions. 

# def create_gprdict(model):   
#     gpr_dict = dict()
#     for rxn in model.reactions:
#         if rxn.gene_reaction_rule:
#             temp = set()
#             for x in [x.strip('() ') for x in rxn.gene_reaction_rule.split(' or ')]:
#                 temp.add(frozenset(y.strip('() ') for y in x.split(' and ')))
#             gpr_dict[rxn.id] = temp
#     return gpr_dict

# """
#     Calculates bound value based on transcriptomics data for reactions in gene reaction rule
    
#     NOTE: 
#     If a reaction R1 has the GPR of 'A and B', it would be parsed to { {A, B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B) ] ) = min(A, B).
#     If a reaction R1 has the GPR of 'A or B', it would be parsed to { {A}, {B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A), min(B) ] ) = sum( [A, B] ).
#     If a reaction R1 has the GPR of '(A and B) or (C and D)', it would be parsed to { {A, B}, {C, D} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B), min(C, D) ] ).
    
#         Parameters
#         ----------
#         model : cobrapy model.
#         Transcriptomics : pandas dataframe with transcriptomics data.Data frame has gene identifiers as index and just one column with transcript values.  
#         rxn : cobrapy model reaction
        
        
#         Returns
#         -------
#         transscript bound value: float.
# """

# def transcript_value_for_rxn(model, transcriptomics_df, rxn):
#     final_transcript_value = 0
#     gene_ids = []
#     for parallel_gene in create_gprdict(model)[rxn.id]:
#         transcript_values = []
#         for gene in parallel_gene:
#             if gene in transcriptomics_df.index:
#                 transcript_values.append(transcriptomics_df.loc[gene].to_numpy()[0])
# #                 print(transcriptomics_df.loc[gene].to_numpy()[0])
#             else:
#                 transcript_values.append(np.inf)
#             min_transcript_val = np.min(transcript_values)
#         final_transcript_value = final_transcript_value + min_transcript_val
# #         if final_transcript_value==newinfbound:
# #             display(rxn.id)
# #             gene_ids.append(rxn.id)
#     return final_transcript_value

# #############################################
# ##################EFLUX2#######################
# ############################################
# """
#     Provides EFLUX2 predictions as explained in Machado et. al (2014) 
    
#         Parameters
#         ----------
#         model : cobrapy model.
#         Transcriptomics : pandas dataframe with transcriptomics data.Data frame has gene identifiers as index and just one column with transcript values.  
        
#         Returns
#         -------
#         eflux2_sol as output from eflux2_model.optimize().
        
# """
# def EFlux2(model, Transcriptomics):
#     # copy model and set tolerance
#     eflux2_model = model.copy()
#     eflux2_model.tolerance = 1e-9
    
#     # set the flux bounds for each reaction using the transcriptomics data    
#     for rxn in eflux2_model.reactions:
#         if 'EX_' not in str(rxn.id):
#             if rxn.gene_reaction_rule:
#                 if rxn.lower_bound < 0.0:
#                     rxn.lower_bound = -transcript_value_for_rxn(model, Transcriptomics, rxn)
#                 else:
#                     pass
#                 if rxn.upper_bound > 0.0:
#                     rxn.upper_bound = transcript_value_for_rxn(model, Transcriptomics, rxn)
#                 else:
#                     pass
#             else:
#                 """When there is no GPR, the arbitrary bounds are removed. 
#                 Common arbitrary bound value of 1000 for E.coli, might be different depending on the model, e.g., 99999.0 for iMM904 yeast model in BiGG"""
#                 if rxn.lower_bound <= -1000:
#                     rxn.lower_bound = -np.Inf
#                 if rxn.upper_bound >= 1000:
#                     rxn.upper_bound = np.Inf 
    
#     # solve FBA problem with transcriptomic bounds
#     fba_solution = eflux2_model.optimize()
#     print('FBA status', fba_solution.status)
#     print('FBA solution', fba_solution.objective_value)
    
#     display(eflux2_model.summary(solution=fba_solution))

#     # constrain the biomass to the optimal value
#     for r in eflux2_model.reactions:
#         if r.objective_coefficient:
#             r.lower_bound = fba_solution.objective_value
#             r.upper_bound = fba_solution.objective_value

#     # Inspect media
# #     display(eflux2_model.medium)
    
#     # Minimize the sum of squared flux values
#     """Note: Because of quadratic objective still have to use cplex objective formulation.
#     Optlang does not support quadratic type of constraints and objectives yet."""
#     eflux2_model.objective = eflux2_model.problem.Objective(add([rxn.flux_expression**2 for rxn in eflux2_model.reactions]), direction='min')
    
#     # solve the minimization of squared fluxes problem
#     EFlux2_solution = eflux2_model.optimize()
    
#     #display(eflux2_model.summary())
#     display(eflux2_model.summary(solution=EFlux2_solution))

    
#     print('E-Flux2 status', EFlux2_solution.status)
#     print('E-Flux2 solution', EFlux2_solution.objective_value)
#     print()
        
#     return EFlux2_solution

# # this function takes in a model, transcript dataframe, a substrate string, 
# # and a substrate uptake rate. It returns a cobrapy solution object
# def get_eflux2_solution(model, transcript_df, substrate, substrate_uptake_rate=100):  
#     with model:
#         # set the biomass reaction and media composition based on carbon source
#         medium = model.medium
#         if substrate == 'phenol':
#             # block glucose biomass reaction and the default CarveMe biomass reactions
#             model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
#             model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
#             model.reactions.get_by_id('Growth').upper_bound = 0
#             model.reactions.get_by_id('Growth').lower_bound = 0
            
#             # make maximizing the phenol biomass reaction the objective function
#             model.objective = 'Growth_Phenol'
            
#             # first, set all media components to surplus levels
#             medium = {key:1000 for (key,value) in model.medium.items()}
            
#             # set the phenol uptake rate the specified value
#             medium["EX_phenol_e"] = substrate_uptake_rate
            
#             # remove all non-phenol carbon sources
#             medium["EX_glc__D_e"] = 0
#             medium['EX_guaiacol_e'] = 0
#             medium['EX_vanlt_e'] = 0
            
#             # allow for excess phenol
#             medium['EX_phenol_e'] = 1000

#         elif substrate == 'glucose':
#             # block phenol biomass reaction and the default CarveMe biomass reactions
#             model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
#             model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
#             model.reactions.get_by_id('Growth').upper_bound = 0
#             model.reactions.get_by_id('Growth').lower_bound = 0
#             # model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
        
#             # make maximizing the glucose biomass reaction the objective function
#             model.objective = 'Growth_Glucose'
            
#             # first, set all media components to surplus levels
#             medium = {key:1000 for (key,value) in model.medium.items()}
            
#             #remove all non-glucose carbon sources:
#             medium["EX_phenol_e"] = 0
#             medium['EX_guaiacol_e'] = 0
#             medium['EX_vanlt_e'] = 0
            
#             # allow for excess glucose
#             medium['EX_glc__D_e'] = 1000
#         else:
#             print('Unknown substrate: Please choose among phenol and glucose')
            
#         # update model medium and solve the E-Flux2 problem
#         model.medium = medium
#         EFlux2_solution = EFlux2(model, transcript_df)
        
#     return EFlux2_solution

    