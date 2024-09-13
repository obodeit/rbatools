import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *

conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03', 'Mean_01']
#conditions = ['Hackett_C01']

Input_Data = pandas.read_csv('../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nlim_Frick_fluxes.csv', sep=';', decimal=',', index_col=0)
#Input_Data = pandas.read_csv('../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nielsen.csv', sep=';', decimal=',', index_col=0)

default_kapps_from_calibration=pandas.read_csv("../default_kapps_refactored_WF.csv",index_col=0)
specific_kapps_from_calibration=pandas.read_csv("../specific_kapps_refactored_WF.csv",index_col=0)
process_efficiencies_from_calibration=pandas.read_csv("../process_efficiencies_refactored_WF.csv",index_col=0)
pg_fractions_from_calibration=pandas.read_csv("../pg_fractions_refactored_WF.csv",index_col=0)
compartment_sizes_from_calibration = pandas.read_csv("../compartment_sizes_refactored_WF.csv",index_col=0)

growth_rates={condition:growth_rate_from_input(input=Input_Data, condition=condition) for condition in conditions}

Simulation = SessionRBA('../Yeast_models/Yeast_iMM904_RBA_model')
#Simulation = SessionRBA('../Yeast_iMM904_RBA_model_no_BMcompo_targets')
Simulation.add_exchange_reactions()

functions_to_include_list=[]

list_of_rxns_to_impose=[]
simulation_results_Js_not_imposed=[]
for condition in conditions:
    flux_bounds_data=flux_bounds_from_input(input=Input_Data,rba_session=Simulation, condition=condition, specific_exchanges=[],specific_directions=[],also_consider_iso_enzmes=False)
    if len(list_of_rxns_to_impose)>0:
        Exchanges_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in list_of_rxns_to_impose if i in list(flux_bounds_data["Reaction_ID"])}
    else:
        Exchanges_to_impose=None
    Simulation.reload_model()
    print("Simulation: {}".format(condition))
    simulation_result=perform_simulations(condition=condition,
                                                  definition_file=Input_Data,
                                                  rba_session=Simulation,
                                                  compartment_sizes=compartment_sizes_from_calibration,
                                                  pg_fractions=pg_fractions_from_calibration,
                                                  process_efficiencies=process_efficiencies_from_calibration,
                                                  #process_efficiencies=regressed_process_efficiencies,
                                                  Default_Kapps=default_kapps_from_calibration,
                                                  #Default_Kapps=regressed_default_kapps,
                                                  Specific_Kapps=specific_kapps_from_calibration,
                                                  #Specific_Kapps=regressed_specific_kapps,
                                                  Exchanges_to_impose=Exchanges_to_impose,
                                                  #sims_to_perform=["Prokaryotic","Eukaryotic",'Eukaryotic_fixed_sizes','Fixed_PG_Eukaryotic','Fixed_PG_Eukaryotic_fixed_sizes'],
                                                  sims_to_perform=["Prokaryotic"],
                                                  #sims_to_perform=["Prokaryotic"],
                                                  #feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                  feasible_stati=["optimal","feasible"],
                                                  try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                                  print_output=True,
                                                  #variability_analysis=["R_EX_glc__D_e","R_EX_etoh_e","R_EX_ac_e","R_EX_glyc_e","R_EX_acald_e","R_EX_lac__D_e","R_EX_succ_e","R_EX_o2_e"],
                                                  variability_analysis=None,
                                                  mu_factor_for_variability=0.99,
                                                  apply_model=False,
                                                  functions_to_include_list=functions_to_include_list,
                                                  transporter_multiplier=1,
                                                  start_val=0,
                                                  #start_val=growth_rates[condition],
                                                  Mu_approx_precision= 0.00001
                                                  )
    Simulation.Results=simulation_result["Simulation_Results"]
    Simulation.write_results(session_name=condition)
    Simulation.SimulationData.export_proteo_map(type='genes')
    Simulation.SimulationData.export_escher_map(type='fluxes')
    Simulation.SimulationData.export_escher_map(type='investment')
