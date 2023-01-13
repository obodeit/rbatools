import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *

conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03']
#conditions = ['Hackett_N005', 'Hackett_N01', 'Hackett_N016', 'Hackett_N03']
#conditions = ['Hackett_P005', 'Hackett_P01', 'Hackett_P016', 'Hackett_P022']


Input_Data = pandas.read_csv(
    '../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nlim.csv', sep=';', decimal=',', index_col=0)

default_kapps_from_calibration=pandas.read_csv("../default_kapps_refactored_WF.csv",index_col=0)
specific_kapps_from_calibration=pandas.read_csv("../specific_kapps_refactored_WF.csv",index_col=0)
process_efficiencies_from_calibration=pandas.read_csv("../process_efficiencies_refactored_WF.csv",index_col=0)
pg_fractions_from_calibration=pandas.read_csv("../pg_fractions_refactored_WF.csv",index_col=0)
compartment_sizes_from_calibration = pandas.read_csv("../compartment_sizes_refactored_WF.csv",index_col=0)

growth_rates={condition:growth_rate_from_input(input=Input_Data, condition=condition) for condition in conditions}

regressed_specific_kapps=regression_on_specific_enzyme_efficiencies(Spec_Kapps=specific_kapps_from_calibration.copy(),min_kapp=360,max_kapp=360000000,conditions=conditions,growth_rates=growth_rates,impose_on_isoenzymes=False,monotonous_quadratic=True,fill_in_missing_conditions=True)
regressed_default_kapps=regression_on_default_enzyme_efficiencies(default_kapps=default_kapps_from_calibration.copy(),min_kapp=360,max_kapp=360000000,conditions=conditions,growth_rates=growth_rates,monotonous_quadratic=True)
regressed_process_efficiencies=regression_on_process_efficiencies(Process_efficiencies=process_efficiencies_from_calibration.copy(),min_efficiency=360,max_efficiency=360000000,conditions=conditions,growth_rates=growth_rates,monotonous_quadratic=True)

regressed_specific_kapps.to_csv("spec_kapp_reg.csv")
regressed_default_kapps.to_csv("def_kapp_reg.csv")
regressed_process_efficiencies.to_csv("proc_eff_reg.csv")

plot_specific_enzyme_efficiencies(point_calibration_kapps=specific_kapps_from_calibration,
                    regressed_kapps=regressed_specific_kapps,
                    conditions=conditions,
                    growth_rates=growth_rates,
                    filename="Specific_Kapp_Plots_sim.pdf")

Simulation = SessionRBA('../Yeast_iMM904_RBA_model')
Simulation.add_exchange_reactions()


functions_to_include_list=["Specific_Kapps","Default_Kapps","Compartment_Sizes","PG_Fractions","Process_Efficiencies"]
functions_to_include_list=["Specific_Kapps","Default_Kapps","Process_Efficiencies"]
functions_to_include_list=["Default_Kapps","Process_Efficiencies"]
functions_to_include_list=["Compartment_Sizes","PG_Fractions","Default_Kapps","Process_Efficiencies"]
functions_to_include_list=["PG_Fractions","Default_Kapps","Process_Efficiencies"]
functions_to_include_list=["Compartment_Sizes","PG_Fractions"]

list_of_rxns_to_impose=[]
simulation_results_Js_not_imposed=[]
for condition in conditions:
    flux_bounds_data=flux_bounds_from_input(input=Input_Data, condition=condition, specific_exchanges=[])
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
                                                  #sims_to_perform=["DefaultKapp","Prokaryotic","Eukaryotic"],
                                                  #sims_to_perform=["Prokaryotic","Eukaryotic"],
                                                  sims_to_perform=["Prokaryotic"],
                                                  feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                  try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                  print_output=True,
                                                  variability_analysis=["R_EX_glc__D_e","R_EX_etoh_e","R_EX_ac_e","R_EX_glyc_e","R_EX_acald_e","R_EX_lac__D_e","R_EX_succ_e","R_EX_o2_e"],
                                                  #variability_analysis=None,
                                                  mu_factor_for_variability=0.99,
                                                  apply_model=False,
                                                  functions_to_include_list=functions_to_include_list,
                                                  transporter_multiplier=1,
                                                  start_val=growth_rates[condition],
                                                  Mu_approx_precision=0.0001
                                                  )
    simulation_results_Js_not_imposed.append(simulation_result)

simulation_results_Js_imposed=[]
for condition in conditions:
    flux_bounds_data=flux_bounds_from_input(input=Input_Data, condition=condition, specific_exchanges=None,specific_directions=[])
    list_of_rxns_to_impose=["R_EX_glc__D_e","R_EX_etoh_e","R_EX_ac_e","R_EX_glyc_e","R_EX_acald_e","R_EX_succ_e"]
    if len(list_of_rxns_to_impose)>0:
        Exchanges_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in list_of_rxns_to_impose if i in list(flux_bounds_data["Reaction_ID"])}
    else:
        Exchanges_to_impose=None
    list_of_rxns_to_impose=[]
    Simulation.reload_model()
    print("Simulation fixed exchanges: {}".format(condition))
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
                                                  #sims_to_perform=["DefaultKapp","Prokaryotic","Eukaryotic"],
                                                  #sims_to_perform=["Prokaryotic","Eukaryotic"],
                                                  sims_to_perform=["Prokaryotic"],
                                                  feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                  try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                  print_output=True,
                                                  variability_analysis=None,
                                                  mu_factor_for_variability=0.99,
                                                  apply_model=False,
                                                  functions_to_include_list=functions_to_include_list,
                                                  transporter_multiplier=1,
                                                  start_val=growth_rates[condition],
                                                  Mu_approx_precision=0.0001
                                                  )
    simulation_results_Js_imposed.append(simulation_result)

Min_VA=pandas.DataFrame()
Max_VA=pandas.DataFrame()
for result in simulation_results_Js_not_imposed:
    try:
        for exchange in list(result['FeasibleRange_prok'].keys()):
            Min_VA.loc[exchange,result["Condition"]]=abs(result['FeasibleRange_prok'][exchange]["Min"])
            Max_VA.loc[exchange,result["Condition"]]=abs(result['FeasibleRange_prok'][exchange]["Max"])
            Min_VA.loc[exchange,"Real_{}".format(result["Condition"])]=result['FeasibleRange_prok'][exchange]["Min"]
            Max_VA.loc[exchange,"Real_{}".format(result["Condition"])]=result['FeasibleRange_prok'][exchange]["Max"]
    except:
        print("")
Min_VA.to_csv("../VAmin_mean.csv")
Max_VA.to_csv("../VAmax_mean.csv")

predicted_growth_rates=pandas.DataFrame(columns=conditions)
predicted_glucose_exchanges=pandas.DataFrame(columns=conditions)
predicted_oxygen_exchanges=pandas.DataFrame(columns=conditions)
predicted_ethanol_exchanges=pandas.DataFrame(columns=conditions)
predicted_acetate_exchanges=pandas.DataFrame(columns=conditions)
predicted_glycerol_exchanges=pandas.DataFrame(columns=conditions)
predicted_acald_exchanges=pandas.DataFrame(columns=conditions)
predicted_lactate_exchanges=pandas.DataFrame(columns=conditions)
predicted_succinate_exchanges=pandas.DataFrame(columns=conditions)

predicted_growth_rates.loc["Mean",:]=extract_predicted_growth_rates(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic')
predicted_glucose_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_glc__D')
predicted_oxygen_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_o2')
predicted_ethanol_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_etoh')
predicted_acetate_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_ac')
predicted_glycerol_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_glyc')
predicted_acald_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_acald')
predicted_lactate_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_lac__D')
predicted_succinate_exchanges.loc["Mean",:]=extract_predicted_exchange_fluxes(inputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_succ')

predicted_growth_rates.to_csv("../Mean_predicted_growth_rates.csv")
predicted_glucose_exchanges.to_csv("../Mean_predicted_glucose_exchanges.csv")
predicted_oxygen_exchanges.to_csv("../Mean_predicted_oxygen_exchanges.csv")
predicted_ethanol_exchanges.to_csv("../Mean_predicted_ethanol_exchanges.csv")
predicted_acetate_exchanges.to_csv("../Mean_predicted_acetate_exchanges.csv")
predicted_glycerol_exchanges.to_csv("../Mean_predicted_glycerol_exchanges.csv")
predicted_acald_exchanges.to_csv("../Mean_predicted_acald_exchanges.csv")
predicted_lactate_exchanges.to_csv("../Mean_predicted_lactate_exchanges.csv")
predicted_succinate_exchanges.to_csv("../Mean_predicted_succinate_exchanges.csv")

pred_prot=extract_proteomes_from_simulation_results(simulation_outputs=simulation_results_Js_not_imposed,type="Prokaryotic")
pred_prot_fixed=extract_proteomes_from_simulation_results(simulation_outputs=simulation_results_Js_imposed,type="Prokaryotic")
pred_prot.to_csv("../Predicted_proteomes.csv")
measured_proteomes=pandas.read_csv("../Corrected_calibration_proteomes.csv",index_col=0)

plot_predicted_fluxes(simulation_outputs=simulation_results_Js_not_imposed,types=["Prokaryotic"],input_definition=Input_Data)
plot_protein_protein_comparison(predicted_proteomes=pred_prot,measured_proteomes=measured_proteomes,conditions=conditions)
plot_predicted_fluxes(simulation_outputs=simulation_results_Js_not_imposed,types=["DefaultKapp","Prokaryotic","Eukaryotic"],input_definition=Input_Data)
plot_predicted_fluxes(simulation_outputs=simulation_results_Js_imposed,types=["DefaultKapp","Prokaryotic","Eukaryotic"],input_definition=Input_Data)
plot_protein_protein_comparison(predicted_proteomes=pred_prot,measured_proteomes=measured_proteomes,conditions=conditions)

#FD=get_flux_distribution(simulation_outputs=simulation_results_Js_not_imposed,result_object='Simulation_Results_Euk', run='Eukaryotic')
FD=get_flux_distribution(simulation_outputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic')
FD.to_csv("../FluxDistrubution.csv")
ExFD=get_exchange_flux_distribution(simulation_outputs=simulation_results_Js_not_imposed,result_object='Simulation_Results', run='Prokaryotic')
ExFD.to_csv("../ExchangeFluxDistrubution.csv")
try:
    out=pandas.DataFrame()
    for i in simulation_results_Js_not_imposed:
        out.loc["Mu",i["Condition"]]=i["Mu_euk"]
        for comp in i['Euk_CompSizes'].keys():
            out.loc[comp,i["Condition"]]=i['Euk_CompSizes'][comp]
    out.to_csv("../CompSizesPredicted_Euk.csv")
except:
    print("")