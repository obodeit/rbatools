import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *
import time
from multiprocessing import Pool , cpu_count
from joblib import Parallel, delayed


def simulate_run(input_dict):
    compartment_sizes=pandas.DataFrame()
    pg_fractions=pandas.DataFrame()
    process_efficiencies=pandas.DataFrame()
    specific_kapps=pandas.DataFrame()
    default_kapps=pandas.DataFrame()
    Simulation = SessionRBA(input_dict["xml_dir"])
    Simulation.add_exchange_reactions()

    for condition in input_dict["conditions"]:
        default_kapps_bootstrapped=pandas.read_csv("{}/DefaultKapps_{}.csv".format(input_dict["folder"],condition),index_col=0)
        specific_kapps_bootstrapped=pandas.read_csv("{}/SpecificKapps_{}.csv".format(input_dict["folder"],condition),index_col=0)
        specific_kapps_directions_bootstrapped=pandas.read_csv("{}/SpecificKappDirections_{}.csv".format(input_dict["folder"],condition),index_col=0)
        compartment_sizes_bootstrapped=pandas.read_csv("{}/CompartmentSizes_{}.csv".format(input_dict["folder"],condition),index_col=0)
        pg_fractions_bootstrapped=pandas.read_csv("{}/PgFractions_{}.csv".format(input_dict["folder"],condition),index_col=0)
        process_efficiencies_bootstrapped=pandas.read_csv("{}/ProcessEfficiencies_{}.csv".format(input_dict["folder"],condition),index_col=0)
        for i in compartment_sizes_bootstrapped.index:
            compartment_sizes.loc[i,condition]=compartment_sizes_bootstrapped.loc[i,input_dict["run"]]
        for i in pg_fractions_bootstrapped.index:
            pg_fractions.loc[i,condition]=pg_fractions_bootstrapped.loc[i,input_dict["run"]]
        for i in process_efficiencies_bootstrapped.index:
            process_efficiencies.loc[i,condition]=process_efficiencies_bootstrapped.loc[i,input_dict["run"]]
        for i in default_kapps_bootstrapped.index:
            default_kapps.loc[condition,"Default Kapp"]=default_kapps_bootstrapped.loc[i,input_dict["run"]]
            default_kapps.loc[condition,"ID"]=condition
        for i in specific_kapps_bootstrapped.index:
            if i in specific_kapps_directions_bootstrapped.index:
                if not i.endswith("_enzyme"):
                    enzyme=i+"_enzyme"
                else:
                    enzyme=i
                specific_kapps.loc[enzyme,"Enzyme_ID"]=enzyme
                specific_kapps.loc[enzyme,condition]=specific_kapps_bootstrapped.loc[i,input_dict["run"]]
                specific_kapps.loc[enzyme,"{}_Flux".format(condition)]=specific_kapps_directions_bootstrapped.loc[i,input_dict["run"]]

    if len(input_dict["list_of_rxns_to_impose"])>0:
        print("Fixed {}".format(input_dict["run"]))
        simulation_outputs_fixed=[]
        for condition in input_dict["conditions"]:
            flux_bounds_data=flux_bounds_from_input(input=input_dict["Input_Data"],rba_session=Simulation, condition=condition, specific_exchanges=None,specific_directions=[],also_consider_iso_enzmes=False)
            Exchanges_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in input_dict["list_of_rxns_to_impose"] if i in list(flux_bounds_data["Reaction_ID"])}
            Simulation.reload_model()
            simulation_result=perform_simulations(condition=condition,
                                                        definition_file=input_dict["Input_Data"],
                                                        rba_session=Simulation,
                                                        compartment_sizes=compartment_sizes,
                                                        pg_fractions=pg_fractions,
                                                        process_efficiencies=process_efficiencies,
                                                        Default_Kapps=default_kapps,
                                                        Specific_Kapps=specific_kapps,
                                                        Exchanges_to_impose=Exchanges_to_impose,
                                                        sims_to_perform=["Prokaryotic"],
                                                        feasible_stati=["optimal","feasible"],
                                                        try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                                        print_output=False,
                                                        variability_analysis=None,
                                                        mu_factor_for_variability=0.99,
                                                        apply_model=False,
                                                        functions_to_include_list=[],
                                                        transporter_multiplier=1,
                                                        start_val=input_dict["growth_rates"][condition]
                                                        )
            simulation_outputs_fixed.append(simulation_result)
        predicted_growth_rates_fixed_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_growth_rates(inputs=simulation_outputs_fixed,result_object='Simulation_Results', run='Prokaryotic'))))
    else:
        predicted_growth_rates_fixed_exchanges=dict(zip(input_dict["conditions"],[numpy.nan]*len(input_dict["conditions"])))

    print("Free {}".format(input_dict["run"]))
    simulation_outputs=[]
    for condition in input_dict["conditions"]:
        Simulation.reload_model()
        simulation_result=perform_simulations(condition=condition,
                                                    definition_file=input_dict["Input_Data"],
                                                    rba_session=Simulation,
                                                    compartment_sizes=compartment_sizes,
                                                    pg_fractions=pg_fractions,
                                                    process_efficiencies=process_efficiencies,
                                                    Default_Kapps=default_kapps,
                                                    Specific_Kapps=specific_kapps,
                                                    Exchanges_to_impose=None,
                                                    sims_to_perform=["Prokaryotic"],
                                                    feasible_stati=["optimal","feasible"],
                                                    try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                                    print_output=False,
                                                    variability_analysis=input_dict["variability_analysis"],
                                                    mu_factor_for_variability=0.99,
                                                    apply_model=False,
                                                    functions_to_include_list=[],
                                                    transporter_multiplier=1,
                                                    start_val=input_dict["growth_rates"][condition])
        simulation_outputs.append(simulation_result)

    predicted_growth_rates=dict(zip(input_dict["conditions"],list(extract_predicted_growth_rates(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic'))))
    predicted_glucose_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_glc__D'))))
    predicted_oxygen_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_o2'))))
    predicted_ethanol_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_etoh'))))
    predicted_acetate_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_ac'))))
    predicted_glycerol_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_glyc'))))
    predicted_acald_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_acald'))))
    predicted_lactate_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_lac__D'))))
    predicted_succinate_exchanges=dict(zip(input_dict["conditions"],list(extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_succ'))))
    out={input_dict["run"]:{"predicted_growth_rates_fixed_exchanges":predicted_growth_rates_fixed_exchanges,
                            "predicted_growth_rates":predicted_growth_rates,
                            "predicted_glucose_exchanges":predicted_glucose_exchanges,
                            "predicted_oxygen_exchanges":predicted_oxygen_exchanges,
                            "predicted_ethanol_exchanges":predicted_ethanol_exchanges,
                            "predicted_acetate_exchanges":predicted_acetate_exchanges,
                            "predicted_glycerol_exchanges":predicted_glycerol_exchanges,
                            "predicted_acald_exchanges":predicted_acald_exchanges,
                            "predicted_lactate_exchanges":predicted_lactate_exchanges,
                            "predicted_succinate_exchanges":predicted_succinate_exchanges,
                            "variability_analysis":{}}}
    if input_dict["variability_analysis"] is not None:
        VA_dict={}
        for result in simulation_outputs:
            VA_dict[result["Condition"]]=result['FeasibleRange_prok']
        out[input_dict["run"]]["variability_analysis"]=VA_dict
    return(out)


def main(folder,xml_dir,conditions,growth_rates,VA_runs=None,VA_variables=None,list_of_rxns_to_impose=[],n_parallel_processes=None):
    Input_Data = pandas.read_csv('../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nlim.csv', sep=';', decimal=',', index_col=0)

    default_kapps_01=pandas.read_csv("{}/DefaultKapps_Hackett_C01.csv".format(folder),index_col=0)
    runs=list(default_kapps_01.columns)
    #runs=[list(default_kapps_01.columns)[0],list(default_kapps_01.columns)[1],list(default_kapps_01.columns)[2],list(default_kapps_01.columns)[3]]
    input_dicts=[]
    for run in runs:
        dict_to_add={"conditions":conditions,
                     "growth_rates":growth_rates,
                     "folder":folder,
                     "run":run,
                     "xml_dir":xml_dir,
                     "list_of_rxns_to_impose":list_of_rxns_to_impose,
                     "variability_analysis":VA_variables}
        if VA_runs is not None:
            if run not in VA_runs:
                dict_to_add["variability_analysis"]=None
        dict_to_add["Input_Data"]=Input_Data
        input_dicts.append(dict_to_add)

    if n_parallel_processes!=1:
        if n_parallel_processes is None:
            num_cores=cpu_count()
            n_jobs=min(num_cores,len(input_dicts))
        else:
            n_jobs=n_parallel_processes

        pool=Pool(n_jobs)

        #simulation_results_pool=pool.map_async(simulate_run,input_dicts)
        #simulation_results=simulation_results_pool.get()

        simulation_results=pool.imap_unordered(simulate_run,input_dicts)

        #simulation_results=Parallel(n_jobs=n_jobs)(delayed(simulate_run)(input_dict) for input_dict in input_dicts)
    else:
        simulation_results=[]
        for input_dict in input_dicts:
            simulation_results.append(simulate_run(input_dict))

    predicted_growth_rates=pandas.DataFrame(columns=conditions)
    predicted_growth_rates_fixed_exchanges=pandas.DataFrame(columns=conditions)
    predicted_glucose_exchanges=pandas.DataFrame(columns=conditions)
    predicted_oxygen_exchanges=pandas.DataFrame(columns=conditions)
    predicted_ethanol_exchanges=pandas.DataFrame(columns=conditions)
    predicted_acetate_exchanges=pandas.DataFrame(columns=conditions)
    predicted_glycerol_exchanges=pandas.DataFrame(columns=conditions)
    predicted_acald_exchanges=pandas.DataFrame(columns=conditions)
    predicted_lactate_exchanges=pandas.DataFrame(columns=conditions)
    predicted_succinate_exchanges=pandas.DataFrame(columns=conditions)

    for result in simulation_results:
        run=list(result.keys())[0]
        for condition in conditions:
            predicted_growth_rates_fixed_exchanges.loc[run,condition]=result[run]["predicted_growth_rates_fixed_exchanges"][condition]
            predicted_growth_rates.loc[run,condition]=result[run]["predicted_growth_rates"][condition]
            predicted_glucose_exchanges.loc[run,condition]=result[run]["predicted_glucose_exchanges"][condition]
            predicted_oxygen_exchanges.loc[run,condition]=result[run]["predicted_oxygen_exchanges"][condition]
            predicted_ethanol_exchanges.loc[run,condition]=result[run]["predicted_ethanol_exchanges"][condition]
            predicted_acetate_exchanges.loc[run,condition]=result[run]["predicted_acetate_exchanges"][condition]
            predicted_glycerol_exchanges.loc[run,condition]=result[run]["predicted_glycerol_exchanges"][condition]
            predicted_acald_exchanges.loc[run,condition]=result[run]["predicted_acald_exchanges"][condition]
            predicted_lactate_exchanges.loc[run,condition]=result[run]["predicted_lactate_exchanges"][condition]
            predicted_succinate_exchanges.loc[run,condition]=result[run]["predicted_succinate_exchanges"][condition]

        if len(list(result[run]["variability_analysis"].keys()))>0:
            VA_min_DF=pandas.DataFrame(columns=conditions)
            VA_max_DF=pandas.DataFrame(columns=conditions)
            for condition in conditions:
                VA_result=result[run]["variability_analysis"][condition]
                for var in VA_result.keys():
                    VA_min_DF.loc[var,condition]=VA_result[var]["Min"]
                    VA_max_DF.loc[var,condition]=VA_result[var]["Max"]
            VA_min_DF.to_csv("{}/VAmin_{}.csv".format(folder,run))
            VA_max_DF.to_csv("{}/VAmax_{}.csv".format(folder,run))

    predicted_growth_rates_fixed_exchanges.to_csv("{}/Bootstrapped_predicted_growth_rates_fixed_exchanges.csv".format(folder))
    predicted_growth_rates.to_csv("{}/Bootstrapped_predicted_growth_rates.csv".format(folder))
    predicted_glucose_exchanges.to_csv("{}/Bootstrapped_predicted_glucose_exchanges.csv".format(folder))
    predicted_oxygen_exchanges.to_csv("{}/Bootstrapped_predicted_oxygen_exchanges.csv".format(folder))
    predicted_ethanol_exchanges.to_csv("{}/Bootstrapped_predicted_ethanol_exchanges.csv".format(folder))
    predicted_acetate_exchanges.to_csv("{}/Bootstrapped_predicted_acetate_exchanges.csv".format(folder))
    predicted_glycerol_exchanges.to_csv("{}/Bootstrapped_predicted_glycerol_exchanges.csv".format(folder))
    predicted_acald_exchanges.to_csv("{}/Bootstrapped_predicted_acald_exchanges.csv".format(folder))
    predicted_lactate_exchanges.to_csv("{}/Bootstrapped_predicted_lactate_exchanges.csv".format(folder))
    predicted_succinate_exchanges.to_csv("{}/Bootstrapped_predicted_succinate_exchanges.csv".format(folder))


if __name__ == "__main__":
    t0=time.time()
    main(folder="../Bootstrapping_Results",
         xml_dir="../Yeast_iMM904_RBA_model",
         conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03'],
         growth_rates={'Hackett_C005':0.05, 'Hackett_C01':0.105, 'Hackett_C016':0.154, 'Hackett_C022':0.214, 'Hackett_C03':0.294},
         VA_runs=["mean_noNoise","origmean","mean_run_param"],
         VA_variables=["R_EX_glc__D_e","R_EX_etoh_e","R_EX_ac_e","R_EX_glyc_e","R_EX_acald_e","R_EX_lac__D_e","R_EX_succ_e","R_EX_o2_e"],
         list_of_rxns_to_impose=["R_EX_glc__D_e","R_EX_etoh_e","R_EX_ac_e","R_EX_glyc_e","R_EX_acald_e","R_EX_succ_e"],
         #list_of_rxns_to_impose=[],
         n_parallel_processes=4)
    print("Total time: {}".format(time.time()-t0))
