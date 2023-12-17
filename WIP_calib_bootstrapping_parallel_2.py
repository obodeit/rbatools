#main()
    ###prep
    ###sampling
    ###generate list of input dicts over runs
        #[{run:{...}}]
    ###calibrate each run
    # for each element in input dict:
        # run_full_calibration
    ###output

import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *
import time
import warnings
from multiprocessing import Pool , cpu_count

def run_calibration_over_conditions_for_bootstrapping_run(input_dict_for_run):
    run_id=list(input_dict_for_run.keys())[0]
    print("Bootstrapping run: {}".format(run_id))
    input_dict=input_dict_for_run[run_id]
    growth_rates={condition:growth_rate_from_input(input=input_dict["definition_file"], condition=condition) for condition in input_dict["conditions"]}
    input_dict["preliminary_run"]=True

    input_dicts_for_conditions=[]
    for condition in input_dict["conditions"]:
        input_dicts_for_conditions.append(input_dict.copy())
        input_dicts_for_conditions[-1]["condition"]=condition
    
    calib_dicts_1=[calibration(i,print_outputs= False) for i in input_dicts_for_conditions]

    calibration_results_1=[]
    for condition in input_dict["conditions"]:
        for i in calib_dicts_1:
            if condition in list(i.keys()):
                calibration_results_1.append(i[condition])

    compartment_sizes_from_calibration_1=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=calibration_results_1)
    pg_fractions_from_calibration_1=extract_pg_fractions_from_calibration_outputs(calibration_outputs=calibration_results_1)
    regressed_compartment_sizes_1=regression_on_compartment_sizes(Comp_sizes=compartment_sizes_from_calibration_1,conditions=input_dict["conditions"],growth_rates=growth_rates,monotonous_quadratic=False)
    regressed_pg_fractions_1=regression_on_pg_fractions(PG_sizes=pg_fractions_from_calibration_1,conditions=input_dict["conditions"],growth_rates=growth_rates,monotonous_quadratic=False)

    for i in input_dicts_for_conditions:
        i.update({"Compartment_sizes":regressed_compartment_sizes_1,
                  "PG_fractions":regressed_pg_fractions_1,
                  "preliminary_run":False})

    calib_dicts_2=[]
    for i in input_dicts_for_conditions:
        calib_dicts_2.append(calibration(i,print_outputs= True))

    out={}
    for condition in input_dict["conditions"]:
        for calib_results in calib_dicts_2:
            if condition in list(calib_results.keys()):
                [condition]
                proteome_corrected={i:calib_results[condition]["Proteome"].loc[i,"copy_number"] for i in calib_results[condition]["Proteome"].index}
                proteome_input={i:input_dict["proteome"].loc[i,condition] for i in input_dict["proteome"].index}
                densities={i:round(calib_results[condition]["Densities_PGs"].loc[i,"Density"],6) for i in calib_results[condition]["Densities_PGs"].index}
                pgs={i:round(calib_results[condition]["Densities_PGs"].loc[i,"PG_fraction"],6) for i in calib_results[condition]["Densities_PGs"].index}
                spec_kapps={i:{"Kapp":calib_results[condition]["Specific_Kapps"].loc[calib_results[condition]["Specific_Kapps"]["Enzyme_ID"]==i,"Kapp"].values[0],"Flux":calib_results[condition]["Specific_Kapps"].loc[calib_results[condition]["Specific_Kapps"]["Enzyme_ID"]==i,"Flux"].values[0]} for i in calib_results[condition]["Specific_Kapps"]["Enzyme_ID"]}
                process_effs={i:calib_results[condition]["Process_Efficiencies"].loc[calib_results[condition]["Process_Efficiencies"]['Process']==i,"Value"].values[0] for i in calib_results[condition]["Process_Efficiencies"]['Process']}
                def_kapp={"default_efficiency":calib_results[condition]["Default_Kapps"]['default_efficiency']}

                out[condition]={"proteome_input":proteome_input,
                                "proteome_corrected":proteome_corrected,
                                "densities":densities,
                                "pgs":pgs,
                                "spec_kapps":spec_kapps,
                                "process_effs":process_effs,
                                "def_kapp":def_kapp}

    return({run_id:out})

def calibration(input_dict,print_outputs=True):
    Simulation = SessionRBA(input_dict["xml_dir"])
    Simulation.add_exchange_reactions()
    calib_results = calibration_workflow(proteome=input_dict["proteome"],
                                         condition=input_dict["condition"],
                                         gene_ID_column='Gene',
                                         definition_file=input_dict["definition_file"],
                                         rba_session=Simulation,
                                         process_efficiency_estimation_input=input_dict["process_efficiency_estimation_input"],
                                         spec_kapps=None,
                                         process_efficiencies=None,
                                         Compartment_sizes=input_dict["Compartment_sizes"],
                                         PG_fractions=input_dict["PG_fractions"],
                                         transporter_multiplier=1,
                                         prelim_run=input_dict["preliminary_run"],
                                         #Mu_approx_precision=  0.000000001,
                                         Mu_approx_precision= 0.000001,
                                         feasible_stati=["optimal","feasible"],
                                         #feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                         min_kapp=None,
                                         print_outputs=print_outputs,
                                         use_mean_enzyme_composition_for_calibration=False,
                                         global_protein_scaling_coeff=1000/6.022e23,
                                         max_kapp_threshold=10000000000)
    return({input_dict["condition"]:calib_results})

def generate_input_proteome(fold_changes,
                            absolute_data,
                            matching_column_in_absolute_data,
                            matching_column_in_fold_change_data,
                            conditions_in_fold_change_data_to_restore,
                            ID_column_in_restored_data,
                            full_annotations):
    restored_relative_Data = infer_copy_numbers_from_reference_copy_numbers(fold_changes=fold_changes,
                                                                        absolute_data=absolute_data,
                                                                        matching_column_in_fold_change_data=matching_column_in_fold_change_data,
                                                                        matching_column_in_absolute_data=matching_column_in_absolute_data,
                                                                        conditions_in_fold_change_data_to_restore=conditions_in_fold_change_data_to_restore)
    restored_relative_Data_with_annotations= add_annotations_to_proteome(input=restored_relative_Data, ID_column=ID_column_in_restored_data, annotations=full_annotations)
    return(restored_relative_Data_with_annotations)

def build_full_annotations(rba_session,
                           model_ribosomal_processes,
                           external_ribosome_annotations,
                           Compartment_Annotations_external,
                           Absolute_Proteome,
                           protein_ID_column,
                           Uniprot,
                           Relative_Proteome):

    ribosomal_proteins = find_ribosomal_proteins(rba_session=rba_session, model_processes=model_ribosomal_processes, external_annotations=external_ribosome_annotations)
    model_protein_compartment_map = build_model_compartment_map(rba_session=rba_session)
    Compartment_Annotations = build_compartment_annotations(Compartment_Annotations_external=Compartment_Annotations_external, model_protein_compartment_map=model_protein_compartment_map)
    annotations_Absolute = build_dataset_annotations(input=Absolute_Proteome, ID_column=protein_ID_column, Uniprot=Uniprot, Compartment_Annotations=Compartment_Annotations, model_protein_compartment_map=model_protein_compartment_map,ribosomal_proteins=ribosomal_proteins)
    annotations_Relative = build_dataset_annotations(input=Relative_Proteome, ID_column=protein_ID_column, Uniprot=Uniprot,Compartment_Annotations=Compartment_Annotations, model_protein_compartment_map=model_protein_compartment_map,ribosomal_proteins=ribosomal_proteins)
    full_annotations = build_full_annotations_from_dataset_annotations(annotations_list=[annotations_Absolute, annotations_Relative])
    return(full_annotations)

def main(conditions,number_samples,n_parallel_processes=None,chunk_size=1,import_initial_parameters=False):
    ###prep
    Input_Data = pandas.read_csv('../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nlim_Frick_fluxes.csv', sep=';', decimal=',', index_col=0)
    Process_Efficiency_Estimation_Input = pandas.read_csv('../DataSetsYeastRBACalibration/Process_Efficiency_Estimation_Input.csv', sep=';', decimal=',')
    Uniprot = pandas.read_csv('../Yeast_iMM904_RBA_model/uniprot.csv', sep='\t')
    Compartment_Annotations_external = pandas.read_csv('../DataSetsYeastRBACalibration/Manually_curated_Protein_Locations_for_Calibration.csv', index_col=None, sep=';')
    Ribosomal_Proteins_Uniprot = pandas.read_csv('../DataSetsYeastRBACalibration/uniprot_ribosomal_proteins.csv', index_col=None, sep=';')
    Hackett_Clim_FCs = pandas.read_csv('../DataSetsYeastRBACalibration/Hacket_ProteinFCs.csv',sep=";")
    Nielsen_01 = pandas.read_csv('../DataSetsYeastRBACalibration/Nielsen01_ProteomicsData.csv',sep=";",index_col=0)
    Simulation = SessionRBA('../Yeast_iMM904_RBA_model')

    if import_initial_parameters:
        default_kapps_imported=pandas.read_csv("../Bootstrapping_Results/Input_parameters/default_kapps_refactored_WF.csv",index_col=0)
        specific_kapps_imported=pandas.read_csv("../Bootstrapping_Results/Input_parameters/specific_kapps_refactored_WF.csv",index_col=0)
        process_efficiencies_imported=pandas.read_csv("../Bootstrapping_Results/Input_parameters/process_efficiencies_refactored_WF.csv",index_col=0)
    else:
        default_kapps_imported=None
        specific_kapps_imported=None
        process_efficiencies_imported=None

    picogram_togram_coefficient = 1e12
    Nielsen_01['Mean_01'] *= picogram_togram_coefficient
    Nielsen_01['Rep1_01'] *= picogram_togram_coefficient
    Nielsen_01['Rep2_01'] *= picogram_togram_coefficient
    Nielsen_01['Rep3_01'] *= picogram_togram_coefficient

    full_annotations=build_full_annotations(rba_session=Simulation,
                        model_ribosomal_processes=['TranslationC', 'TranslationM'],
                        external_ribosome_annotations=Ribosomal_Proteins_Uniprot,
                        Compartment_Annotations_external=Compartment_Annotations_external,
                        Absolute_Proteome=Nielsen_01,
                        protein_ID_column='Gene',
                        Uniprot=Uniprot,
                        Relative_Proteome=Hackett_Clim_FCs)

    ###sampling
    Sampled_Proteomes=sample_copy_numbers_from_residuals_quantiles(Input_data=Nielsen_01.copy(),
                                                                    replicate_cols=["Rep1_01","Rep2_01","Rep3_01"],
                                                                    mean_col='Mean_01',
                                                                    replicate_threshold=2,
                                                                    filter_list=[],
                                                                    #filter_list=list(Hackett_Clim_FCs["Gene"]),
                                                                    target_size=number_samples,
                                                                    number_quantiles=12,
                                                                    transform_residuals=False,
                                                                    regression_type="inverse_lin",
                                                                    start_run_id=0,
                                                                    mean_no_noise=True,
                                                                    sample_mean=True)    
    """
    Sampled_Proteomes=generate_multiple_input_proteomes_from_mean(Input_data=Nielsen_01.copy(),
                                                                  mean_col='Mean_01',
                                                                  n=number_samples)
    """
    print(Sampled_Proteomes)
    sampled_runs=list(Sampled_Proteomes.columns)
    if chunk_size == 1:
        list_of_runs_to_execute=[sampled_runs]
    else:
        list_of_runs_to_execute=[sampled_runs[i * chunk_size:(i + 1) * chunk_size] for i in range((len(sampled_runs) + chunk_size - 1) // chunk_size )]

    print(list_of_runs_to_execute)
    count=0
    for runs_to_execute in list_of_runs_to_execute:
        ###generate list of input dicts over runs
        input_dicts=[]
        for run_ID in runs_to_execute:
            df=pandas.DataFrame(index=list(Sampled_Proteomes.index))
            df["Gene"]=list(df.index)
            df[run_ID]=Sampled_Proteomes[run_ID]
            restored_proteome_for_run=generate_input_proteome(fold_changes=Hackett_Clim_FCs,
                                                        absolute_data=df.loc[pandas.isna(df[run_ID]) == False],
                                                        matching_column_in_absolute_data=run_ID,
                                                        matching_column_in_fold_change_data='Hackett_C01',
                                                        conditions_in_fold_change_data_to_restore=conditions,
                                                        ID_column_in_restored_data='ID',
                                                        full_annotations=full_annotations)
            input_dicts.append({run_ID:{"xml_dir":'../Yeast_iMM904_RBA_model',
                                    "definition_file":Input_Data,
                                    "proteome":restored_proteome_for_run,
                                    "process_efficiency_estimation_input":Process_Efficiency_Estimation_Input,
                                    "Compartment_sizes":None,
                                    "PG_fractions":None,
                                    "conditions":conditions,
                                    "Process_efficiencies":process_efficiencies_imported,
                                    "Specific_kapps":specific_kapps_imported,
                                    "Default_kapps":default_kapps_imported
                                    }})

        ###calibrate each run
        if n_parallel_processes!=1:
            if n_parallel_processes is None:
                num_cores=cpu_count()
                n_jobs=min(num_cores,len(input_dicts))
            else:
                n_jobs=min(n_parallel_processes,len(input_dicts))            
            pool=Pool(n_jobs)
            bootstrapping_runs=pool.imap_unordered(run_calibration_over_conditions_for_bootstrapping_run,input_dicts)
        else:
            bootstrapping_runs=[run_calibration_over_conditions_for_bootstrapping_run(input_dict) for input_dict in input_dicts]
        print(bootstrapping_runs)
        ###generate output
        if count==0:
            Reconstructed_Proteomes={condition:pandas.DataFrame() for condition in conditions}
            Corrected_Proteomes={condition:pandas.DataFrame() for condition in conditions}
            Compartment_sizes={condition:pandas.DataFrame() for condition in conditions}
            Pg_fractions={condition:pandas.DataFrame() for condition in conditions}
            Process_efficiencies={condition:pandas.DataFrame() for condition in conditions}
            Specific_Kapps={condition:pandas.DataFrame() for condition in conditions}
            Specific_Kapps_directions={condition:pandas.DataFrame() for condition in conditions}
            Default_Kapps={condition:pandas.DataFrame() for condition in conditions}

            for run in bootstrapping_runs:
                run_ID=list(run.keys())[0]
                for condition in run[run_ID].keys():
                    for i in run[run_ID][condition]["proteome_input"].keys():
                        Reconstructed_Proteomes[condition].loc[i,run_ID]=run[run_ID][condition]["proteome_input"][i]
                    for i in run[run_ID][condition]["proteome_corrected"].keys():
                        Corrected_Proteomes[condition].loc[i,run_ID]=run[run_ID][condition]["proteome_corrected"][i]
                    for i in run[run_ID][condition]["densities"].keys():
                        Compartment_sizes[condition].loc[i,run_ID]=run[run_ID][condition]["densities"][i]
                    for i in run[run_ID][condition]["pgs"].keys():
                        Pg_fractions[condition].loc[i,run_ID]=run[run_ID][condition]["pgs"][i]
                    for i in run[run_ID][condition]["process_effs"].keys():
                        Process_efficiencies[condition].loc[i,run_ID]=run[run_ID][condition]["process_effs"][i]
                    for i in run[run_ID][condition]["spec_kapps"].keys():
                        Specific_Kapps[condition].loc[i,run_ID]=run[run_ID][condition]["spec_kapps"][i]["Kapp"]
                        Specific_Kapps_directions[condition].loc[i,run_ID]=run[run_ID][condition]["spec_kapps"][i]["Flux"]
                    Default_Kapps[condition].loc["Default_Kapp",run_ID]=run[run_ID][condition]["def_kapp"]["default_efficiency"]

            #if number_samples>1:
            #    for condition in conditions:
            #        for params in [Reconstructed_Proteomes,Corrected_Proteomes,Compartment_sizes,Pg_fractions,Specific_Kapps,Process_efficiencies,Default_Kapps]:
            #            df=params[condition]
            #            df["mean_run_param"]=[numpy.nanmean([df.loc[i,j] for j in df.columns if j.startswith("run_")]) for i in df.index]
            #        Specific_Kapps_directions[condition]["mean_run_param"]=Specific_Kapps_directions[condition]["mean_noNoise"]

        else:
            #load stored thing from csv
            Reconstructed_Proteomes={}
            Corrected_Proteomes={}
            Compartment_sizes={}
            Pg_fractions={}
            Process_efficiencies={}
            Specific_Kapps={}
            Specific_Kapps_directions={}
            Default_Kapps={}
            for condition in conditions:
                Reconstructed_Proteomes[condition]=pandas.read_csv("../Bootstrapping_Results/ReconstructedProteomes_{}.csv".format(condition),index_col=0)
                Corrected_Proteomes[condition]=pandas.read_csv("../Bootstrapping_Results/CorrectedProteomes_{}.csv".format(condition),index_col=0)
                Compartment_sizes[condition]=pandas.read_csv("../Bootstrapping_Results/CompartmentSizes_{}.csv".format(condition),index_col=0)
                Pg_fractions[condition]=pandas.read_csv("../Bootstrapping_Results/PgFractions_{}.csv".format(condition),index_col=0)
                Process_efficiencies[condition]=pandas.read_csv("../Bootstrapping_Results/ProcessEfficiencies_{}.csv".format(condition),index_col=0)
                Specific_Kapps[condition]=pandas.read_csv("../Bootstrapping_Results/SpecificKapps_{}.csv".format(condition),index_col=0)
                Specific_Kapps_directions[condition]=pandas.read_csv("../Bootstrapping_Results/SpecificKappDirections_{}.csv".format(condition),index_col=0)
                Default_Kapps[condition]=pandas.read_csv("../Bootstrapping_Results/DefaultKapps_{}.csv".format(condition),index_col=0)

            #add genereated cols to loaded ones
            for run in bootstrapping_runs:
                run_ID=list(run.keys())[0]
                for condition in run[run_ID].keys():
                    for i in run[run_ID][condition]["proteome_input"].keys():
                        Reconstructed_Proteomes[condition].loc[i,run_ID]=run[run_ID][condition]["proteome_input"][i]
                    for i in run[run_ID][condition]["proteome_corrected"].keys():
                        Corrected_Proteomes[condition].loc[i,run_ID]=run[run_ID][condition]["proteome_corrected"][i]
                    for i in run[run_ID][condition]["densities"].keys():
                        Compartment_sizes[condition].loc[i,run_ID]=run[run_ID][condition]["densities"][i]
                    for i in run[run_ID][condition]["pgs"].keys():
                        Pg_fractions[condition].loc[i,run_ID]=run[run_ID][condition]["pgs"][i]
                    for i in run[run_ID][condition]["process_effs"].keys():
                        Process_efficiencies[condition].loc[i,run_ID]=run[run_ID][condition]["process_effs"][i]
                    for i in run[run_ID][condition]["spec_kapps"].keys():
                        Specific_Kapps[condition].loc[i,run_ID]=run[run_ID][condition]["spec_kapps"][i]["Kapp"]
                        Specific_Kapps_directions[condition].loc[i,run_ID]=run[run_ID][condition]["spec_kapps"][i]["Flux"]
                    Default_Kapps[condition].loc["Default_Kapp",run_ID]=run[run_ID][condition]["def_kapp"]["default_efficiency"]

        #store 
        if (count+1)==len(list_of_runs_to_execute):
            for condition in conditions:
                for params in [Reconstructed_Proteomes,Corrected_Proteomes,Compartment_sizes,Pg_fractions,Specific_Kapps,Process_efficiencies,Default_Kapps]:
                    df=params[condition]
                    df["mean_run_param"]=[numpy.nanmean([df.loc[i,j] for j in df.columns if j.startswith("run_")]) for i in df.index]
                Specific_Kapps_directions[condition]["mean_run_param"]=Specific_Kapps_directions[condition]["mean_noNoise"]
        for condition in conditions:
            Reconstructed_Proteomes[condition].to_csv("../Bootstrapping_Results/ReconstructedProteomes_{}.csv".format(condition))
            Corrected_Proteomes[condition].to_csv("../Bootstrapping_Results/CorrectedProteomes_{}.csv".format(condition))
            Compartment_sizes[condition].to_csv("../Bootstrapping_Results/CompartmentSizes_{}.csv".format(condition))
            Pg_fractions[condition].to_csv("../Bootstrapping_Results/PgFractions_{}.csv".format(condition))
            Specific_Kapps[condition].to_csv("../Bootstrapping_Results/SpecificKapps_{}.csv".format(condition))
            Specific_Kapps_directions[condition].to_csv("../Bootstrapping_Results/SpecificKappDirections_{}.csv".format(condition))
            Process_efficiencies[condition].to_csv("../Bootstrapping_Results/ProcessEfficiencies_{}.csv".format(condition))
            Default_Kapps[condition].to_csv("../Bootstrapping_Results/DefaultKapps_{}.csv".format(condition))

        count+=1
#

if __name__ == "__main__":
    warnings.simplefilter('ignore', UserWarning)
    warnings.simplefilter('ignore', FutureWarning)
    warnings.simplefilter('ignore', RuntimeWarning)
    t0=time.time()
    main(n_parallel_processes=None,
         #conditions = ['Hackett_C01'],
         conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03'],
         number_samples=30,
         chunk_size=8,
         import_initial_parameters=True)
    print("Total time: {}".format(time.time()-t0))
