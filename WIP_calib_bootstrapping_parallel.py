import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *
import time
from multiprocessing import Pool , cpu_count
from joblib import Parallel, delayed

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

def run_full_calib_from_input_dict(input_dict):
    run=list(input_dict.keys())[0]
    Simulation = SessionRBA(input_dict[run]["xml_dir"])
    Simulation.add_exchange_reactions()
    calibration_results_1=[]
    for condition in input_dict[run]["conditions"]:
        Simulation.reload_model()
        calib_results = calibration_workflow(proteome=input_dict[run]["proteome"],
                                            condition=condition,
                                            reference_condition=run,
                                            gene_ID_column='Gene',
                                            definition_file=input_dict[run]["Input"],
                                            rba_session=Simulation,
                                            process_efficiency_estimation_input=input_dict[run]["Process_Input"],
                                            spec_kapps=None,
                                            corrected_spec_kapps=False,
                                            take_lowest_RSS_parameters=True,
                                            process_efficiencies=None,
                                            Compartment_sizes=None,
                                            PG_fractions=None,
                                            prelim_run=True)
        calibration_results_1.append(calib_results)

    pg_fractions_from_calibration_1=extract_pg_fractions_from_calibration_outputs(calibration_outputs=calibration_results_1)
    compartment_sizes_from_calibration_1 = extract_compartment_sizes_from_calibration_outputs(calibration_outputs=calibration_results_1)
    regressed_compartment_sizes_1=regression_on_compartment_sizes(Comp_sizes=compartment_sizes_from_calibration_1,conditions=input_dict[run]["conditions"],growth_rates=input_dict[run]["growth_rates"])
    regressed_pg_fractions_1=regression_on_pg_fractions(PG_sizes=pg_fractions_from_calibration_1,conditions=input_dict[run]["conditions"],growth_rates=input_dict[run]["growth_rates"])

    out={}
    for condition in input_dict[run]["conditions"]:
        Simulation.reload_model()
        print("Run {} - Calib 2: {}".format(run,condition))
        calib_results = calibration_workflow(proteome=input_dict[run]["proteome"],
                                            condition=condition,
                                            reference_condition=run,
                                            gene_ID_column='Gene',
                                            definition_file=input_dict[run]["Input"],
                                            rba_session=Simulation,
                                            process_efficiency_estimation_input=input_dict[run]["Process_Input"],
                                            spec_kapps=None,
                                            corrected_spec_kapps=True,
                                            take_lowest_RSS_parameters=True,
                                            process_efficiencies=None,
                                            Compartment_sizes=regressed_compartment_sizes_1,
                                            PG_fractions=regressed_pg_fractions_1,
                                            transporter_multiplier=1,
                                            prelim_run=False,
                                            final_global_scaling_free_exchanges=True,
                                            Mu_approx_precision=0.00001,
                                            feasible_stati=["optimal","feasible"],
                                            min_kapp=None,
                                            fixed_mu_when_above_target_mu_in_correction=True,
                                            mu_misprediction_tolerance=0.05,
                                            print_outputs=False)

        proteome_input={i:input_dict[run]["proteome"].loc[i,condition] for i in input_dict[run]["proteome"].index}
        proteome_corrected={i:calib_results["Proteome"].loc[i,"copy_number"] for i in calib_results["Proteome"].index}
        densities={i:round(calib_results["Densities_PGs"].loc[i,"Density"],6) for i in calib_results["Densities_PGs"].index}
        pgs={i:round(calib_results["Densities_PGs"].loc[i,"PG_fraction"],6) for i in calib_results["Densities_PGs"].index}
        spec_kapps={i:{"Kapp":calib_results["Specific_Kapps"].loc[calib_results["Specific_Kapps"]["Enzyme_ID"]==i,"Kapp"].values[0],"Flux":calib_results["Specific_Kapps"].loc[calib_results["Specific_Kapps"]["Enzyme_ID"]==i,"Flux"].values[0]} for i in calib_results["Specific_Kapps"]["Enzyme_ID"]}
        process_effs={i:calib_results["Process_Efficiencies"].loc[calib_results["Process_Efficiencies"]['Process']==i,"Value"].values[0] for i in calib_results["Process_Efficiencies"]['Process']}
        def_kapp={"default_efficiency":calib_results["Default_Kapps"]['default_efficiency']}

        out[condition]={"proteome_input":proteome_input,
                        "proteome_corrected":proteome_corrected,
                        "densities":densities,
                        "pgs":pgs,
                        "spec_kapps":spec_kapps,
                        "process_effs":process_effs,
                        "def_kapp":def_kapp}
    return({run:out})
##################################################
##################################################

def main(conditions,growth_rates,number_samples=1,n_parallel_processes=None):
    Input_Data = pandas.read_csv('../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nlim.csv', sep=';', decimal=',', index_col=0)
    Process_Efficiency_Estimation_Input = pandas.read_csv('../DataSetsYeastRBACalibration/Process_Efficiency_Estimation_Input.csv', sep=';', decimal=',')
    Uniprot = pandas.read_csv('../Yeast_iMM904_RBA_model/uniprot.csv', sep='\t')
    Compartment_Annotations_external = pandas.read_csv('../DataSetsYeastRBACalibration/Manually_curated_Protein_Locations_for_Calibration.csv', index_col=None, sep=';')
    Ribosomal_Proteins_Uniprot = pandas.read_csv('../DataSetsYeastRBACalibration/uniprot_ribosomal_proteins.csv', index_col=None, sep=';')
    Hackett_Clim_FCs = pandas.read_csv('../DataSetsYeastRBACalibration/Hacket_ProteinFCs.csv',sep=";")
    Nielsen_01 = pandas.read_csv('../DataSetsYeastRBACalibration/Nielsen01_ProteomicsData.csv',sep=";",index_col=0)
    Simulation = SessionRBA('../Yeast_iMM904_RBA_model')

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

    Sampled_Proteomes=sample_copy_numbers_from_residuals_quantiles(Input_data=Nielsen_01.copy(),
                                                                    replicate_cols=["Rep1_01","Rep2_01","Rep3_01"],
                                                                    mean_col='Mean_01',
                                                                    replicate_threshold=2,
                                                                    filter_list=[],
                                                                    #filter_list=list(Hackett_Clim_FCs["Gene"]),
                                                                    target_size=number_samples,
                                                                    number_quantiles=12,
                                                                    transform_residuals=False,
                                                                    regression_type="inverse_lin")

    #Sampled_Proteomes.to_csv("Sampled_test.csv")

    input_dicts=[]
    for run in Sampled_Proteomes.columns:
        df=pandas.DataFrame(index=list(Sampled_Proteomes.index))
        df["Gene"]=list(df.index)
        df[run]=Sampled_Proteomes[run]
        restored_Hackett_Data=generate_input_proteome(fold_changes=Hackett_Clim_FCs,
                                                    absolute_data=df.loc[pandas.isna(df[run]) == False],
                                                    matching_column_in_absolute_data=run,
                                                    matching_column_in_fold_change_data='Hackett_C01',
                                                    conditions_in_fold_change_data_to_restore=['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03'],
                                                    ID_column_in_restored_data='ID',
                                                    full_annotations=full_annotations)
        input_dicts.append({run:{"proteome":restored_Hackett_Data,
                                    "xml_dir":"../Yeast_iMM904_RBA_model",
                                    "Input":Input_Data,
                                    "Process_Input":Process_Efficiency_Estimation_Input,
                                    "conditions":conditions,
                                    "growth_rates":growth_rates}})
    if n_parallel_processes!=1:
        if n_parallel_processes is None:
            num_cores=cpu_count()
            n_jobs=min(num_cores,len(input_dicts))
        else:
            n_jobs=n_parallel_processes
        pool=Pool(n_jobs)
        bootstrapping_runs=pool.imap_unordered(run_full_calib_from_input_dict,input_dicts)
        #bootstrapping_runs=Parallel(n_jobs=n_jobs)(delayed(run_full_calib_from_input_dict)(input_dict) for input_dict in input_dicts)
    else:
        bootstrapping_runs=[run_full_calib_from_input_dict(input_dict) for input_dict in input_dicts]

    Reconstructed_Proteomes={condition:pandas.DataFrame() for condition in conditions}
    Corrected_Proteomes={condition:pandas.DataFrame() for condition in conditions}
    Compartment_sizes={condition:pandas.DataFrame() for condition in conditions}
    PG_fractions={condition:pandas.DataFrame() for condition in conditions}
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
                PG_fractions[condition].loc[i,run_ID]=run[run_ID][condition]["pgs"][i]
            for i in run[run_ID][condition]["process_effs"].keys():
                Process_efficiencies[condition].loc[i,run_ID]=run[run_ID][condition]["process_effs"][i]
            for i in run[run_ID][condition]["spec_kapps"].keys():
                Specific_Kapps[condition].loc[i,run_ID]=run[run_ID][condition]["spec_kapps"][i]["Kapp"]
                Specific_Kapps_directions[condition].loc[i,run_ID]=run[run_ID][condition]["spec_kapps"][i]["Flux"]
            Default_Kapps[condition].loc["Default_Kapp",run_ID]=run[run_ID][condition]["def_kapp"]["default_efficiency"]

    if number_samples>1:
        for condition in conditions:
            for params in [Reconstructed_Proteomes,Corrected_Proteomes,Compartment_sizes,PG_fractions,Specific_Kapps,Process_efficiencies,Default_Kapps]:
                df=params[condition]
                df["mean_run_param"]=[numpy.nanmean([df.loc[i,j] for j in df.columns if j.startswith("run_")]) for i in df.index]
            Specific_Kapps_directions[condition]["mean_run_param"]=Specific_Kapps_directions[condition]["mean_noNoise"]

    for condition in conditions:
        Reconstructed_Proteomes[condition].to_csv("../Bootstrapping_Results/ReconstructedProteomes_{}.csv".format(condition))
        Corrected_Proteomes[condition].to_csv("../Bootstrapping_Results/CorrectedProteomes_{}.csv".format(condition))
        Compartment_sizes[condition].to_csv("../Bootstrapping_Results/CompartmentSizes_{}.csv".format(condition))
        PG_fractions[condition].to_csv("../Bootstrapping_Results/PgFractions_{}.csv".format(condition))
        Specific_Kapps[condition].to_csv("../Bootstrapping_Results/SpecificKapps_{}.csv".format(condition))
        Specific_Kapps_directions[condition].to_csv("../Bootstrapping_Results/SpecificKappDirections_{}.csv".format(condition))
        Process_efficiencies[condition].to_csv("../Bootstrapping_Results/ProcessEfficiencies_{}.csv".format(condition))
        Default_Kapps[condition].to_csv("../Bootstrapping_Results/DefaultKapps_{}.csv".format(condition))

if __name__ == "__main__":
    t0=time.time()
    main(conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03'],
         growth_rates={'Hackett_C005':0.05, 'Hackett_C01':0.105, 'Hackett_C016':0.154, 'Hackett_C022':0.214, 'Hackett_C03':0.294},
         number_samples=50,
         n_parallel_processes=4)
    print("Total time: {}".format(time.time()-t0))
