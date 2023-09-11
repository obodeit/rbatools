import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *
import time
import warnings
from multiprocessing import Pool , cpu_count

def run_calibration_over_conditions(input_dict,conditions,n_parallel_processes=None,bootstrapping=False):
    growth_rates={condition:growth_rate_from_input(input=input_dict["definition_file"], condition=condition) for condition in conditions}
    input_dict["preliminary_run"]=True
    input_dicts_for_conditions=[]
    for condition in conditions:
        input_dicts_for_conditions.append(input_dict.copy())
        input_dicts_for_conditions[-1]["condition"]=condition
    calib_dicts_1=[calibration(i,print_outputs= not bootstrapping) for i in input_dicts_for_conditions]

    calibration_results_1=[]
    for condition in conditions:
        for i in calib_dicts_1:
            if condition in list(i.keys()):
                calibration_results_1.append(i[condition])

    compartment_sizes_from_calibration_1=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=calibration_results_1)
    pg_fractions_from_calibration_1=extract_pg_fractions_from_calibration_outputs(calibration_outputs=calibration_results_1)
    regressed_compartment_sizes_1=regression_on_compartment_sizes(Comp_sizes=compartment_sizes_from_calibration_1,conditions=conditions,growth_rates=growth_rates,monotonous_quadratic=False)
    regressed_pg_fractions_1=regression_on_pg_fractions(PG_sizes=pg_fractions_from_calibration_1,conditions=conditions,growth_rates=growth_rates,monotonous_quadratic=False)

    for i in input_dicts_for_conditions:
        i.update({"Compartment_sizes":regressed_compartment_sizes_1,
                  "PG_fractions":regressed_pg_fractions_1,
                  "preliminary_run":False})

    if n_parallel_processes is None:
        num_cores=cpu_count()
        n_jobs=min(num_cores,len(conditions))
        pool=Pool(n_jobs)
        calib_dicts_2=pool.map(calibration,input_dicts_for_conditions)
    elif n_parallel_processes==1:
        #calib_dicts_2=[]
        #for i in input_dicts:
        #    print(i["condition"])
        #    calib_dicts_2.append(calibration(i))
        calib_dicts_2=[calibration(i, print_outputs= not bootstrapping) for i in input_dicts_for_conditions]
    else:
        pool=Pool(n_parallel_processes)
        calib_dicts_2=pool.map(calibration,input_dicts_for_conditions)

    calibration_results_2=[]
    corrected_proteomes_DF=pandas.DataFrame()
    for condition in conditions:
        for i in calib_dicts_2:
            if condition in list(i.keys()):
                calibration_results_2.append(i[condition])
                for gene in i[condition]["Proteome"].index:
                    corrected_proteomes_DF.loc[gene,condition]=i[condition]["Proteome"].loc[gene,"copy_number"]

    process_efficiencies_from_calibration=extract_process_capacities_from_calibration_outputs(calibration_outputs=calibration_results_2)
    specific_kapps_from_calibration=extract_specific_kapps_from_calibration_outputs(calibration_outputs=calibration_results_2)
    default_kapps_from_calibration=extract_default_kapps_from_calibration_outputs(calibration_outputs=calibration_results_2)

    if not bootstrapping:
        corrected_proteomes_DF.to_csv("../Corrected_calibration_proteomes.csv",sep=",",decimal=".")
        regressed_pg_fractions_1.to_csv("../pg_fractions_refactored_WF.csv")
        regressed_compartment_sizes_1.to_csv("../compartment_sizes_refactored_WF.csv")
        default_kapps_from_calibration.to_csv("../default_kapps_refactored_WF.csv")
        specific_kapps_from_calibration.to_csv("../specific_kapps_refactored_WF.csv")
        process_efficiencies_from_calibration.to_csv("../process_efficiencies_refactored_WF.csv")

        print("Total Runtime: {}".format(time.time()-initial_time))

        plot_compartment_sizes_and_pg(point_calibration_sizes=compartment_sizes_from_calibration_1,
                                    point_calibration_pg=pg_fractions_from_calibration_1,
                                    regressed_sizes=regressed_compartment_sizes_1,
                                    regressed_pg=regressed_pg_fractions_1,
                                    conditions=conditions,
                                    growth_rates=growth_rates,
                                    filename="../Compartment_sizes_and_PG.pdf")
        try:
            plot_rss_trajectory(calibration_results_2)
        except:
            print("No RSS trajectory")

def calibration(input_dict,print_outputs=True):
    Simulation = SessionRBA(input_dict["xml_dir"])
    Simulation.add_exchange_reactions()
    calib_results = calibration_workflow(proteome=input_dict["proteome"],
                                         condition=input_dict["condition"],
                                         reference_condition=input_dict["reference_condition"],
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

def main(conditions,n_parallel_processes=None):
    Input_Data = pandas.read_csv('../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nlim.csv', sep=';', decimal=',', index_col=0)
    Process_Efficiency_Estimation_Input = pandas.read_csv('../DataSetsYeastRBACalibration/Process_Efficiency_Estimation_Input.csv', sep=';', decimal=',')
    Uniprot = pandas.read_csv('../Yeast_iMM904_RBA_model/uniprot.csv', sep='\t')
    Compartment_Annotations_external = pandas.read_csv('../DataSetsYeastRBACalibration/Manually_curated_Protein_Locations_for_Calibration.csv', index_col=None, sep=';')
    Ribosomal_Proteins_Uniprot = pandas.read_csv('../DataSetsYeastRBACalibration/uniprot_ribosomal_proteins.csv', index_col=None, sep=';')
    Hackett_Clim_FCs = pandas.read_csv('../DataSetsYeastRBACalibration/Hacket_ProteinFCs.csv',sep=";")
    Nielsen_01 = pandas.read_csv('../DataSetsYeastRBACalibration/Nielsen01_ProteomicsData.csv',sep=";",index_col=0)
    Simulation = SessionRBA('../Yeast_iMM904_RBA_model')

    picogram_togram_coefficient = 1e12
    Reference_Condition='Mean_01'

    full_annotations=build_full_annotations(rba_session=Simulation,
                        model_ribosomal_processes=['TranslationC', 'TranslationM'],
                        external_ribosome_annotations=Ribosomal_Proteins_Uniprot,
                        Compartment_Annotations_external=Compartment_Annotations_external,
                        Absolute_Proteome=Nielsen_01,
                        protein_ID_column='Gene',
                        Uniprot=Uniprot,
                        Relative_Proteome=Hackett_Clim_FCs)

    Nielsen_01[Reference_Condition] *= picogram_togram_coefficient
    restored_Hackett_Data=generate_input_proteome(fold_changes=Hackett_Clim_FCs,
                                                absolute_data=Nielsen_01.loc[pandas.isna(Nielsen_01[Reference_Condition]) == False],
                                                matching_column_in_absolute_data=Reference_Condition,
                                                matching_column_in_fold_change_data='Hackett_C01',
                                                conditions_in_fold_change_data_to_restore=conditions,
                                                ID_column_in_restored_data='ID',
                                                full_annotations=full_annotations)

    restored_Hackett_Data.to_csv("../origRestoredProteome.csv")
    initial_time=time.time()

    input_dict={}
    input_dict["xml_dir"]='../Yeast_iMM904_RBA_model'
    input_dict["proteome"]=restored_Hackett_Data
    input_dict["reference_condition"]=Reference_Condition
    input_dict["definition_file"]=Input_Data
    input_dict["process_efficiency_estimation_input"]=Process_Efficiency_Estimation_Input
    input_dict["Compartment_sizes"]=None
    input_dict["PG_fractions"]=None

    run_calibration_over_conditions(input_dict=input_dict,
                                    conditions=conditions,
                                    n_parallel_processes=n_parallel_processes,
                                    bootstrapping=False)

if __name__ == "__main__":
    warnings.simplefilter('ignore', UserWarning)
    warnings.simplefilter('ignore', FutureWarning)
    warnings.simplefilter('ignore', RuntimeWarning)
    #warnings.simplefilter('ignore', SettingWithCopyWarning)
    main(n_parallel_processes=3,
        #conditions = ['Hackett_C03']
        conditions = ['Hackett_C03','Hackett_C005', 'Hackett_C022', 'Hackett_C016', 'Hackett_C01']
        #conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03']
        #conditions = ['Hackett_N005', 'Hackett_N01', 'Hackett_N016', 'Hackett_N03']
        #conditions = ['Hackett_P005', 'Hackett_P01', 'Hackett_P016', 'Hackett_P022']
        )
