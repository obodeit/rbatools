import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *
from rbatools.bootstrapping_utils import *
from rbatools.rba_xml_utils import *
from rbatools.regression_utils import *
from rbatools.plotting_analysis_utils import *
from rbatools.other_utils import *
import time
import warnings
from multiprocessing import Pool , cpu_count


def run_calibration_over_conditions(input_dict,n_parallel_processes=None,regression_on_compartments=False,conditions_for_compartment_regression=None):
    initial_time=time.time()
    growth_rates={condition:growth_rate_from_input(input=input_dict["definition_file"], condition=condition) for condition in input_dict["conditions"]}
    input_dict["preliminary_run"]=True
    input_dicts_for_conditions=[]
    for condition in input_dict["conditions"]:
        input_dicts_for_conditions.append(input_dict.copy())
        input_dicts_for_conditions[-1]["condition"]=condition
    calib_dicts_1=[calibration(i,print_outputs= True) for i in input_dicts_for_conditions]

    calibration_results_1=[]
    for condition in input_dict["conditions"]:
        for i in calib_dicts_1:
            if condition in list(i.keys()):
                calibration_results_1.append(i[condition])

    compartment_sizes_from_calibration_1=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=calibration_results_1)
    pg_fractions_from_calibration_1=extract_pg_fractions_from_calibration_outputs(calibration_outputs=calibration_results_1)
    if conditions_for_compartment_regression is not None:
        regressed_compartment_sizes_1=regression_on_compartment_sizes(filename='../Calibration_stuff/CompSize_Plots_refined3.pdf',Comp_sizes=compartment_sizes_from_calibration_1,conditions=conditions_for_compartment_regression,growth_rates=growth_rates,monotonous_quadratic=False)
        regressed_pg_fractions_1=regression_on_pg_fractions(filename='../Calibration_stuff/PGfraction_Plots_refined3.pdf',PG_sizes=pg_fractions_from_calibration_1,conditions=conditions_for_compartment_regression,growth_rates=growth_rates,monotonous_quadratic=False)
        non_regression_conditions=[i for i in input_dict["conditions"] if i not in conditions_for_compartment_regression]
        if len(non_regression_conditions)!=0:
            for i in non_regression_conditions:
                for comp in regressed_compartment_sizes_1.index:
                    regressed_compartment_sizes_1.loc[comp,i]=compartment_sizes_from_calibration_1.loc[comp,i]
                for comp in regressed_pg_fractions_1.index:
                    regressed_pg_fractions_1.loc[comp,i]=pg_fractions_from_calibration_1.loc[comp,i]
    else:
        regressed_compartment_sizes_1=regression_on_compartment_sizes(filename='../Calibration_stuff/CompSize_Plots_refined3.pdf',Comp_sizes=compartment_sizes_from_calibration_1,conditions=input_dict["conditions"],growth_rates=growth_rates,monotonous_quadratic=False)
        regressed_pg_fractions_1=regression_on_pg_fractions(filename='../Calibration_stuff/PGfraction_Plots_refined3.pdf',PG_sizes=pg_fractions_from_calibration_1,conditions=input_dict["conditions"],growth_rates=growth_rates,monotonous_quadratic=False)

    if regression_on_compartments:
        for i in input_dicts_for_conditions:
            i.update({"Compartment_sizes":regressed_compartment_sizes_1,
                      "PG_fractions":regressed_pg_fractions_1,
                      "preliminary_run":False})
    else:
        for i in input_dicts_for_conditions:
            i.update({"Compartment_sizes":compartment_sizes_from_calibration_1,
                    "PG_fractions":pg_fractions_from_calibration_1,
                    "preliminary_run":False})

    if n_parallel_processes is None:
        num_cores=cpu_count()
        n_jobs=min(num_cores,len(input_dict["conditions"]))
        pool=Pool(n_jobs)
        calib_dicts_2=pool.map(calibration,input_dicts_for_conditions)
    elif n_parallel_processes==1:
        #calib_dicts_2=[]
        #for i in input_dicts:
        #    print(i["condition"])
        #    calib_dicts_2.append(calibration(i))
        calib_dicts_2=[calibration(i, print_outputs=True) for i in input_dicts_for_conditions]
    else:
        pool=Pool(n_parallel_processes)
        calib_dicts_2=pool.map(calibration,input_dicts_for_conditions)

    calibration_results_2=[]
    corrected_proteomes_DF=pandas.DataFrame()
    for condition in input_dict["conditions"]:
        for i in calib_dicts_2:
            if condition in list(i.keys()):
                calibration_results_2.append(i[condition])
                for gene in i[condition]["Proteome"].index:
                    corrected_proteomes_DF.loc[gene,condition]=i[condition]["Proteome"].loc[gene,"copy_number"]

    process_efficiencies_from_calibration=extract_process_capacities_from_calibration_outputs(calibration_outputs=calibration_results_2)
    specific_kapps_from_calibration=extract_specific_kapps_from_calibration_outputs(calibration_outputs=calibration_results_2)
    default_kapps_from_calibration=extract_default_kapps_from_calibration_outputs(calibration_outputs=calibration_results_2)

    corrected_proteomes_DF.to_csv("../Calibration_stuff/Corrected_calibration_proteomes.csv",sep=",",decimal=".")
    regressed_pg_fractions_1.to_csv("../Calibration_stuff/pg_fractions_refactored_WF.csv")
    pg_fractions_from_calibration_1.to_csv("../Calibration_stuff/pg_fractions_refactored_WF_original_not_regressed.csv")
    regressed_compartment_sizes_1.to_csv("../Calibration_stuff/compartment_sizes_refactored_WF.csv")
    compartment_sizes_from_calibration_1.to_csv("../Calibration_stuff/compartment_sizes_refactored_WF_original_not_regressed.csv")
    default_kapps_from_calibration.to_csv("../Calibration_stuff/default_kapps_refactored_WF.csv")
    specific_kapps_from_calibration.to_csv("../Calibration_stuff/specific_kapps_refactored_WF.csv")
    process_efficiencies_from_calibration.to_csv("../Calibration_stuff/process_efficiencies_refactored_WF.csv")

    print("Total Runtime: {}".format(time.time()-initial_time))

    plot_compartment_sizes_and_pg(point_calibration_sizes=compartment_sizes_from_calibration_1,
                                point_calibration_pg=pg_fractions_from_calibration_1,
                                regressed_sizes=regressed_compartment_sizes_1,
                                regressed_pg=regressed_pg_fractions_1,
                                conditions=input_dict["conditions"],
                                growth_rates=growth_rates,
                                filename="../Calibration_stuff/Compartment_sizes_and_PG.pdf")
    try:
        plot_rss_trajectory(filename="../Calibration_stuff/RRS_trajectory.csv",calibration_outputs=calibration_results_2)
    except:
        print("No RSS trajectory")

def calibration(input_dict,print_outputs=True):
    
    Simulation = SessionRBA(input_dict["xml_dir"])

    Simulation.add_exchange_reactions()
    calib_results = calibration_workflow(proteome=input_dict["proteome"],
                                         condition=input_dict["condition"],
                                         gene_ID_column='Gene',
                                         definition_file=input_dict["definition_file"],
                                         rba_session=Simulation,
                                         process_efficiency_estimation_input=input_dict["process_efficiency_estimation_input"],
                                         spec_kapps=input_dict["Specific_kapps"],
                                         default_kapps=input_dict["Default_kapps"],
                                         process_efficiencies=input_dict["Process_efficiencies"],
                                         Compartment_sizes=input_dict["Compartment_sizes"],
                                         PG_fractions=input_dict["PG_fractions"],
                                         transporter_multiplier=1,
                                         prelim_run=input_dict["preliminary_run"],
                                         #Mu_approx_precision=  0.000000001,
                                         Mu_approx_precision= 0.0000001,
                                         feasible_stati=["optimal","feasible"],
                                         #feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                         min_kapp=None,
                                         print_outputs=print_outputs,
                                         use_mean_enzyme_composition_for_calibration=False,
                                         global_protein_scaling_coeff=1000/6.022e23,
                                         #max_kapp_threshold=100000000000
                                         max_kapp_threshold=None,
                                         output_dir="../Calibration_stuff"
                                         )
    return({input_dict["condition"]:calib_results})

def generate_input_proteome(fold_changes,
                            absolute_data,
                            matching_column_in_absolute_data,
                            matching_column_in_fold_change_data,
                            conditions_in_fold_change_data_to_restore,
                            ID_column_in_restored_data,
                            full_annotations,
                            reference_condition=None):
    proteomics_data = infer_copy_numbers_from_reference_copy_numbers(fold_changes=fold_changes,
                                                                        absolute_data=absolute_data,
                                                                        matching_column_in_fold_change_data=matching_column_in_fold_change_data,
                                                                        matching_column_in_absolute_data=matching_column_in_absolute_data,
                                                                        conditions_in_fold_change_data_to_restore=conditions_in_fold_change_data_to_restore)

    if reference_condition is not None:
        if reference_condition in list(absolute_data.columns):
            for gene in absolute_data["Gene"]:
                if gene not in proteomics_data["ID"]:
                    proteomics_data.loc[gene,"ID"]=gene
                xxx=absolute_data.loc[absolute_data["Gene"]==gene,reference_condition]
                proteomics_data.loc[proteomics_data["ID"]==gene,reference_condition]=xxx

    proteomics_data_with_annotations= add_annotations_to_proteome(input=proteomics_data, ID_column=ID_column_in_restored_data, annotations=full_annotations)
    return(proteomics_data_with_annotations)

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

def main(model_xml_dir,conditions,n_parallel_processes=None,regression_on_compartments=False,conditions_for_compartment_regression=None):

    Input_Data = pandas.read_csv('../Calibration_stuff/DataSetsYeastRBACalibration/Calibration_InputDefinition_plusFrick_fluxes_COARSE_Settings.csv', sep=';', decimal=',', index_col=0)
    #Input_Data = pandas.read_csv('../Calibration_stuff/DataSetsYeastRBACalibration/Calibration_InputDefinition_plusFrick_fluxes.csv', sep=';', decimal=',', index_col=0)
    
    Process_Efficiency_Estimation_Input = pandas.read_csv('../Calibration_stuff/DataSetsYeastRBACalibration/Process_Efficiency_Estimation_Input.csv', sep=';', decimal=',')
    Uniprot = pandas.read_csv('{}/uniprot.csv'.format(model_xml_dir), sep='\t')
    Compartment_Annotations_external = pandas.read_csv('../Calibration_stuff/DataSetsYeastRBACalibration/Manually_curated_Protein_Locations_for_Calibration.csv', index_col=None, sep=';')
    Ribosomal_Proteins_Uniprot = pandas.read_csv('../Calibration_stuff/DataSetsYeastRBACalibration/uniprot_ribosomal_proteins.csv', index_col=None, sep=';')

    Hackett_Clim_FCs = pandas.read_csv('../Calibration_stuff/DataSetsYeastRBACalibration/Hacket_ProteinFCs.csv',sep=";")
    Nielsen_01 = pandas.read_csv('../Calibration_stuff/DataSetsYeastRBACalibration/Nielsen01_ProteomicsData.csv',sep=";",index_col=0)

    Simulation = SessionRBA(model_xml_dir)

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
    
    full_annotations.to_csv("../Calibration_stuff/full_annot.csv")

    Nielsen_01[Reference_Condition] *= picogram_togram_coefficient

    if Reference_Condition in conditions:
        reference_condition_for_input_proteome=Reference_Condition
    else:
        reference_condition_for_input_proteome=None

    proteomics_data=generate_input_proteome(fold_changes=Hackett_Clim_FCs,
                                                absolute_data=Nielsen_01.loc[pandas.isna(Nielsen_01[Reference_Condition]) == False],
                                                matching_column_in_absolute_data=Reference_Condition,
                                                matching_column_in_fold_change_data='Hackett_C01',
                                                conditions_in_fold_change_data_to_restore=conditions,
                                                ID_column_in_restored_data='ID',
                                                full_annotations=full_annotations,
                                                reference_condition=reference_condition_for_input_proteome)

    proteomics_data.to_csv("../Calibration_stuff/original_proteomes_for_calibration.csv")

    input_dict={}

    input_dict["conditions"]=conditions
    input_dict["xml_dir"]=model_xml_dir
    input_dict["proteome"]=proteomics_data
    input_dict["definition_file"]=Input_Data
    input_dict["process_efficiency_estimation_input"]=Process_Efficiency_Estimation_Input
    input_dict["Compartment_sizes"]=None
    input_dict["PG_fractions"]=None
    input_dict["Specific_kapps"]=None
    input_dict["Default_kapps"]=None
    input_dict["Process_efficiencies"]=None

    run_calibration_over_conditions(input_dict=input_dict,
                                    n_parallel_processes=n_parallel_processes,
                                    regression_on_compartments=regression_on_compartments,
                                    conditions_for_compartment_regression=conditions_for_compartment_regression)

if __name__ == "__main__":
    warnings.simplefilter('ignore', UserWarning)
    warnings.simplefilter('ignore', FutureWarning)
    warnings.simplefilter('ignore', RuntimeWarning)
    #warnings.simplefilter('ignore', SettingWithCopyWarning)
    main(model_xml_dir='../Yeast_iMM904_RBA_model',
         n_parallel_processes=5,
        regression_on_compartments=True,
        #conditions = ['Hackett_C01','Mean_01'],
        #conditions = ['Hackett_C03','Hackett_C005', 'Hackett_C022', 'Hackett_C016', 'Hackett_C01','Mean_01'],
        conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03'],
        #conditions = ['Hackett_N005', 'Hackett_N01', 'Hackett_N016', 'Hackett_N03'],
        #conditions = ['Hackett_P005', 'Hackett_P01', 'Hackett_P016', 'Hackett_P022'],
        conditions_for_compartment_regression=['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03']
        )
