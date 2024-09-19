import rba
import copy
import pandas
import time
import numpy
#import json
from scipy.stats.mstats import gmean
from rbatools.other_utils import medium_concentrations_from_input , machinery_efficiency_correction_settings_from_input , enzyme_efficiency_estimation_settings_from_input , flux_bounds_from_input , growth_rate_from_input , proteome_fractions_from_input , perform_simulations , perform_simulations_fixed_Mu


def calibration_workflow(proteome,
                         condition,
                         gene_ID_column,
                         definition_file,
                         rba_session,
                         process_efficiency_estimation_input=None,
                         spec_kapps=None,
                         default_kapps=None,
                         process_efficiencies=None,
                         Compartment_sizes=None,
                         PG_fractions=None,
                         transporter_multiplier=3,
                         prelim_run=False,
                         Mu_approx_precision=0.00001,
                         feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                         min_kapp=None,
                         print_outputs=True,
                         global_protein_scaling_coeff=1,
                         use_mean_enzyme_composition_for_calibration=False,
                         max_kapp_threshold=None,
                         output_dir=""):
    """
    _summary_

    Parameters
    ----------
    proteome : pandas.DataFrame
        _description_
    condition : str
        _description_
    gene_ID_column : str
        _description_
    definition_file : pandas.DataFrame
        _description_
    rba_session : rbatools.rba_session.SessionRBA
        _description_
    process_efficiency_estimation_input : pandas.DataFrame, optional
        _description_, by default None
    spec_kapps : _type_, optional
        _description_, by default None
    default_kapps : _type_, optional
        _description_, by default None
    process_efficiencies : _type_, optional
        _description_, by default None
    Compartment_sizes : _type_, optional
        _description_, by default None
    PG_fractions : _type_, optional
        _description_, by default None
    transporter_multiplier : int, optional
        _description_, by default 3
    prelim_run : bool, optional
        _description_, by default False
    Mu_approx_precision : float, optional
        _description_, by default 0.00001
    feasible_stati : list, optional
        _description_, by default ["optimal","feasible","feasible_only_before_unscaling"]
    min_kapp : float, optional
        _description_, by default None
    print_outputs : bool, optional
        _description_, by default True
    global_protein_scaling_coeff : int, optional
        _description_, by default 1
    use_mean_enzyme_composition_for_calibration : bool, optional
        _description_, by default False
    max_kapp_threshold : float, optional
        _description_, by default None
    output_dir : str, optional
        _description_, by default ""
    """
    correction_settings=machinery_efficiency_correction_settings_from_input(input=definition_file, condition=condition)
    enzyme_efficiency_estimation_settings=enzyme_efficiency_estimation_settings_from_input(input=definition_file, condition=condition)

    t0 = time.time()
    correction_results_compartement_sizes = correction_pipeline(input=proteome,
                                             condition=condition,
                                             definition_file=definition_file,
                                             compartments_to_replace={'DEF':"c", 'DEFA':"c", 'Def':"c"},
                                             compartments_no_original_PG=['n', 'Secreted'],
                                             fractions_entirely_replaced_with_expected_value=['Ribosomes'],
                                             imposed_compartment_fractions=proteome_fractions_from_input(input=definition_file, condition=condition),
                                             directly_corrected_compartments=['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'],
                                             merged_compartments={'c': 'Ribosomes'},
                                             min_compartment_fraction=0.00000)
    rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))
    if prelim_run:
        compartment_densities_and_PGs = extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=correction_results_compartement_sizes,rows_to_exclude=["Ribosomes","Total"]+[i for i in correction_results_compartement_sizes.index if i.startswith("pg_")])
        correction_results_compartement_sizes.to_csv(str(output_dir+'/Correction_overview_HackettNielsen_'+condition+'.csv'))
        return({"Densities_PGs":compartment_densities_and_PGs,
                "Condition":condition})

    if Compartment_sizes is not None:
        for i in Compartment_sizes.index:
            correction_results_compartement_sizes.loc[i,"new_protein_fraction"]=Compartment_sizes.loc[i,condition]
            if i in PG_fractions.index:
                correction_results_compartement_sizes.loc[i,"new_PG_fraction"]=PG_fractions.loc[i,condition]

    ### GENERALISE for all merged compartments###
    for i in correction_results_compartement_sizes.index:
        abundance_coeff=1
        if i =="c":
            abundance_coeff=(correction_results_compartement_sizes.loc[i,"new_protein_fraction"]*(1-correction_results_compartement_sizes.loc[i,"new_PG_fraction"])-correction_results_compartement_sizes.loc["Ribosomes","new_protein_fraction"])/(correction_results_compartement_sizes.loc[i,"original_protein_fraction"]*(1-correction_results_compartement_sizes.loc[i,"original_PG_fraction"]))
        else:
            abundance_coeff=(correction_results_compartement_sizes.loc[i,"new_protein_fraction"]*(1-correction_results_compartement_sizes.loc[i,"new_PG_fraction"]))/(correction_results_compartement_sizes.loc[i,"original_protein_fraction"]*(1-correction_results_compartement_sizes.loc[i,"original_PG_fraction"]))
        proteome.loc[proteome["Location"]==i,condition]*=abundance_coeff
        correction_results_compartement_sizes.loc[i,"copy_number_scaling"]=abundance_coeff
    ###
    correction_results_compartement_sizes.to_csv(str(output_dir+'/Correction_overview_HackettNielsen_corrected_'+condition+'.csv'))
    if process_efficiencies is None:
        if process_efficiency_estimation_input is not None:
            process_efficiencies = determine_apparent_process_efficiencies(growth_rate=growth_rate_from_input(input=definition_file,
                                                                           condition=condition),
                                                                           input=process_efficiency_estimation_input,
                                                                           rba_session=rba_session,
                                                                           protein_data=proteome.copy(),
                                                                           proteome_summary=correction_results_compartement_sizes.copy(),
                                                                           condition=condition,
                                                                           gene_id_col=gene_ID_column,
                                                                           fit_nucleotide_assembly_machinery=True)
    else:
        process_efficiencies=import_process_efficiencies(input_data=process_efficiencies,rba_session=rba_session,condition=condition)

    ### define coeff as input ###
    proteome[condition]*=global_protein_scaling_coeff
    ###
    process_efficiencies.to_csv(output_dir+'/ProcEffsOrig_{}.csv'.format(condition))

    compartment_densities_and_PGs = extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=correction_results_compartement_sizes,
                                                                                              rows_to_exclude=["Ribosomes","Total"]+[i for i in correction_results_compartement_sizes.index if i.startswith("pg_")])
    if use_mean_enzyme_composition_for_calibration:
        generate_mean_enzyme_composition_model(rba_session,condition)

        
    flux_bounds_fba=flux_bounds_from_input(input=definition_file,
                                           rba_session=rba_session, 
                                           condition=condition, 
                                           specific_exchanges=None, 
                                           specific_directions=None,
                                           also_consider_iso_enzmes=False)

    Specific_Kapps_Results = estimate_specific_enzyme_efficiencies(rba_session=rba_session, 
                                                                   proteomicsData=build_input_proteome_for_specific_kapp_estimation(proteome, condition), 
                                                                   flux_bounds=flux_bounds_fba, 
                                                                   compartment_densities_and_pg=compartment_densities_and_PGs,
                                                                   mu=growth_rate_from_input(input=definition_file, condition=condition), 
                                                                   biomass_function=enzyme_efficiency_estimation_settings['biomass_function_in_model'], 
                                                                   target_biomass_function=enzyme_efficiency_estimation_settings['use_target_biomass_function'],
                                                                   #rRNA_target_rna_weight_ratio=enzyme_efficiency_estimation_settings['rRNA_target_rna_weight_ratio'],
                                                                   parsimonious_fba=enzyme_efficiency_estimation_settings['parsimonious_fba'], 
                                                                   chose_most_likely_isoreactions=enzyme_efficiency_estimation_settings['chose_most_likely_isoreactions'],
                                                                   impose_on_all_isoreactions=enzyme_efficiency_estimation_settings['impose_on_all_isoreactions'], 
                                                                   zero_on_all_isoreactions=enzyme_efficiency_estimation_settings['zero_on_all_isoreactions'],
                                                                   #node_degree_identical_enzyme_network=0,
                                                                   node_degree_identical_enzyme_network=enzyme_efficiency_estimation_settings['node_degree_identical_enzyme_network'],
                                                                   impose_on_identical_enzymes=enzyme_efficiency_estimation_settings['impose_on_identical_enzymes'],
                                                                   condition=condition, 
                                                                   store_output=True,
                                                                   rxns_to_ignore_when_parsimonious=[],
                                                                   use_bm_flux_of_one=True,
                                                                   output_dir=output_dir)

    Specific_Kapps=Specific_Kapps_Results["Overview"]
    if min_kapp is not None:
        Specific_Kapps.loc[(Specific_Kapps['Kapp']<min_kapp)&(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp']=min_kapp
    Spec_Kapp_estim_FD=Specific_Kapps_Results["Flux_Distribution"]
    fba_flux_directions={}
    if correction_settings['impose_directions_from_fba_during_correction']:
        for fba_rxn in list(Spec_Kapp_estim_FD.index):
            if fba_rxn not in list(flux_bounds_fba.index):
                flux_value=Spec_Kapp_estim_FD.loc[fba_rxn,'FluxValues']
                if fba_rxn in list(rba_session.get_reactions()):
                    if flux_value<0:
                        fba_flux_directions.update({fba_rxn:{"LB":numpy.nan,"UB":0}})
                    elif flux_value>0:
                        fba_flux_directions.update({fba_rxn:{"LB":0,"UB":numpy.nan}})
                    #elif flux_value==0:
                    #    fba_flux_directions.update({fba_rxn:{"LB":0,"UB":0}})
                    for iso_rxn in rba_session.get_reaction_information(fba_rxn)['Twins']:
                        if flux_value<0:
                            fba_flux_directions.update({iso_rxn:{"LB":numpy.nan,"UB":0}})
                        elif flux_value>0:
                            fba_flux_directions.update({iso_rxn:{"LB":0,"UB":numpy.nan}})
                        #elif flux_value==0:
                        #    fba_flux_directions.update({fba_rxn:{"LB":0,"UB":0}})

    Specific_Kapps.to_csv(output_dir+'/Specific_Kapps_Hackett__{}.csv'.format(condition), sep=";", decimal=",")

        ## NEW
        #if enzyme_efficiency_estimation_settings['use_target_biomass_function']:
         #   rba_session.reload_model()
          #  rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

    if spec_kapps is not None:
        Specific_Kapps=import_specific_enzyme_efficiencies(input_data=spec_kapps,rba_session=rba_session,condition=condition)

    if default_kapps is None:
        spec_kapp_median=Specific_Kapps.loc[(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp'].median()
        Default_Kapps={"default_efficiency":spec_kapp_median,"default_transporter_efficiency":transporter_multiplier*spec_kapp_median}
    else:
        Default_Kapps=import_default_enzyme_efficiencies(input_data=default_kapps,condition=condition,default_transporter_kapp_coefficient=transporter_multiplier)

    flux_bounds_data=flux_bounds_from_input(input=definition_file,rba_session=rba_session, condition=condition, specific_exchanges=None, specific_directions=None,also_consider_iso_enzmes=True)
    Exchanges_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in list(flux_bounds_data["Reaction_ID"])}
    if correction_settings['impose_directions_from_fba_during_correction']:
        for rxn in fba_flux_directions.keys():
            if rxn not in Exchanges_to_impose.keys():
                Exchanges_to_impose.update({rxn:fba_flux_directions[rxn]})
        Exchanges_to_impose.update(fba_flux_directions)

    Specific_Kapps_original=Specific_Kapps.copy()
    Default_Kapps_original=Default_Kapps.copy()
    process_efficiencies_original=process_efficiencies.copy()

    condition_to_look_up="Prokaryotic"
    growth_rate_to_look_up="Mu_prok"
    results_to_look_up="Simulation_Results"

    if correction_settings['correct_efficiencies']:
        #rba_session.Problem.set_constraint_types(dict(zip(Simulation.get_enzyme_constraints(),['E']*len(Simulation.get_enzyme_constraints()))))

        #condition_to_look_up="Eukaryotic_fixed_sizes"
        #growth_rate_to_look_up="Mu_euk_fixed"
        #results_to_look_up="Simulation_Results_Euk_fixed"

        steady_count=0
        increasing_RSS_count=0
        iteration_count=0
        continuation_criterion_correction=True
        previous_RSS=numpy.nan
        efficiencies_over_correction_iterations=[]
        rss_trajectory=[]
        minimum_iteration_number=correction_settings['minimum_iteration_number']
        iteration_limit=correction_settings['iteration_limit']
        steady_limit=correction_settings['steady_rss_limit']
        increasing_RSS_limit=correction_settings['increasing_rss_limit']
        rss_tolerance=correction_settings['rss_tolerance']
        increasing_RSS_factor=correction_settings['increasing_rss_factor']

        #proteomics_data_df=build_input_proteome_for_specific_kapp_estimation(proteome, condition)
        #proteomics_data={i:proteomics_data_df.loc[i,'copy_number'] for i in proteomics_data_df.index}
        proteomics_data={}

        while continuation_criterion_correction:
            iteration_count+=1
            ### GLOBAL SCALING
            
            results_global_scaling=global_efficiency_scaling(condition=condition,
                                                                definition_file=definition_file,
                                                                rba_session=rba_session,
                                                                compartment_densities_and_pg=compartment_densities_and_PGs,
                                                                process_efficiencies=process_efficiencies,
                                                                default_kapps=Default_Kapps,
                                                                specific_kapps=Specific_Kapps,
                                                                exchanges_to_impose=Exchanges_to_impose,
                                                                feasible_stati=feasible_stati,
                                                                transporter_multiplier=transporter_multiplier,
                                                                mu_approx_precision=Mu_approx_precision,
                                                                mu_misprediction_tolerance=correction_settings['tolerance_global_scaling'],
                                                                condition_to_look_up=condition_to_look_up,
                                                                growth_rate_to_look_up=growth_rate_to_look_up,
                                                                results_to_look_up=results_to_look_up,
                                                                fixed_mu_when_above_target_mu_in_correction=correction_settings['fixed_growth_rate_global_scaling'],
                                                                n_th_root_mispred=1,
                                                                print_outputs=False,
                                                                adjust_root=correction_settings['abjust_root_of_correction_coeffs_global_scaling'],
                                                                proteomics_constraints_input=proteomics_data)

            if results_global_scaling is None:
                break
            Simulation_results=results_global_scaling["simulation_results"]
            Specific_Kapps=results_global_scaling["specific_kapps"]
            Default_Kapps=results_global_scaling["default_kapps"]
            process_efficiencies=results_global_scaling["process_efficiencies"]

            if len(list(Simulation_results[results_to_look_up].keys()))!=0:
                efficiencies_over_correction_iterations.append({"Specific_Kapps":Specific_Kapps.copy(),"Default_Kapps":Default_Kapps.copy(),"Process_Efficiencies":process_efficiencies.copy()})
                KappCorrectionResults=efficiency_correction(enzyme_efficiencies=Specific_Kapps,
                                                                simulation_results=Simulation_results[results_to_look_up],
                                                                protein_data=build_input_proteome_for_specific_kapp_estimation(proteome, condition),
                                                                rba_session=rba_session,
                                                                condition_to_look_up=condition_to_look_up,
                                                                default_enzyme_efficiencies=Default_Kapps,
                                                                tolerance=None,
                                                                #tolerance=1.05,
                                                                max_kapp=999999999,
                                                                #max_kapp=None,
                                                                n_th_root_mispred=1,
                                                                process_efficiencies=process_efficiencies,
                                                                correct_default_kapp_enzymes=True,
                                                                only_consider_misprediction_for_predicted_nonzero_enzymes=True)
                current_RSS=KappCorrectionResults["Sum_of_squared_residuals"]

                rss_trajectory.append(current_RSS)

                Specific_Kapps=KappCorrectionResults["Kapps"]
                if min_kapp is not None:
                    Specific_Kapps.loc[(Specific_Kapps['Kapp']<min_kapp)&(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp']=min_kapp
                process_efficiencies=KappCorrectionResults["ProcessEfficiencies"]
                ###
                spec_kapp_median=Specific_Kapps.loc[(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp'].median()
                Default_Kapps={"default_efficiency":spec_kapp_median,"default_transporter_efficiency":transporter_multiplier*spec_kapp_median}

                if iteration_count>=minimum_iteration_number:
                    if not pandas.isna(previous_RSS):
                        if not pandas.isna(current_RSS):
                            if (1-rss_tolerance)<=current_RSS/previous_RSS<=(1+rss_tolerance):
                                steady_count+=1
                            else:
                                steady_count=0
                                #
                                if current_RSS>=increasing_RSS_factor*previous_RSS:
                                    increasing_RSS_count+=1
                                else:
                                    increasing_RSS_count=0
                                #
                if print_outputs:
                    print("{} - {} - RSS:{} - Relative change:{} - n_inc:{} - n_steady:{}".format(condition,iteration_count,current_RSS,current_RSS/previous_RSS,increasing_RSS_count,steady_count))
                    #try:
                    #    print([condition,Specific_Kapps.loc["R_ALCD2ir_duplicate_2",'Kapp']])
                    #except:
                    #    pass
                previous_RSS=current_RSS

                if current_RSS>rss_trajectory[0]:
                    continuation_criterion_correction=False
                if steady_count>=steady_limit:
                    continuation_criterion_correction=False
                elif iteration_count>=iteration_limit:
                    continuation_criterion_correction=False
                elif increasing_RSS_count>=increasing_RSS_limit:
                    continuation_criterion_correction=False
            else:
                if iteration_count>=iteration_limit:
                    continuation_criterion_correction=False
        #
        if len(rss_trajectory)>0:
            ### add max condition here ###
            if max_kapp_threshold is not None:
                feasible_iterations=[]
                for i in range(len(rss_trajectory)):
                    max_kapp_in_iteration=efficiencies_over_correction_iterations[i]["Specific_Kapps"]["Kapp"].max(skipna=True)
                    if max_kapp_in_iteration <= max_kapp_threshold:
                        feasible_iterations.append(i)
                filtered_rss_trajectory=[rss_trajectory[i] for i in feasible_iterations]
                lowest_RSS_index=feasible_iterations[filtered_rss_trajectory.index(min(filtered_rss_trajectory))]
            ###
            else:
                lowest_RSS_index=rss_trajectory.index(min(rss_trajectory))
            Default_Kapps_to_return=efficiencies_over_correction_iterations[lowest_RSS_index]["Default_Kapps"]
            Specific_Kapps_to_return=efficiencies_over_correction_iterations[lowest_RSS_index]["Specific_Kapps"]
            process_efficiencies_to_return=efficiencies_over_correction_iterations[lowest_RSS_index]["Process_Efficiencies"]
        else:
            Default_Kapps_to_return=Default_Kapps
            Specific_Kapps_to_return=Specific_Kapps
            process_efficiencies_to_return=process_efficiencies

        if correction_settings['final_global_scaling_after_correction']:
            if correction_settings['final_global_scaling_without_imposed_exchanges']:
                Exchanges_to_impose_here=None
            else:
                Exchanges_to_impose_here=Exchanges_to_impose
            results_global_scaling=global_efficiency_scaling(condition=condition,
                                                             definition_file=definition_file,
                                                             rba_session=rba_session,
                                                             compartment_densities_and_pg=compartment_densities_and_PGs,
                                                             process_efficiencies=process_efficiencies_to_return,
                                                             default_kapps=Default_Kapps_to_return,
                                                             specific_kapps=Specific_Kapps_to_return,
                                                             exchanges_to_impose=Exchanges_to_impose_here,
                                                             feasible_stati=feasible_stati,
                                                             transporter_multiplier=transporter_multiplier,
                                                             mu_approx_precision=Mu_approx_precision,
                                                             mu_misprediction_tolerance=correction_settings['tolerance_global_scaling'],
                                                             condition_to_look_up=condition_to_look_up,
                                                             growth_rate_to_look_up=growth_rate_to_look_up,
                                                             results_to_look_up=results_to_look_up,
                                                             fixed_mu_when_above_target_mu_in_correction=False)
            
            Simulation_results=results_global_scaling["simulation_results"]
            Specific_Kapps_to_return=results_global_scaling["specific_kapps"]
            Default_Kapps_to_return=results_global_scaling["default_kapps"]
            process_efficiencies_to_return=results_global_scaling["process_efficiencies"]
    else:
        if correction_settings['final_global_scaling_without_imposed_exchanges']:
            Exchanges_to_impose_here=None
        else:
            Exchanges_to_impose_here=Exchanges_to_impose
        results_global_scaling=global_efficiency_scaling(condition=condition,
                                                         definition_file=definition_file,
                                                         rba_session=rba_session,
                                                         compartment_densities_and_pg=compartment_densities_and_PGs,
                                                         process_efficiencies=process_efficiencies,
                                                         default_kapps=Default_Kapps,
                                                         specific_kapps=Specific_Kapps,
                                                         exchanges_to_impose=Exchanges_to_impose_here,
                                                         feasible_stati=feasible_stati,
                                                         transporter_multiplier=transporter_multiplier,
                                                         mu_approx_precision=Mu_approx_precision,
                                                         mu_misprediction_tolerance=correction_settings['tolerance_global_scaling'],
                                                         condition_to_look_up=condition_to_look_up,
                                                         growth_rate_to_look_up=growth_rate_to_look_up,
                                                         results_to_look_up=results_to_look_up,
                                                         fixed_mu_when_above_target_mu_in_correction=False)
            
        Simulation_results=results_global_scaling["simulation_results"]
        Specific_Kapps_to_return=results_global_scaling["specific_kapps"]
        Default_Kapps_to_return=results_global_scaling["default_kapps"]
        process_efficiencies_to_return=results_global_scaling["process_efficiencies"]
        rss_trajectory=None
        lowest_RSS_index=None

    if use_mean_enzyme_composition_for_calibration:
        rba_session.reload_model()
        for kapp_reaction in Specific_Kapps_to_return.index:
            respective_enzyme=rba_session.get_reaction_information(kapp_reaction)["Enzyme"]
            for iso_enzyme in rba_session.get_enzyme_information(respective_enzyme)['Isozymes']:
                respective_iso_reaction=rba_session.get_enzyme_information(iso_enzyme)['Reaction']
                if respective_iso_reaction not in Specific_Kapps_to_return.index:
                    Specific_Kapps_to_return.loc[respective_iso_reaction,"Enzyme_ID"]=iso_enzyme
                    Specific_Kapps_to_return.loc[respective_iso_reaction, 'Comment'] = 'isoenzyme'
                    Specific_Kapps_to_return.loc[respective_iso_reaction, 'Taken from'] = kapp_reaction
                    Specific_Kapps_to_return.loc[respective_iso_reaction, 'Flux'] = Specific_Kapps_to_return.loc[kapp_reaction, 'Flux']
                    if enzyme_efficiency_estimation_settings['zero_on_all_isoreactions']:
                        Specific_Kapps_to_return.loc[respective_iso_reaction,'Kapp'] = 0.0
                    elif enzyme_efficiency_estimation_settings['impose_on_all_isoreactions']:
                        Specific_Kapps_to_return.loc[respective_iso_reaction,'Kapp'] = Specific_Kapps_to_return.loc[kapp_reaction, 'Kapp']

    if print_outputs:
        print("")
        print("Runtime - {} : {}".format(condition,time.time() - t0))
        print("")

    return({"RSS_trajectory":rss_trajectory,
            "Chosen_RSS_index":lowest_RSS_index,
            "Densities_PGs":compartment_densities_and_PGs,
            "Condition":condition,
            'Proteome': build_input_proteome_for_specific_kapp_estimation(proteome, condition),
            'correction_results_compartement_sizes': correction_results_compartement_sizes,
            'Default_Kapps': Default_Kapps_to_return,
            'Specific_Kapps': Specific_Kapps_to_return,
            'Process_Efficiencies': process_efficiencies_to_return,
            'Default_Kapps_original': Default_Kapps_original,
            'Specific_Kapps_original': Specific_Kapps_original,
            'Process_Efficiencies_original': process_efficiencies_original})


def correction_pipeline(input,
                          condition,
                          definition_file,
                          compartments_to_replace,
                          compartments_no_original_PG,
                          fractions_entirely_replaced_with_expected_value,
                          imposed_compartment_fractions,
                          directly_corrected_compartments,
                          merged_compartments,
                          min_compartment_fraction):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    condition : _type_
        _description_
    definition_file : _type_
        _description_
    compartments_to_replace : _type_
        _description_
    compartments_no_original_PG : _type_
        _description_
    fractions_entirely_replaced_with_expected_value : _type_
        _description_
    imposed_compartment_fractions : _type_
        _description_
    directly_corrected_compartments : _type_
        _description_
    merged_compartments : _type_
        _description_
    min_compartment_fraction : _type_
        _description_
    """
    out = build_proteome_overview(input=input, condition=condition, compartments_to_replace=compartments_to_replace,
                                  compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=True)
    factor_a=1/(1-sum([out.loc["pg_{}".format(i),"original_protein_fraction"] for i in compartments_no_original_PG]))
    #factor_b=(1-sum([imposed_compartment_fractions[i] for i in fractions_entirely_replaced_with_expected_value]+[out.loc["Ribosomes","original_protein_fraction"]]))/(1-factor_a*out.loc["Ribosomes","original_protein_fraction"])
    factor_b=(1-sum([imposed_compartment_fractions[i] for i in fractions_entirely_replaced_with_expected_value+compartments_no_original_PG]))/(1-factor_a*sum([out.loc[i,"original_protein_fraction"] for i in fractions_entirely_replaced_with_expected_value]))
    factor_c=(1-sum([imposed_compartment_fractions[i] for i in fractions_entirely_replaced_with_expected_value+compartments_no_original_PG]))/(1-sum([imposed_compartment_fractions[i] for i in fractions_entirely_replaced_with_expected_value+compartments_no_original_PG])-factor_a*factor_b*sum([out.loc[i,"original_protein_fraction"] for i in compartments_no_original_PG]))
    new_total=0
    for i in out.index:
        if i not in compartments_to_replace.keys():
            if not i.startswith("pg_"):
                if not i in ["Total"]:
                    out.loc[i,"Factor_A"]=factor_a
                    if not i in fractions_entirely_replaced_with_expected_value:
                        out.loc[i,"Factor_B"]=factor_b
                        if not i in compartments_no_original_PG:
                            out.loc[i,"Factor_C"]=factor_c
                            if i in merged_compartments.keys():
                                out.loc[i,"new_protein_fraction"]=sum([imposed_compartment_fractions[j] for j in [merged_compartments[i]]])+out.loc[i,"original_protein_fraction"]*factor_a*factor_b*factor_c
                                out.loc[i,"new_PG_fraction"]=out.loc[i,"original_PG_fraction"]*out.loc[i,"original_protein_fraction"]*factor_a*factor_b*factor_c/out.loc[i,"new_protein_fraction"]
                                new_total+=out.loc[i,"new_protein_fraction"]
                            else:
                                out.loc[i,"new_protein_fraction"]=out.loc[i,"original_protein_fraction"]*factor_a*factor_b*factor_c
                                out.loc[i,"new_PG_fraction"]=out.loc[i,"original_PG_fraction"]
                                new_total+=out.loc[i,"new_protein_fraction"]
                        else:
                            out.loc[i,"new_protein_fraction"]=imposed_compartment_fractions[i]
                            out.loc[i,"new_PG_fraction"]=1-(out.loc[i,"original_protein_fraction"]*factor_a*factor_b)/out.loc[i,"new_protein_fraction"]
                            new_total+=out.loc[i,"new_protein_fraction"]
                    else:
                        out.loc[i,"new_protein_fraction"]=imposed_compartment_fractions[i]
    out.loc["Total","new_protein_fraction"]=new_total
    out["Location"]=out.index
    return(out)


def build_proteome_overview(input, condition, compartments_to_replace={'DEF':"c", 'DEFA':"c", 'Def':"c"}, compartments_no_original_PG=['n', 'Secreted'], ribosomal_proteins_as_extra_compartment=True):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    condition : _type_
        _description_
    compartments_to_replace : dict, optional
        _description_, by default {'DEF':"c", 'DEFA':"c", 'Def':"c"}
    compartments_no_original_PG : list, optional
        _description_, by default ['n', 'Secreted']
    ribosomal_proteins_as_extra_compartment : bool, optional
        _description_, by default True
    """
    out = determine_compartment_occupation(Data_input=input, Condition=condition, compartments_to_replace=compartments_to_replace,
                                           compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=ribosomal_proteins_as_extra_compartment, only_in_model=False)
    out_in_model = determine_compartment_occupation(Data_input=input, Condition=condition, compartments_to_replace=compartments_to_replace,
                                                    compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=ribosomal_proteins_as_extra_compartment, only_in_model=True)
    for comp in out.index:
        out.loc[comp,'original_protein_fraction'] = out.loc[comp,'original_amino_acid_occupation']/out.loc["Total",'original_amino_acid_occupation']
        out['original_PG_fraction'] = 1-out_in_model['original_amino_acid_occupation'] / out['original_amino_acid_occupation']
    return(out)


def determine_compartment_occupation(Data_input,
                                       Condition,
                                       mass_col='AA_residues',
                                       only_in_model=False,
                                       compartments_to_replace={'DEF':"c"},
                                       compartments_no_original_PG=[],
                                       ribosomal_proteins_as_extra_compartment=True):
    """
    _summary_

    Parameters
    ----------
    Data_input : _type_
        _description_
    Condition : _type_
        _description_
    mass_col : str, optional
        _description_, by default 'AA_residues'
    only_in_model : bool, optional
        _description_, by default False
    compartments_to_replace : dict, optional
        _description_, by default {'DEF':"c"}
    compartments_no_original_PG : list, optional
        _description_, by default []
    ribosomal_proteins_as_extra_compartment : bool, optional
        _description_, by default True
    """
    out=pandas.DataFrame()
    if only_in_model:
        Data_intermediate = Data_input.loc[Data_input['InModel'] >= 1,:]
    else:
        Data_intermediate=Data_input
    Data=pandas.DataFrame()
    for i in Data_intermediate.index:
        for j in Data_intermediate.columns:
            if j != "AA_abundance":
                Data.loc[i,j]=Data_intermediate.loc[i,j]
        Data.loc[i,"AA_abundance"]=Data_intermediate.loc[i,Condition]*Data_intermediate.loc[i,mass_col]

    out.loc["Total","original_amino_acid_occupation"]=sum([i for i in list(Data["AA_abundance"]) if not pandas.isna(i)])
    for i in compartments_to_replace.keys():
        Data.loc[Data['Location'] == i,'Location']=compartments_to_replace[i]
    for i in compartments_no_original_PG:
        intermediate_Data_pg= Data.loc[(Data['Location'] == i) & (Data['InModel'] == 0)].copy()
        out.loc["pg_{}".format(i),"original_amino_acid_occupation"]=sum([i for i in list(intermediate_Data_pg["AA_abundance"]) if not pandas.isna(i)])
        intermediate_Data_non_pg= Data.loc[(Data['Location'] == i) & (Data['InModel'] == 1)].copy()
        out.loc[i,"original_amino_acid_occupation"]=sum([i for i in list(intermediate_Data_non_pg["AA_abundance"]) if not pandas.isna(i)])
    if ribosomal_proteins_as_extra_compartment:
        Data_R = Data.loc[Data['IsRibosomal'] == 1].copy()
        out.loc['Ribosomes', "original_amino_acid_occupation"] =sum([i for i in list(Data_R["AA_abundance"]) if not pandas.isna(i)])
    for comp in list(set(list(Data['Location']))):
        if comp in out.index:
            continue
        if ribosomal_proteins_as_extra_compartment:
            intermediate_Data= Data.loc[(Data['Location'] == comp) & (Data['IsRibosomal'] == 0)].copy()
        else:
            intermediate_Data= Data.loc[Data['Location'] == comp].copy()
        out.loc[comp, "original_amino_acid_occupation"] =sum([i for i in list(intermediate_Data["AA_abundance"]) if not pandas.isna(i)])
    return(out)


def extract_compsizes_and_pgfractions_from_correction_summary(corrsummary,rows_to_exclude):
    """
    _summary_

    Parameters
    ----------
    corrsummary : _type_
        _description_
    rows_to_exclude : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in list(corrsummary.index):
        if i in rows_to_exclude:
            continue
        else:
            out.loc[i,"Compartment_ID"]=i
            out.loc[i,"Density"]=corrsummary.loc[i,"new_protein_fraction"]
            out.loc[i,"PG_fraction"]=corrsummary.loc[i,"new_PG_fraction"]
    return(out)


def determine_apparent_process_efficiencies(growth_rate, input, rba_session, proteome_summary, protein_data, condition, gene_id_col,fit_nucleotide_assembly_machinery=False):
    """
    _summary_

    Parameters
    ----------
    growth_rate : _type_
        _description_
    input : pandas.DataFrame
        _description_
    rba_session : _type_
        _description_
    proteome_summary : _type_
        _description_
    protein_data : _type_
        _description_
    condition : _type_
        _description_
    gene_id_col : _type_
        _description_
    fit_nucleotide_assembly_machinery : bool, optional
        _description_, by default False
    """
    process_efficiencies = pandas.DataFrame()
    for i in input.index:
        process_ID = input.loc[i, 'Process_ID']
        process_name = input.loc[i, 'Process_Name']
        process_client_compartments = input.loc[i, 'Client_Compartments'].split(' , ')
        constituting_proteins = {rba_session.get_protein_information(protein=i)['ProtoID']: rba_session.get_protein_information(protein=i)['AAnumber'] for i in rba_session.get_process_information(process=process_name)['Composition'].keys()}
        Total_client_fraction = sum([proteome_summary.loc[i, 'new_protein_fraction']
                                     for i in process_client_compartments])
        n_AAs_in_machinery = 0
        machinery_size = 0
        for i in constituting_proteins.keys():
            if i in protein_data['ID']:
                n_AAs_in_machinery += protein_data.loc[protein_data['ID'] == i, condition].values[0] * \
                    protein_data.loc[protein_data['ID'] == i, 'AA_residues'].values[0]
                machinery_size += constituting_proteins[i]
        # right reference amounth?
        if n_AAs_in_machinery > 0:
            relative_Protein_fraction_of_machinery = n_AAs_in_machinery / proteome_summary.loc['Total', 'original_amino_acid_occupation']
            specific_capacity = growth_rate*Total_client_fraction/relative_Protein_fraction_of_machinery
            apparent_capacity = specific_capacity*machinery_size
            process_efficiencies.loc[process_name, 'Process'] = process_ID
            process_efficiencies.loc[process_name, 'Parameter'] = str(process_ID+'_apparent_efficiency')
            process_efficiencies.loc[process_name, 'Value'] = apparent_capacity
    median_process_efficiency=numpy.median(numpy.array(process_efficiencies["Value"]))
    for i in input.index:
        process_ID = input.loc[i, 'Process_ID']
        process_name = input.loc[i, 'Process_Name']
        if not process_name in process_efficiencies.index:
            process_efficiencies.loc[process_name, 'Process'] = process_ID
            process_efficiencies.loc[process_name, 'Parameter'] = str(process_ID+'_apparent_efficiency')
            process_efficiencies.loc[process_name, 'Value'] = median_process_efficiency
    if fit_nucleotide_assembly_machinery:
        original_Mu=rba_session.Mu
        rba_session.set_growth_rate(growth_rate)
        machinery_production_fluxes=determine_macromolecule_synthesis_machinery_demand(rba_session)
        for machinery in machinery_production_fluxes.keys():
            process_info=rba_session.get_process_information(process=machinery)
            stoichiometrically_scaled_subunit_concentrations=[]
            for su in process_info["Composition"].keys():
                protoProteinID=rba_session.get_protein_information(protein=su)["ProtoID"]
                if protoProteinID in protein_data['ID']:
                    copy_number=protein_data.loc[protein_data['ID'] == protoProteinID, condition].values[0]
                    scaled_copy_number=copy_number/process_info["Composition"][su]
                    if not pandas.isna(scaled_copy_number):
                        stoichiometrically_scaled_subunit_concentrations.append(scaled_copy_number)
            if len(stoichiometrically_scaled_subunit_concentrations)>0:
                machinery_concentration=gmean(stoichiometrically_scaled_subunit_concentrations)
                apparent_process_efficiency=machinery_production_fluxes[machinery]/machinery_concentration
                process_efficiencies.loc[machinery, 'Process'] = process_info["ID"]
                process_efficiencies.loc[machinery, 'Parameter'] = str( process_info["ID"]+'_apparent_efficiency')
                process_efficiencies.loc[machinery, 'Value'] = apparent_process_efficiency
        rba_session.set_growth_rate(original_Mu)
    return(process_efficiencies)


def determine_macromolecule_synthesis_machinery_demand(rba_session):
    """
    _summary_

    Parameters
    ----------
    rba_session : _type_
        _description_
    """
    production_fluxes_macromolecules={}
    for i in rba_session.model.targets.target_groups._elements_by_id["transcription_targets"].concentrations._elements:
        species=i.species
        if "value" in i.__dict__.keys():
            if not species in production_fluxes_macromolecules.keys():
                production_fluxes_macromolecules[species]=0
            parameter=i.value
            Monomoers_per_molecule=sum(list(rba_session.get_macro_molecule_information(macro_molecule=species)["Composition"].values()))
            current_parameter_value=rba_session.get_current_parameter_value(parameter=parameter)
            production_fluxes_macromolecules[species]+=Monomoers_per_molecule*current_parameter_value*rba_session.Mu
    for i in rba_session.model.targets.target_groups._elements_by_id["transcription_targets"].production_fluxes._elements:
        species=i.species
        if "value" in i.__dict__.keys():
            if not species in production_fluxes_macromolecules.keys():
                production_fluxes_macromolecules[species]=0
            parameter=i.value
            Monomoers_per_molecule=sum(list(rba_session.get_macro_molecule_information(macro_molecule=species)["Composition"].values()))
            current_parameter_value=rba_session.get_current_parameter_value(parameter=parameter)
            production_fluxes_macromolecules[species]+=Monomoers_per_molecule*current_parameter_value
    for i in rba_session.model.targets.target_groups._elements_by_id["replication_targets"].concentrations._elements:
        species=i.species
        if "value" in i.__dict__.keys():
            if not species in production_fluxes_macromolecules.keys():
                production_fluxes_macromolecules[species]=0
            parameter=i.value
            Monomoers_per_molecule=sum(list(rba_session.get_macro_molecule_information(macro_molecule=species)["Composition"].values()))
            current_parameter_value=rba_session.get_current_parameter_value(parameter=parameter)
            production_fluxes_macromolecules[species]+=Monomoers_per_molecule*current_parameter_value*rba_session.Mu
    for i in rba_session.model.targets.target_groups._elements_by_id["replication_targets"].production_fluxes._elements:
        species=i.species
        if "value" in i.__dict__.keys():
            if not species in production_fluxes_macromolecules.keys():
                production_fluxes_macromolecules[species]=0
            parameter=i.value
            Monomoers_per_molecule=sum(list(rba_session.get_macro_molecule_information(macro_molecule=species)["Composition"].values()))
            current_parameter_value=rba_session.get_current_parameter_value(parameter=parameter)
            production_fluxes_macromolecules[species]+=Monomoers_per_molecule*current_parameter_value
    process_demand_monomer_fluxes={}
    for prod_flux in production_fluxes_macromolecules.keys():
        for required_process in rba_session.get_macro_molecule_information(macro_molecule=prod_flux)["ProcessRequirements"].keys():
            if not required_process in process_demand_monomer_fluxes.keys():
                process_demand_monomer_fluxes[required_process]=0
            process_demand_monomer_fluxes[required_process]+=production_fluxes_macromolecules[prod_flux]
    return(process_demand_monomer_fluxes)


def build_input_proteome_for_specific_kapp_estimation(proteomics_data, condition):
    """
    _summary_

    Parameters
    ----------
    proteomics_data : _type_
        _description_
    condition : _type_
        _description_
    """
    out = pandas.DataFrame()
    out['ID'] = proteomics_data['ID']
    out['copy_number'] = proteomics_data[condition]
    return(out)


def estimate_specific_enzyme_efficiencies(rba_session, 
                                          proteomicsData, 
                                          flux_bounds, 
                                          compartment_densities_and_pg,
                                          mu, 
                                          biomass_function=None, 
                                          target_biomass_function=True, 
                                          parsimonious_fba=True, 
                                          chose_most_likely_isoreactions=False,
                                          impose_on_all_isoreactions=True, 
                                          zero_on_all_isoreactions=False,
                                          impose_on_identical_enzymes=True,
                                          node_degree_identical_enzyme_network=1,
                                          condition=None, 
                                          store_output=True,
                                          rxns_to_ignore_when_parsimonious=[],
                                          use_bm_flux_of_one=False,
                                          output_dir=""):
    """
    Parameters
    ----------
    proteomicsData : pandas.DataFrame (in mmol/gDW)
    flux_bounds : pandas.DataFrame  (in mmol/(gDW*h))
    mu : float (in 1/h)
    biomass_function : str
    target_biomass_function : bool
    atp_maintenance_to_biomassfunction : bool
    eukaryotic : bool
    """

    #####
    # 1: Determine Flux Distribution from parsimonious FBA#

    FluxDistribution=determine_calibration_flux_distribution(rba_session=rba_session,
                                                             mu=mu,
                                                             flux_bounds=flux_bounds,
                                                             compartment_densities_and_pg=compartment_densities_and_pg,
                                                             biomass_function=biomass_function,
                                                             target_biomass_function=target_biomass_function,
                                                             parsimonious_fba=parsimonious_fba,
                                                             rxns_to_ignore_when_parsimonious=rxns_to_ignore_when_parsimonious,
                                                             condition=condition,
                                                             use_bm_flux_of_one=use_bm_flux_of_one
                                                             )
    FluxDistribution.to_csv(output_dir+'/Calib_FluxDist_'+condition+'_.csv', sep=';')

    if chose_most_likely_isoreactions:
        ### INCORPORATE ANA´s REMARK ON NOT ELIMINATING ISOENZYMES WITH SUs ONLY PRESENT IN THIS ISOENZYMES ###
        # 2: Determine list of all pre_selected isoenzymes --> "pre_selected_enzymes"#
        ProtoProteinMap = rba_session.ModelStructure.ProteinInfo.return_protein_iso_form_map()
        measured_proteins_isoform_map = {p_ID: ProtoProteinMap[p_ID] for p_ID in list(proteomicsData['ID']) if p_ID in list(ProtoProteinMap.keys())}

        # identify all model reactions, associated with the measured proteins
        measured_proteins_reaction_map = determine_reactions_associated_with_measured_proto_protein(measured_proteins_isoform_map=measured_proteins_isoform_map,
                                                                                                        rba_session=rba_session)

        chosen_isoreactions=pre_select_iso_reactions(measured_proteins_reaction_map=measured_proteins_reaction_map,
                                                     rba_session=rba_session,
                                                     chose_most_quantified=True,
                                                     keep_isorxns_specific_to_quantified_proteins=True)
        pre_selected_enzymes=[]
        for proto_rxn in chosen_isoreactions.keys():
            for iso_rxn in chosen_isoreactions[proto_rxn]:
                if rba_session.get_reaction_information(iso_rxn) is not None:
                    respective_enzyme=rba_session.get_reaction_information(iso_rxn)["Enzyme"]
                    if respective_enzyme not in pre_selected_enzymes:
                        pre_selected_enzymes.append(respective_enzyme)

    else:
        # 2: Determine list of all nonzero/quantified model enzymes from proteomics data --> "pre_selected_enzymes"#
        pre_selected_enzymes=[]
        for model_enzyme in rba_session.get_enzymes():
            if rba_session.get_enzyme_information(model_enzyme) is not None:
                model_enzyme_composition=rba_session.get_enzyme_information(model_enzyme)["Subunits"]
                #model_enzyme_concentration=determine_machinery_concentration(rba_session=rba_session,
                #                                                             machinery_composition=model_enzyme_composition,
                #                                                             proteomicsData=proteomicsData,
                #                                                             proto_proteins=False)
                model_enzyme_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                             machinery_composition=model_enzyme_composition,
                                                                             proteomicsData=proteomicsData,
                                                                             proto_proteins=False)
                
                if model_enzyme_concentration!=0:
                    pre_selected_enzymes.append(model_enzyme)

    # 3: Determine list of all reactions with nonzero/quantified flux in FBA --> "nonzero_reactions"#
    nonzero_flux_reactions=[]
    for reaction in FluxDistribution.index:
        if numpy.isfinite(FluxDistribution.loc[reaction,'FluxValues']):
            if FluxDistribution.loc[reaction,'FluxValues']!=0:
                nonzero_flux_reactions.append(reaction)

    # 4: Generate dictionary with all model_enzymes with nonzero concentration in data 
    # and their associated nonzer fba-flux proto reaction#
    nonzero_concentration_enzymes_with_associated_fba_flux={}
    for enzyme in pre_selected_enzymes:
        associated_reaction=rba_session.get_enzyme_information(enzyme)["Reaction"]
        associated_fba_flux_carrying_reaction=determine_fba_flux_carrying_isoreaction(rba_session=rba_session,
                                                                                      reaction_id=associated_reaction,
                                                                                      flux_distribution=FluxDistribution)
        if associated_fba_flux_carrying_reaction is not None:
            nonzero_concentration_enzymes_with_associated_fba_flux[enzyme]=associated_fba_flux_carrying_reaction

    # 5: ...#
    overview_out = pandas.DataFrame()
    for rxn in list(set(list(nonzero_concentration_enzymes_with_associated_fba_flux.values()))):
        # 5.1: Determine model_enzymes constituting pseudo complex, associated with reaction#
        if rba_session.get_reaction_information(rxn) is None:
            continue
        rxn_enzyme=rba_session.get_reaction_information(rxn)["Enzyme"]
        if rba_session.get_enzyme_information(rxn_enzyme) is None:
            continue

        associated_nonzero_concentration_isoenzymatic_complexes=[i for i in list([rxn_enzyme]+rba_session.get_enzyme_information(rxn_enzyme)["Isozymes"]) if i in list(pre_selected_enzymes)]
        model_enzymes_associated_with_pseudo_complex={}
        for isoenzymatic_complex in associated_nonzero_concentration_isoenzymatic_complexes:
            model_enzymes_associated_with_pseudo_complex[isoenzymatic_complex]=[isoenzymatic_complex]
            if node_degree_identical_enzyme_network==1:
                ### ??? identical enzymes isoform specific? ###
                for identical_composition_enzyme in rba_session.get_enzyme_information(rxn_enzyme)["EnzymesWithIdenticalSubunitComposition"]:
                    if identical_composition_enzyme in nonzero_concentration_enzymes_with_associated_fba_flux.keys():
                        model_enzymes_associated_with_pseudo_complex[isoenzymatic_complex].append(identical_composition_enzyme)
                    ### potentially add secondary edges here
        # 5.2 ...#
        total_flux_dict={}
        individual_constituent_concentrations={i:[] for i in model_enzymes_associated_with_pseudo_complex.keys()}
        for associated_pseudocomplex_enzyme in model_enzymes_associated_with_pseudo_complex.keys():
            respective_composition=rba_session.get_enzyme_information(associated_pseudocomplex_enzyme)["Subunits"]
            respective_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                             machinery_composition=respective_composition,
                                                                             proteomicsData=proteomicsData,
                                                                             proto_proteins=False)
            if (respective_concentration!=0) and (numpy.isfinite(respective_concentration)):
                individual_constituent_concentrations[associated_pseudocomplex_enzyme].append(respective_concentration)
            for identical_composition_enzyme in model_enzymes_associated_with_pseudo_complex[associated_pseudocomplex_enzyme]:
                if identical_composition_enzyme in nonzero_concentration_enzymes_with_associated_fba_flux.keys():
                    associated_fba_rxn=nonzero_concentration_enzymes_with_associated_fba_flux[identical_composition_enzyme]
                    if associated_fba_rxn not in total_flux_dict.keys():
                        total_flux_dict[associated_fba_rxn]=abs(FluxDistribution.loc[associated_fba_rxn,'FluxValues'])

        pseudocomplex_constituent_concentrations=[individual_constituent_concentrations[i][0] for i in individual_constituent_concentrations.keys() if len(individual_constituent_concentrations[i])!=0]
        if len(pseudocomplex_constituent_concentrations)>0:
            concentration_pseudo_complex=len(list(total_flux_dict.keys()))*gmean(numpy.array(pseudocomplex_constituent_concentrations))
        else:
            concentration_pseudo_complex=numpy.nan

        # 5.4 ...#
        total_flux_pseudo_complex=gmean(numpy.array(list(total_flux_dict.values())))
        total_flux_pseudo_complex=sum(total_flux_dict.values())
        # 5.5 ...#
        flux_direction=numpy.nan
        if FluxDistribution.loc[rxn,'FluxValues'] > 0.0:
            flux_direction=1
        elif FluxDistribution.loc[rxn,'FluxValues'] < 0.0:
            flux_direction=-1
        if (total_flux_pseudo_complex!=0) and numpy.isfinite(total_flux_pseudo_complex) and (concentration_pseudo_complex!=0) and numpy.isfinite(concentration_pseudo_complex):
            kapp_pseudo_complex=total_flux_pseudo_complex/concentration_pseudo_complex
            # 5.6 ...#
            for considered_isoenzyme in individual_constituent_concentrations.keys():
                respective_reaction=rba_session.get_enzyme_information(considered_isoenzyme)["Reaction"]
                overview_out.loc[respective_reaction,"Enzyme_ID"]=considered_isoenzyme
                overview_out.loc[respective_reaction, 'Flux_FBA'] = total_flux_pseudo_complex
                overview_out.loc[respective_reaction, 'Flux'] = flux_direction
                overview_out.loc[respective_reaction, 'Concentration'] = concentration_pseudo_complex
                overview_out.loc[respective_reaction,'Kapp'] = kapp_pseudo_complex
                overview_out.loc[respective_reaction, 'Comment'] = 'estimated'
                overview_out.loc[respective_reaction,"Pseudo complex members"]=" , ".join(model_enzymes_associated_with_pseudo_complex)
            # 5.7 ...#
            for isozyme in list(rba_session.get_enzyme_information(rba_session.get_reaction_information(rxn)["Enzyme"])["Isozymes"])+[rba_session.get_reaction_information(rxn)["Enzyme"]]:
                respective_reaction=rba_session.get_enzyme_information(isozyme)["Reaction"]
                if respective_reaction not in overview_out.index:
                    overview_out.loc[respective_reaction,"Enzyme_ID"]=isozyme
                    overview_out.loc[respective_reaction, 'Comment'] = 'isoenzyme'
                    overview_out.loc[respective_reaction, 'Taken from'] = rxn
                    overview_out.loc[respective_reaction, 'Flux'] = flux_direction
                    if zero_on_all_isoreactions:
                        overview_out.loc[respective_reaction,'Kapp'] = 0.0
                    elif impose_on_all_isoreactions:
                        overview_out.loc[respective_reaction,'Kapp'] = kapp_pseudo_complex
                        
    # 6: ...#
    if impose_on_identical_enzymes:
        for rxn in rba_session.get_reactions():
            if rxn not in overview_out.index:
                #
                iso_rxn_already_in_kapp_output=False
                for iso_rxn in rba_session.get_reaction_information(rxn)["Twins"]:
                    if iso_rxn in overview_out.index:
                        iso_rxn_already_in_kapp_output=True
                if iso_rxn_already_in_kapp_output:
                    continue
                #
                respective_enzyme=rba_session.get_reaction_information(rxn)["Enzyme"]
                if respective_enzyme in rba_session.get_enzymes():
                    identical_enzymes=rba_session.get_enzyme_information(respective_enzyme)["EnzymesWithIdenticalSubunitComposition"]
                    if identical_enzymes:
                        kapps_to_average={}
                        for ident_enz in identical_enzymes:
                            respective_reaction=rba_session.get_enzyme_information(ident_enz)["Reaction"]
                            if respective_reaction in overview_out.index:
                                if zero_on_all_isoreactions:
                                    applicable_comments=["estimated"]
                                else:
                                    applicable_comments=["estimated","isoenzyme"]
                                ident_comment=overview_out.loc[respective_reaction,'Comment']
                                if ident_comment in applicable_comments:
                                    ident_kapp=overview_out.loc[respective_reaction,'Kapp']
                                    if not pandas.isna(ident_kapp):
                                        kapps_to_average[ident_enz]=ident_kapp
                        if list(kapps_to_average.values()):
                            overview_out.loc[rxn,"Enzyme_ID"]=respective_enzyme
                            overview_out.loc[rxn,'Kapp'] = numpy.mean(numpy.array(list(kapps_to_average.values())))
                            overview_out.loc[rxn, 'Comment'] = 'twinenzyme'
                            overview_out.loc[rxn, 'Taken from'] = " , ".join(list(kapps_to_average.keys()))
    # 7: ...#
    rba_session.rebuild_from_model()
    rba_session.set_medium(rba_session.Medium)

    overview_out.sort_index(inplace=True)

    # 8: ...#
    if store_output:
        if condition is not None:
            overview_out.to_csv(output_dir+'/SpecKapp_Network_overview_'+condition+'_.csv', sep=';')
        else:
            overview_out.to_csv(output_dir+'/SpecKapp_Network_overview_.csv', sep=';')
    # 9: ...#
    return({"Overview":overview_out,"Flux_Distribution":FluxDistribution})


def determine_calibration_flux_distribution(rba_session,
                                            mu,
                                            flux_bounds,
                                            compartment_densities_and_pg,
                                            biomass_function,
                                            target_biomass_function,
                                            parsimonious_fba,
                                            rxns_to_ignore_when_parsimonious,
                                            condition=None,
                                            use_bm_flux_of_one=False
                                            ):                                        
    """
    _summary_

    Parameters
    ----------
    rba_session : _type_
        _description_
    mu : _type_
        _description_
    flux_bounds : _type_
        _description_
    compartment_densities_and_pg : _type_
        _description_
    biomass_function : _type_
        _description_
    target_biomass_function : _type_
        _description_
    parsimonious_fba : _type_
        _description_
    rxns_to_ignore_when_parsimonious : _type_
        _description_
    condition : _type_, optional
        _description_, by default None
    use_bm_flux_of_one : bool, optional
        _description_, by default False
    """
    for comp in list(compartment_densities_and_pg['Compartment_ID']):
        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_pg.loc[compartment_densities_and_pg['Compartment_ID'] == comp, 'Density'].values[0]
        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = 0.0

    rba_session.rebuild_from_model()
    rba_session.set_medium(rba_session.Medium)
    rba_session.add_exchange_reactions()
    rba_session.set_growth_rate(mu)

    if target_biomass_function:
        derive_bm_from_rbasolution=False
        derive_bm_from_targets=True
        original_medium = copy.deepcopy(rba_session.Medium)
        rba_session.set_medium({i:100.0 for i in original_medium.keys()})
        original_density_constraint_signs=rba_session.Problem.get_constraint_types(constraints=[i for i in rba_session.get_density_constraints() if i in rba_session.Problem.LP.row_names])
        rba_session.Problem.set_constraint_types({i:"E" for i in rba_session.get_density_constraints() if i in rba_session.Problem.LP.row_names})
        solved=rba_session.solve()
        if solved:
            derive_bm_from_rbasolution=True
            derive_bm_from_targets=False
            rba_session.Problem.set_constraint_types(original_density_constraint_signs)
        else:
            print("{} - Solution with equality density not obtained - Status: {}".format(condition,rba_session.Problem.SolutionStatus))
            rba_session.Problem.set_constraint_types(original_density_constraint_signs)
            rba_session.set_growth_rate(mu)
            solved2=rba_session.solve()
            if solved2:
                derive_bm_from_rbasolution=True
                derive_bm_from_targets=False
            else:
                print("{} - Solution with inequality density not obtained - Status: {}".format(condition,rba_session.Problem.SolutionStatus))

        rba_session.set_medium(original_medium)
        rba_session.build_fba_model(rba_derived_biomass_function=True,
                                    from_rba_solution=derive_bm_from_rbasolution,
                                    from_targets=derive_bm_from_targets)
        BMfunction = 'R_BIOMASS_targetsRBA'
    else:
        rba_session.build_fba_model(rba_derived_biomass_function=False)
        BMfunction = biomass_function

    for j in [i for i in rba_session.Medium.keys() if rba_session.Medium[i] == 0]:
        Exrxn = 'R_EX_'+j.split('M_')[-1]+'_e'
        if Exrxn in list(rba_session.FBA.LP.col_names):
            rba_session.FBA.set_ub({Exrxn: 0})

    rxn_LBs = {}
    rxn_UBs = {}
    for rx in flux_bounds['Reaction_ID']:
        lb = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'LB'].values[0]
        ub = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'UB'].values[0]
        if not pandas.isna(lb):
            rxn_LBs.update({rx: lb})
        if not pandas.isna(ub):
            rxn_UBs.update({rx: ub})
    rba_session.FBA.set_ub(rxn_UBs)
    rba_session.FBA.set_lb(rxn_LBs)

    rba_session.FBA.clear_objective()

    rba_session.FBA.set_objective({BMfunction: -1})

    rba_session.FBA.solve_lp()
    if rba_session.FBA.Solved:
        BMfluxOld = rba_session.FBA.SolutionValues[BMfunction]
    else:
        rba_session.FBA.clear_objective()
        rba_session.FBA.set_objective({BMfunction: -1.0})
        rba_session.FBA.set_lb({BMfunction:0.0})
        rba_session.FBA.solve_lp()
        BMfluxOld = rba_session.FBA.SolutionValues[BMfunction]
    if parsimonious_fba:
        rba_session.FBA.parsimonise(rxns_to_ignore_in_objective=rxns_to_ignore_when_parsimonious)
        rba_session.FBA.set_lb(rxn_LBs)
        rba_session.FBA.set_ub(rxn_UBs)
        if use_bm_flux_of_one:
            if BMfluxOld >= 1.0:
                rba_session.FBA.set_lb({BMfunction: 1.0})
                rba_session.FBA.set_ub({BMfunction: 1.0})
                rba_session.FBA.solve_lp()
                if not rba_session.FBA.Solved:
                    rba_session.FBA.set_lb({BMfunction: BMfluxOld})
                    rba_session.FBA.set_ub({BMfunction: BMfluxOld})
                    rba_session.FBA.solve_lp()
            else:
                rba_session.FBA.set_lb({BMfunction: BMfluxOld})
                rba_session.FBA.set_ub({BMfunction: BMfluxOld})
                rba_session.FBA.solve_lp()
        else:
            rba_session.FBA.set_lb({BMfunction: BMfluxOld})
            rba_session.FBA.set_ub({BMfunction: BMfluxOld})
            rba_session.FBA.solve_lp()

    fba_solution=rba_session.FBA.SolutionValues
    FluxDistribution = pandas.DataFrame()
    for rxn in fba_solution.keys():
        if rxn.startswith("Backward_"):
            FluxDistribution.loc[rxn.split("Backward_")[1],'FluxValues']=-rba_session.FBA.SolutionValues[rxn]
        else:
            FluxDistribution.loc[rxn,'FluxValues']=rba_session.FBA.SolutionValues[rxn]
    return(FluxDistribution)


def determine_reactions_associated_with_measured_proto_protein(measured_proteins_isoform_map,rba_session):
    """
    _summary_

    Parameters
    ----------
    measured_proteins_isoform_map : _type_
        _description_
    rba_session : _type_
        _description_
    """
    # identify all model reactions, associated with the measured proteins
    out = {}
    for p_ID in measured_proteins_isoform_map.keys():
        reactions_associated_with_proto_protein = []
        for isoform in measured_proteins_isoform_map[p_ID]:
            for reaction in rba_session.ModelStructure.ProteinInfo.Elements[isoform]['associatedReactions']:
                reactions_associated_with_proto_protein.append(reaction)
        out[p_ID] = list(set(reactions_associated_with_proto_protein))
    return(out)


def pre_select_iso_reactions(measured_proteins_reaction_map,rba_session,chose_most_quantified,keep_isorxns_specific_to_quantified_proteins=False):
    """
    _summary_

    Parameters
    ----------
    measured_proteins_reaction_map : _type_
        _description_
    rba_session : _type_
        _description_
    chose_most_quantified : _type_
        _description_
    keep_isorxns_specific_to_quantified_proteins : bool, optional
        _description_, by default False
    """
    # choose most likely iso-reaction for each measured-protein associated reaction
    proto_rxn_dict = {}
    for p_ID in measured_proteins_reaction_map.keys():
        for rxn in measured_proteins_reaction_map[p_ID]:
            rxn_to_split=str(rxn)
            proto_rxn = rxn_to_split.split('_duplicate')[0]
            if proto_rxn in list(proto_rxn_dict.keys()):
                if rxn in list(proto_rxn_dict[proto_rxn].keys()):
                    proto_rxn_dict[proto_rxn][rxn] += 1
                else:
                    proto_rxn_dict[proto_rxn].update({rxn: 1})
            else:
                proto_rxn_dict[proto_rxn] = {rxn: 1}
    out = {}
    for prx in proto_rxn_dict.keys():
        unique_subunit_dict={}
        for irx in proto_rxn_dict[prx].keys():
            enzyme = rba_session.ModelStructure.ReactionInfo.Elements[irx]['Enzyme']
            unique_subunit_dict[irx] = len(list(rba_session.get_enzyme_information(enzyme)["Subunits"].keys()))
        if chose_most_quantified:
            max_val = max([proto_rxn_dict[prx][i]/unique_subunit_dict[i] for i in proto_rxn_dict[prx].keys()])
            list_isorxns = [i for i in proto_rxn_dict[prx].keys() if proto_rxn_dict[prx][i]/unique_subunit_dict[i] == max_val]
            if len(list_isorxns)>1:
                max_SU_number=max([unique_subunit_dict[i] for i in list_isorxns])
                selected=[i for i in list_isorxns if unique_subunit_dict[i]==max_SU_number]
            else:
                selected=list_isorxns
            ## proposal of Ana in Oberseminar:
            if keep_isorxns_specific_to_quantified_proteins:
                protein_with_only_one_function_isoreactions=[measured_proteins_reaction_map[i][0] for i in measured_proteins_reaction_map.keys() if len(measured_proteins_reaction_map[i])==1]
                for irx in proto_rxn_dict[prx].keys():
                    if not irx in selected:
                        if irx in protein_with_only_one_function_isoreactions:
                            selected.append(irx)
        else:
            selected = [i for i in proto_rxn_dict[prx].keys() if proto_rxn_dict[prx][i] != 0]
        selected.sort()
        out[prx] = selected
    return(out)


def determine_machinery_concentration_by_weighted_geometric_mean(rba_session,machinery_composition,proteomicsData,proto_proteins=False):
    """
    _summary_

    Parameters
    ----------
    rba_session : _type_
        _description_
    machinery_composition : _type_
        _description_
    proteomicsData : _type_
        _description_
    proto_proteins : bool, optional
        _description_, by default False
    """
    subunit_derived_concentrations={}
    for subunit in machinery_composition.keys():
        if proto_proteins:
            if subunit in list(proteomicsData['ID']):
                copy_number_subunit=proteomicsData.loc[proteomicsData['ID']==subunit, 'copy_number'].values[0]
                if numpy.isfinite(copy_number_subunit):
                    if copy_number_subunit!=0:
                        if numpy.isfinite(machinery_composition[subunit]):
                            if machinery_composition[subunit]>0:
                                subunit_derived_concentrations[subunit]={"Concentration":copy_number_subunit,"Stoichiometry":machinery_composition[subunit]}
        else:
            subunit_proto_id=rba_session.get_protein_information(subunit)["ProtoID"]
            if subunit_proto_id in list(proteomicsData['ID']):
                copy_number_subunit=proteomicsData.loc[proteomicsData['ID']==subunit_proto_id, 'copy_number'].values[0]
                if numpy.isfinite(copy_number_subunit):
                    if copy_number_subunit!=0:
                        if numpy.isfinite(machinery_composition[subunit]):
                            if machinery_composition[subunit]>0:
                                subunit_derived_concentrations[subunit_proto_id]={"Concentration":copy_number_subunit,"Stoichiometry":machinery_composition[subunit]}
    if list(subunit_derived_concentrations.keys()):
        concentrations_of_subunits=[]
        stoichiometries_of_subunits=[]
        for i in subunit_derived_concentrations.keys():
            concentrations_of_subunits.append(subunit_derived_concentrations[i]["Concentration"]/subunit_derived_concentrations[i]["Stoichiometry"])
            stoichiometries_of_subunits.append(subunit_derived_concentrations[i]["Stoichiometry"])
        concentration=weighted_geometric_mean(data=concentrations_of_subunits,weights=stoichiometries_of_subunits)
    else:
        concentration=0
    return(concentration)


def determine_machinery_concentration(rba_session,machinery_composition,proteomicsData,proto_proteins=False):
    """
    _summary_

    Parameters
    ----------
    rba_session : _type_
        _description_
    machinery_composition : _type_
        _description_
    proteomicsData : _type_
        _description_
    proto_proteins : bool, optional
        _description_, by default False
    """
    subunit_derived_concentrations={}
    for subunit in machinery_composition.keys():
        if proto_proteins:
            if subunit in list(proteomicsData['ID']):
                copy_number_subunit=proteomicsData.loc[proteomicsData['ID']==subunit, 'copy_number'].values[0]
                if numpy.isfinite(copy_number_subunit):
                    if copy_number_subunit!=0:
                        if numpy.isfinite(machinery_composition[subunit]):
                            if machinery_composition[subunit]>0:
                                if subunit in subunit_derived_concentrations.keys():
                                    subunit_derived_concentrations[subunit]+=copy_number_subunit/machinery_composition[subunit]
                                else:
                                    subunit_derived_concentrations[subunit]=copy_number_subunit/machinery_composition[subunit]
        else:
            subunit_proto_id=rba_session.get_protein_information(subunit)["ProtoID"]
            if subunit_proto_id in list(proteomicsData['ID']):
                copy_number_subunit=proteomicsData.loc[proteomicsData['ID']==subunit_proto_id, 'copy_number'].values[0]
                if numpy.isfinite(copy_number_subunit):
                    if copy_number_subunit!=0:
                        if numpy.isfinite(machinery_composition[subunit]):
                            if machinery_composition[subunit]>0:
                                if subunit_proto_id in subunit_derived_concentrations.keys():
                                    subunit_derived_concentrations[subunit_proto_id]+=copy_number_subunit/machinery_composition[subunit]
                                else:
                                    subunit_derived_concentrations[subunit_proto_id]=copy_number_subunit/machinery_composition[subunit]
    if list(subunit_derived_concentrations.keys()):
        concentration=gmean(numpy.array(list(subunit_derived_concentrations.values())))
    else:
        concentration=0
    return(concentration)


def determine_fba_flux_carrying_isoreaction(rba_session,reaction_id,flux_distribution):
    """
    _summary_

    Parameters
    ----------
    rba_session : _type_
        _description_
    reaction_id : _type_
        _description_
    flux_distribution : _type_
        _description_
    """
    out=None
    if reaction_id in rba_session.get_reactions():
        if reaction_id in flux_distribution.index:
            if numpy.isfinite(flux_distribution.loc[reaction_id,'FluxValues']):
                if flux_distribution.loc[reaction_id,'FluxValues']!=0:
                    out=reaction_id
        else:
            for iso_reaction_id in rba_session.get_reaction_information(reaction_id)["Twins"]:
                if iso_reaction_id in flux_distribution.index:
                    if numpy.isfinite(flux_distribution.loc[iso_reaction_id,'FluxValues']):
                        if flux_distribution.loc[iso_reaction_id,'FluxValues']!=0:
                            out=iso_reaction_id
    return(out)


def global_efficiency_scaling(condition,
                              definition_file,
                              rba_session,
                              compartment_densities_and_pg,
                              process_efficiencies,
                              default_kapps,
                              specific_kapps,
                              exchanges_to_impose,
                              feasible_stati,
                              transporter_multiplier,
                              mu_approx_precision,
                              mu_misprediction_tolerance,
                              condition_to_look_up,
                              growth_rate_to_look_up,
                              results_to_look_up,
                              fixed_mu_when_above_target_mu_in_correction,
                              n_th_root_mispred=1,
                              print_outputs=False,
                              adjust_root=True,
                              proteomics_constraints_input={}):
    """
    _summary_

    Parameters
    ----------
    condition : _type_
        _description_
    definition_file : _type_
        _description_
    rba_session : _type_
        _description_
    compartment_densities_and_pg : _type_
        _description_
    process_efficiencies : _type_
        _description_
    default_kapps : _type_
        _description_
    specific_kapps : _type_
        _description_
    exchanges_to_impose : _type_
        _description_
    feasible_stati : _type_
        _description_
    transporter_multiplier : _type_
        _description_
    mu_approx_precision : _type_
        _description_
    mu_misprediction_tolerance : _type_
        _description_
    condition_to_look_up : _type_
        _description_
    growth_rate_to_look_up : _type_
        _description_
    results_to_look_up : _type_
        _description_
    fixed_mu_when_above_target_mu_in_correction : _type_
        _description_
    n_th_root_mispred : int, optional
        _description_, by default 1
    print_outputs : bool, optional
        _description_, by default False
    adjust_root : bool, optional
        _description_, by default True
    proteomics_constraints_input : dict, optional
        _description_, by default {}
    """
    mu_measured=growth_rate_from_input(input=definition_file, condition=condition)

    default_kapps_for_scaling=copy.deepcopy(default_kapps)
    process_efficiencies_for_scaling=copy.deepcopy(process_efficiencies)
    specific_kapps_for_scaling=copy.deepcopy(specific_kapps)

    total_product_correction_factor=1
    total_product_correction_factors=[total_product_correction_factor]

    simulation_results = perform_simulations(condition=condition,
                                             rba_session=rba_session,
                                             definition_file=definition_file,
                                             compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                             pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                             process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies_for_scaling}]),
                                             Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":default_kapps_for_scaling}]),
                                             Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":specific_kapps_for_scaling}]),
                                             Exchanges_to_impose=exchanges_to_impose,
                                             sims_to_perform=[condition_to_look_up],
                                             feasible_stati=feasible_stati,
                                             try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                             print_output=print_outputs,
                                             apply_model=False,
                                             transporter_multiplier=transporter_multiplier,
                                             start_val=0,
                                             #start_val=mu_measured,
                                             Mu_approx_precision=mu_approx_precision,
                                             max_mu_in_dichotomy=4.0)

    mumax_predicted=simulation_results[growth_rate_to_look_up]
    if print_outputs:
        print("Global start - {}:{}".format(condition,mumax_predicted))
    predicted_growth_rates=[mumax_predicted]

    if mumax_predicted == 0:
        mu_misprediction_factor=10
    else:
        mu_misprediction_factor=mu_measured/mumax_predicted
    if print_outputs:
        print("Global start - {}:{}".format(condition,mu_misprediction_factor))
    mu_iteration_count=0
    runs_of_sign=0
    last_misprediction_direction=0
    if fixed_mu_when_above_target_mu_in_correction:
        while not mu_measured <= mumax_predicted <= (mu_measured+mu_measured*mu_misprediction_tolerance):
            if mu_iteration_count>=10:
                break

            mu_iteration_count+=1

            if runs_of_sign>=2:
                if adjust_root:
                    n_th_root_mispred=2

            default_kapps_for_scaling["default_efficiency"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            default_kapps_for_scaling["default_transporter_efficiency"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            process_efficiencies_for_scaling.loc[:,"Value"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            specific_kapps_for_scaling.loc[:,"Kapp"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            
            total_product_correction_factor*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            total_product_correction_factors.append(total_product_correction_factor)

            simulation_results = perform_simulations(condition=condition,
                                                    rba_session=rba_session,
                                                    definition_file=definition_file,
                                                    compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                    pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                    process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies_for_scaling}]),
                                                    Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":default_kapps_for_scaling}]),
                                                    Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":specific_kapps_for_scaling}]),
                                                    Exchanges_to_impose=exchanges_to_impose,
                                                    sims_to_perform=[condition_to_look_up],
                                                    feasible_stati=feasible_stati,
                                                    try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                                    print_output=print_outputs,
                                                    apply_model=False,
                                                    transporter_multiplier=transporter_multiplier,
                                                    start_val=0,
                                                    #start_val=mu_measured,
                                                    Mu_approx_precision=mu_approx_precision,
                                                    max_mu_in_dichotomy=2*mu_measured)
            mumax_predicted=simulation_results[growth_rate_to_look_up]

            predicted_growth_rates.append(mumax_predicted)

            if mumax_predicted == 0:
                mu_misprediction_factor=10
            else:
                mu_misprediction_factor=mu_measured/mumax_predicted
                if mu_misprediction_factor<1:
                    current_misprediction_direction = 1
                elif mu_misprediction_factor>1:
                    current_misprediction_direction = -1
                else:
                    current_misprediction_direction=0

                if (current_misprediction_direction!=0) and (last_misprediction_direction!=0):
                    if current_misprediction_direction==last_misprediction_direction:
                        runs_of_sign=0
                    else:
                        if current_misprediction_direction==1:
                            runs_of_sign+=1
                else:
                    runs_of_sign=0

                last_misprediction_direction=current_misprediction_direction

            if print_outputs:
                print("Measured: {} - Predicted: {} - mispred coeff: {} - root: {} - runs_of_sign: {}".format(mu_measured,mumax_predicted,mu_misprediction_factor,n_th_root_mispred,runs_of_sign))
        if mumax_predicted==0:
            return(None)
        prediction_residulas_growth_rates=[mu_measured-mu_pred for mu_pred in predicted_growth_rates]
        minimum_prediction_residual=min([abs(i) for i in prediction_residulas_growth_rates])
        list_with_minimal_residuals=[]
        for i in prediction_residulas_growth_rates:
            if abs(i)==minimum_prediction_residual:
                list_with_minimal_residuals.append(i)
        index_of_best_iterration=prediction_residulas_growth_rates.index(min(list_with_minimal_residuals))
        best_cumulative_correction_factor=total_product_correction_factors[index_of_best_iterration]

        default_kapps_out=copy.deepcopy(default_kapps)
        process_efficiencies_out=copy.deepcopy(process_efficiencies)
        specific_kapps_out=copy.deepcopy(specific_kapps)

        default_kapps_out["default_efficiency"]*=best_cumulative_correction_factor
        default_kapps_out["default_transporter_efficiency"]*=best_cumulative_correction_factor
        process_efficiencies_out.loc[:,"Value"]*=best_cumulative_correction_factor
        specific_kapps_out.loc[:,"Kapp"]*=best_cumulative_correction_factor

        simulation_results_fixed = perform_simulations_fixed_Mu(condition=condition,
                                                                rba_session=rba_session,
                                                                definition_file=definition_file,
                                                                compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                                pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                                process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies_out}]),
                                                                Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":default_kapps_out}]),
                                                                Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":specific_kapps_out}]),
                                                                Exchanges_to_impose=exchanges_to_impose,
                                                                sims_to_perform=[condition_to_look_up],
                                                                feasible_stati=feasible_stati,
                                                                try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                                                print_output=print_outputs,
                                                                apply_model=False,
                                                                transporter_multiplier=transporter_multiplier,
                                                                proteomics_constraints_input=proteomics_constraints_input)
        if len(list(simulation_results_fixed[results_to_look_up].keys()))>0:
            simulation_results=simulation_results_fixed

    else:
        while not (mu_measured-mu_measured*mu_misprediction_tolerance) <= mumax_predicted <= (mu_measured+mu_measured*mu_misprediction_tolerance):
            mu_iteration_count+=1
            if mu_iteration_count>=10:
                break
            if runs_of_sign>=2:
                if adjust_root:
                    n_th_root_mispred=2

            default_kapps_for_scaling["default_efficiency"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            default_kapps_for_scaling["default_transporter_efficiency"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            process_efficiencies_for_scaling.loc[:,"Value"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            specific_kapps_for_scaling.loc[:,"Kapp"]*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            
            total_product_correction_factor*=numpy.power(mu_misprediction_factor,1/n_th_root_mispred)
            total_product_correction_factors.append(total_product_correction_factor)

            simulation_results = perform_simulations(condition=condition,
                                                    rba_session=rba_session,
                                                    definition_file=definition_file,
                                                    compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                    pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                    process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies_for_scaling}]),
                                                    Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":default_kapps_for_scaling}]),
                                                    Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":specific_kapps_for_scaling}]),
                                                    Exchanges_to_impose=exchanges_to_impose,
                                                    sims_to_perform=[condition_to_look_up],
                                                    feasible_stati=feasible_stati,
                                                    try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                                    print_output=print_outputs,
                                                    apply_model=False,
                                                    transporter_multiplier=transporter_multiplier,
                                                    start_val=0,
                                                    #start_val=mu_measured,
                                                    Mu_approx_precision=mu_approx_precision,
                                                    max_mu_in_dichotomy=2*mu_measured)
            mumax_predicted=simulation_results[growth_rate_to_look_up]

            predicted_growth_rates.append(mumax_predicted)

            if mumax_predicted == 0:
                mu_misprediction_factor=10
            else:
                mu_misprediction_factor=mu_measured/mumax_predicted
                if mu_misprediction_factor<1:
                    current_misprediction_direction = 1
                elif mu_misprediction_factor>1:
                    current_misprediction_direction = -1
                else:
                    current_misprediction_direction=0

                if (current_misprediction_direction!=0) and (last_misprediction_direction!=0):
                    if current_misprediction_direction==last_misprediction_direction:
                        runs_of_sign=0
                    else:
                        if current_misprediction_direction==1:
                            runs_of_sign+=1
                else:
                    runs_of_sign=0

                last_misprediction_direction=current_misprediction_direction

            if print_outputs:
                print("Measured: {} - Predicted: {} - mispred coeff: {} - root: {} - runs_of_sign: {}".format(mu_measured,
                                                                                                                           mumax_predicted,
                                                                                                                           mu_misprediction_factor,
                                                                                                                           n_th_root_mispred,
                                                                                                                           runs_of_sign))

        prediction_residulas_growth_rates=[mu_measured-mu_pred for mu_pred in predicted_growth_rates]
        minimum_prediction_residual=min([abs(i) for i in prediction_residulas_growth_rates])
        list_with_minimal_residuals=[]
        for i in prediction_residulas_growth_rates:
            if abs(i)==minimum_prediction_residual:
                list_with_minimal_residuals.append(i)
        index_of_best_iterration=prediction_residulas_growth_rates.index(min(list_with_minimal_residuals))
        best_cumulative_correction_factor=total_product_correction_factors[index_of_best_iterration]

        default_kapps_out=copy.deepcopy(default_kapps)
        process_efficiencies_out=copy.deepcopy(process_efficiencies)
        specific_kapps_out=copy.deepcopy(specific_kapps)

        default_kapps_out["default_efficiency"]*=best_cumulative_correction_factor
        default_kapps_out["default_transporter_efficiency"]*=best_cumulative_correction_factor
        process_efficiencies_out.loc[:,"Value"]*=best_cumulative_correction_factor
        specific_kapps_out.loc[:,"Kapp"]*=best_cumulative_correction_factor

        simulation_results = perform_simulations(condition=condition,
                                                rba_session=rba_session,
                                                definition_file=definition_file,
                                                compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_pg}]),
                                                process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies_out}]),
                                                Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":default_kapps_out}]),
                                                Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":specific_kapps_out}]),
                                                Exchanges_to_impose=exchanges_to_impose,
                                                sims_to_perform=[condition_to_look_up],
                                                feasible_stati=feasible_stati,
                                                try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                                print_output=print_outputs,
                                                apply_model=False,
                                                transporter_multiplier=transporter_multiplier,
                                                #start_val=mu_measured,
                                                start_val=0,
                                                Mu_approx_precision=mu_approx_precision,
                                                max_mu_in_dichotomy=2*mu_measured,
                                                proteomics_constraints_input=proteomics_constraints_input)

        if simulation_results[growth_rate_to_look_up] == 0:
            return(None)
    return({"specific_kapps":specific_kapps_out,"default_kapps":default_kapps_out,"process_efficiencies":process_efficiencies_out,"correction_factor":best_cumulative_correction_factor,"simulation_results":simulation_results})


def efficiency_correction(enzyme_efficiencies,
                               simulation_results,
                               protein_data,
                               rba_session,
                               condition_to_look_up,
                               default_enzyme_efficiencies,
                               tolerance=2,
                               n_th_root_mispred=2,
                               process_efficiencies=None,
                               correct_default_kapp_enzymes=False,
                               only_consider_misprediction_for_predicted_nonzero_enzymes=False,
                               max_kapp=None):
    """
    _summary_

    Parameters
    ----------
    enzyme_efficiencies : _type_
        _description_
    simulation_results : _type_
        _description_
    protein_data : _type_
        _description_
    rba_session : _type_
        _description_
    condition_to_look_up : _type_
        _description_
    default_enzyme_efficiencies : _type_
        _description_
    tolerance : int, optional
        _description_, by default 2
    n_th_root_mispred : int, optional
        _description_, by default 2
    process_efficiencies : _type_, optional
        _description_, by default None
    correct_default_kapp_enzymes : bool, optional
        _description_, by default False
    only_consider_misprediction_for_predicted_nonzero_enzymes : bool, optional
        _description_, by default False
    max_kapp : _type_, optional
        _description_, by default None
    """
    proto_protein_quantities={}
    for i in simulation_results["Proteins"].index:
        proto_protein_ID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if proto_protein_ID in proto_protein_quantities:
            proto_protein_quantities[proto_protein_ID]+=simulation_results["Proteins"].loc[i,condition_to_look_up]
        else:
            proto_protein_quantities[proto_protein_ID]=simulation_results["Proteins"].loc[i,condition_to_look_up]
    
    proto_protein_isoform_map={} 
    for i in rba_session.get_proteins():
        proto_protein_ID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if proto_protein_ID in proto_protein_isoform_map:
            proto_protein_isoform_map[proto_protein_ID].append(i)
        else:
            proto_protein_isoform_map[proto_protein_ID]=[i]

    enzyme_efficiencies_out=enzyme_efficiencies.copy() 
    process_efficiencies_out=process_efficiencies.copy()
    subunit_misprediction_factors_enzymes={}
    subunit_stoichiometries_enzymes={}
    subunit_misprediction_factors_default_efficiency_enzymes={}
    subunit_stoichiometries_default_efficiency_enzymes={}
    subunit_misprediction_factors_processes={}
    subunit_stoichiometries_processes={}
    squared_residuals=[]
    for proto_protein_ID in proto_protein_quantities:
        if proto_protein_ID in proto_protein_isoform_map:
            if proto_protein_ID in protein_data["ID"]:
                predicted_protein=proto_protein_quantities[proto_protein_ID]
                measured_protein=protein_data.loc[proto_protein_ID,"copy_number"]
                if (predicted_protein>0) & (measured_protein>0):
                    misprediction_coeff=predicted_protein/measured_protein
                    if (numpy.isfinite(misprediction_coeff)) and (misprediction_coeff!=0):
                        squared_residuals.append((numpy.log(predicted_protein)-numpy.log(measured_protein))**2)
                        for protein in proto_protein_isoform_map[proto_protein_ID]:
                            for enzyme in rba_session.get_protein_information(protein=protein)["associatedEnzymes"]:
                                if only_consider_misprediction_for_predicted_nonzero_enzymes:
                                    if enzyme not in simulation_results["Enzymes"].index:
                                        continue
                                    elif numpy.isfinite(simulation_results["Enzymes"].loc[enzyme,condition_to_look_up])==False:
                                        continue
                                    elif simulation_results["Enzymes"].loc[enzyme,condition_to_look_up]==0:
                                        continue
                                if enzyme in list(enzyme_efficiencies["Enzyme_ID"]):
                                    if enzyme not in subunit_misprediction_factors_enzymes.keys():
                                        subunit_misprediction_factors_enzymes[enzyme]=[misprediction_coeff]
                                        subunit_stoichiometries_enzymes[enzyme]=[rba_session.get_enzyme_information(enzyme)["Subunits"][protein]]
                                    else:
                                        subunit_misprediction_factors_enzymes[enzyme].append(misprediction_coeff)
                                        subunit_stoichiometries_enzymes[enzyme].append(rba_session.get_enzyme_information(enzyme)["Subunits"][protein])
                                else:
                                    if enzyme not in subunit_misprediction_factors_default_efficiency_enzymes.keys():
                                        subunit_misprediction_factors_default_efficiency_enzymes[enzyme]=[misprediction_coeff]
                                        subunit_stoichiometries_default_efficiency_enzymes[enzyme]=[rba_session.get_enzyme_information(enzyme)["Subunits"][protein]]
                                    else:
                                        subunit_misprediction_factors_default_efficiency_enzymes[enzyme].append(misprediction_coeff)
                                        subunit_stoichiometries_default_efficiency_enzymes[enzyme].append(rba_session.get_enzyme_information(enzyme)["Subunits"][protein])
                            for process in rba_session.get_protein_information(protein=protein)["SupportsProcess"]:
                                if process in process_efficiencies.index:
                                    if process not in subunit_misprediction_factors_processes.keys():
                                        subunit_misprediction_factors_processes[process]=[misprediction_coeff]
                                        subunit_stoichiometries_processes[process]=[rba_session.get_process_information(process)['Composition'][protein]]
                                    else:
                                        subunit_misprediction_factors_processes[process].append(misprediction_coeff)
                                        subunit_stoichiometries_processes[process].append(rba_session.get_process_information(process)['Composition'][protein])

    enzyme_correction_coefficients={}
    for enzyme in subunit_misprediction_factors_enzymes.keys():
        enzyme_correction_coefficients[enzyme]=numpy.power(weighted_geometric_mean(data=subunit_misprediction_factors_enzymes[enzyme],weights=subunit_stoichiometries_enzymes[enzyme]),1/n_th_root_mispred)
    enzymes_already_handled=[]
    for enzyme in enzyme_correction_coefficients.keys():
        if enzyme not in enzymes_already_handled:
            iso_enzyme_correction_coeffs={enzyme:enzyme_correction_coefficients[enzyme]}
            for iso_enzyme in rba_session.get_enzyme_information(enzyme)["Isozymes"]:
                if iso_enzyme in enzyme_correction_coefficients.keys():
                    iso_enzyme_correction_coeffs.update({iso_enzyme:enzyme_correction_coefficients[iso_enzyme]})
                else:
                    iso_enzyme_correction_coeffs.update({iso_enzyme:numpy.nan})
            correction_coeff=numpy.median(numpy.array([i for i in list(iso_enzyme_correction_coeffs.values()) if numpy.isfinite(i)]))
            for iso_enzyme_to_consider in iso_enzyme_correction_coeffs.keys():
                if iso_enzyme_to_consider in list(enzyme_efficiencies["Enzyme_ID"]):
                    enzymes_already_handled.append(iso_enzyme_to_consider)
                    old_kapp=enzyme_efficiencies.loc[enzyme_efficiencies["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"].values[0]
                    new_kapp=old_kapp*correction_coeff
                    if (max_kapp is not None) and (new_kapp > max_kapp):
                        continue
                    if tolerance is None:
                        enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"]=new_kapp
                    else:
                        if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                            enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"]=new_kapp

    if correct_default_kapp_enzymes:
        default_efficiency_enzyme_correction_coefficients={}
        for enzyme in subunit_misprediction_factors_default_efficiency_enzymes.keys():
            if enzyme not in list(enzyme_efficiencies["Enzyme_ID"]):
                default_efficiency_enzyme_correction_coefficients[enzyme]=numpy.power(weighted_geometric_mean(data=subunit_misprediction_factors_default_efficiency_enzymes[enzyme],weights=subunit_stoichiometries_default_efficiency_enzymes[enzyme]),1/n_th_root_mispred)
        for enzyme in default_efficiency_enzyme_correction_coefficients.keys():
            if enzyme not in enzymes_already_handled:
                respective_reaction=rba_session.get_enzyme_information(enzyme)["Reaction"].split("_duplicate_")[0]
                if respective_reaction in simulation_results["Reactions"].index:
                    flux=simulation_results["Reactions"].loc[respective_reaction,condition_to_look_up]
                    try:
                        if flux < 0:
                            flux_direction=-1
                            efficiency_parameter= rba_session.model.enzymes.enzymes._elements_by_id[enzyme].backward_efficiency
                        elif flux > 0:
                            flux_direction=1
                            efficiency_parameter= rba_session.model.enzymes.enzymes._elements_by_id[enzyme].forward_efficiency
                        else:
                            continue
                    except:
                        continue
                    iso_enzyme_correction_coeffs={enzyme:default_efficiency_enzyme_correction_coefficients[enzyme]}
                    for iso_enzyme in rba_session.get_enzyme_information(enzyme)["Isozymes"]:
                        if iso_enzyme in default_efficiency_enzyme_correction_coefficients.keys():
                            iso_enzyme_correction_coeffs.update({iso_enzyme:default_efficiency_enzyme_correction_coefficients[iso_enzyme]})
                        else:
                            iso_enzyme_correction_coeffs.update({iso_enzyme:numpy.nan})
                    correction_coeff=numpy.median(numpy.array([i for i in list(iso_enzyme_correction_coeffs.values()) if numpy.isfinite(i)]))
                    for iso_enzyme_to_consider in iso_enzyme_correction_coeffs.keys():
                        enzymes_already_handled.append(iso_enzyme_to_consider)
                        old_kapp=default_enzyme_efficiencies[efficiency_parameter]
                        new_kapp=old_kapp*correction_coeff
                        associated_reaction=rba_session.get_enzyme_information(iso_enzyme_to_consider)["Reaction"]
                        if (max_kapp is not None) and (new_kapp > max_kapp):
                            continue
                        if tolerance is None:
                            enzyme_efficiencies_out.loc[associated_reaction,"Kapp"]=new_kapp
                            enzyme_efficiencies_out.loc[associated_reaction,"Enzyme_ID"]=iso_enzyme_to_consider
                            enzyme_efficiencies_out.loc[associated_reaction,"Flux"]=flux_direction
                            enzyme_efficiencies_out.loc[associated_reaction,"Comment"]="Corrected Default"
                        else:
                            if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                                enzyme_efficiencies_out.loc[associated_reaction,"Kapp"]=new_kapp
                                enzyme_efficiencies_out.loc[associated_reaction,"Enzyme_ID"]=iso_enzyme_to_consider
                                enzyme_efficiencies_out.loc[associated_reaction,"Flux"]=flux_direction
                                enzyme_efficiencies_out.loc[associated_reaction,"Comment"]="Corrected Default"

    process_correction_coefficients={}
    for process in subunit_misprediction_factors_processes.keys():
        if process in list(process_efficiencies.index):
            #process_correction_coefficients[process]=numpy.power(numpy.median(subunit_misprediction_factors_processes[process]),1/n_th_root_mispred)
            process_correction_coefficients[process]=numpy.power(weighted_geometric_mean(data=subunit_misprediction_factors_processes[process],weights=subunit_stoichiometries_processes[process]),1/n_th_root_mispred)
    for process in process_correction_coefficients.keys():
        correction_coeff=process_correction_coefficients[process]
        old_efficiency=process_efficiencies.loc[process,"Value"]
        new_efficiency=old_efficiency*correction_coeff
        if tolerance is None:
            process_efficiencies_out.loc[process,"Value"]=new_efficiency
        else:
            if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                process_efficiencies_out.loc[process,"Value"]=new_efficiency

    return({"Sum_of_squared_residuals":sum(squared_residuals)/len(squared_residuals),"Kapps":enzyme_efficiencies_out,"ProcessEfficiencies":process_efficiencies_out,"Process_MispredictionFactors":process_correction_coefficients,"Enzyme_MispredictionFactors":enzyme_correction_coefficients})


############    
def import_process_efficiencies(input_data,rba_session,condition):
    """
    _summary_

    Parameters
    ----------
    input_data : _type_
        _description_
    rba_session : _type_
        _description_
    condition : _type_
        _description_
    """
    id_name_map={rba_session.get_process_information(i)["ID"]:rba_session.get_process_information(i)["Name"] for i in rba_session.get_processes()}
    out=pandas.DataFrame()
    for i in list(input_data.index):
        out.loc[id_name_map[i],"Value"]=input_data.loc[i,condition]
        out.loc[id_name_map[i],"Process"]=i
        out.loc[id_name_map[i],"Parameter"]=str(i+'_apparent_efficiency')
    return(out)


def import_specific_enzyme_efficiencies(input_data,rba_session,condition):
    """
    _summary_

    Parameters
    ----------
    input_data : _type_
        _description_
    rba_session : _type_
        _description_
    condition : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in list(input_data["Enzyme_ID"]):
        out.loc[rba_session.get_enzyme_information(i)["Reaction"],"Enzyme_ID"]=i
        out.loc[rba_session.get_enzyme_information(i)["Reaction"],"Flux_FBA"]=input_data.loc[i,condition+"_Flux"]
        if numpy.isfinite(input_data.loc[i,condition+"_Flux"]):
            if input_data.loc[i,condition+"_Flux"]<0:
                out.loc[rba_session.get_enzyme_information(i)["Reaction"],"Flux"]=-1.0
            else:
                out.loc[rba_session.get_enzyme_information(i)["Reaction"],"Flux"]=1.0
        else:
            out.loc[rba_session.get_enzyme_information(i)["Reaction"],"Flux"]=numpy.nan
        out.loc[rba_session.get_enzyme_information(i)["Reaction"],"Kapp"]=input_data.loc[i,condition]
    return(out)


def import_default_enzyme_efficiencies(input_data,condition,default_transporter_kapp_coefficient=1):
    """
    _summary_

    Parameters
    ----------
    input_data : _type_
        _description_
    condition : _type_
        _description_
    default_transporter_kapp_coefficient : int, optional
        _description_, by default 1
    """
    return({'default_efficiency': input_data.loc[condition,"Default Kapp"], 'default_transporter_efficiency':input_data.loc[condition,"Default Kapp"]*default_transporter_kapp_coefficient})


def weighted_geometric_mean(data,weights=None):
    """
    _summary_

    Parameters
    ----------
    data : _type_
        _description_
    weights : _type_, optional
        _description_, by default None
    """
    if weights is None:
        value=1
        for i in data:
            value*=i
        out=value**(1/len(data))
    else:
        value=1
        for i in range(len(data)):
            value*=(data[i]**weights[i])
        out=value**(1/(sum(weights)))
    return(out)


def extract_compartment_sizes_from_calibration_outputs(calibration_outputs):
    """
    _summary_

    Parameters
    ----------
    calibration_outputs : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in calibration_outputs:
        for j in i["Densities_PGs"].index:
            out.loc[j,i["Condition"]]=round(i["Densities_PGs"].loc[j,"Density"],6)
    return(out)


def extract_pg_fractions_from_calibration_outputs(calibration_outputs):
    """
    _summary_

    Parameters
    ----------
    calibration_outputs : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in calibration_outputs:
        for j in i["Densities_PGs"].index:
            out.loc[j,i["Condition"]]=round(i["Densities_PGs"].loc[j,"PG_fraction"],6)
    return(out)


def extract_specific_kapps_from_calibration_outputs(calibration_outputs):
    """
    _summary_

    Parameters
    ----------
    calibration_outputs : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in calibration_outputs:
        for j in i['Specific_Kapps']["Enzyme_ID"]:
            out.loc[j,"Enzyme_ID"]=j
            out.loc[j,i["Condition"]]=i['Specific_Kapps'].loc[i['Specific_Kapps']["Enzyme_ID"]==j,"Kapp"].values[0]
            out.loc[j,i["Condition"]+"_Flux"]=i['Specific_Kapps'].loc[i['Specific_Kapps']["Enzyme_ID"]==j,"Flux"].values[0]
    return(out)


def extract_default_kapps_from_calibration_outputs(calibration_outputs):
    """
    _summary_

    Parameters
    ----------
    calibration_outputs : _type_
        _description_
    """
    out = pandas.DataFrame()
    for i in calibration_outputs:
        out.loc[i["Condition"], 'ID'] = i["Condition"]
        out.loc[i["Condition"], 'Default Kapp'] = i['Default_Kapps']["default_efficiency"]
    return(out)


def extract_process_capacities_from_calibration_outputs(calibration_outputs):
    """
    _summary_

    Parameters
    ----------
    calibration_outputs : _type_
        _description_
    """
    out = pandas.DataFrame()
    for i in calibration_outputs:
        for j in list(i['Process_Efficiencies']['Process']):
            val=i['Process_Efficiencies'].loc[i['Process_Efficiencies']['Process']==j,'Value'].values[0]
            out.loc[j,i["Condition"]]=val
    return(out)

############    
def generate_mean_enzyme_composition_model(rba_session,condition):
    enzymes_already_handled=[]
    model_enzymes=rba_session.get_enzymes()
    for i in model_enzymes:
        if i in enzymes_already_handled:
            continue
        all_iso_enzymes=list([i]+rba_session.get_enzyme_information(i)['Isozymes'])
        enzymes_already_handled+=all_iso_enzymes    
        if len(all_iso_enzymes)==1:
            continue
        reactant_dict={}
        product_dict={}
        for j in all_iso_enzymes:
            try:
                for reactant in rba_session.model.enzymes.enzymes.get_by_id(j).machinery_composition.reactants._elements:
                    species = reactant.species
                    if species in reactant_dict.keys():
                        reactant_dict[species]+=reactant.stoichiometry
                    else:
                        reactant_dict[species]=reactant.stoichiometry
                for product in rba_session.model.enzymes.enzymes.get_by_id(j).machinery_composition.products._elements:
                    species = product.species
                    if species in product_dict.keys():
                        product_dict[species]+=product.stoichiometry
                    else:
                        product_dict[species]=product.stoichiometry
            except:
                continue
        mean_reactants={species:reactant_dict[species]/len(all_iso_enzymes) for species in reactant_dict.keys()}
        mean_products={species:product_dict[species]/len(all_iso_enzymes) for species in product_dict.keys()}
        proto_enzyme=all_iso_enzymes[[len(i) for i in all_iso_enzymes].index(min([len(i) for i in all_iso_enzymes]))]

        if len(list(mean_reactants.keys()))!=0:
            rba_session.model.enzymes.enzymes.get_by_id(proto_enzyme).machinery_composition.reactants=rba.xml.ListOfReactants()
            for reactant in mean_reactants.keys():
                species_reference=rba.xml.SpeciesReference(species=reactant, stoichiometry=mean_reactants[reactant])
                rba_session.model.enzymes.enzymes.get_by_id(proto_enzyme).machinery_composition.reactants.append(species_reference)
        if len(list(mean_products.keys()))!=0:
            rba_session.model.enzymes.enzymes.get_by_id(proto_enzyme).machinery_composition.products=rba.xml.ListOfProducts()
            for product in mean_products.keys():
                species_reference=rba.xml.SpeciesReference(species=product, stoichiometry=mean_products[product])
                rba_session.model.enzymes.enzymes.get_by_id(proto_enzyme).machinery_composition.products.append(species_reference)
        for iso_enzyme in all_iso_enzymes:
            if iso_enzyme != proto_enzyme:
                reaction_id=rba_session.model.enzymes.enzymes.get_by_id(iso_enzyme).reaction
                rba_session.model.enzymes.enzymes.remove(rba_session.model.enzymes.enzymes.get_by_id(iso_enzyme))
                rba_session.model.metabolism.reactions.remove(rba_session.model.metabolism.reactions.get_by_id(reaction_id))
            #delete other isoenzymes and isorxns
    rba_session.build_model_structure(file_name='/ModelStructure_meancompo_{}.json'.format(condition),print_warnings=False)
 
