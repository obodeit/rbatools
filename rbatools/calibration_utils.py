import rba
import copy
import pandas
import time
import numpy
#import json
from scipy.stats.mstats import gmean
from rbatools.other_utils import extract_compsizes_and_pgfractions_from_correction_summary , build_proteome_overview , correct_proteome , correct_compartment_fractions , medium_concentrations_from_input , machinery_efficiency_correction_settings_from_input , enzyme_efficiency_estimation_settings_from_input , flux_bounds_from_input , growth_rate_from_input , proteome_fractions_from_input , perform_simulations , perform_simulations_fixed_Mu

def calibration_workflow_2(proteome,
                         condition,
                         definition_file,
                         rba_session,
                         process_efficiency_estimation_input=None,
                         spec_kapps=None,
                         default_kapps=None,
                         process_efficiencies=None,
                         compartment_sizes=None,
                         pg_fractions=None,
                         print_outputs=True,
                         output_dir=""):

    ### DEFINE THESE IN SETTINGS-CSV
    use_mean_enzyme_composition_for_calibration=False
    Mu_approx_precision=0.0000001
    transporter_multiplier=1
    #feasible_stati=["optimal","feasible"]
    feasible_stati=["optimal","feasible","feasible_only_before_unscaling"]
    

    ### SPECIFY THAT ALWAYS THE PROKARYOTIC RBA-ALGORITHM IS USED FOR CALIBRATION
    condition_to_look_up="Prokaryotic"
    growth_rate_to_look_up="Mu_prok"
    results_to_look_up="Simulation_Results"

    # IF comp sizes and pg fractions not provided: estimate.
    # Otherwise don´t
    t0 = time.time()

    rba_session.set_growth_rate(growth_rate_from_input(input=definition_file,condition=condition))

    compartment_occupation_overview = build_proteome_overview(input=proteome, 
                                                              condition=condition, 
                                                              compartments_to_replace={},
                                                              compartments_no_original_PG=[], 
                                                              ribosomal_proteins_as_extra_compartment=False)

    if (compartment_sizes is not None) and (pg_fractions is not None):
        compartment_densities_and_PGs=pandas.DataFrame()
        for comp in list(compartment_sizes.index):
            compartment_densities_and_PGs.loc[comp,"Compartment_ID"]=comp
            compartment_densities_and_PGs.loc[comp,"Density"]=compartment_sizes.loc[comp,condition]
            compartment_densities_and_PGs.loc[comp,"PG_fraction"]=pg_fractions.loc[comp,condition]
    ### MOVE TO script level ###
    else:
        compartment_densities_and_PGs = extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=compartment_occupation_overview,
                                                                                                  rows_to_exclude=["Ribosomes","Total"]+[i for i in compartment_occupation_overview.index if i.startswith("pg_")])

    if process_efficiencies is not None:
        process_efficiencies=import_process_efficiencies(input_data=process_efficiencies,rba_session=rba_session,condition=condition)
    else:
        if process_efficiency_estimation_input is not None:
            process_efficiencies = determine_apparent_process_efficiencies_2(growth_rate=growth_rate_from_input(input=definition_file,condition=condition), 
                                                                             input=process_efficiency_estimation_input,
                                                                             rba_session=rba_session,
                                                                             compartment_densities_and_pg=compartment_densities_and_PGs, 
                                                                             protein_data=proteome.copy(),
                                                                             condition=condition,
                                                                             fit_nucleotide_assembly_machinery=True)
            
            #process_efficiencies_old = determine_apparent_process_efficiencies(growth_rate=growth_rate_from_input(input=definition_file,condition=condition),
            #                                                               input=process_efficiency_estimation_input,
            #                                                               rba_session=rba_session,
            #                                                               protein_data=proteome.copy(),
            #                                                               compartment_densities_and_PGs=compartment_densities_and_PGs,
            #                                                               #total_amino_acid_abundance_in_proteome=correction_results_compartement_sizes.loc['Total', 'original_amino_acid_occupation'],
            #                                                               total_amino_acid_abundance_in_proteome=compartment_occupation_overview.loc['Total', 'original_amino_acid_occupation'],
            #                                                               condition=condition,
            #                                                               fit_nucleotide_assembly_machinery=True)

        #process_efficiencies.to_csv(output_dir+'/ProcEffsOrig_{}.csv'.format(condition))
    process_efficiencies_original=process_efficiencies.copy()

    if use_mean_enzyme_composition_for_calibration:
        generate_mean_enzyme_composition_model(rba_session,condition)
            
    correction_settings=machinery_efficiency_correction_settings_from_input(input=definition_file, condition=condition)
    
    correction_settings['impose_directions_from_fba_during_correction']=False

    if spec_kapps is not None:
        Specific_Kapps=import_specific_enzyme_efficiencies(input_data=spec_kapps,rba_session=rba_session,condition=condition)
        #fba_flux_directions_from_enzyme_efficiency_estimation={}
        correction_settings['impose_directions_from_fba_during_correction']=False
        # or raise error
    else:
        enzyme_efficiency_estimation_settings=enzyme_efficiency_estimation_settings_from_input(input=definition_file, condition=condition)

        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))
        
        flux_bounds_fba=flux_bounds_from_input(input=definition_file,
                                            rba_session=rba_session, 
                                            condition=condition, 
                                            specific_exchanges=None, 
                                            specific_directions=None,
                                            also_consider_iso_enzmes=False)

        Specific_Kapps_Results = estimate_specific_enzyme_efficiencies_new(rba_session=rba_session, 
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
                                                                    impose_on_identical_enzymes=enzyme_efficiency_estimation_settings['impose_on_identical_enzymes'],
                                                                    condition=condition, 
                                                                    store_output=True,
                                                                    rxns_to_ignore_when_parsimonious=[],
                                                                    use_bm_flux_of_one=True, # False
                                                                    pseudocomplex_method='new',
                                                                    output_dir=output_dir)

        Specific_Kapps=Specific_Kapps_Results["Overview"]
        
        Spec_Kapp_estim_FD=Specific_Kapps_Results["Flux_Distribution"]


        fba_flux_directions_from_enzyme_efficiency_estimation={}
        if correction_settings['impose_directions_from_fba_during_correction']:
            for fba_rxn in list(Spec_Kapp_estim_FD.index):
                if fba_rxn not in list(flux_bounds_fba.index):
                    flux_value=Spec_Kapp_estim_FD.loc[fba_rxn,'FluxValues']
                    if fba_rxn in list(rba_session.get_reactions()):
                        if flux_value<0:
                            fba_flux_directions_from_enzyme_efficiency_estimation.update({fba_rxn:{"LB":numpy.nan,"UB":0}})
                        elif flux_value>0:
                            fba_flux_directions_from_enzyme_efficiency_estimation.update({fba_rxn:{"LB":0,"UB":numpy.nan}})
                        for iso_rxn in rba_session.get_reaction_information(fba_rxn)['Twins']:
                            if flux_value<0:
                                fba_flux_directions_from_enzyme_efficiency_estimation.update({iso_rxn:{"LB":numpy.nan,"UB":0}})
                            elif flux_value>0:
                                fba_flux_directions_from_enzyme_efficiency_estimation.update({iso_rxn:{"LB":0,"UB":numpy.nan}})

        #Specific_Kapps.to_csv(output_dir+'/Specific_Kapps_Hackett__{}.csv'.format(condition), sep=";", decimal=",")

        ## NEW
        #if enzyme_efficiency_estimation_settings['use_target_biomass_function']:
         #   rba_session.reload_model()
          #  rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))
    Specific_Kapps_original=Specific_Kapps.copy()

    if default_kapps is not None:
        Default_Kapps=import_default_enzyme_efficiencies(input_data=default_kapps,condition=condition,default_transporter_kapp_coefficient=transporter_multiplier)
    else:
        # DEfine default enzyme effficency as median over estimated specific ones 
        Default_Kapps=calculate_default_enzyme_efficiency_as_median_over_specific_efficiencies(specific_enzyme_efficiencies=Specific_Kapps,
                                                                                               transporter_multiplier=transporter_multiplier)
    Default_Kapps_original=Default_Kapps.copy()

    if correction_settings['correct_efficiencies']:
        flux_bounds_data=flux_bounds_from_input(input=definition_file,rba_session=rba_session, condition=condition, specific_exchanges=None, specific_directions=None,also_consider_iso_enzmes=True)
        flux_bounds_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in list(flux_bounds_data["Reaction_ID"])}
        if correction_settings['impose_directions_from_fba_during_correction']:
            for rxn in fba_flux_directions_from_enzyme_efficiency_estimation.keys():
                if rxn not in flux_bounds_to_impose.keys():
                    flux_bounds_to_impose.update({rxn:fba_flux_directions_from_enzyme_efficiency_estimation[rxn]})
            flux_bounds_to_impose.update(fba_flux_directions_from_enzyme_efficiency_estimation)
        
        #rba_session.Problem.set_constraint_types(dict(zip(Simulation.get_enzyme_constraints(),['E']*len(Simulation.get_enzyme_constraints()))))

        steady_count=0
        increasing_RSS_count=0
        iteration_count=0
        previous_RSS=numpy.nan

        efficiencies_over_correction_iterations=[]
        rss_trajectory=[]

        continuation_criterion_correction=True
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
                                                                exchanges_to_impose=flux_bounds_to_impose,
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
                                                                adjust_root=correction_settings['adjust_root_of_correction_coeffs_global_scaling'])

            if results_global_scaling is None:
                break
            
            Simulation_results=results_global_scaling["simulation_results"]
            if len(list(Simulation_results[results_to_look_up].keys()))!=0:
                current_RSS=calculate_log_protein_prediction_mean_squared_error(protein_data=build_input_proteome_for_specific_kapp_estimation(proteome, condition),
                                                                                simulation_results=Simulation_results[results_to_look_up],
                                                                                rba_session=rba_session,
                                                                                condition_to_look_up=condition_to_look_up)
                rss_trajectory.append(current_RSS)

                Specific_Kapps=results_global_scaling["specific_kapps"]
                Default_Kapps=results_global_scaling["default_kapps"]
                process_efficiencies=results_global_scaling["process_efficiencies"]
                efficiencies_over_correction_iterations.append({"Specific_Kapps":Specific_Kapps.copy(),"Default_Kapps":Default_Kapps.copy(),"Process_Efficiencies":process_efficiencies.copy()})

                #build_input_proteome_for_specific_kapp_estimation(protein_data, condition)
                KappCorrectionResults=efficiency_correction_new(enzyme_efficiencies=Specific_Kapps,
                                                                simulation_results=Simulation_results[results_to_look_up],
                                                                protein_data=build_input_proteome_for_specific_kapp_estimation(proteome, condition),
                                                                rba_session=rba_session,
                                                                condition_to_look_up=condition_to_look_up,
                                                                default_enzyme_efficiencies=Default_Kapps,
                                                                tolerance=None,
                                                                #tolerance=1.05,
                                                                n_th_root_mispred=1,
                                                                process_efficiencies=process_efficiencies,
                                                                correct_default_kapp_enzymes=True,
                                                                only_consider_misprediction_for_predicted_nonzero_enzymes=True)

                Specific_Kapps=KappCorrectionResults["Kapps"]
                process_efficiencies=KappCorrectionResults["ProcessEfficiencies"]
                Default_Kapps=calculate_default_enzyme_efficiency_as_median_over_specific_efficiencies(specific_enzyme_efficiencies=Specific_Kapps,
                                                                                                       transporter_multiplier=transporter_multiplier)

                if iteration_count>=correction_settings['minimum_iteration_number']:
                    if not pandas.isna(previous_RSS):
                        if not pandas.isna(current_RSS):
                            if (1-correction_settings['rss_tolerance'])<=current_RSS/previous_RSS<=(1+correction_settings['rss_tolerance']):
                                steady_count+=1
                            else:
                                steady_count=0
                                #
                                if current_RSS>=correction_settings['increasing_rss_factor']*previous_RSS:
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

                #put in function: maybe also counts on top
                if current_RSS>rss_trajectory[0]:
                    continuation_criterion_correction=False
                if steady_count>=correction_settings['steady_rss_limit']:
                    continuation_criterion_correction=False
                elif iteration_count>=correction_settings['iteration_limit']:
                    continuation_criterion_correction=False
                elif increasing_RSS_count>=correction_settings['increasing_rss_limit']:
                    continuation_criterion_correction=False
            else:
                if iteration_count>=correction_settings['iteration_limit']:
                    continuation_criterion_correction=False
        #
        if len(rss_trajectory)>0:
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
                flux_bounds_to_impose=None
            results_global_scaling=global_efficiency_scaling(condition=condition,
                                                             definition_file=definition_file,
                                                             rba_session=rba_session,
                                                             compartment_densities_and_pg=compartment_densities_and_PGs,
                                                             process_efficiencies=process_efficiencies_to_return,
                                                             default_kapps=Default_Kapps_to_return,
                                                             specific_kapps=Specific_Kapps_to_return,
                                                             exchanges_to_impose=flux_bounds_to_impose,
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
            flux_bounds_to_impose=None
        else:
            flux_bounds_data=flux_bounds_from_input(input=definition_file,rba_session=rba_session, condition=condition, specific_exchanges=None, specific_directions=None,also_consider_iso_enzmes=True)
            flux_bounds_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in list(flux_bounds_data["Reaction_ID"])}
            if correction_settings['impose_directions_from_fba_during_correction']:
                for rxn in fba_flux_directions_from_enzyme_efficiency_estimation.keys():
                    if rxn not in flux_bounds_to_impose.keys():
                        flux_bounds_to_impose.update({rxn:fba_flux_directions_from_enzyme_efficiency_estimation[rxn]})
                flux_bounds_to_impose.update(fba_flux_directions_from_enzyme_efficiency_estimation)
            
        results_global_scaling=global_efficiency_scaling(condition=condition,
                                                         definition_file=definition_file,
                                                         rba_session=rba_session,
                                                         compartment_densities_and_pg=compartment_densities_and_PGs,
                                                         process_efficiencies=process_efficiencies,
                                                         default_kapps=Default_Kapps,
                                                         specific_kapps=Specific_Kapps,
                                                         exchanges_to_impose=flux_bounds_to_impose,
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
            #'correction_results_compartement_sizes': correction_results_compartement_sizes,
            'Default_Kapps': Default_Kapps_to_return,
            'Specific_Kapps': Specific_Kapps_to_return,
            'Process_Efficiencies': process_efficiencies_to_return,
            'Default_Kapps_original': Default_Kapps_original,
            'Specific_Kapps_original': Specific_Kapps_original,
            'Process_Efficiencies_original': process_efficiencies_original})

def calculate_default_enzyme_efficiency_as_median_over_specific_efficiencies(specific_enzyme_efficiencies,transporter_multiplier=1):
    spec_kapp_median=specific_enzyme_efficiencies.loc[(specific_enzyme_efficiencies['Kapp']!=0)&(pandas.isna(specific_enzyme_efficiencies['Kapp'])==False),'Kapp'].median()
    #spec_kapp_median=numpy.median(specific_enzyme_efficiencies.loc[(specific_enzyme_efficiencies['Kapp']!=0)&(pandas.isna(specific_enzyme_efficiencies['Kapp'])==False),'Kapp'].unique())
    return({"default_efficiency":spec_kapp_median,"default_transporter_efficiency":transporter_multiplier*spec_kapp_median})

def determine_apparent_process_efficiencies_2(growth_rate,input, rba_session,compartment_densities_and_pg, protein_data, condition,fit_nucleotide_assembly_machinery=False):

    for comp in list(compartment_densities_and_pg['Compartment_ID']):
        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_pg.loc[compartment_densities_and_pg['Compartment_ID'] == comp, 'Density'].values[0]
        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = 0.0
    rba_session.rebuild_from_model()
    rba_session.set_growth_rate(growth_rate)

#    #{'ProtoID':[isoform_loc_1,isoform_loc_2,...]}

    protoprotein_isoprotein_map = rba_session.ModelStructure.ProteinInfo.return_protein_iso_form_map()

    process_machinery_concentrations={}
    for process in rba_session.get_processes():
        if process in list(input['Process_Name']):
            process_info=rba_session.get_process_information(process)
            complex_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                                machinery_composition=process_info["Composition"],
                                                                                proteomicsData=build_input_proteome_for_specific_kapp_estimation(protein_data, condition),
                                                                                proto_proteins=False)
            if (numpy.isfinite(complex_concentration))and(complex_concentration>0):
                process_machinery_concentrations[process]=complex_concentration

    effective_client_protein_aa_concentrations={}
    for protein in rba_session.get_proteins():
        protein_info=rba_session.get_protein_information(protein=protein)
        proto_id=protein_info['ProtoID']
        if proto_id in list(protein_data['ID']):
            concentration_proto_protein=protein_data.loc[protein_data['ID'] == proto_id, condition].values[0]
            protoprotein_isoform_distribution={}
            #{isoform1:0.1,isoform2:0.06}
            #protoprotein=isoform1+isoform2
            for isoform in protoprotein_isoprotein_map[proto_id]:
                protoprotein_isoform_distribution[isoform]=compartment_densities_and_pg.loc[protein_info['Compartment'],"Density"]
            compartment_fraction_sum_all_isoforms=sum(list(protoprotein_isoform_distribution.values()))
            for isoform in protoprotein_isoform_distribution.keys():
                effective_concentration_isoform=concentration_proto_protein*protoprotein_isoform_distribution[isoform]/compartment_fraction_sum_all_isoforms
                for process_required in protein_info['ProcessRequirements'].keys():
                    if process_required in effective_client_protein_aa_concentrations.keys():
                        effective_client_protein_aa_concentrations[process_required]+=effective_concentration_isoform*protein_info['ProcessRequirements'][process_required]
                    else:
                        effective_client_protein_aa_concentrations[process_required]=effective_concentration_isoform*protein_info['ProcessRequirements'][process_required]
    #{target_species_protein1:its_concentration,}
    #A:0.2 B:1 C:0.1  n_A*0.2+n_B*1+n_C*0.1 n_AA=100 Transl:100 Chap=10
    protein_targets={rba_session.get_target_information(i)["TargetEntity"]:rba_session.get_current_parameter_value(rba_session.get_target_information(i)["TargetParameterID"]) for i in rba_session.get_targets() if rba_session.get_target_information(i)["TargetEntity"] in rba_session.get_proteins()}
    for protein in protein_targets.keys():
        protein_info=rba_session.get_protein_information(protein=protein)
        concentration_target_protein=protein_targets[protein]
        for process_required in protein_info['ProcessRequirements'].keys():
            if process_required in effective_client_protein_aa_concentrations.keys():
                effective_client_protein_aa_concentrations[process_required]+=concentration_target_protein*protein_info['ProcessRequirements'][process_required]
            else:
                effective_client_protein_aa_concentrations[process_required]=concentration_target_protein*protein_info['ProcessRequirements'][process_required]
    
    efficiencies_processing_machineries={i:growth_rate*effective_client_protein_aa_concentrations[i]/process_machinery_concentrations[i] for i in effective_client_protein_aa_concentrations.keys() if i in process_machinery_concentrations.keys()}
    
    median_process_efficiency=numpy.median(numpy.array(list(efficiencies_processing_machineries.values())))

    ### NUCLEOTIDE ASSEMBLY PROCESS KAPPS ###        
    if fit_nucleotide_assembly_machinery:
        machinery_production_fluxes=determine_nucleotide_synthesis_machinery_demand(rba_session)
        for process in machinery_production_fluxes.keys():
            if (numpy.isfinite(machinery_production_fluxes[process]))and(machinery_production_fluxes[process]>0):
                if process not in process_machinery_concentrations.keys():
                    process_info=rba_session.get_process_information(process)
                    complex_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                                        machinery_composition=process_info["Composition"],
                                                                                        proteomicsData=build_input_proteome_for_specific_kapp_estimation(protein_data, condition),
                                                                                        proto_proteins=False)
                    if (numpy.isfinite(complex_concentration))and(complex_concentration>0):
                        efficiencies_processing_machineries[process]=machinery_production_fluxes[process]/complex_concentration
              
    for process in rba_session.get_processes():
        if process in list(input['Process_Name']):
            if process not in efficiencies_processing_machineries.keys():
                efficiencies_processing_machineries[process]=median_process_efficiency
    
    process_efficiencies = pandas.DataFrame()
    for process_name in efficiencies_processing_machineries.keys():
        process_info=rba_session.get_process_information(process_name)
        process_efficiencies.loc[process_name, 'Process'] = process_info["ID"]
        process_efficiencies.loc[process_name, 'Parameter'] = str(process_info["ID"]+'_apparent_efficiency')
        process_efficiencies.loc[process_name, 'Value'] = efficiencies_processing_machineries[process_name]
    return(process_efficiencies)

# Maybe reformulate on protein level#
def determine_apparent_process_efficiencies(growth_rate, input, rba_session,compartment_densities_and_PGs,total_amino_acid_abundance_in_proteome, protein_data, condition,fit_nucleotide_assembly_machinery=False):
    process_efficiencies = pandas.DataFrame()
    for i in input.index:
        process_ID = input.loc[i, 'Process_ID']
        process_name = input.loc[i, 'Process_Name']
        process_client_compartments = input.loc[i, 'Client_Compartments'].split(' , ')
        constituting_proteins = {rba_session.get_protein_information(protein=i)['ProtoID']: rba_session.get_protein_information(protein=i)['AAnumber'] for i in rba_session.get_process_information(process=process_name)['Composition'].keys()}
        Total_client_fraction = sum([compartment_densities_and_PGs.loc[i,"Density"] for i in process_client_compartments])
        
        #MACHINERY COST???
        n_AAs_in_machinery = 0
        machinery_size = 0
        for i in constituting_proteins.keys(): # for protein in complex
            if i in protein_data['ID']:
                n_AAs_in_machinery += protein_data.loc[protein_data['ID'] == i, condition].values[0] * protein_data.loc[protein_data['ID'] == i, 'AA_residues'].values[0]
                machinery_size += constituting_proteins[i] 
        # right reference amount?

        if n_AAs_in_machinery > 0:
            relative_Protein_fraction_of_machinery = n_AAs_in_machinery / total_amino_acid_abundance_in_proteome #fraction of proteome, made up of complex nAA_complex/nAA_total
            #mmol_AA/(mmol_AA*h)
            specific_capacity = growth_rate*Total_client_fraction/relative_Protein_fraction_of_machinery # k_spec * nAA_complex/nAA_total = mu * nAA_clients/nAA_total 
                                                                                                         # k_spec = 
            apparent_capacity = specific_capacity*machinery_size #mmol_AA/(mmol_Mach*h)
            process_efficiencies.loc[process_name, 'Process'] = process_ID
            process_efficiencies.loc[process_name, 'Parameter'] = str(process_ID+'_apparent_efficiency')
            process_efficiencies.loc[process_name, 'Value'] = apparent_capacity
    
    ### DEFAULT PROCESS KAPP ###        
    median_process_efficiency=numpy.median(numpy.array(process_efficiencies["Value"]))
    for i in input.index:
        process_ID = input.loc[i, 'Process_ID']
        process_name = input.loc[i, 'Process_Name']
        if not process_name in process_efficiencies.index:
            process_efficiencies.loc[process_name, 'Process'] = process_ID
            process_efficiencies.loc[process_name, 'Parameter'] = str(process_ID+'_apparent_efficiency')
            process_efficiencies.loc[process_name, 'Value'] = median_process_efficiency

    ### NUCLEOTIDE ASSEMBLY PROCESS KAPPS ###        
    if fit_nucleotide_assembly_machinery:
        machinery_production_fluxes=determine_nucleotide_synthesis_machinery_demand(rba_session)
        #
        for machinery in machinery_production_fluxes.keys():
            process_info=rba_session.get_process_information(process=machinery)
            #### rba_session proteome <-- input ####
            stoichiometrically_scaled_subunit_concentrations=[]
            subunit_stoichiometries=[]
            for su in process_info["Composition"].keys():
                protoProteinID=rba_session.get_protein_information(protein=su)["ProtoID"]
                if protoProteinID in protein_data['ID']:
                    copy_number=protein_data.loc[protein_data['ID'] == protoProteinID, condition].values[0]
                    scaled_copy_number=copy_number/process_info["Composition"][su]
                    if not pandas.isna(scaled_copy_number):
                        stoichiometrically_scaled_subunit_concentrations.append(scaled_copy_number)
                        subunit_stoichiometries.append(process_info["Composition"][su])
            if len(stoichiometrically_scaled_subunit_concentrations)>0:
                machinery_concentration=weighted_geometric_mean(data=stoichiometrically_scaled_subunit_concentrations,weights=subunit_stoichiometries)
            ####    
                #machinery_concentration=gmean(stoichiometrically_scaled_subunit_concentrations)
                apparent_process_efficiency=machinery_production_fluxes[machinery]/machinery_concentration
                process_efficiencies.loc[machinery, 'Process'] = process_info["ID"]
                process_efficiencies.loc[machinery, 'Parameter'] = str( process_info["ID"]+'_apparent_efficiency')
                process_efficiencies.loc[machinery, 'Value'] = apparent_process_efficiency
    return(process_efficiencies)


def determine_nucleotide_synthesis_machinery_demand(rba_session):
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
    #FluxDistribution.to_csv(output_dir+'/Calib_FluxDist_'+condition+'_.csv', sep=';')

    # describe what setting mean
    ### INCORPORATE ANA´s REMARK ON NOT ELIMINATING ISOENZYMES WITH SUs ONLY PRESENT IN THIS ISOENZYMES ###
    # 2: Determine list of all pre_selected isoenzymes --> "pre_selected_enzymes"#
    ProtoProteinMap = rba_session.ModelStructure.ProteinInfo.return_protein_iso_form_map()
    #{'Protein':['Protein__loc__comp']}
    measured_proteins_isoform_map = {p_ID: ProtoProteinMap[p_ID] for p_ID in list(proteomicsData['ID']) if p_ID in list(ProtoProteinMap.keys())}
        
    # identify all model reactions, associated with the measured proteins
    measured_proteins_reaction_map = determine_reactions_associated_with_measured_proto_protein(measured_proteins_isoform_map=measured_proteins_isoform_map,
                                                                                                        rba_session=rba_session)

    # for each protorxn generate a list of chosen isorxns to proceed with
    pre_selected_isoreactions=pre_select_iso_reactions(measured_proteins_reaction_map=measured_proteins_reaction_map,
                                                     rba_session=rba_session,
                                                     chose_most_quantified=chose_most_likely_isoreactions,
                                                     keep_isorxns_specific_to_quantified_proteins=True)

    # generate list of (iso)enzymes, associated to any selected isoreaction
    # pre_selected_isoreactions={R_A:[R_A_duplicate_1,R_A_duplicate_2],R_B:[R_B_duplicate_1],...,R_Z:[R_Z_duplicate_5]}
    pre_selected_enzymes=[] #[R_A_duplicate_1_enzyme,R_A_duplicate_2_enzyme,R_B_duplicate_1_enzyme,...R_Z_duplicate_5_enzyme]
    for proto_rxn in pre_selected_isoreactions.keys():
        for iso_rxn in pre_selected_isoreactions[proto_rxn]:
            if rba_session.get_reaction_information(iso_rxn) is not None:
                respective_enzyme=rba_session.get_reaction_information(iso_rxn)["Enzyme"]
                if respective_enzyme not in pre_selected_enzymes:
                    pre_selected_enzymes.append(respective_enzyme)

    # 3: Determine list of all reactions with nonzero/quantified flux in FBA --> "nonzero_reactions"#
    nonzero_flux_reactions=[]
    for reaction in FluxDistribution.index:
        if numpy.isfinite(FluxDistribution.loc[reaction,'FluxValues']):
            if FluxDistribution.loc[reaction,'FluxValues']!=0:
                nonzero_flux_reactions.append(reaction)

    # 4: Generate dictionary with all model_enzymes with nonzero concentration in data 
    # and their associated nonzero fba-flux proto reaction#
    nonzero_concentration_enzymes_with_associated_fba_flux={} #{R_A_duplicate_1_enzyme:R_A,R_A_duplicate_2_enzyme:R_A,...,R_Z_duplicate_2_enzyme:R_Z}
    for enzyme in pre_selected_enzymes:
        associated_reaction=rba_session.get_enzyme_information(enzyme)["Reaction"]
        associated_fba_flux_carrying_reaction=determine_fba_flux_carrying_isoreaction(rba_session=rba_session,
                                                                                      reaction_id=associated_reaction,
                                                                                      flux_distribution=FluxDistribution)
        if associated_fba_flux_carrying_reaction is not None:
            nonzero_concentration_enzymes_with_associated_fba_flux[enzyme]=associated_fba_flux_carrying_reaction

    # 5: ...#
    overview_out = pandas.DataFrame()
    for rxn in list(set(list(nonzero_concentration_enzymes_with_associated_fba_flux.values()))): # for unique reaction
        
        # 5.1: Determine model_enzymes constituting pseudo complex, associated with reaction#
        if rba_session.get_reaction_information(rxn) is None:
            continue
        rxn_enzyme=rba_session.get_reaction_information(rxn)["Enzyme"]
        if rba_session.get_enzyme_information(rxn_enzyme) is None:
            continue

        #generate list of all isoenzymes, which were preselected catalyzing rxn:
        associated_nonzero_concentration_isoenzymatic_complexes=[i for i in list([rxn_enzyme]+rba_session.get_enzyme_information(rxn_enzyme)["Isozymes"]) if i in list(pre_selected_enzymes)]
        
        model_enzymes_associated_with_pseudo_complex={} 
        # {R_A_duplicate_2_enzyme:[R_A_duplicate_2_enzyme,R_X_enzyme]} 
        #  R_X_enzyme and R_A_duplicate_2_enzyme identical complex (if both are nonzero measured abundant and carry flux in FBA)
        for isoenzymatic_complex in associated_nonzero_concentration_isoenzymatic_complexes:
            model_enzymes_associated_with_pseudo_complex[isoenzymatic_complex]=[isoenzymatic_complex] #  scaling by flux proportion would go here
            ### ??? identical enzymes isoform specific? ###
            for identical_composition_enzyme in rba_session.get_enzyme_information(rxn_enzyme)["EnzymesWithIdenticalSubunitComposition"]:
                if identical_composition_enzyme in nonzero_concentration_enzymes_with_associated_fba_flux.keys():
                    model_enzymes_associated_with_pseudo_complex[isoenzymatic_complex].append(identical_composition_enzyme)
        
        # 5.2 ...#
        total_flux_dict={}
        individual_constituent_concentrations={i:[] for i in model_enzymes_associated_with_pseudo_complex.keys()}
        #individual_constituent_concentrations={R_A_duplicate_2_enzyme:[concentration_R_A_duplicate_2_enzyme]} 
        #total_flux_dict={R_A:J_A,R_X:J_X}
        for associated_pseudocomplex_enzyme in model_enzymes_associated_with_pseudo_complex.keys():
            #Calculate complex concentration
            respective_composition=rba_session.get_enzyme_information(associated_pseudocomplex_enzyme)["Subunits"]
            respective_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                             machinery_composition=respective_composition,
                                                                             proteomicsData=proteomicsData,
                                                                             proto_proteins=False)
            # if finite and nonzero append complex concentration
            if (respective_concentration!=0) and (numpy.isfinite(respective_concentration)):
                individual_constituent_concentrations[associated_pseudocomplex_enzyme].append(respective_concentration)
            # iterate over other catalytic activities of complex
            for identical_composition_enzyme in model_enzymes_associated_with_pseudo_complex[associated_pseudocomplex_enzyme]:
                if identical_composition_enzyme in nonzero_concentration_enzymes_with_associated_fba_flux.keys():
                # if protein measured and fba flux of other catalytic activity:
                    # identify rxn-id in flux-distribution
                    associated_fba_rxn=nonzero_concentration_enzymes_with_associated_fba_flux[identical_composition_enzyme]
                    if associated_fba_rxn not in total_flux_dict.keys():
                        # if flux not already considered record absolute value of flux
                        total_flux_dict[associated_fba_rxn]=abs(FluxDistribution.loc[associated_fba_rxn,'FluxValues'])

        #calculate pseudocomplex concentration, associated to rxn in question as geometric mean of constituting complex abundances
        pseudocomplex_constituent_concentrations=[individual_constituent_concentrations[i][0] for i in individual_constituent_concentrations.keys() if len(individual_constituent_concentrations[i])!=0]
        if len(pseudocomplex_constituent_concentrations)>0:
            # maybe sum #
            concentration_pseudo_complex=gmean(numpy.array(pseudocomplex_constituent_concentrations))
        else:
            concentration_pseudo_complex=numpy.nan

        # 5.4 ...#
        #calculate pseudocomplex flux, associated to rxn in question as mean of constituting fluxes
        #flux_pseudo_complex=gmean(numpy.array(list(total_flux_dict.values())))
        flux_pseudo_complex=sum(total_flux_dict.values())/len(list(total_flux_dict.keys()))
        
        # 5.5 ...#
        # decide if backward or forward kapp, based on fba flux direction
        flux_direction=numpy.nan
        if FluxDistribution.loc[rxn,'FluxValues'] > 0.0:
            flux_direction=1
        elif FluxDistribution.loc[rxn,'FluxValues'] < 0.0:
            flux_direction=-1
        if (flux_pseudo_complex!=0) and numpy.isfinite(flux_pseudo_complex) and (concentration_pseudo_complex!=0) and numpy.isfinite(concentration_pseudo_complex):
            # if flux and concentration of pseudo complex finite and nonzero:
            kapp_pseudo_complex=flux_pseudo_complex/concentration_pseudo_complex
            # 5.6 write down results
            for considered_isoenzyme in individual_constituent_concentrations.keys():
                respective_reaction=rba_session.get_enzyme_information(considered_isoenzyme)["Reaction"]
                overview_out.loc[respective_reaction,"Enzyme_ID"]=considered_isoenzyme
                overview_out.loc[respective_reaction, 'Flux_FBA'] = flux_pseudo_complex
                overview_out.loc[respective_reaction, 'Flux'] = flux_direction
                overview_out.loc[respective_reaction, 'Concentration'] = concentration_pseudo_complex
                overview_out.loc[respective_reaction,'Kapp'] = kapp_pseudo_complex
                overview_out.loc[respective_reaction, 'Comment'] = 'estimated'
                overview_out.loc[respective_reaction,"Pseudo complex members"]=" , ".join(model_enzymes_associated_with_pseudo_complex)
            
            # 5.7 deal with (eliminated) isoenzymes#
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
                        
    # 6: impose on enzymes with identical composition, which do not have a flux predicted in FBA#
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
    #if store_output:
    #    if condition is not None:
    #        overview_out.to_csv(output_dir+'/SpecKapp_Network_overview_'+condition+'_.csv', sep=';')
    #    else:
    #        overview_out.to_csv(output_dir+'/SpecKapp_Network_overview_.csv', sep=';')
    # 9: ...#
    return({"Overview":overview_out,"Flux_Distribution":FluxDistribution})


def estimate_specific_enzyme_efficiencies_new(rba_session, 
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
                                          condition=None, 
                                          store_output=True,
                                          rxns_to_ignore_when_parsimonious=[],
                                          use_bm_flux_of_one=False,
                                          output_dir="",
                                          pseudocomplex_method='old'):
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
    #FluxDistribution.to_csv(output_dir+'/Calib_FluxDist_'+condition+'_.csv', sep=';')

    # describe what setting mean
    ### INCORPORATE ANA´s REMARK ON NOT ELIMINATING ISOENZYMES WITH SUs ONLY PRESENT IN THIS ISOENZYMES ###
    # 2: Determine list of all pre_selected isoenzymes --> "pre_selected_enzymes"#
    ProtoProteinMap = rba_session.ModelStructure.ProteinInfo.return_protein_iso_form_map()
    #{'Protein':['Protein__loc__comp']}
    measured_proteins_isoform_map = {p_ID: ProtoProteinMap[p_ID] for p_ID in list(proteomicsData['ID']) if p_ID in list(ProtoProteinMap.keys())}
        
    # identify all model reactions, associated with the measured proteins
    measured_proteins_reaction_map = determine_reactions_associated_with_measured_proto_protein(measured_proteins_isoform_map=measured_proteins_isoform_map,
                                                                                                        rba_session=rba_session)

    # for each protorxn generate a list of chosen isorxns to proceed with
    pre_selected_isoreactions=pre_select_iso_reactions(measured_proteins_reaction_map=measured_proteins_reaction_map,
                                                     rba_session=rba_session,
                                                     chose_most_quantified=chose_most_likely_isoreactions,
                                                     keep_isorxns_specific_to_quantified_proteins=True)

    # generate list of (iso)enzymes, associated to any selected isoreaction
    # pre_selected_isoreactions={R_A:[R_A_duplicate_1,R_A_duplicate_2],R_B:[R_B_duplicate_1],...,R_Z:[R_Z_duplicate_5]}
    pre_selected_enzymes=[] #[R_A_duplicate_1_enzyme,R_A_duplicate_2_enzyme,R_B_duplicate_1_enzyme,...R_Z_duplicate_5_enzyme]
    for proto_rxn in pre_selected_isoreactions.keys():
        for iso_rxn in pre_selected_isoreactions[proto_rxn]:
            if rba_session.get_reaction_information(iso_rxn) is not None:
                respective_enzyme=rba_session.get_reaction_information(iso_rxn)["Enzyme"]
                if respective_enzyme not in pre_selected_enzymes:
                    pre_selected_enzymes.append(respective_enzyme)

    # 3: Determine list of all reactions with nonzero/quantified flux in FBA --> "nonzero_reactions"#
    nonzero_flux_reactions=[]
    for reaction in FluxDistribution.index:
        if numpy.isfinite(FluxDistribution.loc[reaction,'FluxValues']):
            if FluxDistribution.loc[reaction,'FluxValues']!=0:
                nonzero_flux_reactions.append(reaction)

    #############################
    overview_out = pandas.DataFrame()
    for reaction in nonzero_flux_reactions:
        if reaction not in rba_session.get_reactions():
            continue
        reaction_enzyme=rba_session.get_reaction_information(reaction)['Enzyme']
        if reaction_enzyme not in rba_session.get_enzymes():
            continue
        all_preselected_reaction_isoenzymes=[i for i in list([reaction_enzyme]+rba_session.get_enzyme_information(reaction_enzyme)['Isozymes']) if i in pre_selected_enzymes]
        if len(all_preselected_reaction_isoenzymes)<1:
            continue
        
        reaction_flux=FluxDistribution.loc[reaction,'FluxValues']

        all_isoenzyme_concentrations_specific_network_scaled_by_flux={} #isoenzymes of rxn and their allocated share of concentration{}

        all_fluxes_specific_network={reaction:abs(reaction_flux)} # all rxns in network and their fluxes (needed for old way)
        all_isoenzyme_concentrations_specific_network={} #isoenzymes of rxn and their concentration (needed for old way)
        all_catalytic_activities_in_pseudo_complex=[] #name of all (iso)enzymes in network (needed only for output info)

        for isoenzyme in all_preselected_reaction_isoenzymes: # all (preselected) isoenzymes of rxn 
            all_catalytic_activities_in_pseudo_complex.append(isoenzyme)
            total_flux_isoenzyme_all_catalytic_activities=abs(reaction_flux) # for concentration partitioning, other fluxes added later
            # iterate over other catalytic activities of isoenzyme:
            for other_catalytic_activity_enzyme in rba_session.get_enzyme_information(isoenzyme)['IdenticalEnzymes']:
                all_catalytic_activities_in_pseudo_complex.append(other_catalytic_activity_enzyme)
                associated_reaction=rba_session.get_enzyme_information(other_catalytic_activity_enzyme)["Reaction"]
                # iterate over all isorxns of it:
                for isoreaction in list([associated_reaction]+rba_session.get_reaction_information(associated_reaction)['Twins']):
                    if isoreaction in nonzero_flux_reactions: # does it carry flux?
                        isoreaction_flux=FluxDistribution.loc[isoreaction,'FluxValues']
                        all_fluxes_specific_network[isoreaction]=abs(isoreaction_flux)
                        total_flux_isoenzyme_all_catalytic_activities+=abs(isoreaction_flux) # added for concentration partitioning
            #determine isoenzyme-concentration
            isoenzyme_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                             machinery_composition=rba_session.get_enzyme_information(isoenzyme)["Subunits"],
                                                                             proteomicsData=proteomicsData,
                                                                             proto_proteins=False)
            all_isoenzyme_concentrations_specific_network[isoenzyme]=isoenzyme_concentration
            #record allocated concentration partition by flux partitioning:
            all_isoenzyme_concentrations_specific_network_scaled_by_flux[isoenzyme]=isoenzyme_concentration*abs(reaction_flux)/(total_flux_isoenzyme_all_catalytic_activities)

        if pseudocomplex_method=='old':
            concentration_pseudo_complex=gmean(numpy.array(list(all_isoenzyme_concentrations_specific_network.values())))
            flux_pseudo_complex=numpy.mean(list(all_fluxes_specific_network.values()))

        if pseudocomplex_method=='new':
            concentration_pseudo_complex=sum(list(all_isoenzyme_concentrations_specific_network_scaled_by_flux.values()))
            flux_pseudo_complex=abs(reaction_flux)

        # 5.5 ...#
        if (flux_pseudo_complex!=0) and numpy.isfinite(flux_pseudo_complex) and (concentration_pseudo_complex!=0) and numpy.isfinite(concentration_pseudo_complex):
            # decide if backward or forward kapp, based on fba flux direction
            flux_direction=numpy.nan
            if reaction_flux > 0.0:
                flux_direction=1
            elif reaction_flux < 0.0:
                flux_direction=-1

            kapp_pseudo_complex=flux_pseudo_complex/concentration_pseudo_complex
            # 5.6 write down results
            for isoenzyme in all_preselected_reaction_isoenzymes:
                respective_reaction=rba_session.get_enzyme_information(isoenzyme)["Reaction"]
                overview_out.loc[respective_reaction,"Enzyme_ID"]=isoenzyme
                overview_out.loc[respective_reaction, 'Flux_FBA'] = flux_pseudo_complex
                overview_out.loc[respective_reaction, 'Flux'] = flux_direction
                overview_out.loc[respective_reaction, 'Concentration'] = concentration_pseudo_complex
                overview_out.loc[respective_reaction,'Kapp'] = kapp_pseudo_complex
                overview_out.loc[respective_reaction, 'Comment'] = 'estimated'
                overview_out.loc[respective_reaction,"Pseudo complex members"]=" , ".join(all_catalytic_activities_in_pseudo_complex)

            not_preselected_reaction_isoenzymes=[i for i in list([reaction_enzyme]+rba_session.get_enzyme_information(reaction_enzyme)['Isozymes']) if i not in all_preselected_reaction_isoenzymes]
            for not_selected_isozyme in not_preselected_reaction_isoenzymes:
                respective_reaction=rba_session.get_enzyme_information(not_selected_isozyme)["Reaction"]
                if respective_reaction not in overview_out.index:
                    overview_out.loc[respective_reaction,"Enzyme_ID"]=not_selected_isozyme
                    overview_out.loc[respective_reaction, 'Comment'] = 'eliminated isoenzyme'
                    overview_out.loc[respective_reaction, 'Flux'] = flux_direction
                    if zero_on_all_isoreactions:
                        overview_out.loc[respective_reaction,'Kapp'] = 0.0
                    elif impose_on_all_isoreactions:
                        overview_out.loc[respective_reaction,'Kapp'] = kapp_pseudo_complex

    # 6: impose on enzymes with identical composition, which do not have a flux predicted in FBA#
    if impose_on_identical_enzymes:
        for reaction in rba_session.get_reactions():
            reaction_already_considered_in_previous_estimation=False
            enzyme_kapps_to_infer_kapp_from={}
            for isoreaction in list([reaction]+rba_session.get_reaction_information(reaction)["Twins"]):
                if isoreaction in list(overview_out.index):
                    reaction_already_considered_in_previous_estimation=True
                    break
                isoenzyme=rba_session.get_reaction_information(isoreaction)["Enzyme"]
                if isoenzyme not in rba_session.get_enzymes():
                    continue
                for other_catalytic_activity_enzyme in rba_session.get_enzyme_information(isoenzyme)['IdenticalEnzymes']:
                    if other_catalytic_activity_enzyme in list(overview_out["Enzyme_ID"]):
                        respective_kapp=float(overview_out.loc[overview_out["Enzyme_ID"]==other_catalytic_activity_enzyme,'Kapp'].values[0])
                        if numpy.isfinite(respective_kapp):
                            if respective_kapp != 0:
                                enzyme_kapps_to_infer_kapp_from[other_catalytic_activity_enzyme]=respective_kapp                        
            if reaction_already_considered_in_previous_estimation:
                continue
            if len(list(enzyme_kapps_to_infer_kapp_from.keys()))>0:
                for isoreaction in list([reaction]+rba_session.get_reaction_information(reaction)["Twins"]):
                    overview_out.loc[isoreaction,"Enzyme_ID"]=rba_session.get_reaction_information(isoreaction)["Enzyme"]
                    overview_out.loc[isoreaction,'Kapp'] = numpy.mean(list(enzyme_kapps_to_infer_kapp_from.values()))
                    overview_out.loc[isoreaction, 'Comment'] = 'inferred from other enzymes'
                    overview_out.loc[isoreaction,"Enzymes inferred from"]=" , ".join(list(enzyme_kapps_to_infer_kapp_from.keys()))

    #############################
 
    # 7: ...#
    rba_session.rebuild_from_model()
    rba_session.set_medium(rba_session.Medium)

    overview_out.sort_index(inplace=True)

    # 8: ...#
    #if store_output:
    #    if condition is not None:
    #        overview_out.to_csv(output_dir+'/SpecKapp_Network_overview_'+condition+'_.csv', sep=';')
    #    else:
    #        overview_out.to_csv(output_dir+'/SpecKapp_Network_overview_.csv', sep=';')
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
        original_density_constraint_signs=rba_session.Problem.get_constraint_types(constraints=[i for i in rba_session.get_density_constraints() if i in rba_session.Problem.LP.row_names])
        rba_session.Problem.set_constraint_types({i:"E" for i in rba_session.get_density_constraints() if i in rba_session.Problem.LP.row_names})

        solved1=rba_session.solve()
        if solved1:
            derive_bm_from_rbasolution=True
            derive_bm_from_targets=False
            rba_session.Problem.set_constraint_types(original_density_constraint_signs)
        else:
            print("{} - Solution with equality density  and original medium not obtained - Status: {}".format(condition,rba_session.Problem.SolutionStatus))
            rba_session.set_medium({i:100.0 for i in original_medium.keys()})
            solved2=rba_session.solve()
            rba_session.Problem.set_constraint_types(original_density_constraint_signs)
            if solved2:
                derive_bm_from_rbasolution=True
                derive_bm_from_targets=False
            else:
                print("{} - Solution with equality density not obtained - Status: {}".format(condition,rba_session.Problem.SolutionStatus))
                rba_session.set_growth_rate(mu)
                solved3=rba_session.solve()
                if solved3:
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
        # explain
        #set use_bm_flux_of_one to True by defualt (remove option at all)?
        if use_bm_flux_of_one:
            # if 1: mu predicted = mu measured (reformulate to have perfect correspondence when BM-flux = mu_mrasured)
            # 1.0 --> mu measured
            # if FBA overpredicts imposed growth-rate:
            if BMfluxOld >= 1.0:
                #impose imposed growth rate
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
    proto_rxn_dict = {} #{R_A: {R_A_duplicate_1: 3 , R_A_duplicate_2: 2}} 
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
    for prx in proto_rxn_dict.keys(): # iterate over proto rxns in {R_A: {R_A_duplicate_1: 3 , R_A_duplicate_2: 2 , R_A_duplicate_3: 0}} (R_A in this case)
        unique_subunit_dict={} # how many unique (ignoring subunit-stoichiometry) subunits are in complex
        for irx in proto_rxn_dict[prx].keys(): # iterate over isorxns, associated to protorxn
            enzyme = rba_session.ModelStructure.ReactionInfo.Elements[irx]['Enzyme']
            unique_subunit_dict[irx] = len(list(rba_session.get_enzyme_information(enzyme)["Subunits"].keys()))
        if chose_most_quantified:
            max_coverage = max([proto_rxn_dict[prx][i]/unique_subunit_dict[i] for i in proto_rxn_dict[prx].keys()]) # maximum coverage
            list_isorxns = [i for i in proto_rxn_dict[prx].keys() if proto_rxn_dict[prx][i]/unique_subunit_dict[i] == max_coverage]

            #if len(list_isorxns)>1:
            #    max_SU_number=max([unique_subunit_dict[i] for i in list_isorxns])
            #    selected=[i for i in list_isorxns if unique_subunit_dict[i]==max_SU_number]
            #else: #do this by defaul omit if condition
            #    selected=list_isorxns

            selected=list_isorxns
            
            ## proposal of Ana in Oberseminar:
            if keep_isorxns_specific_to_quantified_proteins:
                # generate list of (iso)reactions, where one subunit is the only "function" of a measured protein
                # which are the only activity of a measured protein, 
                # in line below
                isoreactions_based_on_measured_proteins_with_only_this_activity=[measured_proteins_reaction_map[i][0] for i in measured_proteins_reaction_map.keys() if len(measured_proteins_reaction_map[i])==1]
                for irx in proto_rxn_dict[prx].keys():
                    if not irx in selected:
                        if irx in isoreactions_based_on_measured_proteins_with_only_this_activity:
                            selected.append(irx)
        else: #only exlude isenzymes with no measured subunit
            selected = [i for i in proto_rxn_dict[prx].keys() if proto_rxn_dict[prx][i] != 0]
        selected.sort()
        out[prx] = selected #{R_A: [R_A_duplicate_1: 3 , R_A_duplicate_2: 2 ]} example for whole process (all eventualities of elimination)
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
            if subunit in rba_session.get_proteins():
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

def efficiency_correction_new(enzyme_efficiencies,
                               simulation_results,
                               protein_data,
                               rba_session,
                               condition_to_look_up,
                               default_enzyme_efficiencies,
                               tolerance=2,
                               n_th_root_mispred=2,
                               process_efficiencies=None,
                               correct_default_kapp_enzymes=False,
                               only_consider_misprediction_for_predicted_nonzero_enzymes=False):

    enzyme_efficiencies_out=enzyme_efficiencies.copy() 
    process_efficiencies_out=process_efficiencies.copy()

    # for each protoprotein, determine predicted anbundance as sum over all predicted isoforms
    # {Protein1:Abundance1,Protein2:Abundance2}
    predicted_proto_protein_quantities={} 
    for i in simulation_results["Proteins"].index:
        proto_protein_ID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if proto_protein_ID in predicted_proto_protein_quantities:
            predicted_proto_protein_quantities[proto_protein_ID]+=simulation_results["Proteins"].loc[i,condition_to_look_up]
        else:
            predicted_proto_protein_quantities[proto_protein_ID]=simulation_results["Proteins"].loc[i,condition_to_look_up]

    predicted_proteome=pandas.DataFrame()
    for gene in predicted_proto_protein_quantities.keys():
        predicted_proteome.loc[gene,"ID"]=gene
        predicted_proteome.loc[gene,"copy_number"]=predicted_proto_protein_quantities[gene]
    
    process_correction_coefficients={}
    for process in rba_session.get_processes():
        machinery_composition=rba_session.get_process_information(process)['Composition']
        measured_machinery_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                                machinery_composition=machinery_composition,
                                                                                proteomicsData=protein_data,
                                                                                proto_proteins=False)

        predicted_machinery_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                                machinery_composition=machinery_composition,
                                                                                proteomicsData=predicted_proteome,
                                                                                proto_proteins=False)
        machinery_misprediction_coefficient=1
        if (measured_machinery_concentration != 0) and (predicted_machinery_concentration != 0):
            if (numpy.isfinite(measured_machinery_concentration)) and (numpy.isfinite(predicted_machinery_concentration)):
                #machinery_misprediction_coefficient=predicted_machinery_concentration/measured_machinery_concentration
                machinery_misprediction_coefficient=numpy.power(predicted_machinery_concentration/measured_machinery_concentration,1/n_th_root_mispred)

        if process in list(process_efficiencies.index):
            if numpy.isfinite(machinery_misprediction_coefficient):
                if machinery_misprediction_coefficient!=1:
                    current_efficiency=process_efficiencies.loc[process,"Value"]
                    updated_efficiency=current_efficiency*machinery_misprediction_coefficient
                    if tolerance is None:
                        process_efficiencies_out.loc[process,"Value"]=updated_efficiency
                        process_correction_coefficients[process]=machinery_misprediction_coefficient
                    else:
                        if abs(numpy.log(tolerance)) <= abs(numpy.log(machinery_misprediction_coefficient)):
                            process_efficiencies_out.loc[process,"Value"]=updated_efficiency
                            process_correction_coefficients[process]=machinery_misprediction_coefficient

#    numpy.power(,1/n_th_root_mispred)
     # {(Iso)Enzyme: correction_coeff} 
    enzyme_correction_coefficients_original={}
    for enzyme in rba_session.get_enzymes():
        associated_reaction=rba_session.get_enzyme_information(enzyme)["Reaction"]
        # find out wheter it is forward or backward specific kapp to add :
        if associated_reaction in simulation_results["Reactions"].index:
            flux=simulation_results["Reactions"].loc[associated_reaction,condition_to_look_up]
            if flux!=0:
                enzyme_composition=rba_session.get_enzyme_information(enzyme)['Subunits']
                measured_machinery_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                                        machinery_composition=enzyme_composition,
                                                                                        proteomicsData=protein_data,
                                                                                        proto_proteins=False)

                predicted_machinery_concentration=determine_machinery_concentration_by_weighted_geometric_mean(rba_session=rba_session,
                                                                                        machinery_composition=enzyme_composition,
                                                                                        proteomicsData=predicted_proteome,
                                                                                        proto_proteins=False)
                machinery_misprediction_coefficient=1
                if (measured_machinery_concentration != 0) and (predicted_machinery_concentration != 0):
                    if (numpy.isfinite(measured_machinery_concentration)) and (numpy.isfinite(predicted_machinery_concentration)):
                        machinery_misprediction_coefficient=predicted_machinery_concentration/measured_machinery_concentration
                        #machinery_misprediction_coefficient=numpy.power(predicted_machinery_concentration/measured_machinery_concentration,1/n_th_root_mispred)
                enzyme_correction_coefficients_original[enzyme]=machinery_misprediction_coefficient

    enzymes_already_corrected=[]
    enzyme_correction_coefficients={}
    for enzyme_to_correct in enzyme_correction_coefficients_original.keys():
        if enzyme_to_correct in enzymes_already_corrected:
            continue
        correction_coeffs=[]
        #isozymes=list([enzyme_to_correct]+rba_session.get_enzyme_information(enzyme_to_correct)["Isozymes"])
        isozymes=[enzyme_to_correct]
        enzymes_already_corrected+=isozymes
        for isozyme in isozymes:
            if isozyme in enzyme_correction_coefficients_original.keys():
                correction_coeffs.append(enzyme_correction_coefficients_original[isozyme]) 
        #mean_correction_coeff=gmean(numpy.array(correction_coeffs))
        mean_correction_coeff=numpy.power(gmean(numpy.array(correction_coeffs)),1/n_th_root_mispred)
        if numpy.isfinite(mean_correction_coeff):
            for isozyme in isozymes:
                if isozyme in list(enzyme_efficiencies["Enzyme_ID"]):
                    current_efficiency=enzyme_efficiencies.loc[enzyme_efficiencies["Enzyme_ID"]==isozyme,"Kapp"].values[0]
                    updated_efficiency=current_efficiency*mean_correction_coeff
                    if tolerance is None:
                        enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==isozyme,"Kapp"]=updated_efficiency
                        enzyme_correction_coefficients[enzyme]=mean_correction_coeff
                    else:
                        if abs(numpy.log(tolerance)) <= abs(numpy.log(mean_correction_coeff)):
                            enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==isozyme,"Kapp"]=updated_efficiency
                            enzyme_correction_coefficients[enzyme]=mean_correction_coeff
                else:
                    if correct_default_kapp_enzymes:
                        #respective_reaction=rba_session.get_enzyme_information(isozyme)["Reaction"].split("_duplicate_")[0]
                        respective_reaction=rba_session.get_enzyme_information(isozyme)["Reaction"]
                        # find out wheter it is forward or backward specific kapp to add :
                        if respective_reaction in simulation_results["Reactions"].index:
                            flux=simulation_results["Reactions"].loc[respective_reaction,condition_to_look_up]
                            try:
                                if flux < 0:
                                    flux_direction=-1
                                    efficiency_parameter= rba_session.model.enzymes.enzymes._elements_by_id[isozyme].backward_efficiency
                                elif flux > 0:
                                    flux_direction=1
                                    efficiency_parameter= rba_session.model.enzymes.enzymes._elements_by_id[isozyme].forward_efficiency
                                else:
                                    continue
                            except:
                                continue
                            if efficiency_parameter in default_enzyme_efficiencies.keys():
                                current_efficiency=default_enzyme_efficiencies[efficiency_parameter]
                                updated_efficiency=current_efficiency*mean_correction_coeff
                                if tolerance is None:
                                    enzyme_correction_coefficients[enzyme]=mean_correction_coeff
                                    enzyme_efficiencies_out.loc[respective_reaction,"Kapp"]=updated_efficiency
                                    enzyme_efficiencies_out.loc[respective_reaction,"Enzyme_ID"]=isozyme
                                    enzyme_efficiencies_out.loc[respective_reaction,"Flux"]=flux_direction
                                    enzyme_efficiencies_out.loc[respective_reaction,"Comment"]="Corrected Default"
                                else:
                                    if abs(numpy.log(tolerance)) <= abs(numpy.log(mean_correction_coeff)):
                                        enzyme_correction_coefficients[enzyme]=mean_correction_coeff
                                        enzyme_efficiencies_out.loc[respective_reaction,"Kapp"]=updated_efficiency
                                        enzyme_efficiencies_out.loc[respective_reaction,"Enzyme_ID"]=isozyme
                                        enzyme_efficiencies_out.loc[respective_reaction,"Flux"]=flux_direction
                                        enzyme_efficiencies_out.loc[respective_reaction,"Comment"]="Corrected Default"

    return({"Kapps":enzyme_efficiencies_out,
            "ProcessEfficiencies":process_efficiencies_out,
            "Process_MispredictionFactors":process_correction_coefficients,
            "Enzyme_MispredictionFactors":enzyme_correction_coefficients})

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
                               only_consider_misprediction_for_predicted_nonzero_enzymes=False):
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
    """
    # for each protoprotein, determine predicted anbundance as sum over all predicted isoforms
    # {Protein1:Abundance1,Protein2:Abundance2}
    predicted_proto_protein_quantities={} 
    for i in simulation_results["Proteins"].index:
        proto_protein_ID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if proto_protein_ID in predicted_proto_protein_quantities:
            predicted_proto_protein_quantities[proto_protein_ID]+=simulation_results["Proteins"].loc[i,condition_to_look_up]
        else:
            predicted_proto_protein_quantities[proto_protein_ID]=simulation_results["Proteins"].loc[i,condition_to_look_up]
    
    # for each protoprotein generate list of isoforms
    # {Protein1:[Protein1_c,Protein1_m],Protein2:[Protein2_n]}
    proto_protein_isoform_map=rba_session.ModelStructure.ProteinInfo.return_protein_iso_form_map()
    #proto_protein_isoform_map={} 
    #for i in rba_session.get_proteins():
    #    proto_protein_ID=rba_session.get_protein_information(protein=i)["ProtoID"]
    #    if proto_protein_ID in proto_protein_isoform_map:
    #        proto_protein_isoform_map[proto_protein_ID].append(i)
    #    else:
    #        proto_protein_isoform_map[proto_protein_ID]=[i]
    

    enzyme_efficiencies_out=enzyme_efficiencies.copy() 
    process_efficiencies_out=process_efficiencies.copy()
    subunit_misprediction_factors_enzymes={} # {(Iso)Enzyme:[misprediction_SU1,misprediction_SU2,misprediction_SU3,...]}
    subunit_misprediction_factors_default_efficiency_enzymes={} # {(Iso)Enzyme:[misprediction_SU1,misprediction_SU2,misprediction_SU3,...]}
    subunit_misprediction_factors_processes={} # {Process:[misprediction_SU1,misprediction_SU2,misprediction_SU3,...]}
    subunit_stoichiometries_enzymes={} # {(Iso)Enzyme:[stoichiometry_SU1,stoichiometry_SU2,stoichiometry_SU3,...]}
    subunit_stoichiometries_default_efficiency_enzymes={} # {(Iso)Enzyme:[stoichiometry_SU1,stoichiometry_SU2,stoichiometry_SU3,...]}
    subunit_stoichiometries_processes={} 

    # for machineries determine subunit misprediction- and stoichiometry coeffs 
    # {MAchinery:[misprediction_SU1,misprediction_SU2,misprediction_SU3,...]}
    # {MAchinery:[stoichiometry_SU1,stoichiometry_SU2,stoichiometry_SU3,...]}
    for proto_protein_ID in predicted_proto_protein_quantities:
        if proto_protein_ID in proto_protein_isoform_map: # can proto protein be mapped to model protein-isofrms?
            if proto_protein_ID in protein_data["ID"]: # is protein quantified in measured data
                # measured and predicted concentrations:
                predicted_protein=predicted_proto_protein_quantities[proto_protein_ID]
                measured_protein=protein_data.loc[proto_protein_ID,"copy_number"]
                # check if predicted and measured abundance is nonzero and finite:
                if (predicted_protein>0) & (measured_protein>0) & (numpy.isfinite(predicted_protein)) & (numpy.isfinite(measured_protein)):
                    misprediction_coeff=predicted_protein/measured_protein #determine misprediction ratio
                    # iterate over protein isoforms:
                    for protein in proto_protein_isoform_map[proto_protein_ID]:
                        # iterate over associated enzymes
                        for enzyme in rba_session.get_protein_information(protein=protein)["associatedEnzymes"]:
                            # only continue with enzyme if enzyme (catalytic activity) 
                            # is predicted finite and nonzero in RBA-solution:
                            if only_consider_misprediction_for_predicted_nonzero_enzymes:
                                if enzyme not in simulation_results["Enzymes"].index:
                                    continue
                                elif numpy.isfinite(simulation_results["Enzymes"].loc[enzyme,condition_to_look_up])==False:
                                    continue
                                elif simulation_results["Enzymes"].loc[enzyme,condition_to_look_up]==0:
                                    continue
                            # if enzyme was previously assigned a specific efficiency:
                            if enzyme in list(enzyme_efficiencies["Enzyme_ID"]):
                                # {(Iso)Enzyme:[misprediction_SU1,misprediction_SU2,misprediction_SU3,...]}
                                # {(Iso)Enzyme:[stoichiometry_SU1,stoichiometry_SU2,stoichiometry_SU3,...]}
                                if enzyme not in subunit_misprediction_factors_enzymes.keys():
                                    subunit_misprediction_factors_enzymes[enzyme]=[misprediction_coeff]
                                    subunit_stoichiometries_enzymes[enzyme]=[rba_session.get_enzyme_information(enzyme)["Subunits"][protein]]
                                else:
                                    subunit_misprediction_factors_enzymes[enzyme].append(misprediction_coeff)
                                    subunit_stoichiometries_enzymes[enzyme].append(rba_session.get_enzyme_information(enzyme)["Subunits"][protein])
                            # if enzyme was previously assigned the default efficiency:
                            else:
                                # {(Iso)Enzyme:[misprediction_SU1,misprediction_SU2,misprediction_SU3,...]}
                                # {(Iso)Enzyme:[stoichiometry_SU1,stoichiometry_SU2,stoichiometry_SU3,...]}
                                if enzyme not in subunit_misprediction_factors_default_efficiency_enzymes.keys():
                                    subunit_misprediction_factors_default_efficiency_enzymes[enzyme]=[misprediction_coeff]
                                    subunit_stoichiometries_default_efficiency_enzymes[enzyme]=[rba_session.get_enzyme_information(enzyme)["Subunits"][protein]]
                                else:
                                    subunit_misprediction_factors_default_efficiency_enzymes[enzyme].append(misprediction_coeff)
                                    subunit_stoichiometries_default_efficiency_enzymes[enzyme].append(rba_session.get_enzyme_information(enzyme)["Subunits"][protein])
                        for process in rba_session.get_protein_information(protein=protein)["SupportsProcess"]:
                            if process in process_efficiencies.index:
                                # {Process:[misprediction_SU1,misprediction_SU2,misprediction_SU3,...]}
                                # {Process:[stoichiometry_SU1,stoichiometry_SU2,stoichiometry_SU3,...]}
                                if process not in subunit_misprediction_factors_processes.keys():
                                    subunit_misprediction_factors_processes[process]=[misprediction_coeff]
                                    subunit_stoichiometries_processes[process]=[rba_session.get_process_information(process)['Composition'][protein]]
                                else:
                                    subunit_misprediction_factors_processes[process].append(misprediction_coeff)
                                    subunit_stoichiometries_processes[process].append(rba_session.get_process_information(process)['Composition'][protein])

    enzyme_correction_coefficients={} # {(Iso)Enzyme: correction_coeff} 
    # iterate over enzymes with determined SU misprediction factors
    for enzyme in subunit_misprediction_factors_enzymes.keys():
        #calculate complex misprediction factor as kapp-correction coeff
        enzyme_correction_coefficients[enzyme]=numpy.power(weighted_geometric_mean(data=subunit_misprediction_factors_enzymes[enzyme],weights=subunit_stoichiometries_enzymes[enzyme]),1/n_th_root_mispred)
        #enzyme_correction_coefficients[enzyme]=numpy.power(weighted_harmonic_mean(data=subunit_misprediction_factors_enzymes[enzyme],weights=subunit_stoichiometries_enzymes[enzyme]),1/n_th_root_mispred)
    
    enzymes_already_handled=[]
    for enzyme in enzyme_correction_coefficients.keys():
        if enzyme not in enzymes_already_handled:
            iso_enzyme_correction_coeffs={enzyme:enzyme_correction_coefficients[enzyme]}
            for iso_enzyme in rba_session.get_enzyme_information(enzyme)["Isozymes"]:
                if iso_enzyme in enzyme_correction_coefficients.keys():
                    iso_enzyme_correction_coeffs.update({iso_enzyme:enzyme_correction_coefficients[iso_enzyme]})
                else:
                    iso_enzyme_correction_coeffs.update({iso_enzyme:numpy.nan})
            #iso_enzyme_correction_coeffs={R_A_enzyme:corr_coeff_R_A_enzyme,R_A_duplicate_2_enzyme:corr_coeff_R_A_duplicate_2_enzyme,...}
            #common correction coeff for all isoenzymes as median over individual ones
            correction_coeff=numpy.median(numpy.array([i for i in list(iso_enzyme_correction_coeffs.values()) if numpy.isfinite(i)]))
            #correct each isoenzyme-kapp with common correction factor for proteoenzyme
            #if change larger than tolerance 
            for iso_enzyme_to_consider in iso_enzyme_correction_coeffs.keys():
                if iso_enzyme_to_consider in list(enzyme_efficiencies["Enzyme_ID"]):
                    enzymes_already_handled.append(iso_enzyme_to_consider)
                    old_kapp=enzyme_efficiencies.loc[enzyme_efficiencies["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"].values[0]
                    new_kapp=old_kapp*correction_coeff
                    if tolerance is None:
                        enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"]=new_kapp
                    else:
                        if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                            enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"]=new_kapp

    # Same for default kapp enzymes
    if correct_default_kapp_enzymes:
        default_efficiency_enzyme_correction_coefficients={} # {(Iso)Enzyme: correction_coeff} 
        # iterate over enzymes with determined SU misprediction factors
        for enzyme in subunit_misprediction_factors_default_efficiency_enzymes.keys():
            if enzyme not in list(enzyme_efficiencies["Enzyme_ID"]):
                #calculate complex misprediction factor as kapp-correction coeff
                default_efficiency_enzyme_correction_coefficients[enzyme]=numpy.power(weighted_geometric_mean(data=subunit_misprediction_factors_default_efficiency_enzymes[enzyme],weights=subunit_stoichiometries_default_efficiency_enzymes[enzyme]),1/n_th_root_mispred)
                #default_efficiency_enzyme_correction_coefficients[enzyme]=numpy.power(weighted_harmonic_mean(data=subunit_misprediction_factors_default_efficiency_enzymes[enzyme],weights=subunit_stoichiometries_default_efficiency_enzymes[enzyme]),1/n_th_root_mispred)

        for enzyme in default_efficiency_enzyme_correction_coefficients.keys():
            if enzyme not in enzymes_already_handled:
                respective_reaction=rba_session.get_enzyme_information(enzyme)["Reaction"].split("_duplicate_")[0]
                # find out wheter it is forward or backward specific kapp to add :
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
            #process_correction_coefficients[process]=numpy.power(weighted_harmonic_mean(data=subunit_misprediction_factors_processes[process],weights=subunit_stoichiometries_processes[process]),1/n_th_root_mispred)

    for process in process_correction_coefficients.keys():
        correction_coeff=process_correction_coefficients[process]
        old_efficiency=process_efficiencies.loc[process,"Value"]
        new_efficiency=old_efficiency*correction_coeff
        if tolerance is None:
            process_efficiencies_out.loc[process,"Value"]=new_efficiency
        else:
            if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                process_efficiencies_out.loc[process,"Value"]=new_efficiency

    return({"Kapps":enzyme_efficiencies_out,
            "ProcessEfficiencies":process_efficiencies_out,
            "Process_MispredictionFactors":process_correction_coefficients,
            "Enzyme_MispredictionFactors":enzyme_correction_coefficients})


def calculate_log_protein_prediction_mean_squared_error(protein_data,
                                                        simulation_results,
                                                        rba_session,
                                                        condition_to_look_up):
    predicted_proto_protein_quantities={}
    for i in simulation_results["Proteins"].index:
        proto_protein_ID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if proto_protein_ID in predicted_proto_protein_quantities:
            predicted_proto_protein_quantities[proto_protein_ID]+=simulation_results["Proteins"].loc[i,condition_to_look_up]
        else:
            predicted_proto_protein_quantities[proto_protein_ID]=simulation_results["Proteins"].loc[i,condition_to_look_up]
    
    proto_protein_isoform_map={} 
    for i in rba_session.get_proteins():
        proto_protein_ID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if proto_protein_ID in proto_protein_isoform_map:
            proto_protein_isoform_map[proto_protein_ID].append(i)
        else:
            proto_protein_isoform_map[proto_protein_ID]=[i]

    squared_residuals=[]
    for proto_protein_ID in predicted_proto_protein_quantities:
        if proto_protein_ID in proto_protein_isoform_map:
            if proto_protein_ID in protein_data["ID"]:
                predicted_protein=predicted_proto_protein_quantities[proto_protein_ID]
                measured_protein=protein_data.loc[proto_protein_ID,"copy_number"]
                if (predicted_protein>0) & (measured_protein>0) & (numpy.isfinite(predicted_protein))& (numpy.isfinite(measured_protein)):
                    squared_residuals.append((numpy.log(predicted_protein)-numpy.log(measured_protein))**2)
    return(sum(squared_residuals)/len(squared_residuals))
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

def weighted_harmonic_mean(data,weights=None):
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
        denominator=sum([1/i for i in data])
        enumerator=len(data)
    else:
        denominator=sum(weights[i]/data[i] for i in range(len(data)))
        enumerator=sum(weights)
    out=enumerator/denominator
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

    #rba_session.build_model_structure(file_name='/ModelStructure_meancompo_{}.json'.format(condition),print_warnings=False)
 
