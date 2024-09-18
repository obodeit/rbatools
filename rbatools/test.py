def correct_compartment_sizes_and_pg_fractions():
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

def correct_proteomics_data():
    ### define coeff as input ###
    proteome[condition]*=global_protein_scaling_coeff
    correction_results_compartement_sizes['original_amino_acid_occupation']*=global_protein_scaling_coeff

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
    proteome : _type_
        _description_
    condition : _type_
        _description_
    gene_ID_column : _type_
        _description_
    definition_file : _type_
        _description_
    rba_session : _type_
        _description_
    process_efficiency_estimation_input : _type_, optional
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
    min_kapp : _type_, optional
        _description_, by default None
    print_outputs : bool, optional
        _description_, by default True
    global_protein_scaling_coeff : int, optional
        _description_, by default 1
    use_mean_enzyme_composition_for_calibration : bool, optional
        _description_, by default False
    max_kapp_threshold : _type_, optional
        _description_, by default None
    output_dir : str, optional
        _description_, by default ""
    """
    correction_settings=machinery_efficiency_correction_settings_from_input(input=definition_file, condition=condition)
    enzyme_efficiency_estimation_settings=enzyme_efficiency_estimation_settings_from_input(input=definition_file, condition=condition)

    t0 = time.time()
    
    #########################################
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


def determine_apparent_process_efficiencies(growth_rate, input, rba_session, proteome_summary, protein_data, condition, gene_id_col,fit_nucleotide_assembly_machinery=False):
    """
    _summary_

    Parameters
    ----------
    growth_rate : _type_
        _description_
    input : _type_
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
        #constituting_proteins = {rba_session.get_protein_information(protein=i)['ProtoID']: rba_session.get_protein_information(protein=i)['AAnumber'] for i in rba_session.get_process_information(process=process_name)['Composition'].keys()}
        constituting_proteins={} 
        for i in rba_session.get_process_information(process=process_name)['Composition'].keys():
            proto_ID=rba_session.get_protein_information(protein=i)['ProtoID']
            if proto_ID in constituting_proteins.keys():
                constituting_proteins[proto_ID]+=rba_session.get_protein_information(protein=i)['AAnumber']*rba_session.get_process_information(process=process_name)['Composition'][i]
            else:
                constituting_proteins[proto_ID]=rba_session.get_protein_information(protein=i)['AAnumber']*rba_session.get_process_information(process=process_name)['Composition'][i]
            
        Total_client_fraction = sum([proteome_summary.loc[i, 'new_protein_fraction']
                                     for i in process_client_compartments])
        n_AAs_in_machinery = 0
        machinery_size = 0 ## 
        # iterate over all proteins, composing process machinery and sum-up 
        for i in constituting_proteins.keys():
            if i in protein_data['ID']:
                n_AAs_in_machinery += protein_data.loc[protein_data['ID'] == i, condition].values[0] * protein_data.loc[protein_data['ID'] == i, 'AA_residues'].values[0]
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


###################

def calibration_workflow_original(proteome,
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
    proteome : _type_
        _description_
    condition : _type_
        _description_
    gene_ID_column : _type_
        _description_
    definition_file : _type_
        _description_
    rba_session : _type_
        _description_
    process_efficiency_estimation_input : _type_, optional
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
    min_kapp : _type_, optional
        _description_, by default None
    print_outputs : bool, optional
        _description_, by default True
    global_protein_scaling_coeff : int, optional
        _description_, by default 1
    use_mean_enzyme_composition_for_calibration : bool, optional
        _description_, by default False
    max_kapp_threshold : _type_, optional
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

