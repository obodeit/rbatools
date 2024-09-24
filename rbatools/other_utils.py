import pandas
import numpy
import copy
from rbatools.rba_xml_utils import inject_estimated_efficiencies_into_model , inject_estimated_efficiencies_as_functions_into_model

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


def correct_compartment_fractions(proteome,condition,definition_file,
                                  compartments_to_replace,
                                  compartments_no_original_PG,
                                  fractions_entirely_replaced_with_expected_value,
                                  directly_corrected_compartments,
                                  imposed_compartment_fractions,
                                  merged_compartments,
                                  min_compartment_fraction):
    correction_results_compartement_sizes = correction_pipeline(input=proteome,
                                                condition=condition,
                                                definition_file=definition_file,
                                                compartments_to_replace=compartments_to_replace,
                                                compartments_no_original_PG=compartments_no_original_PG,
                                                fractions_entirely_replaced_with_expected_value=fractions_entirely_replaced_with_expected_value,
                                                imposed_compartment_fractions=imposed_compartment_fractions,
                                                directly_corrected_compartments=directly_corrected_compartments,
                                                merged_compartments=merged_compartments,
                                                min_compartment_fraction=min_compartment_fraction)

    #correction_results_compartement_sizes.to_csv(str(output_dir+'/Correction_overview_HackettNielsen_'+condition+'.csv'))
    return({"Densities_PGs":extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=correction_results_compartement_sizes,rows_to_exclude=["Ribosomes","Total"]+[i for i in correction_results_compartement_sizes.index if i.startswith("pg_")]),
            "Proteome":proteome,
            "Condition":condition,
            "Proteome_summary":correction_results_compartement_sizes})


def correct_proteome(correction_results_compartement_sizes,proteome,condition,Compartment_sizes,PG_fractions):
    for i in Compartment_sizes.index:
        correction_results_compartement_sizes.loc[i,"new_protein_fraction"]=Compartment_sizes.loc[i,condition]
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
    #correction_results_compartement_sizes.to_csv(str(output_dir+'/Correction_overview_HackettNielsen_corrected_'+condition+'.csv'))
    return({"Proteome":proteome,
            "Condition":condition,
            "Proteome_summary":correction_results_compartement_sizes})


def machinery_efficiency_correction_settings_from_input(input, condition):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    condition : _type_
        _description_
    """
    out={}
    out['tolerance_global_scaling']=0.1
    if 'Global_scaling_tolerance_growth_rate_approximation' in list(input['Type']):
        value=input.loc[input['Type'] == 'Global_scaling_tolerance_growth_rate_approximation', condition].values[0]
        if numpy.isfinite(value):
            out['tolerance_global_scaling']=value

    out['fixed_growth_rate_global_scaling']=False
    if 'Solve_at_measured_growth_rate_after_global_scaling' in list(input['Type']):
        value=input.loc[input['Type'] == 'Solve_at_measured_growth_rate_after_global_scaling', condition].values[0]
        if numpy.isfinite(value):
            out['fixed_growth_rate_global_scaling']=bool(value)
    
    out['abjust_root_of_correction_coeffs_global_scaling']=False
    if 'Global_scaling_adjust_roots_of_scaling_coefficient' in list(input['Type']):
        value=input.loc[input['Type'] == 'Global_scaling_adjust_roots_of_scaling_coefficient', condition].values[0]
        if numpy.isfinite(value):
            out['abjust_root_of_correction_coeffs_global_scaling']=bool(value)

    out['minimum_iteration_number']=2
    if 'Minimum_iteration_number_efficiency_correction' in list(input['Type']):
        value=input.loc[input['Type'] == 'Minimum_iteration_number_efficiency_correction', condition].values[0]
        if numpy.isfinite(value):
            out['minimum_iteration_number']=value

    out['iteration_limit']=10
    if 'Maximum_iteration_number_efficiency_correction' in list(input['Type']):
        value=input.loc[input['Type'] == 'Maximum_iteration_number_efficiency_correction', condition].values[0]
        if numpy.isfinite(value):
            out['iteration_limit']=value

    out['steady_rss_limit']=5
    if 'Limit_for_consecutive_iterations_without_change_in_RSS' in list(input['Type']):
        value=input.loc[input['Type'] == 'Limit_for_consecutive_iterations_without_change_in_RSS', condition].values[0]
        if numpy.isfinite(value):
            out['steady_rss_limit']=value

    out['increasing_rss_limit']=3
    if 'Limit_for_consecutive_iterations_with_increasing_change_in_RSS' in list(input['Type']):
        value=input.loc[input['Type'] == 'Limit_for_consecutive_iterations_with_increasing_change_in_RSS', condition].values[0]
        if numpy.isfinite(value):
            out['increasing_rss_limit']=value

    out['rss_tolerance']=0.1
    if 'Range_for_steady_RSS' in list(input['Type']):
        value=input.loc[input['Type'] == 'Range_for_steady_RSS', condition].values[0]
        if numpy.isfinite(value):
            out['rss_tolerance']=value

    out['increasing_rss_factor']=1.1
    if 'Tolerance_for_increasing_RSS' in list(input['Type']):
        value=input.loc[input['Type'] == 'Tolerance_for_increasing_RSS', condition].values[0]
        if numpy.isfinite(value):
            out['increasing_rss_factor']=value

    out['correct_efficiencies']=True
    if 'Correct_machinery_efficiencies' in list(input['Type']):
        value=input.loc[input['Type'] == 'Correct_machinery_efficiencies', condition].values[0]
        if numpy.isfinite(value):
            out['correct_efficiencies']=bool(value)
    ####
    out['impose_directions_from_fba_during_correction']=False
    if 'Impose_directions_from_enzyme_efficiency_estimation_flux_distribution' in list(input['Type']):
        value=input.loc[input['Type'] == 'Impose_directions_from_enzyme_efficiency_estimation_flux_distribution', condition].values[0]
        if numpy.isfinite(value):
            out['impose_directions_from_fba_during_correction']=bool(value)

    out['final_global_scaling_after_correction']=False
    if 'Global_scaling_after_machinery_efficiency_correction' in list(input['Type']):
        value=input.loc[input['Type'] == 'Global_scaling_after_machinery_efficiency_correction', condition].values[0]
        if numpy.isfinite(value):
            out['final_global_scaling_after_correction']=bool(value)

    out['final_global_scaling_without_imposed_exchanges']=False
    if 'Final_global_scaling_without_imposed_exchanges' in list(input['Type']):
        value=input.loc[input['Type'] == 'Final_global_scaling_without_imposed_exchanges', condition].values[0]
        if numpy.isfinite(value):
            out['final_global_scaling_without_imposed_exchanges']=bool(value)
    return(out)


def enzyme_efficiency_estimation_settings_from_input(input, condition):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    condition : _type_
        _description_
    """
    out={}
    out['rRNA_target_rna_weight_ratio']=4.8824
    if 'rRNA_target_rna_weight_ratio' in list(input['Type']):
        value=input.loc[input['Type'] == 'rRNA_target_rna_weight_ratio', condition].values[0]
        if numpy.isfinite(value):
            out['rRNA_target_rna_weight_ratio']=float(value)

    out['use_target_biomass_function']=True
    if 'RBA_derived_biomass_function' in list(input['Type']):
        value=input.loc[input['Type'] == 'RBA_derived_biomass_function', condition].values[0]
        if numpy.isfinite(value):
            out['use_target_biomass_function']=bool(value)

    out['biomass_function_in_model']=None
    if 'Biomass_function_in_model' in list(input['Type']):
        value=input.loc[input['Type'] == 'Biomass_function_in_model', condition].values[0]
        if type(value) is str:
            out['biomass_function_in_model']=bool(value)

    out['parsimonious_fba']=True
    if 'Use_parsimonious_FBA' in list(input['Type']):
        value=input.loc[input['Type'] == 'Use_parsimonious_FBA', condition].values[0]
        if numpy.isfinite(value):
            out['parsimonious_fba']=bool(value)

    out['chose_most_likely_isoreactions']=True
    if 'Eliminate_unlikely_isoenzymes' in list(input['Type']):
        value=input.loc[input['Type'] == 'Eliminate_unlikely_isoenzymes', condition].values[0]
        if numpy.isfinite(value):
            out['chose_most_likely_isoreactions']=bool(value)

    out['impose_on_all_isoreactions']=False
    if 'Impose_efficiency_onto_all_other_isoenzymes' in list(input['Type']):
        value=input.loc[input['Type'] == 'Impose_efficiency_onto_all_other_isoenzymes', condition].values[0]
        if numpy.isfinite(value):
            out['impose_on_all_isoreactions']=bool(value)

    out['zero_on_all_isoreactions']=True
    if 'Impose_zero_onto_all_other_isoenzymes' in list(input['Type']):
        value=input.loc[input['Type'] == 'Impose_zero_onto_all_other_isoenzymes', condition].values[0]
        if numpy.isfinite(value):
            out['zero_on_all_isoreactions']=bool(value)

    out['node_degree_identical_enzyme_network']=1
    if 'Consider_other_catalytic_activities_of_isoenzymes' in list(input['Type']):
        value=input.loc[input['Type'] == 'Consider_other_catalytic_activities_of_isoenzymes', condition].values[0]
        if numpy.isfinite(value):
            out['node_degree_identical_enzyme_network']=value

    out['impose_on_identical_enzymes']=True
    if 'Use_efficiencies_for_enzymes_with_identical_composition' in list(input['Type']):
        value=input.loc[input['Type'] == 'Use_efficiencies_for_enzymes_with_identical_composition', condition].values[0]
        if numpy.isfinite(value):
            out['impose_on_identical_enzymes']=bool(value)

    return(out)


def flux_bounds_from_input(input,rba_session, condition, specific_exchanges=None, specific_directions=None,also_consider_iso_enzmes=True):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    rba_session : _type_
        _description_
    condition : _type_
        _description_
    specific_exchanges : _type_, optional
        _description_, by default None
    specific_directions : _type_, optional
        _description_, by default None
    also_consider_iso_enzmes : bool, optional
        _description_, by default True
    """
    out = pandas.DataFrame(columns=['Reaction_ID', 'LB', 'UB'])

    flux_mean_df = input.loc[input['Type'] == 'Flux_Range_Mean', :]
    flux_mean_SE = input.loc[input['Type'] == 'Flux_Range_StandardError', :]
    if specific_exchanges is None:
        exchanges_to_set = list(flux_mean_df['ID'])
    else:
        exchanges_to_set = specific_exchanges
    for rx in exchanges_to_set:
        reactions_to_consider=[rx]
        for reaction_to_consider in reactions_to_consider:
            mean_val = flux_mean_df.loc[flux_mean_df['ID'] == rx, condition].values[0]
            if not pandas.isna(mean_val):
                out.loc[reaction_to_consider, 'Reaction_ID'] = reaction_to_consider
                out.loc[reaction_to_consider, 'Mean'] = mean_val
                if str(rx+'_SE') in flux_mean_SE['ID']:
                    SE_val = flux_mean_SE.loc[flux_mean_SE['ID'] == str(rx+'_SE'), condition].values[0]
                else:
                    SE_val=numpy.nan
                if not pandas.isna(SE_val):
                    lb = mean_val-SE_val
                    ub = mean_val+SE_val
                    if mean_val < 0:
                        out.loc[reaction_to_consider, 'LB'] = lb
                        if ub > 0:
                            out.loc[reaction_to_consider, 'UB'] = 0
                        else:
                            out.loc[reaction_to_consider, 'UB'] = ub
                    elif mean_val > 0:
                        out.loc[reaction_to_consider, 'UB'] = ub
                        if lb < 0:
                            out.loc[reaction_to_consider, 'LB'] = 0
                        else:
                            out.loc[reaction_to_consider, 'LB'] = lb
                    else:
                        out.loc[reaction_to_consider, 'LB'] = lb
                        out.loc[reaction_to_consider, 'UB'] = ub
                else:
                    out.loc[reaction_to_consider, 'LB'] = mean_val
                    out.loc[reaction_to_consider, 'UB'] = mean_val

    flux_dir_df = input.loc[input['Type'] == 'Flux_Direction', :]
    if specific_directions is None:
        directions_to_set = list(flux_dir_df['ID'])
    else:
        directions_to_set = specific_directions
    for rx in directions_to_set:
        reactions_to_consider=[rx]
        if also_consider_iso_enzmes:
            if rx in rba_session.get_reactions():
                reactions_to_consider+=list(rba_session.get_reaction_information(rx)['Twins'])
        for reaction_to_consider in reactions_to_consider:
            out.loc[reaction_to_consider, 'Reaction_ID'] = reaction_to_consider
            if flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == 1:
                out.loc[reaction_to_consider, 'LB'] = 0
            elif flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == -1:
                out.loc[reaction_to_consider, 'UB'] = 0
            elif flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == 0:
                out.loc[reaction_to_consider, 'LB'] = 0
                out.loc[reaction_to_consider, 'UB'] = 0

    flux_upper_df = input.loc[input['Type'] == 'Flux_Upper_Bound', :]
    for rx in list(flux_upper_df['ID']):
        reactions_to_consider=[rx]
        if also_consider_iso_enzmes:
            if rx in rba_session.get_reactions():
                reactions_to_consider+=list(rba_session.get_reaction_information(rx)['Twins'])
        for reaction_to_consider in reactions_to_consider:
            out.loc[reaction_to_consider, 'Reaction_ID'] = reaction_to_consider
            out.loc[reaction_to_consider, 'UB'] = flux_upper_df.loc[flux_upper_df['ID'] == rx, condition].values[0]

    flux_lower_df = input.loc[input['Type'] == 'Flux_Lower_Bound', :]
    for rx in list(flux_lower_df['ID']):
        reactions_to_consider=[rx]
        if also_consider_iso_enzmes:
            if rx in rba_session.get_reactions():
                reactions_to_consider+=list(rba_session.get_reaction_information(rx)['Twins'])
        for reaction_to_consider in reactions_to_consider:
            out.loc[rx, 'Reaction_ID'] = reaction_to_consider
            out.loc[rx, 'LB'] = flux_lower_df.loc[flux_lower_df['ID'] == rx, condition].values[0]
    return(out)


def growth_rate_from_input(input, condition):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    condition : _type_
        _description_
    """
    return(input.loc[input['Type'] == 'Growth_Rate', condition].values[0])


def proteome_fractions_from_input(input, condition):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    condition : _type_
        _description_
    """
    df = input.loc[input['Type'] == 'Expected_ProteomeFraction', :]
    return(dict(zip(list(df['ID']), list(df[condition]))))


def medium_concentrations_from_input(input, condition):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    condition : _type_
        _description_
    """
    df = input.loc[input['Type'] == 'Medium_Concentration', :]
    return(dict(zip(list(df['ID']), list(df[condition]))))


def generate_compartment_size_and_pg_input(compartment_sizes,pg_fractions,condition):
    """
    _summary_

    Parameters
    ----------
    compartment_sizes : _type_
        _description_
    pg_fractions : _type_
        _description_
    condition : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in compartment_sizes.index:
        out.loc[i,"Compartment_ID"]=i
        out.loc[i,"Density"]=round(compartment_sizes.loc[i,condition],3)
    for i in pg_fractions.index:
        out.loc[i,"Compartment_ID"]=i
        out.loc[i,"PG_fraction"]=round(pg_fractions.loc[i,condition],3)
    return(out)


def generate_process_efficiency_input(process_efficiencies,condition,parameter_name_suffix):
    """
    _summary_

    Parameters
    ----------
    process_efficiencies : _type_
        _description_
    condition : _type_
        _description_
    parameter_name_suffix : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in process_efficiencies.index:
        out.loc[i,"Process"]=i
        out.loc[i,"Parameter"]=i+parameter_name_suffix
        out.loc[i,"Value"]=process_efficiencies.loc[i,condition]
    return(out)


def generate_default_kapp_input(default_kapps,condition,transporter_multiplier):
    """
    _summary_

    Parameters
    ----------
    default_kapps : _type_
        _description_
    condition : _type_
        _description_
    transporter_multiplier : _type_
        _description_
    """
    out={}
    out["default_efficiency"]=default_kapps.loc[condition,"Default Kapp"]
    out["default_transporter_efficiency"]=transporter_multiplier*default_kapps.loc[condition,"Default Kapp"]
    return(out)


def generate_specific_kapp_input(specific_kapps,condition):
    """
    _summary_

    Parameters
    ----------
    specific_kapps : _type_
        _description_
    condition : _type_
        _description_
    """
    out=pandas.DataFrame()
    for i in specific_kapps.index:
        if not pandas.isna(specific_kapps.loc[i,condition]):
            out.loc[i,"Enzyme_ID"]=i
            out.loc[i,"Kapp"]=specific_kapps.loc[i,condition]
            out.loc[i,"Flux"]=specific_kapps.loc[i,condition+"_Flux"]
    return(out)


def perform_simulations(condition,
                        rba_session,
                        definition_file,
                        compartment_sizes,
                        pg_fractions,
                        process_efficiencies,
                        Default_Kapps,
                        Specific_Kapps,
                        Exchanges_to_impose=None,
                        sims_to_perform=['DefaultKapp','Prokaryotic','Eukaryotic'],
                        feasible_stati=['optimal','feasible'],
                        try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                        print_output=True,
                        variability_analysis=None,
                        mu_factor_for_variability=1,
                        apply_model=False,
                        functions_to_include_list=[],
                        transporter_multiplier=3,
                        start_val=numpy.nan,
                        Mu_approx_precision=0.000001,
                        max_mu_in_dichotomy=1,
                        proteomics_constraints_input={}):
    """
    _summary_

    Parameters
    ----------
    condition : _type_
        _description_
    rba_session : _type_
        _description_
    definition_file : _type_
        _description_
    compartment_sizes : _type_
        _description_
    pg_fractions : _type_
        _description_
    process_efficiencies : _type_
        _description_
    Default_Kapps : _type_
        _description_
    Specific_Kapps : _type_
        _description_
    Exchanges_to_impose : _type_, optional
        _description_, by default None
    sims_to_perform : list, optional
        _description_, by default ['DefaultKapp','Prokaryotic','Eukaryotic']
    feasible_stati : list, optional
        _description_, by default ['optimal','feasible']
    try_unscaling_if_sol_status_is_feasible_only_before_unscaling : bool, optional
        _description_, by default True
    print_output : bool, optional
        _description_, by default True
    variability_analysis : _type_, optional
        _description_, by default None
    mu_factor_for_variability : int, optional
        _description_, by default 1
    apply_model : bool, optional
        _description_, by default False
    functions_to_include_list : list, optional
        _description_, by default []
    transporter_multiplier : int, optional
        _description_, by default 3
    start_val : _type_, optional
        _description_, by default numpy.nan
    Mu_approx_precision : float, optional
        _description_, by default 0.000001
    max_mu_in_dichotomy : int, optional
        _description_, by default 1
    proteomics_constraints_input : dict, optional
        _description_, by default {}
    """
    out={'SolutionStatus_def':None,
        'SolutionStatus_prok':None,
        'SolutionStatus_euk':None,
        'SolutionStatus_euk_fixed':None,
        'SolutionStatus_fixed_pg_euk':None,
        'SolutionStatus_fixed_pg_euk_fixed':None,
        'FeasibleRange_def':{},
        'FeasibleRange_prok':{},
        'FeasibleRange_euk':{},
        'FeasibleRange_euk_fixed':{},
        'FeasibleRange_fixed_pg_euk':{},
        'FeasibleRange_fixed_pg_euk_fixed':{},
        'Mu_def':numpy.nan,
        'Mu_prok':numpy.nan,
        'Mu_euk':numpy.nan,
        'Mu_euk_fixed':numpy.nan,
        'Mu_fixed_pg_euk':numpy.nan,
        'Mu_euk_fixed_pg_fixed':numpy.nan,
        'Simulation_Results': {}, 
        'Simulation_Results_Euk': {}, 
        'Simulation_Results_Euk_fixed': {}, 
        'Simulation_Results_DefKapp': {}, 
        'Simulation_Results_fixed_pg_Euk': {}, 
        'Simulation_Results_fixed_pg_Euk_fixed': {}, 
        'Euk_CompSizes': {},
        'Euk_fixed_CompSizes': {},
        'Fixed_pg_Euk_CompSizes': {},
        'Fixed_pg_Euk_fixed_CompSizes': {},
        'Condition':condition,
        'Performed_simulations':sims_to_perform}

    if 'DefaultKapp' in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if 'Default_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if 'Compartment_Sizes' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if 'PG_Fractions' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if 'Process_Efficiencies' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]['LB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['LB'])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]['UB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['UB'])})
        out['Mu_def']=rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy,start_value=start_val, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        out['SolutionStatus_def']=rba_session.Problem.SolutionStatus
        try:
            rba_session.clear_results_and_parameters()
            rba_session.record_results('DefaultKapp')
            if print_output:
                print('Mu Def: {}'.format(out['Mu_def']))
            out['Simulation_Results_DefKapp'] = copy.deepcopy(rba_session.Results)
        except:
            if print_output:
                print('Mu Def: failed')
        if variability_analysis is not None:
            if len(list(out['Simulation_Results_DefKapp'].keys()))!=0:
                rba_session.set_growth_rate(out['Mu_def']*mu_factor_for_variability)
                out['FeasibleRange_def']=rba_session.get_feasible_range(variability_analysis)

    if 'Prokaryotic' in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if 'Specific_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if 'Default_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if 'Compartment_Sizes' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if 'PG_Fractions' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if 'Process_Efficiencies' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]['LB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['LB'])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]['UB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['UB'])})
        if len(list(proteomics_constraints_input.keys()))!=0:
            rba_session.build_protein_usage_constraints(input_proteome=proteomics_constraints_input)
        out['Mu_prok'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy,start_value=start_val, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        out['SolutionStatus_prok']=rba_session.Problem.SolutionStatus
                
        try:
            rba_session.clear_results_and_parameters()
            rba_session.record_results('Prokaryotic')
            if print_output:
                print('Mu Prok: {}'.format(out['Mu_prok']))
            out['Simulation_Results'] = copy.deepcopy(rba_session.Results)
        except:
            if print_output:
                print('Mu Prok: failed')
                print('Status: {} - Mu: {}'.format(out['SolutionStatus_prok'],out['Mu_prok']))
        if variability_analysis is not None:
            if len(list(out['Simulation_Results'].keys()))!=0:
                rba_session.set_growth_rate(out['Mu_prok']*mu_factor_for_variability)
                out['FeasibleRange_prok']=rba_session.get_feasible_range(variability_analysis)

    if 'Eukaryotic' in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if 'Specific_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if 'Default_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if 'Compartment_Sizes' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if 'PG_Fractions' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if 'Process_Efficiencies' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        #rba_session.eukaryotic_densities_calibration(CompartmentRelationships=False)
        #rba_session.eukaryotic_densities_pg_fraction(fixed_size_compartments=[],compartment_fraction_prefix='fraction_protein_')
        pg_fractions_for_euk={comp: str('fraction_non_enzymatic_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        compartment_fractions_for_euk={comp:str('fraction_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        rba_session.make_eukaryotic(amino_acid_concentration_total='amino_acid_concentration',
                                pg_fractions=pg_fractions_for_euk,
                                compartment_fractions=compartment_fractions_for_euk,
                                    compartments_with_imposed_sizes=[],
                                normalise_global_fraction=True,
                                compartment_bound_tolerance=0.0)

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]['LB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['LB'])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]['UB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['UB'])})

        out['Mu_euk'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        #out['Mu_euk'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy,start_value=max_mu_in_dichotomy/2, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        out['SolutionStatus_euk']=rba_session.Problem.SolutionStatus
        try:
            rba_session.clear_results_and_parameters()
            rba_session.record_results('Eukaryotic')
            if print_output:
                print('Mu Euk: {}'.format(out['Mu_euk']))
            out['Simulation_Results_Euk'] = copy.deepcopy(rba_session.Results)
            for comp in list(compartment_fractions_for_euk.keys()):
                out['Euk_CompSizes'][comp] = rba_session.Problem.SolutionValues[str('f_'+comp)]
        except:
            if print_output:
                print('Mu Euk: failed')
                print('Status: {} - Mu: {}'.format(out['SolutionStatus_euk'],out['Mu_euk']))
        if variability_analysis is not None:
            if len(list(out['Simulation_Results_Euk'].keys()))!=0:
                rba_session.set_growth_rate(out['Mu_euk']*mu_factor_for_variability)
                out['FeasibleRange_euk']=rba_session.get_feasible_range(variability_analysis)

    if 'Eukaryotic_fixed_sizes' in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if 'Specific_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if 'Default_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if 'Compartment_Sizes' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if 'PG_Fractions' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if 'Process_Efficiencies' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        #rba_session.eukaryotic_densities_calibration(CompartmentRelationships=False)
        #rba_session.eukaryotic_densities_pg_fraction(fixed_size_compartments=[],compartment_fraction_prefix='fraction_protein_')
        pg_fractions_for_euk={comp: str('fraction_non_enzymatic_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        compartment_fractions_for_euk={comp:str('fraction_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        rba_session.make_eukaryotic(amino_acid_concentration_total='amino_acid_concentration',
                                pg_fractions=pg_fractions_for_euk,
                                compartment_fractions=compartment_fractions_for_euk,
                                compartments_with_imposed_sizes=list(compartment_fractions_for_euk.keys()),
                                normalise_global_fraction=True,
                                compartment_bound_tolerance=0.1,
                                imposed_compartments_without_tolerance=['Secreted','n'])

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]['LB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['LB'])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]['UB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['UB'])})

        out['Mu_euk_fixed'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        #out['Mu_euk_fixed'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy,start_value=max_mu_in_dichotomy/2, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        out['SolutionStatus_euk_fixed']=rba_session.Problem.SolutionStatus
        try:
            rba_session.clear_results_and_parameters()
            rba_session.record_results('Eukaryotic_fixed_sizes')
            if print_output:
                print('Mu Euk fixed: {}'.format(out['Mu_euk_fixed']))
            out['Simulation_Results_Euk_fixed'] = copy.deepcopy(rba_session.Results)
            for comp in list(compartment_fractions_for_euk.keys()):
                out['Euk_fixed_CompSizes'][comp] = rba_session.Problem.SolutionValues[str('f_'+comp)]
        except:
            if print_output:
                print('Mu Euk fixed: failed')
                print('Status: {} - Mu: {}'.format(out['SolutionStatus_euk_fixed'],out['Mu_euk_fixed']))
        if variability_analysis is not None:
            if len(list(out['Simulation_Results_Euk_fixed'].keys()))!=0:
                rba_session.set_growth_rate(out['Mu_euk_fixed']*mu_factor_for_variability)
                out['FeasibleRange_euk_fixed']=rba_session.get_feasible_range(variability_analysis)

    if 'Fixed_PG_Eukaryotic' in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if 'Specific_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if 'Default_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if 'Compartment_Sizes' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if 'PG_Fractions' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if 'Process_Efficiencies' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        #rba_session.eukaryotic_densities_calibration(CompartmentRelationships=False)
        #rba_session.eukaryotic_densities_pg_fraction(fixed_size_compartments=[],compartment_fraction_prefix='fraction_protein_')
        pg_fractions_for_euk={comp: str('fraction_non_enzymatic_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        compartment_fractions_for_euk={comp:str('fraction_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        rba_session.make_eukaryotic_fixed_pg_content(amino_acid_concentration_total='amino_acid_concentration',
                                                    external_compartment_fractions=[],                                         
                                                    pg_fractions=pg_fractions_for_euk,
                                                    compartment_fractions=compartment_fractions_for_euk,
                                                    compartments_with_imposed_sizes=[])

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]['LB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['LB'])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]['UB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['UB'])})

        out['Mu_fixed_pg_euk'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        #out['Mu_fixed_pg_euk'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy,start_value=max_mu_in_dichotomy/2, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        out['SolutionStatus_fixed_pg_euk']=rba_session.Problem.SolutionStatus
        try:
            rba_session.clear_results_and_parameters()
            rba_session.record_results('Fixed_PG_Eukaryotic')
            if print_output:
                print('Mu fixed PG Euk: {}'.format(out['Mu_fixed_pg_euk']))
            out['Simulation_Results_fixed_pg_Euk'] = copy.deepcopy(rba_session.Results)
            for comp in list(compartment_fractions_for_euk.keys()):
                out['Fixed_pg_Euk_CompSizes'][comp] = rba_session.Problem.SolutionValues[str('f_'+comp)]
        except:
            if print_output:
                print('Mu fixed PG Euk: failed')
                print('Status: {} - Mu: {}'.format(out['SolutionStatus_fixed_pg_euk'],out['Mu_fixed_pg_euk']))
        if variability_analysis is not None:
            if len(list(out['Simulation_Results_fixed_pg_Euk'].keys()))!=0:
                rba_session.set_growth_rate(out['Mu_fixed_pg_euk']*mu_factor_for_variability)
                out['FeasibleRange_fixed_pg_euk']=rba_session.get_feasible_range(variability_analysis)

    if 'Fixed_PG_Eukaryotic_fixed_sizes' in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if 'Specific_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if 'Default_Kapps' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if 'Compartment_Sizes' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if 'PG_Fractions' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if 'Process_Efficiencies' in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix='_apparent_efficiency')
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        #rba_session.eukaryotic_densities_calibration(CompartmentRelationships=False)
        #rba_session.eukaryotic_densities_pg_fraction(fixed_size_compartments=[],compartment_fraction_prefix='fraction_protein_')
        pg_fractions_for_euk={comp: str('fraction_non_enzymatic_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        compartment_fractions_for_euk={comp:str('fraction_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        rba_session.make_eukaryotic_fixed_pg_content(amino_acid_concentration_total='amino_acid_concentration',
                                                    external_compartment_fractions=[],                                         
                                                    pg_fractions=pg_fractions_for_euk,
                                                    compartment_fractions=compartment_fractions_for_euk,
                                                    compartments_with_imposed_sizes=list(compartment_fractions_for_euk.keys()))

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]['LB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['LB'])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]['UB'] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]['UB'])})

        out['Mu_fixed_pg_euk_fixed'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        #out['Mu_fixed_pg_euk_fixed'] = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max_value=max_mu_in_dichotomy,start_value=max_mu_in_dichotomy/2, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling,verbose=False)
        out['SolutionStatus_fixed_pg_euk_fixed']=rba_session.Problem.SolutionStatus
        try:
            rba_session.clear_results_and_parameters()
            rba_session.record_results('Fixed_PG_Eukaryotic_fixed_sizes')
            if print_output:
                print('Mu fixed PG Euk fixed: {}'.format(out['Mu_fixed_pg_euk_fixed']))
            out['Simulation_Results_fixed_pg_Euk_fixed'] = copy.deepcopy(rba_session.Results)
            for comp in list(compartment_fractions_for_euk.keys()):
                out['Fixed_pg_Euk_fixed_CompSizes'][comp] = rba_session.Problem.SolutionValues[str('f_'+comp)]
        except:
            if print_output:
                print('Mu fixed PG Euk fixed: failed')
                print('Status: {} - Mu: {}'.format(out['SolutionStatus_fixed_pg_euk_fixed'],out['Mu_fixed_pg_euk_fixed']))
        if variability_analysis is not None:
            if len(list(out['Simulation_Results_fixed_pg_Euk_fixed'].keys()))!=0:
                rba_session.set_growth_rate(out['Mu_fixed_pg_euk_fixed']*mu_factor_for_variability)
                out['FeasibleRange_fixed_pg_euk_fixed']=rba_session.get_feasible_range(variability_analysis)

    rba_session.clear_results_and_parameters()
    return(out)


def perform_simulations_fixed_Mu(condition,
                                 rba_session,
                                 definition_file,
                                 compartment_sizes,
                                 pg_fractions,
                                 process_efficiencies,
                                 Default_Kapps,
                                 Specific_Kapps,
                                 Exchanges_to_impose=None,
                                 metabolite_to_minimize=None,
                                 sims_to_perform=["DefaultKapp","Prokaryotic","Eukaryotic"],
                                 feasible_stati=["optimal","feasible"],
                                 try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
                                 Mu_scaling_coeff=1,
                                 total_protein_scaling_coeff=1,
                                 print_output=False,
                                 variability_analysis=None,
                                 mu_factor_for_variability=1,
                                 apply_model=False,
                                 functions_to_include_list=[],
                                 transporter_multiplier=3,
                                 proteomics_constraints_input={}):
    """
    _summary_

    Parameters
    ----------
    condition : _type_
        _description_
    rba_session : _type_
        _description_
    definition_file : _type_
        _description_
    compartment_sizes : _type_
        _description_
    pg_fractions : _type_
        _description_
    process_efficiencies : _type_
        _description_
    Default_Kapps : _type_
        _description_
    Specific_Kapps : _type_
        _description_
    Exchanges_to_impose : _type_, optional
        _description_, by default None
    metabolite_to_minimize : _type_, optional
        _description_, by default None
    sims_to_perform : list, optional
        _description_, by default ["DefaultKapp","Prokaryotic","Eukaryotic"]
    feasible_stati : list, optional
        _description_, by default ["optimal","feasible"]
    try_unscaling_if_sol_status_is_feasible_only_before_unscaling : bool, optional
        _description_, by default True
    Mu_scaling_coeff : int, optional
        _description_, by default 1
    total_protein_scaling_coeff : int, optional
        _description_, by default 1
    print_output : bool, optional
        _description_, by default False
    variability_analysis : _type_, optional
        _description_, by default None
    mu_factor_for_variability : int, optional
        _description_, by default 1
    apply_model : bool, optional
        _description_, by default False
    functions_to_include_list : list, optional
        _description_, by default []
    transporter_multiplier : int, optional
        _description_, by default 3
    proteomics_constraints_input : dict, optional
        _description_, by default {}
    """

    def_Feasible_Ranges={}
    prok_Feasible_Ranges={}
    euk_Feasible_Ranges={}
    euk_fixed_Feasible_Ranges={}
    def_results={}
    prok_results={}
    euk_results={}
    euk_fixed_results={}
    compartment_fractions_euk={}
    compartment_fractions_euk_fixed={}
    if "DefaultKapp" in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if "Default_Kapps" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if "Compartment_Sizes" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if "PG_Fractions" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if "Process_Efficiencies" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]["LB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["LB"])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]["UB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["UB"])})

        rba_session.set_growth_rate(Mu=Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
        if total_protein_scaling_coeff != 1:
            rba_session.Problem.set_right_hand_side({"{}_density".format(comp):rba_session.Problem.get_right_hand_side("{}_density".format(comp))["{}_density".format(comp)]*total_protein_scaling_coeff for comp in list(compartment_densities_and_PGs['Compartment_ID'])})
        if metabolite_to_minimize is None:
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        else:
            met_conc=rba_session.find_min_substrate_concentration(metabolite=metabolite_to_minimize, precision=0.001, max=1000, recording=False)
            if print_output:
                print("Min {} : {}".format(metabolite_to_minimize,met_conc))
            rba_session.set_medium({metabolite_to_minimize:met_conc})
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        if rba_session.Problem.Solved:
            rba_session.record_results('DefaultKapp')
            if print_output:
                print('Mu Def: {}'.format(growth_rate_from_input(input=definition_file, condition=condition)))
            def_results = copy.deepcopy(rba_session.Results)
            rba_session.clear_results_and_parameters()
            if variability_analysis is not None:
                if len(list(def_results.keys())):
                    rba_session.set_growth_rate(Mu=mu_factor_for_variability*Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
                    def_Feasible_Ranges=rba_session.get_feasible_range(variability_analysis)

    if "Prokaryotic" in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if "Specific_Kapps" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if "Default_Kapps" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if "Compartment_Sizes" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if "PG_Fractions" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if "Process_Efficiencies" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        if len(list(proteomics_constraints_input.keys()))!=0:
            rba_session.build_protein_usage_constraints(input_proteome=proteomics_constraints_input)

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]["LB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["LB"])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]["UB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["UB"])})

        rba_session.set_growth_rate(Mu=Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
        if total_protein_scaling_coeff != 1:
            rba_session.Problem.set_right_hand_side({"{}_density".format(comp):rba_session.Problem.get_right_hand_side("{}_density".format(comp))["{}_density".format(comp)]*total_protein_scaling_coeff for comp in list(compartment_densities_and_PGs['Compartment_ID'])})
        if metabolite_to_minimize is None:
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        else:
            met_conc=rba_session.find_min_substrate_concentration(metabolite=metabolite_to_minimize, precision=0.001, max=1000, recording=False)
            if print_output:
                print("Min {} : {}".format(metabolite_to_minimize,met_conc))
            rba_session.set_medium({metabolite_to_minimize:met_conc})
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        if rba_session.Problem.Solved:
            rba_session.record_results('Prokaryotic')
            if print_output:
                print('Mu Prok: {}'.format(growth_rate_from_input(input=definition_file, condition=condition)))
            prok_results = copy.deepcopy(rba_session.Results)
            rba_session.clear_results_and_parameters()
            if variability_analysis is not None:
                if len(list(prok_results.keys())):
                    rba_session.set_growth_rate(Mu=mu_factor_for_variability*Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
                    prok_Feasible_Ranges=rba_session.get_feasible_range(variability_analysis)

    if "Eukaryotic" in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if "Specific_Kapps" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if "Default_Kapps" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if "Compartment_Sizes" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if "PG_Fractions" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if "Process_Efficiencies" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        pg_fractions_for_euk={comp: str('fraction_non_enzymatic_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        compartment_fractions_for_euk={comp:str('fraction_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        rba_session.make_eukaryotic(amino_acid_concentration_total='amino_acid_concentration',
                                pg_fractions=pg_fractions_for_euk,
                                compartment_fractions=compartment_fractions_for_euk,
                                compartments_with_imposed_sizes=[],
                                normalise_global_fraction=True,
                                compartment_bound_tolerance=0.0)

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]["LB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["LB"])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]["UB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["UB"])})

        rba_session.set_growth_rate(Mu=Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
        if metabolite_to_minimize is None:
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        else:
            met_conc=rba_session.find_min_substrate_concentration(metabolite=metabolite_to_minimize, precision=0.001, max=1000, recording=False)
            if print_output:
                print("Min {} : {}".format(metabolite_to_minimize,met_conc))
            rba_session.set_medium({metabolite_to_minimize:met_conc})
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        if rba_session.Problem.Solved:
            try:
                rba_session.record_results('Eukaryotic')
                if print_output:
                    print('Mu Euk: {}'.format(growth_rate_from_input(input=definition_file, condition=condition)))
                euk_results = copy.deepcopy(rba_session.Results)

                compartment_fractions_euk = {}
                for comp in list(compartment_fractions_for_euk.keys()):
                    compartment_fractions_euk[comp] = rba_session.Problem.SolutionValues[str('f_'+comp)]
                rba_session.clear_results_and_parameters()
            except:
                compartment_fractions_euk={}
                euk_results = {}
            if variability_analysis is not None:
                if len(list(euk_results.keys()))!=0:
                    rba_session.set_growth_rate(Mu=mu_factor_for_variability*Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
                    euk_Feasible_Ranges=rba_session.get_feasible_range(variability_analysis)

    if "Eukaryotic_fixed_sizes" in sims_to_perform:
        rba_session.reload_model()

        if not apply_model:
            #Densities & PG
            if compartment_sizes is not None:
                if pg_fractions is not None:
                    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
                        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
            # Process efficiencies & Def/Spec Kapps
            process_efficiencies_to_inject=None
            Default_Kapps_to_inject=None
            Specific_Kapps_to_inject=None
            if process_efficiencies is not None:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
            if Default_Kapps is not None:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
            # Spec Kapps
            if Specific_Kapps is not None:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
            inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
        else:
            if "Specific_Kapps" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=Specific_Kapps,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=Specific_Kapps,condition=condition)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=None, process_efficiencies=None)
            if "Default_Kapps" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=Default_Kapps,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=Default_Kapps,condition=condition,transporter_multiplier=transporter_multiplier)
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=Default_Kapps_to_inject, process_efficiencies=None)

            if "Compartment_Sizes" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=compartment_sizes,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]

            if "PG_Fractions" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=None,
                                                                    compartment_densities=None,
                                                                    pg_fractions=pg_fractions,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                if compartment_sizes is not None:
                    if pg_fractions is not None:
                        compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes,pg_fractions=pg_fractions,condition=condition)
                        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
                            rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]

            if "Process_Efficiencies" in functions_to_include_list:
                inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                                    specific_kapps=None,
                                                                    default_kapps=None,
                                                                    process_efficiencies=process_efficiencies,
                                                                    compartment_densities=None,
                                                                    pg_fractions=None,
                                                                    round_to_digits=2,
                                                                    transporter_coeff=3)
            else:
                process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies,condition=condition,parameter_name_suffix="_apparent_efficiency")
                inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies_to_inject)

        rba_session.rebuild_from_model()
        # Medium
        rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))

        pg_fractions_for_euk={comp: str('fraction_non_enzymatic_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        compartment_fractions_for_euk={comp:str('fraction_protein_'+comp) for comp in list(compartment_densities_and_PGs['Compartment_ID'])}
        rba_session.make_eukaryotic(amino_acid_concentration_total='amino_acid_concentration',
                                pg_fractions=pg_fractions_for_euk,
                                compartment_fractions=compartment_fractions_for_euk,
                                compartments_with_imposed_sizes=list(compartment_fractions_for_euk.keys()),
                                normalise_global_fraction=True,
                                compartment_bound_tolerance=0.1,
                                imposed_compartments_without_tolerance=['Secreted','n'])

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]["LB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["LB"])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]["UB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["UB"])})

        rba_session.set_growth_rate(Mu=Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
        if metabolite_to_minimize is None:
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        else:
            met_conc=rba_session.find_min_substrate_concentration(metabolite=metabolite_to_minimize, precision=0.001, max=1000, recording=False)
            if print_output:
                print("Min {} : {}".format(metabolite_to_minimize,met_conc))
            rba_session.set_medium({metabolite_to_minimize:met_conc})
            rba_session.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        if rba_session.Problem.Solved:
            try:
                rba_session.record_results('Eukaryotic_fixed_sizes')
                if print_output:
                    print('Mu Euk fixed: {}'.format(growth_rate_from_input(input=definition_file, condition=condition)))
                euk_fixed_results = copy.deepcopy(rba_session.Results)

                compartment_fractions_euk_fixed = {}
                for comp in list(compartment_fractions_for_euk.keys()):
                    compartment_fractions_euk_fixed[comp] = rba_session.Problem.SolutionValues[str('f_'+comp)]
                rba_session.clear_results_and_parameters()
            except:
                compartment_fractions_euk_fixed={}
                euk_fixed_results = {}
            if variability_analysis is not None:
                if len(list(euk_results.keys()))!=0:
                    rba_session.set_growth_rate(Mu=mu_factor_for_variability*Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
                    euk_fixed_Feasible_Ranges=rba_session.get_feasible_range(variability_analysis)

    return({"FeasibleRange_def":def_Feasible_Ranges,
            "FeasibleRange_prok":prok_Feasible_Ranges,
            "FeasibleRange_euk":euk_Feasible_Ranges,
            "FeasibleRange_euk_fixed":euk_fixed_Feasible_Ranges,
            'Simulation_Results': prok_results, 
            'Simulation_Results_Euk': euk_results, 
            'Simulation_Results_Euk_fixed': euk_fixed_results, 
            'Simulation_Results_DefKapp': def_results, 
            'Euk_CompSizes': compartment_fractions_euk,
            'Euk_fixed_CompSizes': compartment_fractions_euk_fixed,
            "Condition":condition})


def find_ribosomal_proteins(rba_session, model_processes=['TranslationC', 'TranslationM'], external_annotations=None):
    """
    _summary_

    Parameters
    ----------
    rba_session : _type_
        _description_
    model_processes : list, optional
        _description_, by default ['TranslationC', 'TranslationM']
    external_annotations : _type_, optional
        _description_, by default None
    """
    out = []
    for i in model_processes:
        out += [rba_session.get_protein_information(protein=j)['ProtoID']
                for j in list(rba_session.get_process_information(process=i)['Composition'].keys()) if j in rba_session.get_proteins()]
    if external_annotations is not None:
        out += list(external_annotations['ID'])
    return(list(set(out)))


def build_model_compartment_map(rba_session):
    """
    _summary_

    Parameters
    ----------
    rba_session : _type_
        _description_
    """
    out = {rba_session.get_protein_information(protein=i)['ProtoID']: rba_session.get_protein_information(protein=i)['Compartment'] for i in list(
        rba_session.get_proteins())}
    return(out)


def build_compartment_annotations(Compartment_Annotations_external, model_protein_compartment_map):
    """
    _summary_

    Parameters
    ----------
    Compartment_Annotations_external : _type_
        _description_
    model_protein_compartment_map : _type_
        _description_
    """
    for i in Compartment_Annotations_external.index:
        if Compartment_Annotations_external.loc[i, 'ID'] in list(model_protein_compartment_map.keys()):
            Compartment_Annotations_external.loc[i, 'modelproteinannotation'] = 1
        else:
            Compartment_Annotations_external.loc[i, 'modelproteinannotation'] = 0
    Compartment_Annotations_internal = pandas.DataFrame()
    Compartment_Annotations_internal['ID'] = list(model_protein_compartment_map.keys())
    Compartment_Annotations_internal['ModelComp'] = list(model_protein_compartment_map.values())
    Compartment_Annotations = pandas.concat(
        [Compartment_Annotations_internal, Compartment_Annotations_external.loc[Compartment_Annotations_external['modelproteinannotation'] == 0, ['ID', 'ModelComp']]], axis=0)
    return(Compartment_Annotations)


def build_dataset_annotations(input, ID_column, Uniprot, Compartment_Annotations, model_protein_compartment_map,ribosomal_proteins=[]):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    ID_column : _type_
        _description_
    Uniprot : _type_
        _description_
    Compartment_Annotations : _type_
        _description_
    model_protein_compartment_map : _type_
        _description_
    ribosomal_proteins : list, optional
        _description_, by default []
    """
    out = pandas.DataFrame()
    for g in list(input[ID_column]):
        out.loc[g, 'ID'] = g
        matches = [i for i in list(Uniprot.loc[pandas.isna(
            Uniprot['Gene names']) == False, 'Gene names']) if g in i]
        mass_prot = numpy.nan
        if len(matches) > 0:
            mass_prot = len(Uniprot.loc[Uniprot['Gene names'] == matches[0], 'Sequence'].values[0])
        out.loc[g, 'AA_residues'] = mass_prot
        if g in list(Compartment_Annotations['ID']):
            out.loc[g, 'Location'] = Compartment_Annotations.loc[Compartment_Annotations['ID']
                                                                 == g, 'ModelComp'].values[0]
        in_model = 0
        if g in model_protein_compartment_map.keys():
            in_model = 1
        is_ribosomal = 0
        if g in ribosomal_proteins:
            is_ribosomal = 1
        out.loc[g, 'InModel'] = in_model
        out.loc[g, 'IsRibosomal'] = is_ribosomal
    return(out)


def build_full_annotations_from_dataset_annotations(annotations_list):
    """
    _summary_

    Parameters
    ----------
    annotations_list : _type_
        _description_
    """
    out = pandas.concat(annotations_list, axis=0)
    index = out.index
    is_duplicate = index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    out = out[not_duplicate]
    return(out)


def infer_copy_numbers_from_reference_copy_numbers(fold_changes, absolute_data, matching_column_in_fold_change_data, matching_column_in_absolute_data, conditions_in_fold_change_data_to_restore):
    """
    _summary_

    Parameters
    ----------
    fold_changes : _type_
        _description_
    absolute_data : _type_
        _description_
    matching_column_in_fold_change_data : _type_
        _description_
    matching_column_in_absolute_data : _type_
        _description_
    conditions_in_fold_change_data_to_restore : _type_
        _description_
    """
    out = pandas.DataFrame()
    for i in list(absolute_data['Gene']):
        if i in list(fold_changes['Gene']):
            FoldChange_match = fold_changes.loc[fold_changes['Gene']== i, matching_column_in_fold_change_data].values[0]
            CopyNumber_match = absolute_data.loc[absolute_data['Gene']== i, matching_column_in_absolute_data].values[0]
            if not pandas.isna(FoldChange_match):
                if not pandas.isna(CopyNumber_match):
                    out.loc[i, 'ID'] = i
                    out.loc[i, 'Absolute_Reference_restored'] = CopyNumber_match/(2**FoldChange_match)
    for gene in list(out['ID']):
        Abs_Ref = out.loc[gene, 'Absolute_Reference_restored']
        for condition in conditions_in_fold_change_data_to_restore:
            if condition in list(fold_changes.columns):
                out.loc[gene, condition] = Abs_Ref * (2**fold_changes.loc[fold_changes['Gene'] == gene, condition].values[0])
    return(out)


def add_annotations_to_proteome(input, ID_column, annotations):
    """
    _summary_

    Parameters
    ----------
    input : _type_
        _description_
    ID_column : _type_
        _description_
    annotations : _type_
        _description_
    """
    for i in input.index:
        if input.loc[i, ID_column] in annotations.index:
            input.loc[i, 'AA_residues'] = annotations.loc[input.loc[i, ID_column], 'AA_residues']
            input.loc[i, 'Location'] = annotations.loc[input.loc[i, ID_column], 'Location']
            input.loc[i, 'InModel'] = annotations.loc[input.loc[i, ID_column], 'InModel']
            input.loc[i, 'IsRibosomal'] = annotations.loc[input.loc[i, ID_column], 'IsRibosomal']
    return(input)


def extract_proteomes_from_calibration_results(calib_results):
    """
    _summary_

    Parameters
    ----------
    calib_results : _type_
        _description_
    """
    proteomes=pandas.DataFrame()
    for calib_result in calib_results:
        condition=calib_result["Condition"]
        proteome=calib_result["Proteome"]
        for i in list(proteome.index):
            proteomes.loc[i,condition]=proteome.loc[i,"copy_number"]
    return(proteomes)


def get_flux_distribution(simulation_outputs,result_object='Simulation_Results', run='Prokaryotic'):
    """
    _summary_

    Parameters
    ----------
    simulation_outputs : _type_
        _description_
    result_object : str, optional
        _description_, by default 'Simulation_Results'
    run : str, optional
        _description_, by default 'Prokaryotic'
    """
    out=pandas.DataFrame(columns=[sim_result["Condition"] for sim_result in simulation_outputs])
    for sim_result in simulation_outputs:
        if len(list(sim_result[result_object].keys()))>0:
            for rx in sim_result[result_object]['uniqueReactions'].index:
                out.loc[rx,sim_result["Condition"]]=sim_result[result_object]['uniqueReactions'].loc[rx, run]
    return(out)


def get_exchange_flux_distribution(simulation_outputs,result_object='Simulation_Results', run='Prokaryotic'):
    """
    _summary_

    Parameters
    ----------
    simulation_outputs : _type_
        _description_
    result_object : str, optional
        _description_, by default 'Simulation_Results'
    run : str, optional
        _description_, by default 'Prokaryotic'
    """
    out=pandas.DataFrame(columns=[sim_result["Condition"] for sim_result in simulation_outputs])
    for sim_result in simulation_outputs:
        if len(list(sim_result[result_object].keys()))>0:
            for rx in sim_result[result_object]['ExchangeFluxes'].index:
                out.loc[rx,sim_result["Condition"]]=sim_result[result_object]['ExchangeFluxes'].loc[rx, run]
    return(out)


def extract_proteomes_from_simulation_results(simulation_outputs,type="Prokaryotic"):
    """
    _summary_

    Parameters
    ----------
    simulation_outputs : _type_
        _description_
    type : str, optional
        _description_, by default "Prokaryotic"
    """
    if type=="Prokaryotic":
        results_object="Simulation_Results"
    elif type=="Eukaryotic":
        results_object="Simulation_Results_Euk"
    elif type=="Eukaryotic_fixed_sizes":
        results_object="Simulation_Results_Euk_fixed"
    elif type=="DefaultKapp":
        results_object="Simulation_Results_DefKapp"
    out=pandas.DataFrame()
    for sim_result in simulation_outputs:
        try:
            condition=sim_result["Condition"]
            proto_proteome=sim_result[results_object]["ProtoProteins"]
            df_col=list(proto_proteome.columns)[0]
            for i in proto_proteome.index:
                out.loc[i,condition]=proto_proteome.loc[i,df_col]
        except:
            pass
    return(out)



##???##
def build_input_for_default_kapp_estimation(input, type='new',max_PG_fraction=0.9,comps_minPG=["Secreted","n","vM","gM"]):
    out = pandas.DataFrame(columns=['Compartment_ID', 'Density', 'PG_fraction'])
    for i in input.index:
        if i not in ['Total', 'Ribosomes']:
            out.loc[i, 'Compartment_ID'] = i
            if type == 'new':
                out.loc[i, 'Density'] = input.loc[i, 'new_protein_fraction']
                out.loc[i, 'PG_fraction'] = input.loc[i, 'new_PG_fraction']
            elif type == 'original':
                out.loc[i, 'Density'] = input.loc[i, 'original_protein_fraction']
                out.loc[i, 'PG_fraction'] = input.loc[i, 'original_PG_fraction']
            if i in comps_minPG:
                if out.loc[i, 'PG_fraction']>max_PG_fraction:
                    out.loc[i, 'PG_fraction']=max_PG_fraction
    total_size=sum(list(out["Density"]))
    out["Density"]/=total_size
    return(out)
