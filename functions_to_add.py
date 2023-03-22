"""
    #global
    fix_mus={'Hackett_C03':False,'Hackett_C005':True, 'Hackett_C022':True, 'Hackett_C01':True, 'Hackett_C016':True}
    adjust_roots={'Hackett_C03':False,'Hackett_C005':True, 'Hackett_C022':True, 'Hackett_C01':True, 'Hackett_C016':True}
    mu_misprediction_tolerances={'Hackett_C03':0.1,'Hackett_C005':0.05, 'Hackett_C022':0.05, 'Hackett_C01':0.05, 'Hackett_C016':0.05}

    #overall
    minimum_iteration_numbers={'Hackett_C03':2,'Hackett_C005':2, 'Hackett_C022':2, 'Hackett_C01':2, 'Hackett_C016':2}
    iteration_limits={'Hackett_C03':10,'Hackett_C005':10, 'Hackett_C022':10, 'Hackett_C01':10, 'Hackett_C016':10}
    steady_RSS_limits={'Hackett_C03':5,'Hackett_C005':5, 'Hackett_C022':5, 'Hackett_C01':5, 'Hackett_C016':5}
    increasing_RSS_limits={'Hackett_C03':3,'Hackett_C005':3, 'Hackett_C022':3, 'Hackett_C01':3, 'Hackett_C016':3}
    rss_tolerances={'Hackett_C03':0.1,'Hackett_C005':0.1, 'Hackett_C022':0.1, 'Hackett_C01':0.1, 'Hackett_C016':0.1}
    increasing_rss_factors={'Hackett_C03':1.1,'Hackett_C005':1.1, 'Hackett_C022':1.1, 'Hackett_C01':1.1, 'Hackett_C016':1.1}


def machinery_efficiency_correction_settings_from_input(input, condition):
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
    out={}
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


def compartment_correction_settings_from_input(input, condition):
    out={}
    out['compartments_to_replace']=bool(input.loc[input['Type'] == '', condition].values[0])
    out['compartments_no_original_PG']=input.loc[input['Type'] == '', condition].values[0]
    out['fractions_entirely_replaced_with_expected_value']=bool(input.loc[input['Type'] == '', condition].values[0])
    out['imposed_compartment_fractions']=bool(input.loc[input['Type'] == '', condition].values[0])
    out['directly_corrected_compartments']=bool(input.loc[input['Type'] == '', condition].values[0])
    out['merged_compartments']=bool(input.loc[input['Type'] == '', condition].values[0])
    out['min_compartment_fraction']=input.loc[input['Type'] == '', condition].values[0]
    return(out)


def general_calibration_settings_from_input(input, condition):

#transporter_multiplier
#Mu_approx_precision
#proteome scaling

"""