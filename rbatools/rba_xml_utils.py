import rba
import pandas
import numpy
import json

def inject_default_kapps(rba_session, default_kapp, default_transporter_kapp):
    if numpy.isfinite(default_kapp):
        rba_session.model.parameters.functions._elements_by_id['default_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_kapp
    if numpy.isfinite(default_transporter_kapp):
        rba_session.model.parameters.functions._elements_by_id['default_transporter_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_transporter_kapp


def inject_process_capacities(rba_session, process_efficiencies, round_to_digits=0, min_value=100):
    """
    Parameters
    ----------
    process_efficiencies : pandas.DataFrame(columns=['Process','Parameter','Value'])
    """
    for i in process_efficiencies.index:
        if numpy.isfinite(process_efficiencies.loc[i, 'Value']):
            val = round(process_efficiencies.loc[i, 'Value'], round_to_digits)
            if val >= min_value:
                if process_efficiencies.loc[i, 'Process'] in rba_session.model.processes.processes._elements_by_id.keys():
                    rba_session.model.processes.processes._elements_by_id[process_efficiencies.loc[i,'Process']].machinery.capacity.value = process_efficiencies.loc[i, 'Parameter']
                    const = rba.xml.parameters.Function(process_efficiencies.loc[i, 'Parameter'], 'constant', parameters={'CONSTANT': val}, variable=None)
                    if process_efficiencies.loc[i, 'Parameter'] not in rba_session.model.parameters.functions._elements_by_id.keys():
                        rba_session.model.parameters.functions.append(const)
                    else:
                        rba_session.model.parameters.functions._elements_by_id[const.id].parameters._elements_by_id['CONSTANT'].value = val


def inject_specific_kapps(rba_session, specific_kapps, round_to_digits=0, min_value=1000):
    """
    Parameters
    ----------
    specific_kapps : pandas.DataFrame
    """
    if 'Enzyme_ID' in list(specific_kapps.columns):
        for enz in list(specific_kapps['Enzyme_ID']):
            if numpy.isfinite(specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Kapp'].values[0]):
                val = round(specific_kapps.loc[specific_kapps['Enzyme_ID']== enz, 'Kapp'].values[0], round_to_digits)
                const = rba.xml.parameters.Function(str(enz + '_kapp__constant'), 'constant', parameters={'CONSTANT': val}, variable=None)
                if str(enz + '_kapp__constant') not in rba_session.model.parameters.functions._elements_by_id.keys():
                    rba_session.model.parameters.functions.append(const)
                else:
                    rba_session.model.parameters.functions._elements_by_id[const.id].parameters._elements_by_id['CONSTANT'].value = val
                for e in rba_session.model.enzymes.enzymes:
                    if e.id == enz:
                        if 'Flux' in list(specific_kapps.columns):
                            if specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Flux'].values[0] > 0:
                                e.forward_efficiency = str(enz + '_kapp__constant')
                            elif specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Flux'].values[0] < 0:
                                e.backward_efficiency = str(enz + '_kapp__constant')
                            else:
                                e.forward_efficiency = str(enz + '_kapp__constant')
                                e.backward_efficiency = str(enz + '_kapp__constant')
                        else:
                            e.forward_efficiency = str(enz + '_kapp__constant')
                            e.backward_efficiency = str(enz + '_kapp__constant')


def inject_default_kapps_as_function(rba_session, default_kapps,transporter_coeff, round_to_digits,x_min,x_max,indep_variable="growth_rate"):
    model_dict=json.loads(str(default_kapps.loc["Regression","Model"]))
    respective_function_default=None
    respective_function_default_transport=None
    parameter_ID_default="default_efficiency"
    parameter_ID_default_transport="default_transporter_efficiency"

    if list(model_dict.keys())[0]=="constant":
        if numpy.isfinite(model_dict["constant"]["CONSTANT"]):
            respective_function_default=rba.xml.parameters.Function(parameter_ID_default, 'constant', parameters={'CONSTANT': model_dict["constant"]["CONSTANT"]}, variable=indep_variable)
            respective_function_default_transport=rba.xml.parameters.Function(parameter_ID_default_transport, 'constant', parameters={'CONSTANT': transporter_coeff*model_dict["constant"]["CONSTANT"]}, variable=indep_variable)
    elif list(model_dict.keys())[0]=="linear":
        if numpy.isfinite(model_dict["linear"]["A"]):
            if numpy.isfinite(model_dict["linear"]["B"]):
                respective_function_default=rba.xml.parameters.Function(parameter_ID_default, 'linear', parameters={'LINEAR_COEF': model_dict["linear"]["A"],'LINEAR_CONSTANT': model_dict["linear"]["B"],'X_MIN': x_min,'X_MAX': x_max,'Y_MIN': model_dict["linear"]["Y_min"],'Y_MAX': model_dict["linear"]["Y_max"],}, variable=indep_variable)
                respective_function_default_transport=rba.xml.parameters.Function(parameter_ID_default_transport, 'linear', parameters={'LINEAR_COEF': transporter_coeff*model_dict["linear"]["A"],'LINEAR_CONSTANT': transporter_coeff*model_dict["linear"]["B"],'X_MIN': x_min,'X_MAX': x_max,'Y_MIN': transporter_coeff*model_dict["linear"]["Y_min"],'Y_MAX': transporter_coeff*model_dict["linear"]["Y_max"],}, variable=indep_variable)
    if respective_function_default is not None:
        if parameter_ID_default in rba_session.model.parameters.functions._elements_by_id.keys():
            rba_session.model.parameters.functions.remove(rba_session.model.parameters.functions._elements_by_id[parameter_ID_default])
        rba_session.model.parameters.functions.append(respective_function_default)
    if respective_function_default_transport is not None:
        if parameter_ID_default_transport in rba_session.model.parameters.functions._elements_by_id.keys():
            rba_session.model.parameters.functions.remove(rba_session.model.parameters.functions._elements_by_id[parameter_ID_default_transport])
        rba_session.model.parameters.functions.append(respective_function_default_transport)


def inject_process_capacities_as_function(rba_session, process_efficiencies, round_to_digits,x_min,x_max,indep_variable="growth_rate"):
    """
    Parameters
    ----------
    process_efficiencies : pandas.DataFrame(columns=['Process','Parameter','Value'])
    """
    for i in process_efficiencies.index:
        parameter_ID=str(i+"_capacity")
        model_dict=json.loads(str(process_efficiencies.loc[i,"Model"]))
        respective_function=None
        if list(model_dict.keys())[0]=="constant":
            if numpy.isfinite(model_dict["constant"]["CONSTANT"]):
                respective_function=rba.xml.parameters.Function(parameter_ID, 'constant', parameters={'CONSTANT': model_dict["constant"]["CONSTANT"]}, variable=indep_variable)
        elif list(model_dict.keys())[0]=="linear":
            if numpy.isfinite(model_dict["linear"]["A"]):
                if numpy.isfinite(model_dict["linear"]["B"]):
                    respective_function=rba.xml.parameters.Function(parameter_ID, 'linear', parameters={'LINEAR_COEF': model_dict["linear"]["A"],'LINEAR_CONSTANT': model_dict["linear"]["B"],'X_MIN': x_min,'X_MAX': x_max,'Y_MIN': model_dict["linear"]["Y_min"],'Y_MAX': model_dict["linear"]["Y_max"],}, variable=indep_variable)
        if respective_function is not None:
            if i in rba_session.model.processes.processes._elements_by_id.keys():
                rba_session.model.processes.processes._elements_by_id[i].machinery.capacity.value = parameter_ID
                if parameter_ID in rba_session.model.parameters.functions._elements_by_id.keys():
                    rba_session.model.parameters.functions.remove(rba_session.model.parameters.functions._elements_by_id[parameter_ID])
                rba_session.model.parameters.functions.append(respective_function)


def inject_compartment_densities_as_function(rba_session, compartment_densities, round_to_digits,x_min,x_max,indep_variable="growth_rate"):
    comp_den_param_map={den.compartment:str("fraction_protein_"+str(den.compartment)) for den in rba_session.model.density.target_densities._elements}
    for i in compartment_densities.index:
        if i in comp_den_param_map.keys():
            parameter_ID=comp_den_param_map[i]
            model_dict=json.loads(str(compartment_densities.loc[i,"Model"]))
            respective_function=None
            if list(model_dict.keys())[0]=="constant":
                if numpy.isfinite(model_dict["constant"]["CONSTANT"]):
                    respective_function=rba.xml.parameters.Function(parameter_ID, 'constant', parameters={'CONSTANT': model_dict["constant"]["CONSTANT"]}, variable=indep_variable)
            elif list(model_dict.keys())[0]=="linear":
                if numpy.isfinite(model_dict["linear"]["A"]):
                    if numpy.isfinite(model_dict["linear"]["B"]):
                        respective_function=rba.xml.parameters.Function(parameter_ID, 'linear', parameters={'LINEAR_COEF': model_dict["linear"]["A"],'LINEAR_CONSTANT': model_dict["linear"]["B"],'X_MIN': x_min,'X_MAX': x_max,'Y_MIN': 0.0,'Y_MAX': 1.0,}, variable=indep_variable)
            if respective_function is not None:
                if parameter_ID in rba_session.model.parameters.functions._elements_by_id.keys():
                    rba_session.model.parameters.functions.remove(rba_session.model.parameters.functions._elements_by_id[parameter_ID])
                rba_session.model.parameters.functions.append(respective_function)
        else:
            if str("fraction_protein_"+i) in rba_session.model.parameters.functions._elements_by_id.keys():
                parameter_ID=str("fraction_protein_"+i)
                model_dict=json.loads(str(compartment_densities.loc[i,"Model"]))
                respective_function=None
                if list(model_dict.keys())[0]=="constant":
                    if numpy.isfinite(model_dict["constant"]["CONSTANT"]):
                        respective_function=rba.xml.parameters.Function(parameter_ID, 'constant', parameters={'CONSTANT': model_dict["constant"]["CONSTANT"]}, variable=indep_variable)
                elif list(model_dict.keys())[0]=="linear":
                    if numpy.isfinite(model_dict["linear"]["A"]):
                        if numpy.isfinite(model_dict["linear"]["B"]):
                            respective_function=rba.xml.parameters.Function(parameter_ID, 'linear', parameters={'LINEAR_COEF': model_dict["linear"]["A"],'LINEAR_CONSTANT': model_dict["linear"]["B"],'X_MIN': x_min,'X_MAX': x_max,'Y_MIN': 0.0,'Y_MAX': 1.0,}, variable=indep_variable)
                if respective_function is not None:
                    rba_session.model.parameters.functions.remove(rba_session.model.parameters.functions._elements_by_id[parameter_ID])
                    rba_session.model.parameters.functions.append(respective_function)


def inject_pg_fractions_as_function(rba_session, pg_fractions, round_to_digits,x_min,x_max,indep_variable="growth_rate"):
    comp_pg_param_map={str(tar.value.split("nonenzymatic_proteins_")[1]):str("fraction_non_enzymatic_protein_"+tar.value.split("nonenzymatic_proteins_")[1]) for tar in rba_session.model.targets.target_groups._elements_by_id["translation_targets"].concentrations._elements}
    for i in pg_fractions.index:
        if i in comp_pg_param_map.keys():
            parameter_ID=comp_pg_param_map[i]
            model_dict=json.loads(str(pg_fractions.loc[i,"Model"]))
            respective_function=None
            if list(model_dict.keys())[0]=="constant":
                if numpy.isfinite(model_dict["constant"]["CONSTANT"]):
                    respective_function=rba.xml.parameters.Function(parameter_ID, 'constant', parameters={'CONSTANT': model_dict["constant"]["CONSTANT"]}, variable=indep_variable)
            elif list(model_dict.keys())[0]=="linear":
                if numpy.isfinite(model_dict["linear"]["A"]):
                    if numpy.isfinite(model_dict["linear"]["B"]):
                        respective_function=rba.xml.parameters.Function(parameter_ID, 'linear', parameters={'LINEAR_COEF': model_dict["linear"]["A"],'LINEAR_CONSTANT': model_dict["linear"]["B"],'X_MIN': x_min,'X_MAX': x_max,'Y_MIN': 0.0,'Y_MAX': 1.0,}, variable=indep_variable)
            if respective_function is not None:
                if parameter_ID in rba_session.model.parameters.functions._elements_by_id.keys():
                    rba_session.model.parameters.functions.remove(rba_session.model.parameters.functions._elements_by_id[parameter_ID])
                rba_session.model.parameters.functions.append(respective_function)


def inject_specific_kapps_as_function(rba_session, specific_kapps, round_to_digits,x_min,x_max,indep_variable="growth_rate"):
    """
    Parameters
    ----------
    specific_kapps : pandas.DataFrame
    """
    if 'Enzyme_ID' in list(specific_kapps.columns):
        for i in specific_kapps['Enzyme_ID']:
            parameter_ID=str(i + '_kapp')
            model_dict=json.loads(str(specific_kapps.loc[i,"Model"]))
            fluxes=[specific_kapps.loc[i,col] for col in list(specific_kapps.columns) if col.endswith("_Flux")]
            fluxes_neg=[j for j in fluxes if j<0]
            fluxes_pos=[j for j in fluxes if j>0]
            respective_function=None
            if list(model_dict.keys())[0]=="constant":
                if numpy.isfinite(model_dict["constant"]["CONSTANT"]):
                    respective_function=rba.xml.parameters.Function(parameter_ID, 'constant', parameters={'CONSTANT': model_dict["constant"]["CONSTANT"]}, variable=indep_variable)
            elif list(model_dict.keys())[0]=="linear":
                if numpy.isfinite(model_dict["linear"]["A"]):
                    if numpy.isfinite(model_dict["linear"]["B"]):
                        respective_function=rba.xml.parameters.Function(parameter_ID, 'linear', parameters={'LINEAR_COEF': model_dict["linear"]["A"],'LINEAR_CONSTANT': model_dict["linear"]["B"],'X_MIN': x_min,'X_MAX': x_max,'Y_MIN': model_dict["linear"]["Y_min"],'Y_MAX': model_dict["linear"]["Y_max"],}, variable=indep_variable)
            if respective_function is not None:
                if parameter_ID in rba_session.model.parameters.functions._elements_by_id.keys():
                    rba_session.model.parameters.functions.remove(rba_session.model.parameters.functions._elements_by_id[parameter_ID])
                rba_session.model.parameters.functions.append(respective_function)
                for enz in rba_session.model.enzymes.enzymes:
                    if enz.id == i:
                        if len(fluxes_pos)>0:
                            enz.forward_efficiency = parameter_ID
                        if len(fluxes_neg)>0:
                            enz.backward_efficiency = parameter_ID
                        if len(fluxes_pos)==0:
                            if len(fluxes_neg)==0:
                                enz.forward_efficiency = parameter_ID
                                enz.backward_efficiency = parameter_ID


def inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=None, round_to_digits=2):
    """
    Parameters
    ----------
    specific_kapps : pandas.DataFrame(columns=['Enzyme_ID','Kapp'])
    default_kapps : {'default_efficiency':value,'default_transporter_efficiency':value}
    process_efficiencies : pandas.DataFrame(columns=['Process','Parameter','Value'])
    """
    if specific_kapps is not None:
        inject_specific_kapps(rba_session=rba_session,specific_kapps=specific_kapps, round_to_digits=round_to_digits)

    if default_kapps is not None:
        if type(default_kapps) is dict:
            inject_default_kapps(rba_session=rba_session,default_kapp=default_kapps['default_efficiency'], default_transporter_kapp=default_kapps['default_transporter_efficiency'])

    if process_efficiencies is not None:
        inject_process_capacities(rba_session=rba_session,process_efficiencies=process_efficiencies)


def inject_estimated_efficiencies_as_functions_into_model(rba_session,
                                                          specific_kapps=None,
                                                          default_kapps=None,
                                                          process_efficiencies=None,
                                                          compartment_densities=None,
                                                          pg_fractions=None,
                                                          round_to_digits=2,
                                                          transporter_coeff=3):
    """
    """
    if compartment_densities is not None:
        inject_compartment_densities_as_function(rba_session=rba_session,compartment_densities=compartment_densities, round_to_digits=round_to_digits,x_min=0.05,x_max=0.3,indep_variable="growth_rate")

    if pg_fractions is not None:
        inject_pg_fractions_as_function(rba_session=rba_session,pg_fractions=pg_fractions, round_to_digits=round_to_digits,x_min=0.05,x_max=0.3,indep_variable="growth_rate")

    if specific_kapps is not None:
        inject_specific_kapps_as_function(rba_session=rba_session,specific_kapps=specific_kapps, round_to_digits=round_to_digits,x_min=0.05,x_max=0.3,indep_variable="growth_rate")

    if default_kapps is not None:
        inject_default_kapps_as_function(rba_session=rba_session,default_kapps=default_kapps,transporter_coeff=transporter_coeff, round_to_digits=round_to_digits,x_min=0.05,x_max=0.3,indep_variable="growth_rate")

    if process_efficiencies is not None:
        inject_process_capacities_as_function(rba_session=rba_session,process_efficiencies=process_efficiencies, round_to_digits=round_to_digits,x_min=0.05,x_max=0.3,indep_variable="growth_rate")

