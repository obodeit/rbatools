import rba
import copy
import pandas
import time
import numpy
import seaborn
import json
import math
import matplotlib.pyplot as plt
import scipy.signal
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats.mstats import gmean
# import matplotlib.pyplot as plt

def generate_proto_protein_map(rba_session):
    protomap={}
    for i in rba_session.get_proteins():
        protoID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if protoID in protomap.keys():
            protomap[protoID].append(i)
        else:
            protomap[protoID]=[i]
    return(protomap)


def determine_kapp_type_classification(protein_DF,rbaSession,specific_kapps,protomap):
    comments=[]
    comp_comments=[]
    for i in protein_DF.index:
        comment=""
        comp_comment=""
        if i in protomap.keys():
            if len(protomap[i])>1:
                comp_comment="several"
            elif len(protomap[i])==1:
                comp_comment="one"
            for j in protomap[i]:
                if j in rbaSession.get_proteins():
                    if len(rbaSession.get_protein_information(protein=j)['associatedEnzymes']) > 0:
                        for k in rbaSession.get_protein_information(protein=j)['associatedEnzymes']:
                            if k in list(specific_kapps["Enzyme_ID"]):
                                if comment in ["","Specific"]:
                                    comment="Specific"
                                elif comment in ["Default","Both"]:
                                    comment="Both"
                                elif comment in ["NonMetabolic"]:
                                    comment="All"
                            else:
                                if comment in ["","Default"]:
                                    comment="Default"
                                elif comment in ["Specific","Both"]:
                                    comment="Both"
                                elif comment in ["NonMetabolic"]:
                                    comment="All"
                    else:
                        if comment in ["Specific","Default","Both","All"]:
                            comment="All"
                        elif comment in ["NonMetabolic",""]:
                            comment="NonMetabolic"
        comments.append(comment)
        comp_comments.append(comp_comment)
    protein_DF["ProteinComment"]=comments
    protein_DF["CompartmentComment"]=comp_comments
    return(protein_DF)


def extract_predicted_growth_rates(inputs=[],result_object='Simulation_Results', run='Prokaryotic'):
    out = []
    for i in inputs:
        try:
            out.append(i[result_object]['Mu'].loc['Mu', run])
        except:
            out.append(numpy.nan)
    return(out)


def find_ribosomal_proteins(rba_session, model_processes=['TranslationC', 'TranslationM'], external_annotations=None):
    out = []
    for i in model_processes:
        out += [rba_session.get_protein_information(protein=j)['ProtoID']
                for j in list(rba_session.get_process_information(process=i)['Composition'].keys()) if j in rba_session.get_proteins()]
    if external_annotations is not None:
        out += list(external_annotations['ID'])
    return(list(set(out)))


def build_model_compartment_map(rba_session):
    out = {rba_session.get_protein_information(protein=i)['ProtoID']: rba_session.get_protein_information(protein=i)['Compartment'] for i in list(
        rba_session.get_proteins())}
    return(out)


def build_compartment_annotations(Compartment_Annotations_external, model_protein_compartment_map):
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
    out = pandas.concat(annotations_list, axis=0)
    index = out.index
    is_duplicate = index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    out = out[not_duplicate]
    return(out)


def infer_copy_numbers_from_reference_copy_numbers(fold_changes, absolute_data, matching_column_in_fold_change_data, matching_column_in_absolute_data, conditions_in_fold_change_data_to_restore):
    out = pandas.DataFrame()
    for i in list(absolute_data['Gene']):
        if i in list(fold_changes['Gene']):
            FoldChange_match = fold_changes.loc[fold_changes['Gene']== i, matching_column_in_fold_change_data].values[0]
            CopyNumber_match = absolute_data.loc[absolute_data['Gene']== i, matching_column_in_absolute_data].values[0]
            if not pandas.isna(FoldChange_match):
                if not pandas.isna(CopyNumber_match):
                    out.loc[i, 'ID'] = i
                    out.loc[i, 'Absolute_Reference'] = CopyNumber_match/(2**FoldChange_match)
    for gene in list(out['ID']):
        Abs_Ref = out.loc[gene, 'Absolute_Reference']
        for condition in conditions_in_fold_change_data_to_restore:
            out.loc[gene, condition] = Abs_Ref * \
                (2**fold_changes.loc[fold_changes['Gene'] == gene, condition].values[0])
    return(out)


def add_annotations_to_proteome(input, ID_column, annotations):
    for i in input.index:
        if input.loc[i, ID_column] in annotations.index:
            input.loc[i, 'AA_residues'] = annotations.loc[input.loc[i, ID_column], 'AA_residues']
            input.loc[i, 'Location'] = annotations.loc[input.loc[i, ID_column], 'Location']
            input.loc[i, 'InModel'] = annotations.loc[input.loc[i, ID_column], 'InModel']
            input.loc[i, 'IsRibosomal'] = annotations.loc[input.loc[i, ID_column], 'IsRibosomal']
    return(input)


def determine_compartment_occupation(Data,
                                       Condition,
                                       mass_col='AA_residues',
                                       only_in_model=False,
                                       compartments_to_replace={'DEF':"c"},
                                       compartments_no_original_PG=[],
                                       ribosomal_proteins_as_extra_compartment=True):
    out=pandas.DataFrame()
    if only_in_model:
        Data = Data.loc[Data['InModel'] >= 1]
    Data["AA_abundance"]=Data[Condition]*Data[mass_col]
#    Data.loc[:,"AA_abundance"]=Data.loc[:,Condition]*Data.loc[:,mass_col]

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


def build_proteome_overview(input, condition, compartments_to_replace={'DEF':"c", 'DEFA':"c", 'Def':"c"}, compartments_no_original_PG=['n', 'Secreted'], ribosomal_proteins_as_extra_compartment=True):
    out = determine_compartment_occupation(Data=input, Condition=condition, compartments_to_replace=compartments_to_replace,
                                           compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=ribosomal_proteins_as_extra_compartment, only_in_model=False)
    out_in_model = determine_compartment_occupation(Data=input, Condition=condition, compartments_to_replace=compartments_to_replace,
                                                    compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=ribosomal_proteins_as_extra_compartment, only_in_model=True)

    for comp in out.index:
        out.loc[comp,'original_protein_fraction'] = out.loc[comp,'original_amino_acid_occupation']/out.loc["Total",'original_amino_acid_occupation']
        out['original_PG_fraction'] = 1-out_in_model['original_amino_acid_occupation'] / out['original_amino_acid_occupation']
    return(out)


def determine_macromolecule_synthesis_machinery_demand(rba_session):
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


def determine_apparent_process_efficiencies(growth_rate, input, rba_session, proteome_summary, protein_data, condition, gene_id_col,fit_nucleotide_assembly_machinery=False):
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
                protein_data.loc[protein_data['ID'] == i, ]
                n_AAs_in_machinery += protein_data.loc[protein_data['ID'] == i, condition].values[0] * \
                    protein_data.loc[protein_data['ID'] == i, 'AA_residues'].values[0]
                machinery_size += constituting_proteins[i]
        # right reference amounth?
        if n_AAs_in_machinery > 0:
            relative_Protein_fraction_of_machinery = n_AAs_in_machinery / proteome_summary.loc['Total', 'original_amino_acid_occupation']
            specific_capacity = growth_rate*Total_client_fraction/relative_Protein_fraction_of_machinery
            apparent_capacity = specific_capacity*machinery_size
            # process_ID[process_name] = apparent_capacity
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


def correction_pipeline(input,
                          condition,
                          definition_file,
                          reference_condition,
                          compartments_to_replace,
                          compartments_no_original_PG,
                          fractions_entirely_replaced_with_expected_value,
                          imposed_compartment_fractions,
                          directly_corrected_compartments,
                          merged_compartments,
                          min_compartment_fraction):
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

def flux_bounds_from_input(input, condition, specific_exchanges=None, specific_directions=None):
    flux_mean_df = input.loc[input['Type'] == 'ExchangeFlux_Mean', :]
    flux_mean_SE = input.loc[input['Type'] == 'ExchangeFlux_StandardError', :]
    out = pandas.DataFrame(columns=['Reaction_ID', 'LB', 'UB'])
    if specific_exchanges is None:
        exchanges_to_set = list(flux_mean_df['ID'])
    else:
        exchanges_to_set = specific_exchanges
    for rx in exchanges_to_set:
        mean_val = flux_mean_df.loc[flux_mean_df['ID'] == rx, condition].values[0]
        if not pandas.isna(mean_val):
            out.loc[rx, 'Reaction_ID'] = rx
            out.loc[rx, 'Mean'] = mean_val
            if str(rx+'_SE') in flux_mean_SE['ID']:
                SE_val = flux_mean_SE.loc[flux_mean_SE['ID'] == str(rx+'_SE'), condition].values[0]
            else:
                SE_val=numpy.nan
            if not pandas.isna(SE_val):
                lb = mean_val-SE_val
                ub = mean_val+SE_val
                if mean_val < 0:
                    out.loc[rx, 'LB'] = lb
                    if ub > 0:
                        out.loc[rx, 'UB'] = 0
                    else:
                        out.loc[rx, 'UB'] = ub
                elif mean_val > 0:
                    out.loc[rx, 'UB'] = ub
                    if lb < 0:
                        out.loc[rx, 'LB'] = 0
                    else:
                        out.loc[rx, 'LB'] = lb
                else:
                    out.loc[rx, 'LB'] = lb
                    out.loc[rx, 'UB'] = ub
            else:
                out.loc[rx, 'LB'] = mean_val
                out.loc[rx, 'UB'] = mean_val
    flux_dir_df = input.loc[input['Type'] == 'Flux_Direction', :]
    if specific_directions is None:
        directions_to_set = list(flux_dir_df['ID'])
    else:
        directions_to_set = specific_directions
    for rx in directions_to_set:
        out.loc[rx, 'Reaction_ID'] = rx
        if flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == 1:
            out.loc[rx, 'LB'] = 0
        elif flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == -1:
            out.loc[rx, 'UB'] = 0
        elif flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == 0:
            out.loc[rx, 'LB'] = 0
            out.loc[rx, 'UB'] = 0
    flux_upper_df = input.loc[input['Type'] == 'Flux_Upper_Bound', :]
    for rx in list(flux_upper_df['ID']):
        out.loc[rx, 'Reaction_ID'] = rx
        out.loc[rx, 'UB'] = flux_upper_df.loc[flux_upper_df['ID'] == rx, condition].values[0]
    flux_lower_df = input.loc[input['Type'] == 'Flux_Lower_Bound', :]
    for rx in list(flux_lower_df['ID']):
        out.loc[rx, 'Reaction_ID'] = rx
        out.loc[rx, 'LB'] = flux_lower_df.loc[flux_lower_df['ID'] == rx, condition].values[0]
    return(out)


def growth_rate_from_input(input, condition):
    return(input.loc[input['Type'] == 'Growth_Rate', condition].values[0])


def proteome_fractions_from_input(input, condition):
    df = input.loc[input['Type'] == 'Expected_ProteomeFraction', :]
    return(dict(zip(list(df['ID']), list(df[condition]))))


def medium_concentrations_from_input(input, condition):
    df = input.loc[input['Type'] == 'Medium_Concentration', :]
    return(dict(zip(list(df['ID']), list(df[condition]))))


def build_input_proteome_for_specific_kapp_estimation(proteomics_data, condition):
    out = pandas.DataFrame()
    out['ID'] = proteomics_data['ID']
    out['copy_number'] = proteomics_data[condition]
    return(out)


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


def perform_simulations(condition,
                        rba_session,
                        definition_file,
                        compartment_sizes,
                        pg_fractions,
                        process_efficiencies,
                        Default_Kapps,
                        Specific_Kapps,
                        Exchanges_to_impose=None,
                        sims_to_perform=["DefaultKapp","Prokaryotic","Eukaryotic"],
                        feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                        try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                        print_output=True,
                        variability_analysis=None,
                        mu_factor_for_variability=1,
                        apply_model=False,
                        functions_to_include_list=[],
                        transporter_multiplier=3,
                        start_val=numpy.nan,
                        Mu_approx_precision=0.000001,
                        max_mu_in_dichotomy=1):
    mumax_def=numpy.nan
    mumax_prok=numpy.nan
    mumax_euk=numpy.nan
    def_Feasible_Ranges={}
    prok_Feasible_Ranges={}
    euk_Feasible_Ranges={}
    def_results={}
    prok_results={}
    euk_results={}
    compartment_fractions={}

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
            # Process efficiencies & Def Kapps
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
        mumax_def = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max=max_mu_in_dichotomy,start_value=start_val, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        try:
            rba_session.record_results('DefaultKapp')
            if print_output:
                print('Mu Def: {}'.format(mumax_def))
            def_results = copy.deepcopy(rba_session.Results)
            rba_session.clear_results_and_parameters()
        except:
            def_results = {}
        if variability_analysis is not None:
            if len(list(def_results.keys()))!=0:
                rba_session.set_growth_rate(mumax_def*mu_factor_for_variability)
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

        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]["LB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["LB"])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]["UB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["UB"])})
        mumax_prok = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max=max_mu_in_dichotomy,start_value=start_val, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        try:
            rba_session.record_results('Prokaryotic')
            if print_output:
                print('Mu Prok: {}'.format(mumax_prok))
            prok_results = copy.deepcopy(rba_session.Results)
            rba_session.clear_results_and_parameters()
        except:
            prok_results = {}
        if variability_analysis is not None:
            if len(list(prok_results.keys()))!=0:
                rba_session.set_growth_rate(mumax_prok*mu_factor_for_variability)
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

        rba_session.eukaryotic_densities_calibration(CompartmentRelationships=False)
        #rba_session.eukaryotic_densities_pg_fraction(fixed_size_compartments=[],compartment_fraction_prefix="fraction_protein_")
        if Exchanges_to_impose is not None:
            rba_session.Problem.set_lb({exrx: Exchanges_to_impose[exrx]["LB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["LB"])})
            rba_session.Problem.set_ub({exrx: Exchanges_to_impose[exrx]["UB"] for exrx in list(Exchanges_to_impose.keys()) if not pandas.isna(Exchanges_to_impose[exrx]["UB"])})
        mumax_euk = rba_session.find_max_growth_rate(precision=Mu_approx_precision,max=max_mu_in_dichotomy,start_value=start_val, feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        try:
            rba_session.record_results('Eukaryotic')
            if print_output:
                print('Mu Euk: {}'.format(mumax_euk))
            euk_results = copy.deepcopy(rba_session.Results)

            compartment_fractions = {}
            totAA = rba_session.Problem.get_right_hand_side('O_total')['O_total']
            for occ in ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']:
                compartment_fractions[occ] = rba_session.Problem.SolutionValues[str('O_'+occ)]/totAA
            rba_session.clear_results_and_parameters()
        except:
            compartment_fractions={}
            euk_results = {}
        if variability_analysis is not None:
            if len(list(euk_results.keys()))!=0:
                rba_session.set_growth_rate(mumax_euk*mu_factor_for_variability)
                euk_Feasible_Ranges=rba_session.get_feasible_range(variability_analysis)

    #rba_session.model.write(output_dir="Yeast_model_test")
    return({"FeasibleRange_def":def_Feasible_Ranges,"FeasibleRange_prok":prok_Feasible_Ranges,"FeasibleRange_euk":euk_Feasible_Ranges,"Mu_def":mumax_def,"Mu_prok":mumax_prok,"Mu_euk":mumax_euk,'Simulation_Results': prok_results, 'Simulation_Results_Euk': euk_results, 'Simulation_Results_DefKapp': def_results, 'Euk_CompSizes': compartment_fractions,"Condition":condition})


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
                                 feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                 try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                 Mu_scaling_coeff=1,
                                 total_protein_scaling_coeff=1,
                                 print_output=False,
                                 variability_analysis=None,
                                 mu_factor_for_variability=1,
                                 apply_model=False,
                                 functions_to_include_list=[],
                                 transporter_multiplier=3):

    def_Feasible_Ranges={}
    prok_Feasible_Ranges={}
    euk_Feasible_Ranges={}
    def_results={}
    prok_results={}
    euk_results={}
    compartment_fractions={}
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

        rba_session.eukaryotic_densities_calibration(CompartmentRelationships=False)
        #rba_session.eukaryotic_densities_pg_fraction(amino_acid_concentration_total='amino_acid_concentration_total',PG_fraction_prefix="fraction_non_enzymatic_protein_")

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
            rba_session.record_results('Eukaryotic')
            if print_output:
                print('Mu Euk: {}'.format(growth_rate_from_input(input=definition_file, condition=condition)))
            euk_results = copy.deepcopy(rba_session.Results)

            compartment_fractions = {}
            totAA = rba_session.Problem.get_right_hand_side('O_total')['O_total']
            for occ in ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']:
                compartment_fractions[occ] = rba_session.Problem.SolutionValues[str('O_'+occ)]/totAA
            rba_session.clear_results_and_parameters()
            if variability_analysis is not None:
                if len(list(euk_results.keys())):
                    rba_session.set_growth_rate(Mu=mu_factor_for_variability*Mu_scaling_coeff*growth_rate_from_input(input=definition_file, condition=condition))
                    euk_Feasible_Ranges=rba_session.get_feasible_range(variability_analysis)

    return({"FeasibleRange_def":def_Feasible_Ranges,"FeasibleRange_prok":prok_Feasible_Ranges,"FeasibleRange_euk":euk_Feasible_Ranges,'Simulation_Results': prok_results, 'Simulation_Results_Euk': euk_results, 'Simulation_Results_DefKapp': def_results, 'Euk_CompSizes': compartment_fractions,"Condition":condition})


def extract_predicted_exchange_fluxes(inputs=[],result_object='Simulation_Results', run='Prokaryotic', metabolite=''):
    out = []
    for i in inputs:
        try:
            out.append(abs(i[result_object]['ExchangeFluxes'].loc[metabolite, run]))
        except:
            out.append(numpy.nan)
    return(out)


def extract_compartment_sizes_from_calibration_outputs(calibration_outputs):
    out=pandas.DataFrame()
    for i in calibration_outputs:
        for j in i["Densities_PGs"].index:
            out.loc[j,i["Condition"]]=round(i["Densities_PGs"].loc[j,"Density"],6)
    return(out)


def extract_pg_fractions_from_calibration_outputs(calibration_outputs):
    out=pandas.DataFrame()
    for i in calibration_outputs:
        for j in i["Densities_PGs"].index:
            out.loc[j,i["Condition"]]=round(i["Densities_PGs"].loc[j,"PG_fraction"],6)
    return(out)


def extract_specific_kapps_from_calibration_outputs(calibration_outputs):
    out=pandas.DataFrame()
    for i in calibration_outputs:
        for j in i['Specific_Kapps']["Enzyme_ID"]:
            out.loc[j,"Enzyme_ID"]=j
            out.loc[j,i["Condition"]]=i['Specific_Kapps'].loc[i['Specific_Kapps']["Enzyme_ID"]==j,"Kapp"].values[0]
            out.loc[j,i["Condition"]+"_Flux"]=i['Specific_Kapps'].loc[i['Specific_Kapps']["Enzyme_ID"]==j,"Flux"].values[0]
    return(out)


def extract_default_kapps_from_calibration_outputs(calibration_outputs):
    out = pandas.DataFrame()
    for i in calibration_outputs:
        out.loc[i["Condition"], 'ID'] = i["Condition"]
        out.loc[i["Condition"], 'Default Kapp'] = i['Default_Kapps']["default_efficiency"]
    return(out)


def extract_process_capacities_from_calibration_outputs(calibration_outputs):
    out = pandas.DataFrame()
    for i in calibration_outputs:
        for j in list(i['Process_Efficiencies']['Process']):
            val=i['Process_Efficiencies'].loc[i['Process_Efficiencies']['Process']==j,'Value'].values[0]
            out.loc[j,i["Condition"]]=val
    return(out)


def generate_compartment_size_and_pg_input(compartment_sizes,pg_fractions,condition):
    out=pandas.DataFrame()
    for i in compartment_sizes.index:
        out.loc[i,"Compartment_ID"]=i
        out.loc[i,"Density"]=compartment_sizes.loc[i,condition]
    for i in pg_fractions.index:
        out.loc[i,"Compartment_ID"]=i
        out.loc[i,"PG_fraction"]=pg_fractions.loc[i,condition]
    return(out)


def generate_process_efficiency_input(process_efficiencies,condition,parameter_name_suffix):
    out=pandas.DataFrame()
    for i in process_efficiencies.index:
        out.loc[i,"Process"]=i
        out.loc[i,"Parameter"]=i+parameter_name_suffix
        out.loc[i,"Value"]=process_efficiencies.loc[i,condition]
    return(out)


def generate_default_kapp_input(default_kapps,condition,transporter_multiplier):
    out={}
    out["default_efficiency"]=default_kapps.loc[condition,"Default Kapp"]
    out["default_transporter_efficiency"]=transporter_multiplier*default_kapps.loc[condition,"Default Kapp"]
    return(out)


def generate_specific_kapp_input(specific_kapps,condition):
    out=pandas.DataFrame()
    for i in specific_kapps.index:
        if not pandas.isna(specific_kapps.loc[i,condition]):
            out.loc[i,"Enzyme_ID"]=i
            out.loc[i,"Kapp"]=specific_kapps.loc[i,condition]
            out.loc[i,"Flux"]=specific_kapps.loc[i,condition+"_Flux"]
    return(out)


def linear_function(x, a, b):
    return((a*x)+b)

def quadratic_function(x, a, b,c):
    return((a*x**2)+(b*x)+c)

def eval_linear_function(x_in, a, b):
    return([(a*x)+b for x in x_in])

def eval_quadratic_function(x_in, a, b , c):
    return([(a*x**2)+(b*x)+c for x in x_in])

def logistic_function(x, y_max ,x0, k, y_min):
    y = y_max / (1 + numpy.exp(-k*(x-x0)))+y_min
    return (y)

def logistic_function_1st_derivative(x, y_max ,x_0, k):
    return((y_max*k*numpy.exp(-k*(x - x_0)))/((1 + numpy.exp(-k*(x-x_0)))**2))

def logistic_function_2nd_derivative(x, y_max ,x_0, k):
    return(((numpy.exp(k*x_0)-numpy.exp(k*x))*numpy.exp(k*(x+x_0))*y_max*k**2)/((numpy.exp(k*x)+numpy.exp(k*x_0))**3))

def eval_logistic_function(x_in,y_max,x0,k,y_min):
    return([logistic_function(x=x,y_max=y_max,x0=x0,k=k,y_min=y_min) for x in x_in])

def eval_logistic_function_2nd_derivative(x_in,y_max,x0,k):
    return([logistic_function_2nd_derivative(x=x,y_max=y_max,x_0=x0,k=k) for x in x_in])

def eval_logistic_function_1st_derivative(x_in,y_max,x0,k):
    return([logistic_function_1st_derivative(x=x,y_max=y_max,x_0=x0,k=k) for x in x_in])

def do_lin_regression(x_to_fit,y_to_fit,min_val):
    if len(x_to_fit)>2:
        popt_lin, pcov_lin = curve_fit(linear_function, xdata=x_to_fit, ydata=y_to_fit)
        a=popt_lin[0]
        b=popt_lin[1]
    elif len(x_to_fit)==2:
        a=(abs(y_to_fit[1])-abs(y_to_fit[0]))/(abs(x_to_fit[1])-abs(x_to_fit[0]))
        b=y_to_fit[0]-(a*x_to_fit[0])
    elif len(x_to_fit)==1:
        a=0
        b=y_to_fit[0]
    elif len(x_to_fit)==0:
        a=0
        b=min_val
    return({"A":round(a,5),"B":round(b,5)})

def do_quadratic_regression(x_to_fit,y_to_fit,min_val):
    if len(x_to_fit)>2:
        popt_quad, pcov_quad = curve_fit(quadratic_function, xdata=x_to_fit, ydata=y_to_fit)
        a=popt_quad[0]
        b=popt_quad[1]
        c=popt_quad[2]
    elif len(x_to_fit)==2:
        a=0
        b=(abs(y_to_fit[1])-abs(y_to_fit[0]))/(abs(x_to_fit[1])-abs(x_to_fit[0]))
        c=y_to_fit[0]-(a*x_to_fit[0])
    elif len(x_to_fit)==1:
        a=0
        b=0
        c=y_to_fit[0]
    elif len(x_to_fit)==0:
        a=0
        b=0
        c=min_val
    return({"A":round(a,5),"B":round(b,5),"C":round(c,5)})

def do_log_regression(x_to_fit,y_to_fit,x_to_plot,max_val,min_val):
    try:
        p0 = [max(y_to_fit), numpy.median(x_to_fit),1,min(y_to_fit)]
        popt_log, pcov_lin = curve_fit(logistic_function, xdata=x_to_fit, ydata=y_to_fit,p0=p0)
        popt_lin, pcov_lin = curve_fit(linear_function, xdata=x_to_fit, ydata=y_to_fit)

        ymax=popt_log[0]
        ymin=popt_log[3]
        xmin=popt_log[1]
        k=popt_log[2]

        if ymin<(popt_log[0]+popt_log[3]):
            ymax=popt_log[0]+popt_log[3]
            ymin=popt_log[3]
        else:
            ymax=popt_log[3]
            ymin=popt_log[0]+popt_log[3]

        y_pred=eval_logistic_function(x_in=x_to_plot,y_max=popt_log[0],x0=popt_log[1],k=popt_log[2],y_min=popt_log[3])
        y_slope=eval_logistic_function_1st_derivative(x_in=x_to_plot,y_max=popt_log[0],x0=popt_log[1],k=popt_log[2])
        y_slope_2=eval_logistic_function_2nd_derivative(x_in=x_to_plot,y_max=popt_log[0],x0=popt_log[1],k=popt_log[2])

        first_derivative_peak_index=list(scipy.signal.find_peaks(x=y_slope)[0])
        first_derivative_valley_index=list(scipy.signal.find_peaks(x=[-i for i in y_slope])[0])
        second_derivative_peak_index=list(scipy.signal.find_peaks(x=y_slope_2)[0])
        second_derivative_valley_index=list(scipy.signal.find_peaks(x=[-i for i in y_slope_2])[0])

        if len(second_derivative_peak_index)>0:
            if len(second_derivative_valley_index)>0:
                #both inflection points on 1st derivative
                xdiff=abs(x_to_plot[second_derivative_valley_index[0]]-x_to_plot[second_derivative_peak_index[0]])
                ydiff=abs(y_pred[second_derivative_valley_index[0]]-y_pred[second_derivative_peak_index[0]])
                if second_derivative_peak_index[0]<second_derivative_valley_index[0]:
                    #positive shape
                    a_lin=ydiff/xdiff
                    constant_for_linear_function=eval_logistic_function([x_to_plot[second_derivative_peak_index[0]]], *popt_log)[0]-(a_lin*x_to_plot[second_derivative_peak_index[0]])
                    if max(y_to_fit)>y_to_fit[-1]:
                        y_max=min([ymax,max_val])
                    else:
                        y_max=min([max(y_to_fit),max_val])
                    if min(y_to_fit)<y_to_fit[0]:
                        y_min=max([ymin,min_val])
                    else:
                        y_min=max([min(y_to_fit),min_val])
                    params={"Y_max":y_max,"Y_min":y_min,"A":a_lin,"B":constant_for_linear_function}
                else:
                    #negative shape
                    a_lin=-ydiff/xdiff
                    constant_for_linear_function=eval_logistic_function([x_to_plot[second_derivative_peak_index[0]]], *popt_log)[0]-(a_lin*x_to_plot[second_derivative_peak_index[0]])
                    if max(y_to_fit)>y_to_fit[0]:
                        y_max=min([ymax,max_val])
                    else:
                        y_max=min([max(y_to_fit),max_val])
                    if min(y_to_fit)<y_to_fit[-1]:
                        y_min=max([ymin,min_val])
                    else:
                        y_min=max([min(y_to_fit),min_val])

                    params={"Y_max":y_max,"Y_min":y_min,"A":a_lin,"B":constant_for_linear_function}
            else:
                #non-complete positive curve
                if len(first_derivative_peak_index)>0:
                    a_lin=y_slope[first_derivative_peak_index[0]]
                    constant_for_linear_function=eval_logistic_function([x_to_plot[first_derivative_peak_index[0]]], *popt_log)[0]-(a_lin*x_to_plot[first_derivative_peak_index[0]])
                    if max(y_to_fit)>y_to_fit[-1]:
                        y_max=min([ymax,max_val])
                    else:
                        y_max=min([max(y_to_fit),max_val])
                    if min(y_to_fit)<y_to_fit[0]:
                        y_min=max([ymin,min_val])
                    else:
                        y_min=max([min(y_to_fit),min_val])

                    params={"Y_max":y_max,"Y_min":y_min,"A":a_lin,"B":constant_for_linear_function}
                else:
                    a_lin=y_slope[second_derivative_peak_index[0]]
                    constant_for_linear_function=eval_logistic_function([x_to_plot[second_derivative_peak_index[0]]], *popt_log)[0]-(a_lin*x_to_plot[second_derivative_peak_index[0]])
                    if max(y_to_fit)>y_to_fit[-1]:
                        y_max=min([ymax,max_val])
                    else:
                        y_max=min([max(y_to_fit),max_val])
                    if min(y_to_fit)<y_to_fit[0]:
                        y_min=max([ymin,min_val])
                    else:
                        y_min=max([min(y_to_fit),min_val])

                    params={"Y_max":y_max,"Y_min":y_min,"A":a_lin,"B":constant_for_linear_function}
        else:
            if len(second_derivative_valley_index)>0:
                #non-complete negative curve
                if len(first_derivative_valley_index)>0:
                    a_lin=y_slope[first_derivative_valley_index[0]]
                    constant_for_linear_function=eval_logistic_function([x_to_plot[first_derivative_valley_index[0]]], *popt_log)[0]-(a_lin*x_to_plot[first_derivative_valley_index[0]])
                    if max(y_to_fit)>y_to_fit[0]:
                        y_max=min([ymax,max_val])
                    else:
                        y_max=min([max(y_to_fit),max_val])
                    if min(y_to_fit)<y_to_fit[-1]:
                        y_min=max([ymin,min_val])
                    else:
                        y_min=max([min(y_to_fit),min_val])

                    params={"Y_max":y_max,"Y_min":y_min,"A":a_lin,"B":constant_for_linear_function}

                else:
                    a_lin=y_slope[second_derivative_valley_index[0]]
                    constant_for_linear_function=eval_logistic_function([x_to_plot[second_derivative_valley_index[0]]], *popt_log)[0]-(a_lin*x_to_plot[second_derivative_valley_index[0]])
                    if max(y_to_fit)>y_to_fit[0]:
                        y_max=min([ymax,max_val])
                    else:
                        y_max=min([max(y_to_fit),max_val])
                    if min(y_to_fit)<y_to_fit[-1]:
                        y_min=max([ymin,min_val])
                    else:
                        y_min=max([min(y_to_fit),min_val])

                    params={"Y_max":y_max,"Y_min":y_min,"A":a_lin,"B":constant_for_linear_function}
            else:
                if len(first_derivative_valley_index)>0:
                    a_lin=y_slope[first_derivative_valley_index[0]]
                    constant_for_linear_function=eval_logistic_function([x_to_plot[first_derivative_valley_index[0]]], *popt_log)[0]-(a_lin*x_to_plot[first_derivative_valley_index[0]])
                    if max(y_to_fit)>y_to_fit[0]:
                        y_max=min([ymax,max_val])
                    else:
                        y_max=min([max(y_to_fit),max_val])
                    if min(y_to_fit)<y_to_fit[-1]:
                        y_min=max([ymin,min_val])
                    else:
                        y_min=max([min(y_to_fit),min_val])

                    params={"Y_max":y_max,"Y_min":y_min,"A":a_lin,"B":constant_for_linear_function}

                else:
                    try:
                        popt_lin, pcov_lin = curve_fit(linear_function, xdata=x_to_fit, ydata=y_to_fit)
                        y_max=min([max(y_to_fit),max_val])
                        y_min=max([min(y_to_fit),min_val])
                        params={"Y_max":y_max,"Y_min":y_min,"A":popt_lin[0],"B":popt_lin[1]}

                    except:
                        y_max=min([max(y_to_fit),max_val])
                        y_min=max([min(y_to_fit),min_val])
                        params={"Y_max":y_max,"Y_min":y_min,"A":0,"B":numpy.mean(y_to_fit)}
    except:
        try:
            lin_pars=do_lin_regression(x_to_fit=x_to_fit,y_to_fit=y_to_fit,min_val=min_val)
            y_max=min([max(y_to_fit),max_val])
            y_min=max([min(y_to_fit),min_val])
            params={"Y_max":y_max,"Y_min":y_min,"A":lin_pars["A"],"B":lin_pars["B"]}
        except:
            y_max=min([max(y_to_fit),max_val])
            y_min=max([min(y_to_fit),min_val])
            params={"Y_max":y_max,"Y_min":y_min,"A":0,"B":numpy.mean(y_to_fit)}
    return({i:round(params[i],5)for i in params.keys()})

def quad_predictions(params,x_to_fit):
    x_to_use=[]
    for x in x_to_fit:
        prelim_x=x
        if params["X_min"] is not None:
            if x < params["X_min"]:
                prelim_x=params["X_min"]
        if params["X_max"] is not None:
            if x > params["X_max"]:
                prelim_x=params["X_max"]
        x_to_use.append(prelim_x)
    quad_prediction_values=eval_quadratic_function(x_in=x_to_use, a=params["A"], b=params["B"] , c=params["C"])
    y_pred=[]
    for i in range(len(x_to_fit)):
        if quad_prediction_values[i]<params["Y_min"]:
            y_pred.append(params["Y_min"])
        elif quad_prediction_values[i]>params["Y_max"]:
            y_pred.append(params["Y_max"])
        else:
            y_pred.append(quad_prediction_values[i])
    return(y_pred)

def lin_predictions(params,x_to_fit):
    x_to_use=[]
    for x in x_to_fit:
        prelim_x=x
        if params["X_min"] is not None:
            if x < params["X_min"]:
                prelim_x=params["X_min"]
        if params["X_max"] is not None:
            if x > params["X_max"]:
                prelim_x=params["X_max"]
        x_to_use.append(prelim_x)
    lin_prediction_values=eval_linear_function(x_in=x_to_use, a=params["A"], b=params["B"])
    y_pred=[]
    for i in range(len(x_to_fit)):
        if lin_prediction_values[i]<params["Y_min"]:
            y_pred.append(params["Y_min"])
        elif lin_prediction_values[i]>params["Y_max"]:
            y_pred.append(params["Y_max"])
        else:
            y_pred.append(lin_prediction_values[i])
    return(y_pred)

def calculate_rss(y_predicted,y_measured):
    if len(y_predicted)==len(y_measured):
        RSS=0
        for i in range(len(y_measured)):
            RSS+=(y_measured[i]-y_predicted[i])**2
        return(RSS)
    else:
        return(None)

def do_regression(x_to_fit,y_to_fit,x_to_plot,max_val,min_val,monotonous_quadratic=False,total_x_range=(0,1)):
    #y_max=min([max(y_to_fit),max_val])
    #y_min=max([min(y_to_fit),min_val])
    y_max=max_val
    y_min=min_val
    x_max=total_x_range[1]
    x_min=total_x_range[0]
    log_regression_results=do_log_regression(x_to_fit=x_to_fit,y_to_fit=y_to_fit,x_to_plot=x_to_plot,max_val=max_val,min_val=min_val)
    log_regression_results.update({"X_max":None,"X_min":None})
    try:
        popt_lin, pcov_lin = curve_fit(linear_function, xdata=x_to_fit, ydata=y_to_fit)
        lin_regression_results={"X_max":None,"X_min":None,"Y_max":round(y_max,5),"Y_min":round(y_min,5),"A":round(popt_lin[0],5),"B":round(popt_lin[1],5)}
    except:
        lin_regression_results={"X_max":None,"X_min":None,"Y_max":round(y_max,5),"Y_min":round(y_min,5),"A":0,"B":round(numpy.mean(y_to_fit),5)}
    try:
        popt_quad, pcov_quad = curve_fit(quadratic_function, xdata=x_to_fit, ydata=y_to_fit)
        quad_regression_results={"X_max":None,"X_min":None,"Y_max":round(y_max,5),"Y_min":round(y_min,5),"A":round(popt_quad[0],5),"B":round(popt_quad[1],5),"C":round(popt_quad[2],5)}
        if quad_regression_results["A"]!=0:
            if monotonous_quadratic:
                extremum_x= -0.5*quad_regression_results["B"]/quad_regression_results["A"]
                if extremum_x > x_min:
                    if extremum_x < x_max:
                        test_params_extremum_is_xmax=quad_regression_results.copy()
                        test_params_extremum_is_xmax["X_max"]=extremum_x
                        test_params_extremum_is_xmin=quad_regression_results.copy()
                        test_params_extremum_is_xmin["X_min"]=extremum_x
                        try:
                            predictions_extremum_is_xmax=quad_predictions(params=test_params_extremum_is_xmax,x_to_fit=x_to_fit)
                            RSS_extremum_is_xmax=calculate_rss(y_predicted=predictions_extremum_is_xmax,y_measured=y_to_fit)
                        except:
                            RSS_extremum_is_xmax=None
                        try:
                            predictions_extremum_is_xmin=quad_predictions(params=test_params_extremum_is_xmin,x_to_fit=x_to_fit)
                            RSS_extremum_is_xmin=calculate_rss(y_predicted=predictions_extremum_is_xmin,y_measured=y_to_fit)
                        except:
                            RSS_extremum_is_xmin=None
                        if RSS_extremum_is_xmax is not None:
                            if RSS_extremum_is_xmin is not None:
                                if RSS_extremum_is_xmin > RSS_extremum_is_xmax:
                                    quad_regression_results=test_params_extremum_is_xmax
                                else:
                                    quad_regression_results=test_params_extremum_is_xmin
                            else:
                                quad_regression_results=test_params_extremum_is_xmax
                        else:
                            if RSS_extremum_is_xmin is not None:
                                quad_regression_results=test_params_extremum_is_xmin
        else:
            quad_regression_results=None
    except:
        quad_regression_results=None

    predictions_on_log_model=lin_predictions(params=log_regression_results,x_to_fit=x_to_fit)
    RSS_log=calculate_rss(y_predicted=predictions_on_log_model,y_measured=y_to_fit)

    predictions_on_lin_model=lin_predictions(params=lin_regression_results,x_to_fit=x_to_fit)
    RSS_lin=calculate_rss(y_predicted=predictions_on_lin_model,y_measured=y_to_fit)

    if quad_regression_results is not None:
        predictions_on_quad_model=quad_predictions(params=quad_regression_results,x_to_fit=x_to_fit)
        RSS_quad=calculate_rss(y_predicted=predictions_on_quad_model,y_measured=y_to_fit)

    out={}
    out["LinParams"]=lin_regression_results
    out["LogParams"]=log_regression_results
    out["QuadParams"]=None
    if RSS_log<=RSS_lin:
        if quad_regression_results is not None:
            out["QuadParams"]=quad_regression_results
            if RSS_log<=RSS_quad:
                out["Type"]="Log"
                out["Parameters"]=log_regression_results
            else:
                out["Type"]="Quad"
                out["Parameters"]=quad_regression_results
        else:
            out["Type"]="Log"
            out["Parameters"]=log_regression_results
    else:
        if quad_regression_results is not None:
            out["QuadParams"]=quad_regression_results
            if RSS_lin<=RSS_quad:
                out["Type"]="Lin"
                out["Parameters"]=lin_regression_results
            else:
                out["Type"]="Quad"
                out["Parameters"]=quad_regression_results
        else:
            out["Type"]="Lin"
            out["Parameters"]=lin_regression_results
    return(out)

def regression_on_specific_enzyme_efficiencies(Spec_Kapps,min_kapp,max_kapp,conditions,growth_rates,only_lin=False,impose_on_isoenzymes=True,monotonous_quadratic=False,fill_in_missing_conditions=False):
    out=pandas.DataFrame(columns=Spec_Kapps.columns)
    for enzyme in Spec_Kapps.index:
        out.loc[enzyme,list(Spec_Kapps.columns)]=Spec_Kapps.loc[enzyme,list(Spec_Kapps.columns)]
        for condition in conditions:
            out.loc[enzyme,condition]=None
        x=[growth_rates[i] for i in conditions]
        y=list(Spec_Kapps.loc[enzyme,conditions])
        x_s=[]
        y_s=[]
        for i in range(len(x)):
            if not pandas.isna(y[i]):
                if y[i]>=max_kapp:
                    y_s.append(max_kapp)
                    x_s.append(x[i])
                elif y[i]>=min_kapp:
                    y_s.append(y[i])
                    x_s.append(x[i])
        if len(x_s)>1:
            x_to_plot=(numpy.linspace(x_s[0], x_s[-1], 1000))
            regression_results=do_regression(x_to_fit=x_s,y_to_fit=y_s,x_to_plot=x_to_plot,max_val=max_kapp,min_val=min_kapp,monotonous_quadratic=monotonous_quadratic,total_x_range=(min(x),max(x)))
            if regression_results["Type"]=="Quad":
                if fill_in_missing_conditions:
                    for condition in conditions:
                        out.loc[enzyme,condition]=quad_predictions(params=regression_results["Parameters"],x_to_fit=[growth_rates[condition]])
                else:
                    for condition in conditions:
                        if Spec_Kapps.loc[enzyme,condition] is not None:
                            out.loc[enzyme,condition]=quad_predictions(params=regression_results["Parameters"],x_to_fit=[growth_rates[condition]])
            else:
                if fill_in_missing_conditions:
                    for condition in conditions:
                        out.loc[enzyme,condition]=lin_predictions(params=regression_results["Parameters"],x_to_fit=[growth_rates[condition]])
                else:
                    for condition in conditions:
                        if Spec_Kapps.loc[enzyme,condition] is not None:
                            out.loc[enzyme,condition]=lin_predictions(params=regression_results["Parameters"],x_to_fit=[growth_rates[condition]])
        elif len(x_s)==1:
            regression_results={"Parameters":{"Y_max":y_s[0],"Y_min":y_s[0],"A":0,"B":y_s[0]}}
            if fill_in_missing_conditions:
                for condition in conditions:
                    out.loc[enzyme,condition]=y_s[0]
            else:
                for condition in conditions:
                    if Spec_Kapps.loc[enzyme,condition] is not None:
                        out.loc[enzyme,condition]=y_s[0]
        elif len(x_s)==0:
            regression_results={"Parameters":{"Y_max":min_kapp,"Y_min":min_kapp,"A":0,"B":min_kapp}}
            for condition in conditions:
                out.loc[enzyme,condition]=min_kapp
        if regression_results["Parameters"]["A"]!=0:
            if "C" in regression_results["Parameters"].keys():
                out.loc[enzyme,"Model"]=json.dumps({"quadratic":regression_results["Parameters"]})
            else:
                out.loc[enzyme,"Model"]=json.dumps({"linear":regression_results["Parameters"]})
        else:
            out.loc[enzyme,"Model"]=json.dumps({"constant":{'CONSTANT':regression_results["Parameters"]["B"]}})
    return(out)

def regression_on_compartment_sizes(Comp_sizes,conditions,growth_rates,monotonous_quadratic=False):
    pp = PdfPages('CompSize_Plots_refined3.pdf')
    out=pandas.DataFrame()
    for comp in Comp_sizes.index:
        x=[growth_rates[i] for i in conditions]
        y=list(Comp_sizes.loc[comp,conditions])
        x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
        regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=1,min_val=0,monotonous_quadratic=monotonous_quadratic,total_x_range=(0,1))
        if regression_results["Type"]=="Quad":
            color="orange"
            predictions_for_plot=quad_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
            predictions_for_file=quad_predictions(params=regression_results["Parameters"],x_to_fit=x)
        else:
            if regression_results["Type"]=="Log":
                color="red"
            elif regression_results["Type"]=="Lin":
                color="green"
            predictions_for_plot=lin_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
            predictions_for_file=lin_predictions(params=regression_results["Parameters"],x_to_fit=x)
        fig = plt.figure(figsize=(10,5))
        plt.scatter(x,y)
        plt.plot(x_to_plot,predictions_for_plot,color=color)
        #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_min"]]*len(x_to_plot),color="black")
        #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_max"]]*len(x_to_plot),color="black")
        #plt.ylim(ymin=0*min_kapp)
        plt.title(comp)
        plt.xlabel(json.dumps(regression_results["Parameters"]))
        pp.savefig(fig)
        plt.close()
        out.loc[comp,conditions]=predictions_for_file
        if regression_results["Parameters"]["A"]!=0:
            if "C" in regression_results["Parameters"].keys():
                out.loc[comp,"Model"]=json.dumps({"quadratic":regression_results["Parameters"]})
            else:
                out.loc[comp,"Model"]=json.dumps({"linear":regression_results["Parameters"]})
        else:
            out.loc[comp,"Model"]=json.dumps({"constant":{'CONSTANT':regression_results["Parameters"]["B"]}})
    pp.close()
    return(out)

def regression_on_pg_fractions(PG_sizes,conditions,growth_rates,monotonous_quadratic=False):
    pp = PdfPages('PGfraction_Plots_refined3.pdf')
    out=pandas.DataFrame()
    for comp in PG_sizes.index:
        x=[growth_rates[i]for i in conditions]
        y=list(PG_sizes.loc[comp,conditions])

        x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
        regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=1,min_val=0,monotonous_quadratic=monotonous_quadratic,total_x_range=(0,1))
        if regression_results["Type"]=="Quad":
            color="orange"
            predictions_for_plot=quad_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
            predictions_for_file=quad_predictions(params=regression_results["Parameters"],x_to_fit=x)
        else:
            if regression_results["Type"]=="Log":
                color="red"
            elif regression_results["Type"]=="Lin":
                color="green"
            predictions_for_plot=lin_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
            predictions_for_file=lin_predictions(params=regression_results["Parameters"],x_to_fit=x)
        fig = plt.figure(figsize=(10,5))
        plt.scatter(x,y)
        plt.plot(x_to_plot,predictions_for_plot,color=color)
        #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_min"]]*len(x_to_plot),color="black")
        #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_max"]]*len(x_to_plot),color="black")
        #plt.ylim(ymin=0*min_kapp)
        plt.title(comp)
        plt.xlabel(json.dumps(regression_results["Parameters"]))
        pp.savefig(fig)
        plt.close()
        out.loc[comp,conditions]=predictions_for_file
        if regression_results["Parameters"]["A"]!=0:
            if "C" in regression_results["Parameters"].keys():
                out.loc[comp,"Model"]=json.dumps({"quadratic":regression_results["Parameters"]})
            else:
                out.loc[comp,"Model"]=json.dumps({"linear":regression_results["Parameters"]})
        else:
            out.loc[comp,"Model"]=json.dumps({"constant":{'CONSTANT':regression_results["Parameters"]["B"]}})
    pp.close()
    return(out)

def regression_on_process_efficiencies(Process_efficiencies,min_efficiency,max_efficiency,conditions,growth_rates,monotonous_quadratic=False):
    pp = PdfPages('ProcessEfficiencies_Plots_refined3.pdf')
    out=pandas.DataFrame()
    for process in Process_efficiencies.index:
        x=[growth_rates[i]for i in conditions]
        y=list(Process_efficiencies.loc[process,conditions])

        x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
        regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=max_efficiency,min_val=min_efficiency,monotonous_quadratic=monotonous_quadratic,total_x_range=(min(x),max(x)))
        if regression_results["Type"]=="Quad":
            color="orange"
            predictions_for_plot=quad_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
            predictions_for_file=quad_predictions(params=regression_results["Parameters"],x_to_fit=x)
        else:
            if regression_results["Type"]=="Log":
                color="red"
            elif regression_results["Type"]=="Lin":
                color="green"
            predictions_for_plot=lin_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
            predictions_for_file=lin_predictions(params=regression_results["Parameters"],x_to_fit=x)
        fig = plt.figure(figsize=(10,5))
        plt.scatter(x,y)
        plt.plot(x_to_plot,predictions_for_plot,color=color)
        #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_min"]]*len(x_to_plot),color="black")
        #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_max"]]*len(x_to_plot),color="black")
        #plt.ylim(ymin=0*min_kapp)
        plt.title(process)
        plt.xlabel(json.dumps(regression_results["Parameters"]))
        pp.savefig(fig)
        plt.close()
        out.loc[process,conditions]=predictions_for_file
        if regression_results["Parameters"]["A"]!=0:
            if "C" in regression_results["Parameters"].keys():
                out.loc[process,"Model"]=json.dumps({"quadratic":regression_results["Parameters"]})
            else:
                out.loc[process,"Model"]=json.dumps({"linear":regression_results["Parameters"]})
        else:
            out.loc[process,"Model"]=json.dumps({"constant":{'CONSTANT':regression_results["Parameters"]["B"]}})
    pp.close()
    return(out)

def regression_on_default_enzyme_efficiencies(default_kapps,min_kapp,max_kapp,conditions,growth_rates,monotonous_quadratic=False):
    pp = PdfPages('DefKapp_Plots_refined3.pdf')
    out=pandas.DataFrame()
    x=[growth_rates[i]for i in conditions]
    y=list(default_kapps.loc[conditions,"Default Kapp"])

    x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
    regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=max_kapp,min_val=min_kapp,monotonous_quadratic=monotonous_quadratic,total_x_range=(min(x),max(x)))
    if regression_results["Type"]=="Quad":
        color="orange"
        predictions_for_plot=quad_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
        predictions_for_file=quad_predictions(params=regression_results["Parameters"],x_to_fit=x)
    else:
        if regression_results["Type"]=="Log":
            color="red"
        elif regression_results["Type"]=="Lin":
            color="green"
        predictions_for_plot=lin_predictions(params=regression_results["Parameters"],x_to_fit=x_to_plot)
        predictions_for_file=lin_predictions(params=regression_results["Parameters"],x_to_fit=x)
    fig = plt.figure(figsize=(10,5))
    plt.scatter(x,y)
    plt.plot(x_to_plot,predictions_for_plot,color=color)
    #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_min"]]*len(x_to_plot),color="black")
    #plt.plot(x_to_plot,[regression_results["Parameters"]["Y_max"]]*len(x_to_plot),color="black")
    plt.title("Default Kapp")
    plt.xlabel(json.dumps(regression_results["Parameters"]))
    pp.savefig(fig)
    plt.close()
    out.loc[:,"Default Kapp"]=predictions_for_file
    for cond in conditions:
        #out.loc[cond,"Default Kapp"]=list(lin_predictions(params=regression_results["Parameters"],x_to_fit=[growth_rates[cond]]))[0]
        out.loc[cond,"ID"]= cond
        out.loc[cond,"Mu"]= growth_rates[cond]
    if regression_results["Parameters"]["A"]!=0:
        if "C" in regression_results["Parameters"].keys():
            out.loc["Regression","Model"]=json.dumps({"quadratic":regression_results["Parameters"]})
        else:
            out.loc["Regression","Model"]=json.dumps({"linear":regression_results["Parameters"]})
    else:
        out.loc["Regression","Model"]=json.dumps({"constant":{'CONSTANT':regression_results["Parameters"]["B"]}})
    pp.close()
    return(out)

def plot_predicted_fluxes(simulation_outputs,types=["DefaultKapp","Prokaryotic","Eukaryotic"],input_definition=None):

    ########
    Mus_o2 = [0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.28, 0.3, 0.35, 0.4]
    Glc_J = [0.3, 0.6, 1.1, 1.7, 2.3, 2.8, 3.4, 4.5, 8.6, 11.1]
    EtOH_J = [0, 0, 0, 0, 0, 0, 0.11, 2.3, 9.5, 13.9]
    O2_J = [0.8,  1.3, 2.5, 3.9, 5.3, 7, 7.4, 6.1, 5.1, 3.7]
    CO2_J = [0.8,  1.4, 2.7, 4.2, 5.7, 7.5, 8, 8.8, 14.9, 18.9]
    Ac_J = [0, 0, 0, 0, 0, 0, 0.08, 0.41, 0.62, 0.6]
    Glyc_J = [0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.15]

    Mu_Hackett_list=[]
    Glc_Hackett_list=[]
    EtOH_Hackett_list=[]
    Ac_Hackett_list=[]
    Acald_Hackett_list=[]
    Glyc_Hackett_list=[]
    Lac_Hackett_list=[]
    Succ_Hackett_list=[]

    Glc_Hackett_SE_list=[]
    EtOH_Hackett_SE_list=[]
    Ac_Hackett_SE_list=[]
    Acald_Hackett_SE_list=[]
    Glyc_Hackett_SE_list=[]
    Lac_Hackett_SE_list=[]
    Succ_Hackett_SE_list=[]
    if input_definition is not None:
        for result in simulation_outputs:
            condition=result["Condition"]
            try:
                Mu_Hackett_list.append(growth_rate_from_input(input=input_definition,condition=condition))
            except:
                Mu_Hackett_list.append(numpy.nan)
            try:
                Glc_Hackett_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_Mean')&(input_definition['ID'] == 'R_EX_glc__D_e'), condition].values[0])
            except:
                Glc_Hackett_list.append(numpy.nan)
            try:
                EtOH_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_Mean')&(input_definition['ID'] == 'R_EX_etoh_e'), condition].values[0])
            except:
                EtOH_Hackett_list.append(numpy.nan)
            try:
                Ac_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_Mean')&(input_definition['ID'] == 'R_EX_ac_e'), condition].values[0])
            except:
                Ac_Hackett_list.append(numpy.nan)
            try:
                Acald_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_Mean')&(input_definition['ID'] == 'R_EX_acald_e'), condition].values[0])
            except:
                Acald_Hackett_list.append(numpy.nan)
            try:
                Glyc_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_Mean')&(input_definition['ID'] == 'R_EX_glyc_e'), condition].values[0])
            except:
                Glyc_Hackett_list.append(numpy.nan)
            try:
                Lac_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_Mean')&(input_definition['ID'] == 'R_EX_lac__D_e'), condition].values[0])
            except:
                Lac_Hackett_list.append(numpy.nan)
            try:
                Succ_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_Mean')&(input_definition['ID'] == 'R_EX_succ_e'), condition].values[0])
            except:
                Succ_Hackett_list.append(numpy.nan)

            try:
                Glc_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_StandardError')&(input_definition['ID'] == 'R_EX_glc__D_e_SE'), condition].values[0])
            except:
                Glc_Hackett_SE_list.append(numpy.nan)
            try:
                EtOH_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_StandardError')&(input_definition['ID'] == 'R_EX_etoh_e_SE'), condition].values[0])
            except:
                EtOH_Hackett_SE_list.append(numpy.nan)
            try:
                Ac_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_StandardError')&(input_definition['ID'] == 'R_EX_ac_e_SE'), condition].values[0])
            except:
                Ac_Hackett_SE_list.append(numpy.nan)
            try:
                Acald_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_StandardError')&(input_definition['ID'] == 'R_EX_acald_e_SE'), condition].values[0])
            except:
                Acald_Hackett_SE_list.append(numpy.nan)
            try:
                Glyc_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_StandardError')&(input_definition['ID'] == 'R_EX_glyc_e_SE'), condition].values[0])
            except:
                Glyc_Hackett_SE_list.append(numpy.nan)
            try:
                Lac_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_StandardError')&(input_definition['ID'] == 'R_EX_lac__D_e_SE'), condition].values[0])
            except:
                Lac_Hackett_SE_list.append(numpy.nan)
            try:
                Succ_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'ExchangeFlux_StandardError')&(input_definition['ID'] == 'R_EX_succ_e_SE'), condition].values[0])
            except:
                Succ_Hackett_SE_list.append(numpy.nan)

        Glc_Hackett_SE=numpy.array(Glc_Hackett_SE_list)
        EtOH_Hackett_SE=numpy.array(EtOH_Hackett_SE_list)
        Ac_Hackett_SE=numpy.array(Ac_Hackett_SE_list)
        Acald_Hackett_SE=numpy.array(Acald_Hackett_SE_list)
        Glyc_Hackett_SE=numpy.array(Glyc_Hackett_SE_list)
        Lac_Hackett_SE=numpy.array(Lac_Hackett_SE_list)
        Succ_Hackett_SE=numpy.array(Succ_Hackett_SE_list)

        Mu_Hackett=numpy.array(Mu_Hackett_list)
        Glc_Hackett=numpy.array(Glc_Hackett_list)
        EtOH_Hackett=numpy.array(EtOH_Hackett_list)
        Ac_Hackett=numpy.array(Ac_Hackett_list)
        Acald_Hackett=numpy.array(Acald_Hackett_list)
        Glyc_Hackett=numpy.array(Glyc_Hackett_list)
        Lac_Hackett=numpy.array(Lac_Hackett_list)
        Succ_Hackett=numpy.array(Succ_Hackett_list)

    ## Hackett#
    #Mu_Hackett = numpy.array([0.0498630244, 0.1054314572, 0.154377453333333, 0.2126503108, 0.293841410333333])
    #Glc_Hackett = numpy.array([0.7367, 1.5462, 2.1722, 5.1571, 9.5962])
    #EtOH_Hackett = numpy.array([0.0127, 0.0529, 0.1084, 4.6066, 14.0672])
    #Ac_Hackett = numpy.array([0.0017, 0.0031, 0.0052, 0.4433, 0.8851])
    #Acald_Hackett = numpy.array([0.0073,0.0167,0.0231,0.0798,0.3282])
    #Glyc_Hackett = numpy.array([0.0035, 0.0077, 0.0065, 0.0579, 0.1699])
    #Lac_Hackett  = numpy.array([0.0063,0.0145,0.0222,0.0572,0.2041])
    #Succ_Hackett = numpy.array([0.0004,0.001,0.0048,0.0193,0.0448])

    #Glc_Hackett_SE = numpy.array([0.0068,0.0257,0.0343,0.0902,1.7616])
    #EtOH_Hackett_SE = numpy.array([0.0052,0.0136,0.0094,0.3675,1.247])
    #Ac_Hackett_SE = numpy.array([0.0001,0.0005,0.0002,0.0223,0.0657])
    #Acald_Hackett_SE = numpy.array([0.0004,0.0011,0.0014,0.0332,0.0633])
    #Glyc_Hackett_SE = numpy.array([0.0001,0.0007,0.0014,0.0165,0.0202])
    #Lac_Hackett_SE  = numpy.array([0.0003,0.0015,0.0016,0.0098,0.0169])
    #Succ_Hackett_SE = numpy.array([0.0001,0.0001,0.0001,0.0006,0.0017])

    Mus_predicted_def = extract_predicted_growth_rates(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp')
    Mus_predicted = extract_predicted_growth_rates(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic')
    Mus_predicted_euk = extract_predicted_growth_rates(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic')

    Glc_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_glc__D')
    EtOH_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_etoh')
    O2_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_o2')
    CO2_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_co2')
    Ac_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_ac')
    Glycerol_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_glyc')
    Acald_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_acald')
    Lac_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_lac__D')
    Succ_Exchange_predicted_def = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_DefKapp', run='DefaultKapp', metabolite='M_succ')

    Glc_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_glc__D')
    EtOH_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_etoh')
    O2_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_o2')
    CO2_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_co2')
    Ac_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_ac')
    Glycerol_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_glyc')
    Acald_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_acald')
    Lac_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_lac__D')
    Succ_Exchange_predicted = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results', run='Prokaryotic', metabolite='M_succ')

    Glc_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_glc__D')
    EtOH_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_etoh')
    O2_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_o2')
    CO2_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_co2')
    Ac_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_ac')
    Glycerol_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_glyc')
    Acald_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_acald')
    Lac_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_lac__D')
    Succ_Exchange_predicted_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk', run='Eukaryotic', metabolite='M_succ')


    Glc_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_glc__D_e',bound_type="Min")
    EtOH_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_etoh_e',bound_type="Min")
    O2_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_o2_e',bound_type="Min")
    Ac_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_ac_e',bound_type="Min")
    Glycerol_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_glyc_e',bound_type="Min")
    Acald_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_acald_e',bound_type="Min")
    Lac_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_lac__D_e',bound_type="Min")
    Succ_VAmin_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_succ_e',bound_type="Min")
    Glc_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_glc__D_e',bound_type="Max")
    EtOH_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_etoh_e',bound_type="Max")
    O2_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_o2_e',bound_type="Max")
    Ac_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_ac_e',bound_type="Max")
    Glycerol_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_glyc_e',bound_type="Max")
    Acald_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_acald_e',bound_type="Max")
    Lac_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_lac__D_e',bound_type="Max")
    Succ_VAmax_def=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_def', variable='R_EX_succ_e',bound_type="Max")

    Glc_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_glc__D_e',bound_type="Min")
    EtOH_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_etoh_e',bound_type="Min")
    O2_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_o2_e',bound_type="Min")
    Ac_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_ac_e',bound_type="Min")
    Glycerol_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_glyc_e',bound_type="Min")
    Acald_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_acald_e',bound_type="Min")
    Lac_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_lac__D_e',bound_type="Min")
    Succ_VAmin_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_succ_e',bound_type="Min")
    Glc_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_glc__D_e',bound_type="Max")
    EtOH_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_etoh_e',bound_type="Max")
    O2_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_o2_e',bound_type="Max")
    Ac_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_ac_e',bound_type="Max")
    Glycerol_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_glyc_e',bound_type="Max")
    Acald_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_acald_e',bound_type="Max")
    Lac_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_lac__D_e',bound_type="Max")
    Succ_VAmax_prok=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_prok', variable='R_EX_succ_e',bound_type="Max")

    Glc_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_glc__D_e',bound_type="Min")
    EtOH_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_etoh_e',bound_type="Min")
    O2_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_o2_e',bound_type="Min")
    Ac_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_ac_e',bound_type="Min")
    Glycerol_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_glyc_e',bound_type="Min")
    Acald_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_acald_e',bound_type="Min")
    Lac_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_lac__D_e',bound_type="Min")
    Succ_VAmin_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_succ_e',bound_type="Min")
    Glc_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_glc__D_e',bound_type="Max")
    EtOH_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_etoh_e',bound_type="Max")
    O2_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_o2_e',bound_type="Max")
    Ac_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_ac_e',bound_type="Max")
    Glycerol_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_glyc_e',bound_type="Max")
    Acald_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_acald_e',bound_type="Max")
    Lac_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_lac__D_e',bound_type="Max")
    Succ_VAmax_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk', variable='R_EX_succ_e',bound_type="Max")

    ###
    fig, axs = plt.subplots(3, 3, figsize=(18, 11), sharex=True)
    axs[0, 0].plot(Mu_Hackett, Mu_Hackett, color='lightgreen')
    legendlist=['Hackett']
    if "DefaultKapp" in types:
        axs[0, 0].scatter(Mu_Hackett, Mus_predicted_def, color='orange')
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[0, 0].scatter(Mu_Hackett, Mus_predicted, color='black')
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[0, 0].scatter(Mu_Hackett, Mus_predicted_euk, color='red')
        legendlist.append("Euk.")
    axs[0, 0].legend(legendlist)
    axs[0, 0].set_title('Predicted vs measured growth-rate')
    axs[0, 0].set_ylabel('$\mu$ [$h^{-1}$]')
    axs[0, 0].set_xlabel('$\mu$ [$h^{-1}$]')

    axs[0, 1].plot(Mus_o2, Glc_J, color='lightblue')
    axs[0, 1].plot(Mu_Hackett, Glc_Hackett, color='lightgreen')
    axs[0, 1].fill_between(Mu_Hackett, Glc_Hackett-Glc_Hackett_SE, Glc_Hackett+Glc_Hackett_SE, color='lightgreen',alpha=0.7)
    legendlist=['van Hoek', 'Hackett','SE - Hackett']
    if "DefaultKapp" in types:
        axs[0, 1].scatter(Mus_predicted_def, Glc_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[0, 1].scatter(Mus_predicted, Glc_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[0, 1].scatter(Mus_predicted_euk, Glc_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[0, 1].legend(legendlist)
    #axs[0, 1].scatter(Mus_predicted_def, Glc_VAmin_def, color='orange',marker=6, alpha=0.5)
    #axs[0, 1].scatter(Mus_predicted_def, Glc_VAmax_def, color='orange',marker=7, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[0, 1].plot([Mus_predicted[i]]*2,[Glc_VAmin_prok[i],Glc_VAmax_prok[i]],color="black",alpha=0.1)
        axs[0, 1].scatter(Mus_predicted, Glc_VAmin_prok, color='black',marker=6, alpha=0.5)
        axs[0, 1].scatter(Mus_predicted, Glc_VAmax_prok, color='black',marker=7, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[0, 1].plot([Mus_predicted_euk[i]]*2,[Glc_VAmin_euk[i],Glc_VAmax_euk[i]],color="red",alpha=0.1)
        axs[0, 1].scatter(Mus_predicted_euk, Glc_VAmin_euk, color='red',marker=6, alpha=0.5)
        axs[0, 1].scatter(Mus_predicted_euk, Glc_VAmax_euk, color='red',marker=7, alpha=0.5)
    axs[0, 1].set_title('Glucose-uptake rate')
    axs[0, 1].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[0, 1].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    axs[0, 2].plot(Mus_o2, O2_J, color='lightblue')
    legendlist=['van Hoek']
    if "DefaultKapp" in types:
        axs[0, 2].scatter(Mus_predicted_def, O2_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[0, 2].scatter(Mus_predicted, O2_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[0, 2].scatter(Mus_predicted_euk, O2_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[0, 2].legend(legendlist)
    #axs[0, 2].scatter(Mus_predicted_def, O2_VAmin_def, color='orange',marker=6, alpha=0.5)
    #axs[0, 2].scatter(Mus_predicted_def, O2_VAmax_def, color='orange',marker=7, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[0, 2].plot([Mus_predicted[i]]*2,[O2_VAmin_prok[i],O2_VAmax_prok[i]],color="black",alpha=0.1)
        axs[0, 2].scatter(Mus_predicted, O2_VAmin_prok, color='black',marker=6, alpha=0.5)
        axs[0, 2].scatter(Mus_predicted, O2_VAmax_prok, color='black',marker=7, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[0, 2].plot([Mus_predicted_euk[i]]*2,[O2_VAmin_euk[i],O2_VAmax_euk[i]],color="red",alpha=0.1)
        axs[0, 2].scatter(Mus_predicted_euk, O2_VAmin_euk, color='red',marker=6, alpha=0.5)
        axs[0, 2].scatter(Mus_predicted_euk, O2_VAmax_euk, color='red',marker=7, alpha=0.5)
    axs[0, 2].set_title('Oxygen-uptake rate')
    axs[0, 2].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[0, 2].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    axs[1, 0].plot(Mus_o2, EtOH_J, color='lightblue')
    axs[1, 0].plot(Mu_Hackett, EtOH_Hackett, color='lightgreen')
    axs[1, 0].fill_between(Mu_Hackett, EtOH_Hackett-EtOH_Hackett_SE, EtOH_Hackett+EtOH_Hackett_SE, color='lightgreen',alpha=0.7)
    legendlist=['van Hoek', 'Hackett','SE - Hackett']
    if "DefaultKapp" in types:
        axs[1, 0].scatter(Mus_predicted_def, EtOH_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[1, 0].scatter(Mus_predicted, EtOH_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[1, 0].scatter(Mus_predicted_euk, EtOH_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[1, 0].legend(legendlist)
    #axs[1, 0].scatter(Mus_predicted_def, EtOH_VAmin_def, color='orange',marker=7, alpha=0.5)
    #axs[1, 0].scatter(Mus_predicted_def, EtOH_VAmax_def, color='orange',marker=6, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[1, 0].plot([Mus_predicted[i]]*2,[EtOH_VAmin_prok[i],EtOH_VAmax_prok[i]],color="black",alpha=0.1)
        axs[1, 0].scatter(Mus_predicted, EtOH_VAmin_prok, color='black',marker=7, alpha=0.5)
        axs[1, 0].scatter(Mus_predicted, EtOH_VAmax_prok, color='black',marker=6, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[1,0].plot([Mus_predicted_euk[i]]*2,[EtOH_VAmin_euk[i],EtOH_VAmax_euk[i]],color="red",alpha=0.1)
        axs[1, 0].scatter(Mus_predicted_euk, EtOH_VAmin_euk, color='red',marker=7, alpha=0.5)
        axs[1, 0].scatter(Mus_predicted_euk, EtOH_VAmax_euk, color='red',marker=6, alpha=0.5)
    axs[1, 0].set_title('Ethanol-excretion rate')
    axs[1, 0].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[1, 0].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    axs[1, 1].plot(Mus_o2, Ac_J, color='lightblue')
    axs[1, 1].plot(Mu_Hackett, Ac_Hackett, color='lightgreen')
    axs[1, 1].fill_between(Mu_Hackett, Ac_Hackett-Ac_Hackett_SE, Ac_Hackett+Ac_Hackett_SE, color='lightgreen',alpha=0.7)
    legendlist=['van Hoek', 'Hackett','SE - Hackett']
    if "DefaultKapp" in types:
        axs[1, 1].scatter(Mus_predicted_def, Ac_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[1, 1].scatter(Mus_predicted, Ac_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[1, 1].scatter(Mus_predicted_euk, Ac_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[1, 1].legend(legendlist)
    #axs[1, 1].scatter(Mus_predicted_def, Ac_VAmin_def, color='orange',marker=7, alpha=0.5)
    #axs[1, 1].scatter(Mus_predicted_def, Ac_VAmax_def, color='orange',marker=6, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[1, 1].plot([Mus_predicted[i]]*2,[Ac_VAmin_prok[i],Ac_VAmax_prok[i]],color="black",alpha=0.1)
        axs[1, 1].scatter(Mus_predicted, Ac_VAmin_prok, color='black',marker=7, alpha=0.5)
        axs[1, 1].scatter(Mus_predicted, Ac_VAmax_prok, color='black',marker=6, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[1, 1].plot([Mus_predicted_euk[i]]*2,[Ac_VAmin_euk[i],Ac_VAmax_euk[i]],color="red",alpha=0.1)
        axs[1, 1].scatter(Mus_predicted_euk, Ac_VAmin_euk, color='red',marker=7, alpha=0.5)
        axs[1, 1].scatter(Mus_predicted_euk, Ac_VAmax_euk, color='red',marker=6, alpha=0.5)
    axs[1, 1].set_title('Acetate-excretion rate')
    axs[1, 1].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[1, 1].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    axs[1, 2].plot(Mus_o2, Glyc_J, color='lightblue')
    axs[1, 2].plot(Mu_Hackett, Glyc_Hackett, color='lightgreen')
    axs[1, 2].fill_between(Mu_Hackett, Glyc_Hackett-Glyc_Hackett_SE, Glyc_Hackett+Glyc_Hackett_SE, color='lightgreen',alpha=0.7)
    legendlist=['van Hoek', 'Hackett','SE - Hackett']
    if "DefaultKapp" in types:
        axs[1, 2].scatter(Mus_predicted_def, Glycerol_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[1, 2].scatter(Mus_predicted, Glycerol_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[1, 2].scatter(Mus_predicted_euk, Glycerol_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[1, 2].legend(legendlist)
    #axs[1, 2].scatter(Mus_predicted_def, Glycerol_VAmin_def, color='orange',marker=7, alpha=0.5)
    #axs[1, 2].scatter(Mus_predicted_def, Glycerol_VAmax_def, color='orange',marker=6, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[1, 2].plot([Mus_predicted[i]]*2,[Glycerol_VAmin_prok[i],Glycerol_VAmax_prok[i]],color="black",alpha=0.1)
        axs[1, 2].scatter(Mus_predicted, Glycerol_VAmin_prok, color='black',marker=7, alpha=0.5)
        axs[1, 2].scatter(Mus_predicted, Glycerol_VAmax_prok, color='black',marker=6, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[1, 2].plot([Mus_predicted_euk[i]]*2,[Glycerol_VAmin_euk[i],Glycerol_VAmax_euk[i]],color="red",alpha=0.1)
        axs[1, 2].scatter(Mus_predicted_euk, Glycerol_VAmin_euk, color='red',marker=7, alpha=0.5)
        axs[1, 2].scatter(Mus_predicted_euk, Glycerol_VAmax_euk, color='red',marker=6, alpha=0.5)
    axs[1, 2].set_title('Glycerol-excretion rate')
    axs[1, 2].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[1, 2].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    axs[2, 0].plot(Mu_Hackett, Acald_Hackett, color='lightgreen')
    axs[2, 0].fill_between(Mu_Hackett, Acald_Hackett-Acald_Hackett_SE, Acald_Hackett+Acald_Hackett_SE, color='lightgreen',alpha=0.7)
    legendlist=['Hackett','SE - Hackett']
    if "DefaultKapp" in types:
        axs[2, 0].scatter(Mus_predicted_def, Acald_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[2, 0].scatter(Mus_predicted, Acald_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[2, 0].scatter(Mus_predicted_euk, Acald_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[2, 0].legend(legendlist)
    #axs[2, 0].scatter(Mus_predicted_def, Acald_VAmin_def, color='orange',marker=7, alpha=0.5)
    #axs[2, 0].scatter(Mus_predicted_def, Acald_VAmax_def, color='orange',marker=6, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[2, 0].plot([Mus_predicted[i]]*2,[Acald_VAmin_prok[i],Acald_VAmax_prok[i]],color="black",alpha=0.1)
        axs[2, 0].scatter(Mus_predicted, Acald_VAmin_prok, color='black',marker=7, alpha=0.5)
        axs[2, 0].scatter(Mus_predicted, Acald_VAmax_prok, color='black',marker=6, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[2, 0].plot([Mus_predicted_euk[i]]*2,[Acald_VAmin_euk[i],Acald_VAmax_euk[i]],color="red",alpha=0.1)
        axs[2, 0].scatter(Mus_predicted_euk, Acald_VAmin_euk, color='red',marker=7, alpha=0.5)
        axs[2, 0].scatter(Mus_predicted_euk, Acald_VAmax_euk, color='red',marker=6, alpha=0.5)
    axs[2, 0].set_title('Acetaldehyde-excretion rate')
    axs[2, 0].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[2, 0].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    axs[2, 1].plot(Mu_Hackett, Lac_Hackett, color='lightgreen')
    axs[2, 1].fill_between(Mu_Hackett, Lac_Hackett-Lac_Hackett_SE, Lac_Hackett+Lac_Hackett_SE, color='lightgreen',alpha=0.7)
    legendlist=['Hackett','SE - Hackett']
    if "DefaultKapp" in types:
        axs[2, 1].scatter(Mus_predicted_def, Lac_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[2, 1].scatter(Mus_predicted, Lac_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[2, 1].scatter(Mus_predicted_euk, Lac_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[2, 1].legend(legendlist)
    #axs[2, 1].scatter(Mus_predicted_def, Lac_VAmin_def, color='orange',marker=7, alpha=0.5)
    #axs[2, 1].scatter(Mus_predicted_def, Lac_VAmax_def, color='orange',marker=6, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[2,1].plot([Mus_predicted[i]]*2,[Lac_VAmin_prok[i],Lac_VAmax_prok[i]],color="black",alpha=0.1)
        axs[2, 1].scatter(Mus_predicted, Lac_VAmin_prok, color='black',marker=7, alpha=0.5)
        axs[2, 1].scatter(Mus_predicted, Lac_VAmax_prok, color='black',marker=6, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[2, 1].plot([Mus_predicted_euk[i]]*2,[Lac_VAmin_euk[i],Lac_VAmax_euk[i]],color="red",alpha=0.1)
        axs[2, 1].scatter(Mus_predicted_euk, Lac_VAmin_euk, color='red',marker=7, alpha=0.5)
        axs[2, 1].scatter(Mus_predicted_euk, Lac_VAmax_euk, color='red',marker=6, alpha=0.5)
    axs[2, 1].set_title('Lactate-excretion rate')
    axs[2, 1].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[2, 1].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    axs[2, 2].plot(Mu_Hackett, Succ_Hackett, color='lightgreen')
    axs[2, 2].fill_between(Mu_Hackett, Succ_Hackett-Succ_Hackett_SE, Succ_Hackett+Succ_Hackett_SE, color='lightgreen',alpha=0.7)
    legendlist=['Hackett','SE - Hackett']
    if "DefaultKapp" in types:
        axs[2, 2].scatter(Mus_predicted_def, Succ_Exchange_predicted_def, color='orange', alpha=0.8)
        legendlist.append("Def.")
    if "Prokaryotic" in types:
        axs[2, 2].scatter(Mus_predicted, Succ_Exchange_predicted, color='black', alpha=0.8)
        legendlist.append("Prok.")
    if "Eukaryotic" in types:
        axs[2, 2].scatter(Mus_predicted_euk, Succ_Exchange_predicted_euk, color='red', alpha=0.8)
        legendlist.append("Euk.")
    axs[2, 2].legend(legendlist)
    #axs[2, 2].scatter(Mus_predicted_def, Succ_VAmin_def, color='orange',marker=7, alpha=0.5)
    #axs[2, 2].scatter(Mus_predicted_def, Succ_VAmax_def, color='orange',marker=6, alpha=0.5)
    if "Prokaryotic" in types:
        for i in range(len(Mus_predicted)):
            axs[2, 2].plot([Mus_predicted[i]]*2,[Succ_VAmin_prok[i],Succ_VAmax_prok[i]],color="black",alpha=0.1)
        axs[2, 2].scatter(Mus_predicted, Succ_VAmin_prok, color='black',marker=7, alpha=0.5)
        axs[2, 2].scatter(Mus_predicted, Succ_VAmax_prok, color='black',marker=6, alpha=0.5)
    if "Eukaryotic" in types:
        for i in range(len(Mus_predicted_euk)):
            axs[2, 2].plot([Mus_predicted_euk[i]]*2,[Succ_VAmin_euk[i],Succ_VAmax_euk[i]],color="red",alpha=0.1)
        axs[2, 2].scatter(Mus_predicted_euk, Succ_VAmin_euk, color='red',marker=7, alpha=0.5)
        axs[2, 2].scatter(Mus_predicted_euk, Succ_VAmax_euk, color='red',marker=6, alpha=0.5)
    axs[2, 2].set_title('Succinate-excretion rate')
    axs[2, 2].set_xlabel('$\mu$ [$h^{-1}$]')
    axs[2, 2].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

    plt.show()

def plot_specific_enzyme_efficiencies(point_calibration_kapps,regressed_kapps,conditions,growth_rates,filename,min_y=0,max_y=10**12):
    out_pdf=PdfPages(filename)
    for enzyme in point_calibration_kapps.index:
        if enzyme in regressed_kapps.index:
            fig=plt.figure(figsize=(6,4))
            plt.scatter([growth_rates[i] for i in conditions],list(point_calibration_kapps.loc[enzyme,conditions]))
            x_to_plot=(numpy.linspace([growth_rates[i] for i in conditions][0], [growth_rates[i] for i in conditions][-1], 1000))
            if list(json.loads(regressed_kapps.loc[enzyme,"Model"]).keys())[0]=="linear":
                plt.scatter([growth_rates[i] for i in conditions],list(regressed_kapps.loc[enzyme,conditions]),color="red",marker="x")
                plt.plot(x_to_plot,lin_predictions(params=json.loads(regressed_kapps.loc[enzyme,"Model"])["linear"],x_to_fit=x_to_plot),color="red")
            elif list(json.loads(regressed_kapps.loc[enzyme,"Model"]).keys())[0]=="constant":
                plt.scatter([growth_rates[i] for i in conditions],list(regressed_kapps.loc[enzyme,conditions]),color="green",marker="x")
                plt.plot(x_to_plot,[json.loads(regressed_kapps.loc[enzyme,"Model"])["constant"]["CONSTANT"]]*len(list(x_to_plot)),color="green")
            elif list(json.loads(regressed_kapps.loc[enzyme,"Model"]).keys())[0]=="quadratic":
                plt.scatter([growth_rates[i] for i in conditions],list(regressed_kapps.loc[enzyme,conditions]),color="orange",marker="x")
                plt.plot(x_to_plot,quad_predictions(params=json.loads(regressed_kapps.loc[enzyme,"Model"])["quadratic"],x_to_fit=x_to_plot),color="orange")
            plt.title(enzyme)
            #plt.ylim(0.9*min_y,1.1*min([max_y,max(list(point_calibration_kapps.loc[enzyme,conditions])),max(list(regressed_kapps.loc[enzyme,conditions]))]))
            plt.xlabel("Growth rate $h^{-1}$")
            plt.ylabel("Kapp $h^{-1}$")
            out_pdf.savefig(fig)
            plt.close()
    out_pdf.close()

def plot_compartment_sizes_and_pg(point_calibration_sizes,point_calibration_pg,regressed_sizes,regressed_pg,conditions,growth_rates,filename):
    out_pdf=PdfPages(filename)
    for compartment in point_calibration_sizes.index:
        fig, ax = plt.subplots(1, 2, figsize=(9, 4))
        ax[0].scatter([growth_rates[i] for i in conditions],list(point_calibration_sizes.loc[compartment,conditions]))
        if compartment in regressed_sizes.index:
            ax[0].plot([growth_rates[i] for i in conditions],list(regressed_sizes.loc[compartment,conditions]),color="red")
        ax[0].set_title("Size: {}".format(compartment))
        ax[0].set_xlabel("Growth rate $h^{-1}$")
        ax[0].set_ylabel("Fraction")
        if compartment in point_calibration_pg.index:
            ax[1].scatter([growth_rates[i] for i in conditions],list(point_calibration_pg.loc[compartment,conditions]))
            if compartment in regressed_pg.index:
                ax[1].plot([growth_rates[i] for i in conditions],list(regressed_pg.loc[compartment,conditions]),color="red")
            ax[1].set_title("PG: {}".format(compartment))
            ax[1].set_xlabel("Growth rate $h^{-1}$")
            ax[1].set_ylabel("Fraction")
        out_pdf.savefig(fig)
        plt.close()
    out_pdf.close()

def plot_rss_trajectory(calibration_outputs):
    out=pandas.DataFrame()
    for cal_input in calibration_outputs:
        for i in range(len(cal_input["RSS_trajectory"])):
            out.loc[i,cal_input["Condition"]]=cal_input["RSS_trajectory"][i]
    out.to_csv("RRS_trajectory.csv")
    n_figs=len(calibration_outputs)
    if n_figs<=3:
        n_rows=1
        n_cols=3
        ax_indices=[(0,0),(0,1),(0,2)]
    elif n_figs<=6:
        n_rows=2
        n_cols=3
        ax_indices=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
    elif n_figs<=9:
        n_rows=3
        n_cols=3
        ax_indices=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(28, 11), sharex=True)
    count=0
    for cal_input in calibration_outputs:
        axs[ax_indices[count][0],ax_indices[count][1]].plot([i+1 for i in range(len(cal_input["RSS_trajectory"]))],cal_input["RSS_trajectory"])
        axs[ax_indices[count][0],ax_indices[count][1]].scatter([cal_input["RSS_trajectory"].index(min(cal_input["RSS_trajectory"]))+1],min(cal_input["RSS_trajectory"]),color="red")
        axs[ax_indices[count][0],ax_indices[count][1]].set_title("RSS: {}".format(cal_input["Condition"]))
        axs[ax_indices[count][0],ax_indices[count][1]].set_xlabel("Iteration")
        axs[ax_indices[count][0],ax_indices[count][1]].set_ylabel("Residual Sum of squares")
        count+=1
    plt.show()

def perform_protein_protein_comparison(rba_session,condition,calibration_object,simulation_object,simulation_type,scaling_coeff=6.023e20):
    ProtoProteinMap=generate_proto_protein_map(rba_session=rba_session)
    out = pandas.DataFrame()
    for i in list(set(list(simulation_object[simulation_type]['ProtoProteins'].index)+list(calibration_object['Proteome']['ID']))):
        if i in list(simulation_object[simulation_type]['ProtoProteins'].index):
            if simulation_object[simulation_type]['ProtoProteins'].loc[i].values[0] !=0:
                out.loc[i, 'ID'] = i
                if i in list(simulation_object[simulation_type]['ProtoProteins'].index):
                    out.loc[i, 'Predicted'] = scaling_coeff * simulation_object[simulation_type]['ProtoProteins'].loc[i].values[0]
                if i in list(calibration_object['Proteome']['ID']):
                    out.loc[i, 'Measured'] = scaling_coeff * calibration_object['Proteome'].loc[calibration_object['Proteome']['ID'] == i, 'copy_number'].values[0]

    out = determine_kapp_type_classification(protein_DF=out,rbaSession=rba_session,specific_kapps=calibration_object['Specific_Kapps'],protomap=ProtoProteinMap)
    return(out)

def plot_predicted_versus_measured_proteomes(rba_session,calibration_outputs,simulation_outputs,type):
    if type=="Prokaryotic":
        results_object="Simulation_Results"
    elif type=="Eukaryotic":
        results_object="Simulation_Results_Euk"
    elif type=="DefaultKapp":
        results_object="Simulation_Results_DefKapp"

def efficiency_correction(specific_kapps,
                               simulation_results,
                               proteinData,
                               rba_session,
                               condition,
                               also_correct_default,
                               default_kapps,
                               tolerance=2,
                               previous_misprediction_factors={},
                               process_efficiencies=None,
                               correct_process_efficiencies=False,
                               correct_default_kapp_enzymes=False):
    protoProtein_results={}
    for i in simulation_results["Proteins"].index:
        protoID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if protoID in protoProtein_results:
            protoProtein_results[protoID]+=simulation_results["Proteins"].loc[i,condition]
        else:
            protoProtein_results[protoID]=simulation_results["Proteins"].loc[i,condition]

    Model_ProtoProtein_Map={}
    for i in rba_session.get_proteins():
        protoID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if protoID in Model_ProtoProtein_Map:
            Model_ProtoProtein_Map[protoID].append(i)
        else:
            Model_ProtoProtein_Map[protoID]=[i]
    kapps_out=specific_kapps.copy()
    process_efficiencies_out=process_efficiencies.copy()
    misprediction_factors_enzymes={}
    misprediction_factors_defaultKapp_enzymes={}
    misprediction_factors_processes={}
    squared_residuals=[]
    for protoID in protoProtein_results:
        if protoID in Model_ProtoProtein_Map:
            if protoID in proteinData["ID"]:
                predicted_protein=protoProtein_results[protoID]
                measured_protein=proteinData.loc[protoID,"copy_number"]
                if (predicted_protein>0) & (measured_protein>0):
                    #misprediction_coeff=predicted_protein/measured_protein
                    misprediction_coeff=numpy.sqrt(predicted_protein/measured_protein)
                    squared_residuals.append((numpy.log(predicted_protein)-numpy.log(measured_protein))**2)
                    #squared_residuals.append((predicted_protein-measured_protein)**2)
                    for protein in Model_ProtoProtein_Map[protoID]:
                        for enzyme in rba_session.get_protein_information(protein=protein)["associatedEnzymes"]:
                            if enzyme in list(specific_kapps["Enzyme_ID"]):
                                if enzyme not in misprediction_factors_enzymes.keys():
                                    misprediction_factors_enzymes[enzyme]=[misprediction_coeff]
                                else:
                                    misprediction_factors_enzymes[enzyme].append(misprediction_coeff)
                            else:
                                if enzyme not in misprediction_factors_defaultKapp_enzymes.keys():
                                    misprediction_factors_defaultKapp_enzymes[enzyme]=[misprediction_coeff]
                                else:
                                    misprediction_factors_defaultKapp_enzymes[enzyme].append(misprediction_coeff)
                        if correct_process_efficiencies:
                            for process in rba_session.get_protein_information(protein=protein)["SupportsProcess"]:
                                if process in process_efficiencies.index:
                                    if process not in misprediction_factors_processes.keys():
                                        misprediction_factors_processes[process]=[misprediction_coeff]
                                    else:
                                        misprediction_factors_processes[process].append(misprediction_coeff)

    enzyme_correction_coefficients={}
    no_change_in_enzyme_efficiencies=True
    for enzyme in misprediction_factors_enzymes.keys():
        if enzyme in list(kapps_out["Enzyme_ID"]):
            old_kapp=specific_kapps.loc[specific_kapps["Enzyme_ID"]==enzyme,"Kapp"].values[0]
            if len(misprediction_factors_enzymes[enzyme])==1:
                correction_coeff=misprediction_factors_enzymes[enzyme][0]
            else:
                correction_coeff=numpy.median(misprediction_factors_enzymes[enzyme])
            enzyme_correction_coefficients[enzyme]=correction_coeff
            if previous_misprediction_factors is not None:
                if enzyme in previous_misprediction_factors.keys():
                    if correction_coeff<1:
                        coeffs_below_1=[i for i in previous_misprediction_factors[enzyme] if i < 1]
                        if len(coeffs_below_1)>=1:
                            if correction_coeff <= max(coeffs_below_1):
                                previous_misprediction_factors[enzyme].append(correction_coeff)
                                continue
                    if correction_coeff>1:
                        coeffs_above_1=[i for i in previous_misprediction_factors[enzyme] if i > 1]
                        if len(coeffs_above_1)>=1:
                            if  correction_coeff >= min(coeffs_above_1):
                                previous_misprediction_factors[enzyme].append(correction_coeff)
                                continue
                    else:
                        previous_misprediction_factors[enzyme].append(correction_coeff)
                else:
                    previous_misprediction_factors[enzyme]=[correction_coeff]
        if tolerance is None:
            new_kapp=old_kapp*correction_coeff
            kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Kapp"]=new_kapp
            no_change_in_enzyme_efficiencies=False
        else:
            if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                new_kapp=old_kapp*correction_coeff
                kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Kapp"]=new_kapp
                no_change_in_enzyme_efficiencies=False

    if correct_default_kapp_enzymes:
        for enzyme in misprediction_factors_defaultKapp_enzymes.keys():
            if enzyme in list(kapps_out["Enzyme_ID"]):
                old_kapp=specific_kapps.loc[specific_kapps["Enzyme_ID"]==enzyme,"Kapp"].values[0]
                if len(misprediction_factors_defaultKapp_enzymes[enzyme])==1:
                    correction_coeff=misprediction_factors_defaultKapp_enzymes[enzyme][0]
                else:
                    correction_coeff=numpy.median(misprediction_factors_defaultKapp_enzymes[enzyme])
                enzyme_correction_coefficients[enzyme]=correction_coeff
                if previous_misprediction_factors is not None:
                    if enzyme in previous_misprediction_factors.keys():
                        if correction_coeff<1:
                            coeffs_below_1=[i for i in previous_misprediction_factors[enzyme] if i < 1]
                            if len(coeffs_below_1)>=1:
                                if correction_coeff <= max(coeffs_below_1):
                                    previous_misprediction_factors[enzyme].append(correction_coeff)
                                    continue
                        if correction_coeff>1:
                            coeffs_above_1=[i for i in previous_misprediction_factors[enzyme] if i > 1]
                            if len(coeffs_above_1)>=1:
                                if  correction_coeff >= min(coeffs_above_1):
                                    previous_misprediction_factors[enzyme].append(correction_coeff)
                                    continue
                        else:
                            previous_misprediction_factors[enzyme].append(correction_coeff)
                    else:
                        previous_misprediction_factors[enzyme]=[correction_coeff]
                if tolerance is None:
                    new_kapp=old_kapp*correction_coeff
                    kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Kapp"]=new_kapp
                    no_change_in_enzyme_efficiencies=False
                else:
                    if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                        new_kapp=old_kapp*correction_coeff
                        kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Kapp"]=new_kapp
                        no_change_in_enzyme_efficiencies=False
            else:
                rxn=rba_session.get_enzyme_information(enzyme=enzyme)["Reaction"].split("_duplicate_")[0]
                if rxn in simulation_results["Reactions"].index:
                    flux=simulation_results["Reactions"].loc[rxn,condition]
                    if flux < 0:
                        direction=-1
                        parameter= rba_session.model.enzymes.enzymes._elements_by_id[enzyme].backward_efficiency
                    elif flux > 0:
                        direction=1
                        parameter= rba_session.model.enzymes.enzymes._elements_by_id[enzyme].forward_efficiency
                    else:
                        continue
                    old_kapp=default_kapps[parameter]
                    if len(misprediction_factors_defaultKapp_enzymes[enzyme])==1:
                        correction_coeff=misprediction_factors_defaultKapp_enzymes[enzyme][0]
                    else:
                        correction_coeff=numpy.median(misprediction_factors_defaultKapp_enzymes[enzyme])
                    enzyme_correction_coefficients[enzyme]=correction_coeff
                    if previous_misprediction_factors is not None:
                        if enzyme in previous_misprediction_factors.keys():
                            if correction_coeff<1:
                                coeffs_below_1=[i for i in previous_misprediction_factors[enzyme] if i < 1]
                                if len(coeffs_below_1)>=1:
                                    if correction_coeff <= max(coeffs_below_1):
                                        previous_misprediction_factors[enzyme].append(correction_coeff)
                                        continue
                            if correction_coeff>1:
                                coeffs_above_1=[i for i in previous_misprediction_factors[enzyme] if i > 1]
                                if len(coeffs_above_1)>=1:
                                    if  correction_coeff >= min(coeffs_above_1):
                                        previous_misprediction_factors[enzyme].append(correction_coeff)
                                        continue
                            else:
                                previous_misprediction_factors[enzyme].append(correction_coeff)
                        else:
                            previous_misprediction_factors[enzyme]=[correction_coeff]
                    if tolerance is None:
                        new_kapp=old_kapp*correction_coeff
                        kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Kapp"]=new_kapp
                        kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Flux"]=direction
                        kapps_out.loc[kapps_out["Comment"]==enzyme,"Flux"]="Corrected Default"
                        no_change_in_enzyme_efficiencies=False
                    else:
                        if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                            new_kapp=old_kapp*correction_coeff
                            kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Kapp"]=new_kapp
                            kapps_out.loc[kapps_out["Enzyme_ID"]==enzyme,"Flux"]=direction
                            kapps_out.loc[kapps_out["Comment"]==enzyme,"Flux"]="Corrected Default"
                            no_change_in_enzyme_efficiencies=False

    no_change_in_process_efficiencies=True
    if correct_process_efficiencies:
        process_correction_coefficients={}
        for process in misprediction_factors_processes.keys():
            if process in list(process_efficiencies.index):
                old_efficiency=process_efficiencies.loc[process,"Value"]
                if len(misprediction_factors_processes[process])==1:
                    correction_coeff=misprediction_factors_processes[process][0]
                else:
                    correction_coeff=numpy.median(misprediction_factors_processes[process])
            process_correction_coefficients[process]=correction_coeff
            if tolerance is None:
                new_efficiency=old_efficiency*correction_coeff
                process_efficiencies_out.loc[process,"Value"]=new_efficiency
                no_change_in_process_efficiencies=False
            else:
                if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                    new_efficiency=old_efficiency*correction_coeff
                    process_efficiencies_out.loc[process,"Value"]=new_efficiency
                    no_change_in_process_efficiencies=False
    changes_applied=True
    if no_change_in_process_efficiencies:
        if no_change_in_enzyme_efficiencies:
            changes_applied=False
    return({"Sum_of_squared_residuals":sum(squared_residuals),"Changes":changes_applied,"Kapps":kapps_out,"ProcessEfficiencies":process_efficiencies_out,"Process_MispredictionFactors":process_correction_coefficients,"Enzyme_MispredictionFactors":enzyme_correction_coefficients,"Misprediction_factor_trajectory":previous_misprediction_factors})

def extract_proteomes_from_calibration_results(calib_results):
    proteomes=pandas.DataFrame()
    for calib_result in calib_results:
        condition=calib_result["Condition"]
        proteome=calib_result["Proteome"]
        for i in list(proteome.index):
            proteomes.loc[i,condition]=proteome.loc[i,"copy_number"]
    return(proteomes)

def extract_compsizes_and_pgfractions_from_correction_summary(corrsummary,rows_to_exclude):
    out=pandas.DataFrame()
    for i in list(corrsummary.index):
        if i in rows_to_exclude:
            continue
        else:
            out.loc[i,"Compartment_ID"]=i
            out.loc[i,"Density"]=corrsummary.loc[i,"new_protein_fraction"]
            out.loc[i,"PG_fraction"]=corrsummary.loc[i,"new_PG_fraction"]
    return(out)

def extract_feasible_bounds(inputs=[],feasible_range_object='FeasibleRange_prok', variable='',bound_type="Min"):
    out = []
    for i in inputs:
        try:
            val=abs(i[feasible_range_object][variable][bound_type])
            if val!=0:
                out.append(val)
            else:
                out.append(0)
        except:
            out.append(0)
    return(out)

def calibration_workflow(proteome,
                         condition,
                         reference_condition,
                         gene_ID_column,
                         definition_file,
                         rba_session,
                         process_efficiency_estimation_input=None,
                         spec_kapps=None,
                         corrected_spec_kapps=False,
                         take_lowest_RSS_parameters=False,
                         pre_optimise_def_kapp=True,
                         process_efficiencies=None,
                         Compartment_sizes=None,
                         PG_fractions=None,
                         transporter_multiplier=3,
                         prelim_run=False,
                         final_global_scaling_free_exchanges=False,
                         Mu_approx_precision=0.00001,
                         min_kapp=None,
                         fixed_mu_when_above_target_mu_in_correction=True,
                         mu_misprediction_tolerance=0.05,
                         print_outputs=True):
    t0 = time.time()
    correction_results = correction_pipeline(input=proteome,
                                             condition=condition,
                                             definition_file=definition_file,
                                             reference_condition=reference_condition,
                                             compartments_to_replace={'DEF':"c", 'DEFA':"c", 'Def':"c"},
                                             compartments_no_original_PG=['n', 'Secreted'],
                                             fractions_entirely_replaced_with_expected_value=['Ribosomes'],
                                             imposed_compartment_fractions=proteome_fractions_from_input(input=definition_file, condition=condition),
                                             directly_corrected_compartments=['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'],
                                             merged_compartments={'c': 'Ribosomes'},
                                             min_compartment_fraction=0.00000)
    rba_session.set_medium(medium_concentrations_from_input(input=definition_file, condition=condition))
    if prelim_run:
        compartment_densities_and_PGs = extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=correction_results,rows_to_exclude=["Ribosomes","Total"]+[i for i in correction_results.index if i.startswith("pg_")])
        correction_results.to_csv(str('Correction_overview_HackettNielsen_'+condition+'.csv'))
        return({"Densities_PGs":compartment_densities_and_PGs,
                "Condition":condition})
    if Compartment_sizes is not None:
        for i in Compartment_sizes.index:
            correction_results.loc[i,"new_protein_fraction"]=Compartment_sizes.loc[i,condition]
            if i in PG_fractions.index:
                correction_results.loc[i,"new_PG_fraction"]=PG_fractions.loc[i,condition]

    for i in correction_results.index:
        abundance_coeff=1
        if i =="c":
            abundance_coeff=(correction_results.loc[i,"new_protein_fraction"]*(1-correction_results.loc[i,"new_PG_fraction"])-correction_results.loc["Ribosomes","new_protein_fraction"])/(correction_results.loc[i,"original_protein_fraction"]*(1-correction_results.loc[i,"original_PG_fraction"]))
        else:
            abundance_coeff=(correction_results.loc[i,"new_protein_fraction"]*(1-correction_results.loc[i,"new_PG_fraction"]))/(correction_results.loc[i,"original_protein_fraction"]*(1-correction_results.loc[i,"original_PG_fraction"]))
        proteome.loc[proteome["Location"]==i,condition]*=abundance_coeff
        correction_results.loc[i,"copy_number_scaling"]=abundance_coeff

    correction_results.to_csv(str('Correction_overview_HackettNielsen_corrected_'+condition+'.csv'))
    if process_efficiencies is None:
        if process_efficiency_estimation_input is not None:
            process_efficiencies = determine_apparent_process_efficiencies(growth_rate=growth_rate_from_input(input=definition_file,
                                                                           condition=condition),
                                                                           input=process_efficiency_estimation_input,
                                                                           rba_session=rba_session,
                                                                           protein_data=proteome.copy(),
                                                                           proteome_summary=correction_results.copy(),
                                                                           condition=condition,
                                                                           gene_id_col=gene_ID_column,
                                                                           fit_nucleotide_assembly_machinery=True)
    proteome[condition]*=1000/6.022e23
    process_efficiencies.to_csv("ProcEffsOrig_{}.csv".format(condition))
    if spec_kapps is None:
        mets_to_ignore=['M_octa','M_acACP','M_hdca','M_palmACP','M_hdcoa','M_ocACP','M_ocdcea','M_ttdca','M_pmtcoa','M_ddca','M_nadp','M_tdcoa','M_hdeACP','M_tdeACP','M_occoa','M_dca','M_tdecoa','M_ddcacoa','M_ttdcea','M_ocdcyaACP','M_myrsACP','M_hdcea','M_ocdca','M_nadph','M_ocdcaACP','M_octeACP','M_dcaACP','M_ddcaACP','M_ocdcya','M_nadh','M_dcacoa','M_nad','M_stcoa']
        Specific_Kapps_Results = estimate_specific_enzyme_efficiencies(rba_session=rba_session,
                                                                        proteomicsData=build_input_proteome_for_specific_kapp_estimation(proteome, condition),
                                                                        flux_bounds=flux_bounds_from_input(input=definition_file, condition=condition, specific_exchanges=None, specific_directions=None),
                                                                        mu=growth_rate_from_input(input=definition_file, condition=condition),
                                                                        biomass_function=None,
                                                                        target_biomass_function=True,
                                                                        parsimonious_fba=True,
                                                                        only_non_ambigous_proteins=False,
                                                                        chose_isoreaction=True,
                                                                        equalize_identical_enzymes=True,
                                                                        only_identical_reactions_with_twin_enzyme=False,
                                                                        equalize_identical_reactions=False,
                                                                        metabolites_to_ignore=mets_to_ignore,
                                                                        impose_on_all_isoreactions=False,
                                                                        zero_on_all_isoreactions=True,
                                                                        condition=condition,
                                                                        store_output=False,
                                                                        rxns_to_ignore_when_parsimonious=[]
                                                                        )
        #Specific_Kapps=Specific_Kapps_Results["Overview"].loc[Specific_Kapps_Results["Overview"]["Kapp"]>=360,:]
        Specific_Kapps=Specific_Kapps_Results["Overview"]
        if min_kapp is not None:
            Specific_Kapps.loc[(Specific_Kapps['Kapp']<min_kapp)&(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp']=min_kapp
        Spec_Kapp_estim_FD=Specific_Kapps_Results["Flux_Distribution"]
        Specific_Kapps.to_csv("Specific_Kapps_Hackett__{}.csv".format(condition), sep=";", decimal=",")
    else:
        print('importing spec kapps')
        Specific_Kapps = pandas.DataFrame()
        Specific_Kapps['Enzyme_ID'] = spec_kapps['ID']
        Specific_Kapps['Kapp'] = spec_kapps[condition]
        Specific_Kapps['Flux'] = spec_kapps[str(condition+'_Flux')]

    #Specific_Kapps['Kapp']*=100
    spec_kapp_median=Specific_Kapps.loc[(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp'].median()
    #spec_kapp_median=numpy.median(list(set(list(Specific_Kapps.loc[(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp']))))
    Default_Kapps={"default_efficiency":spec_kapp_median,"default_transporter_efficiency":transporter_multiplier*spec_kapp_median}
    ##
    zero_flux_rxns=[]
    ##

    flux_bounds_data=flux_bounds_from_input(input=definition_file, condition=condition, specific_exchanges=None, specific_directions=[])
    Exchanges_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in list(flux_bounds_data["Reaction_ID"])}
    compartment_densities_and_PGs = extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=correction_results,rows_to_exclude=["Ribosomes","Total"]+[i for i in correction_results.index if i.startswith("pg_")])
    mu_measured=growth_rate_from_input(input=definition_file, condition=condition)

    Results_to_look_up="Simulation_Results"
    #Results_to_look_up="Simulation_Results_Euk"
    condition_to_look_up="Prokaryotic"
    #condition_to_look_up="Eukaryotic"
    Growth_rate_to_look_up="Mu_prok"
    #Growth_rate_to_look_up="Mu_euk"

    if pre_optimise_def_kapp:
        print("start def")
        #####
        for zero_rxn in zero_flux_rxns:
            Exchanges_to_impose[zero_rxn]={"LB":0,"UB":0}
        #####
        Simulation_results = perform_simulations(condition=condition,
                                                 rba_session=rba_session,
                                                 definition_file=definition_file,
                                                 compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                 pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                 process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies}]),
                                                 Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps}]),
                                                 Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps}]),
                                                 #Specific_Kapps=None,
                                                 Exchanges_to_impose=Exchanges_to_impose,
                                                 sims_to_perform=["DefaultKapp"],
                                                 feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                 try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                 print_output=print_outputs,
                                                 apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
        mumax_predicted=Simulation_results["Mu_def"]
        mu_misprediction_factor=mu_measured/mumax_predicted
        mu_iteration_count=0
        mu_iteration_limit=10
        while not (mu_measured-mu_measured*mu_misprediction_tolerance) <= mumax_predicted <= (mu_measured+mu_measured*mu_misprediction_tolerance):
            mu_iteration_count+=1
            if mu_iteration_count>=mu_iteration_limit:
                break
            Default_Kapps["default_efficiency"]*=mu_misprediction_factor
            Default_Kapps["default_transporter_efficiency"]*=mu_misprediction_factor
            Simulation_results = perform_simulations(condition=condition,
                                                     rba_session=rba_session,
                                                     definition_file=definition_file,
                                                     compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                     pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                     process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies}]),
                                                     Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps}]),
                                                     Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps}]),
                                                     Exchanges_to_impose=Exchanges_to_impose,
                                                     sims_to_perform=["DefaultKapp"],
                                                     feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                     try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                     print_output=print_outputs,
                                                     apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
            mumax_predicted=Simulation_results["Mu_def"]
            mu_misprediction_factor=mu_measured/mumax_predicted

    misprediction_coeffs=pandas.DataFrame()
    #if condition=="Hackett_C03":
    #    corrected_spec_kapps=False
    if corrected_spec_kapps:
        misprediction_trajectory={}
        for zero_rxn in zero_flux_rxns:
            Exchanges_to_impose[zero_rxn]={"LB":0,"UB":0}
        #####
        #Exchanges_to_impose=None
        steady_count=0
        iteration_count=0
        iteration_limit=10
        minimum_iteration_number=2
        steady_limit=2
        continuation_criterion=True
        previous_RSS=numpy.nan
        lowest_RSS_soFar=numpy.nan
        RSS_tolerance=0.01
        efficiencies_over_correction_iterations=[]
        RSS_trajectory=[]
        increasing_RSS_count=0
        increasing_RSS_factor=1
        initial_RSS=numpy.nan
        increasing_RSS_limit=2
        while continuation_criterion:
            iteration_count+=1
            no_change_in_global_scaling=True
            Simulation_results = perform_simulations(condition=condition,
                                                     rba_session=rba_session,
                                                     definition_file=definition_file,
                                                     compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                     pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                     process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies}]),
                                                     Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps}]),
                                                     Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps}]),
                                                     Exchanges_to_impose=Exchanges_to_impose,
                                                     sims_to_perform=[condition_to_look_up],
                                                     feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                     try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                     print_output=print_outputs,
                                                     apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
            mumax_predicted=Simulation_results[Growth_rate_to_look_up]
            mu_misprediction_factor=mu_measured/mumax_predicted
            #print("Measured: {} - Predicted: {} - mispred coeff: {}".format(mu_measured,mumax_predicted,mu_misprediction_factor))
            product_misprediction_factors=1
            mu_iteration_count=0
            while not (mu_measured-mu_measured*mu_misprediction_tolerance) <= mumax_predicted <= (mu_measured+mu_measured*mu_misprediction_tolerance):
                mu_iteration_count+=1
                if mu_iteration_count>=10:
                    break
                no_change_in_global_scaling=False
                product_misprediction_factors*=mu_misprediction_factor
                Default_Kapps["default_efficiency"]*=mu_misprediction_factor
                Default_Kapps["default_transporter_efficiency"]*=mu_misprediction_factor
                process_efficiencies.loc[:,"Value"]*=mu_misprediction_factor
                Specific_Kapps.loc[:,"Kapp"]*=mu_misprediction_factor
                if min_kapp is not None:
                    Specific_Kapps.loc[(Specific_Kapps['Kapp']<min_kapp)&(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp']=min_kapp

                Simulation_results = perform_simulations(condition=condition,
                                                        rba_session=rba_session,
                                                        definition_file=definition_file,
                                                        compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                        pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                        process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies}]),
                                                        Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps}]),
                                                        Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps}]),
                                                        Exchanges_to_impose=Exchanges_to_impose,
                                                        sims_to_perform=[condition_to_look_up],
                                                        feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                        try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                        print_output=print_outputs,
                                                        apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
                mumax_predicted=Simulation_results[Growth_rate_to_look_up]
                mu_misprediction_factor=mu_measured/mumax_predicted
                #print("Measured: {} - Predicted: {} - mispred coeff: {}".format(mu_measured,mumax_predicted,mu_misprediction_factor))

            ### New ###
            if fixed_mu_when_above_target_mu_in_correction:
                if mumax_predicted > mu_measured:
                    Simulation_results_fixed = perform_simulations_fixed_Mu(condition=condition,
                                                                    rba_session=rba_session,
                                                                    definition_file=definition_file,
                                                                    compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                                    pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                                    process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies}]),
                                                                    Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps}]),
                                                                    Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps}]),
                                                                    Exchanges_to_impose=Exchanges_to_impose,
                                                                    sims_to_perform=[condition_to_look_up],
                                                                    feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                                    try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                                    print_output=print_outputs,
                                                                    apply_model=False,transporter_multiplier=transporter_multiplier)
                    if len(list(Simulation_results_fixed[Results_to_look_up].keys()))>0:
                        Simulation_results=Simulation_results_fixed
            ######
            if len(list(Simulation_results[Results_to_look_up].keys()))!=0:
                efficiencies_over_correction_iterations.append({"Specific_Kapps":Specific_Kapps.copy(),"Default_Kapps":Default_Kapps.copy(),"Process_Efficiencies":process_efficiencies.copy()})
                KappCorrectionResults=efficiency_correction(specific_kapps=Specific_Kapps,
                                                                 simulation_results=Simulation_results[Results_to_look_up],
                                                                 proteinData=build_input_proteome_for_specific_kapp_estimation(proteome, condition),
                                                                 rba_session=rba_session,
                                                                 condition=condition_to_look_up,
                                                                 also_correct_default=False,
                                                                 default_kapps=Default_Kapps,
                                                                 tolerance=None,
                                                                 previous_misprediction_factors=None,
                                                                 process_efficiencies=process_efficiencies,
                                                                 correct_process_efficiencies=True,
                                                                 correct_default_kapp_enzymes=True)
                current_RSS=KappCorrectionResults["Sum_of_squared_residuals"]
                #if numpy.isfinite(initial_RSS):
                #    if current_RSS>=increasing_RSS_factor*initial_RSS:
                #        increasing_RSS_count+=1
                #else:
                #    initial_RSS=current_RSS
                if numpy.isfinite(lowest_RSS_soFar):
                    if current_RSS>=increasing_RSS_factor*lowest_RSS_soFar:
                        increasing_RSS_count+=1
                    else:
                        increasing_RSS_count=0

                RSS_trajectory.append(current_RSS)
                lowest_RSS_soFar=min(RSS_trajectory)

                Specific_Kapps=KappCorrectionResults["Kapps"]
                if min_kapp is not None:
                    Specific_Kapps.loc[(Specific_Kapps['Kapp']<min_kapp)&(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp']=min_kapp

                process_efficiencies=KappCorrectionResults["ProcessEfficiencies"]
                ###
                spec_kapp_median=Specific_Kapps.loc[(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp'].median()
                Default_Kapps={"default_efficiency":spec_kapp_median,"default_transporter_efficiency":transporter_multiplier*spec_kapp_median}
                ###
                mispred_enzymes=KappCorrectionResults["Enzyme_MispredictionFactors"]
                mispred_processes=KappCorrectionResults["Process_MispredictionFactors"]
                misprediction_trajectory=KappCorrectionResults["Misprediction_factor_trajectory"]
                for mispred_enzyme in mispred_enzymes.keys():
                    misprediction_coeffs.loc[mispred_enzyme,str(iteration_count)]=mispred_enzymes[mispred_enzyme]*product_misprediction_factors
                for mispred_process in mispred_processes.keys():
                    misprediction_coeffs.loc[str("Process_"+mispred_process),str(iteration_count)]=mispred_processes[mispred_process]*product_misprediction_factors
                misprediction_coeffs.loc["Global Misprediction",str(iteration_count)]=product_misprediction_factors

                if iteration_count>=minimum_iteration_number:
                    if not pandas.isna(previous_RSS):
                        if not pandas.isna(current_RSS):
                            if (1-RSS_tolerance)<=current_RSS/previous_RSS<=(1+RSS_tolerance):
                                steady_count+=1
                                #steady_count=0
                            else:
                                steady_count=0
                if print_outputs:
                    print("{} : {} : {}".format(iteration_count,current_RSS,current_RSS/previous_RSS))
                previous_RSS=current_RSS
            if steady_count>=steady_limit:
                continuation_criterion=False
            elif iteration_count>=iteration_limit:
                continuation_criterion=False
            elif increasing_RSS_count>=increasing_RSS_limit:
                continuation_criterion=False

        if len(RSS_trajectory)>0:
            if take_lowest_RSS_parameters:
                lowest_RSS_index=RSS_trajectory.index(min(RSS_trajectory))
                Default_Kapps_to_return=efficiencies_over_correction_iterations[lowest_RSS_index]["Default_Kapps"]
                Specific_Kapps_to_return=efficiencies_over_correction_iterations[lowest_RSS_index]["Specific_Kapps"]
                process_efficiencies_to_return=efficiencies_over_correction_iterations[lowest_RSS_index]["Process_Efficiencies"]
            else:
                Default_Kapps_to_return=efficiencies_over_correction_iterations[-1]["Default_Kapps"]
                Specific_Kapps_to_return=efficiencies_over_correction_iterations[-1]["Specific_Kapps"]
                process_efficiencies_to_return=efficiencies_over_correction_iterations[-1]["Process_Efficiencies"]
        else:
            Default_Kapps_to_return=Default_Kapps
            Specific_Kapps_to_return=Specific_Kapps
            process_efficiencies_to_return=process_efficiencies

        if final_global_scaling_free_exchanges:
            Exchanges_to_impose_here=None
        else:
            Exchanges_to_impose_here=Exchanges_to_impose
        Simulation_results = perform_simulations(condition=condition,
                                                rba_session=rba_session,
                                                definition_file=definition_file,
                                                compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies_to_return}]),
                                                Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps_to_return}]),
                                                Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps_to_return}]),
                                                Exchanges_to_impose=Exchanges_to_impose_here,
                                                sims_to_perform=["Prokaryotic"],
                                                feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                print_output=print_outputs,
                                                apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
        mumax_predicted=Simulation_results["Mu_prok"]
        mu_misprediction_factor=mu_measured/mumax_predicted
        product_misprediction_factors=1
        mu_iteration_count=0
        while not (mu_measured-mu_measured*mu_misprediction_tolerance) <= mumax_predicted <= (mu_measured+mu_measured*mu_misprediction_tolerance):
            mu_iteration_count+=1
            if mu_iteration_count>=10:
                break
            no_change_in_global_scaling=False
            product_misprediction_factors*=mu_misprediction_factor
            Default_Kapps_to_return["default_efficiency"]*=mu_misprediction_factor
            Default_Kapps_to_return["default_transporter_efficiency"]*=mu_misprediction_factor
            process_efficiencies_to_return["Value"]*=mu_misprediction_factor
            Specific_Kapps_to_return["Kapp"]*=mu_misprediction_factor
            if min_kapp is not None:
                Specific_Kapps_to_return.loc[(Specific_Kapps['Kapp']<min_kapp)&(Specific_Kapps_to_return['Kapp']!=0)&(pandas.isna(Specific_Kapps_to_return['Kapp'])==False),'Kapp']=min_kapp
            Simulation_results = perform_simulations(condition=condition,
                                                    rba_session=rba_session,
                                                    definition_file=definition_file,
                                                    compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                    pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                    process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies_to_return}]),
                                                    Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps_to_return}]),
                                                    Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps_to_return}]),
                                                    Exchanges_to_impose=None,
                                                    sims_to_perform=["Prokaryotic"],
                                                    feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                    try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                    print_output=print_outputs,
                                                    apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
            mumax_predicted=Simulation_results["Mu_prok"]
            mu_misprediction_factor=mu_measured/mumax_predicted
    else:
        if final_global_scaling_free_exchanges:
            Exchanges_to_impose_here=None
        else:
            Exchanges_to_impose_here=Exchanges_to_impose
        Simulation_results = perform_simulations(condition=condition,
                                                 rba_session=rba_session,
                                                 definition_file=definition_file,
                                                 compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                 pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                 process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies}]),
                                                 Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps}]),
                                                 Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps}]),
                                                 Exchanges_to_impose=Exchanges_to_impose_here,
                                                 sims_to_perform=["Prokaryotic"],
                                                 feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                 try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                 print_output=print_outputs,
                                                 apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
        mumax_predicted=Simulation_results["Mu_prok"]
        mu_misprediction_factor=mu_measured/mumax_predicted
        product_misprediction_factors=1
        mu_iteration_count=0
        while not (mu_measured-mu_measured*mu_misprediction_tolerance) <= mumax_predicted <= (mu_measured+mu_measured*mu_misprediction_tolerance):
            mu_iteration_count+=1
            if mu_iteration_count>=10:
                break
            no_change_in_global_scaling=False
            product_misprediction_factors*=mu_misprediction_factor
            Default_Kapps["default_efficiency"]*=mu_misprediction_factor
            Default_Kapps["default_transporter_efficiency"]*=mu_misprediction_factor
            process_efficiencies.loc[:,"Value"]*=mu_misprediction_factor
            Specific_Kapps.loc[:,"Kapp"]*=mu_misprediction_factor
            if min_kapp is not None:
                Specific_Kapps.loc[(Specific_Kapps['Kapp']<min_kapp)&(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp']=min_kapp
            Simulation_results = perform_simulations(condition=condition,
                                                    rba_session=rba_session,
                                                    definition_file=definition_file,
                                                    compartment_sizes=extract_compartment_sizes_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                    pg_fractions=extract_pg_fractions_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Densities_PGs":compartment_densities_and_PGs}]),
                                                    process_efficiencies=extract_process_capacities_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Process_Efficiencies":process_efficiencies}]),
                                                    Default_Kapps=extract_default_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Default_Kapps":Default_Kapps}]),
                                                    Specific_Kapps=extract_specific_kapps_from_calibration_outputs(calibration_outputs=[{"Condition":condition,"Specific_Kapps":Specific_Kapps}]),
                                                    Exchanges_to_impose=Exchanges_to_impose_here,
                                                    sims_to_perform=[condition_to_look_up],
                                                    feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],
                                                    try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False,
                                                    print_output=print_outputs,
                                                    apply_model=False,transporter_multiplier=transporter_multiplier,start_val=mu_measured,Mu_approx_precision=Mu_approx_precision,max_mu_in_dichotomy=2*mu_measured)
            mumax_predicted=Simulation_results["Mu_prok"]
            mu_misprediction_factor=mu_measured/mumax_predicted
        RSS_trajectory=[]
        Default_Kapps_to_return=Default_Kapps
        Specific_Kapps_to_return=Specific_Kapps
        process_efficiencies_to_return=process_efficiencies


    if print_outputs:
        print(time.time() - t0)
    return({"RSS_trajectory":RSS_trajectory,
            "Densities_PGs":compartment_densities_and_PGs,
            "Misprediction_Coeffs":misprediction_coeffs,
            "Condition":condition,
            'Proteome': build_input_proteome_for_specific_kapp_estimation(proteome, condition),
            'Correction_Results': correction_results,
            'Default_Kapps': Default_Kapps_to_return,
            'Specific_Kapps': Specific_Kapps_to_return,
            'Process_Efficiencies': process_efficiencies_to_return})

def get_flux_distribution(simulation_outputs,result_object='Simulation_Results', run='Prokaryotic'):
    out=pandas.DataFrame(columns=[sim_result["Condition"] for sim_result in simulation_outputs])
    for sim_result in simulation_outputs:
        if len(list(sim_result[result_object].keys()))>0:
            for rx in sim_result[result_object]['uniqueReactions'].index:
                out.loc[rx,sim_result["Condition"]]=sim_result[result_object]['uniqueReactions'].loc[rx, run]
    return(out)

def get_exchange_flux_distribution(simulation_outputs,result_object='Simulation_Results', run='Prokaryotic'):
    out=pandas.DataFrame(columns=[sim_result["Condition"] for sim_result in simulation_outputs])
    for sim_result in simulation_outputs:
        if len(list(sim_result[result_object].keys()))>0:
            for rx in sim_result[result_object]['ExchangeFluxes'].index:
                out.loc[rx,sim_result["Condition"]]=sim_result[result_object]['ExchangeFluxes'].loc[rx, run]
    return(out)

def extract_proteomes_from_simulation_results(simulation_outputs,type="Prokaryotic"):
    if type=="Prokaryotic":
        results_object="Simulation_Results"
    elif type=="Eukaryotic":
        results_object="Simulation_Results_Euk"
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
            x=1
    return(out)

def plot_protein_protein_comparison(predicted_proteomes,measured_proteomes,conditions):
    from sklearn.linear_model import LinearRegression
    plot_dimensions={1:(1,1),2:(1,2),3:(1,3),4:(2,3),5:(2,3),6:(2,3),7:(3,3),8:(3,3),9:(3,3)}
    plot_indices={1:(1,1),2:(1,2),3:(1,3),4:(2,1),5:(2,2),6:(2,3),7:(3,1),8:(3,2),9:(3,3)}
    number_conditions=len(conditions)
    fig, axs = plt.subplots(plot_dimensions[number_conditions][0], plot_dimensions[number_conditions][1], figsize=(28, 7), sharex=True)
    condition_count=0
    for condition in conditions:
        condition_count+=1
        fig_row=plot_indices[condition_count][0]-1
        fig_col=plot_indices[condition_count][1]-1
        protein_comparison=pandas.DataFrame()
        if condition in list(predicted_proteomes.columns):
            if condition in list(measured_proteomes.columns):
                for protein in predicted_proteomes.index:
                    protein_comparison.loc[protein,"Predicted"]=6.023e20 *predicted_proteomes.loc[protein,condition]
                for protein in measured_proteomes.index:
                    protein_comparison.loc[protein,"Measured"]=6.023e20 *measured_proteomes.loc[protein,condition]
                protein_comparison["Ratio"]=protein_comparison["Predicted"]/protein_comparison["Measured"]
                protein_comparison_to_proceed=protein_comparison.loc[(numpy.isfinite(protein_comparison["Ratio"]))&(protein_comparison["Ratio"]!=0)]
                #y_reg = numpy.reshape(numpy.array([numpy.log10(i) for i in list(protein_comparison_to_proceed["Predicted"])]), (len(list(protein_comparison_to_proceed["Predicted"])), 1))
                #x_reg = numpy.reshape(numpy.array([numpy.log10(i) for i in list(protein_comparison_to_proceed["Measured"])]), (len(list(protein_comparison_to_proceed["Measured"])), 1))
                x_reg = numpy.reshape(numpy.array([i for i in list(protein_comparison_to_proceed["Predicted"])]), (len(list(protein_comparison_to_proceed["Predicted"])), 1))
                y_reg = numpy.reshape(numpy.array([i for i in list(protein_comparison_to_proceed["Measured"])]), (len(list(protein_comparison_to_proceed["Measured"])), 1))
                regressor = LinearRegression(fit_intercept=False)
                regressor.fit(x_reg, y_reg)
                predictions = regressor.predict(x_reg)
                total_max=max([max(list(protein_comparison_to_proceed["Predicted"])),max(list(protein_comparison_to_proceed["Measured"]))])
                total_min=min([min(list(protein_comparison_to_proceed["Predicted"])),min(list(protein_comparison_to_proceed["Measured"]))])
                #[numpy.log10(i) for i in list(protein_comparison_to_prodeed["Predicted"])]
                axs[fig_row, fig_col].plot([numpy.log10(total_min), numpy.log10(total_max)], [numpy.log10(total_min), numpy.log10(total_max)], color='green', linewidth=3,alpha=0.6)
                axs[fig_row, fig_col].scatter([numpy.log10(i) for i in list(protein_comparison_to_proceed["Predicted"])],[numpy.log10(i) for i in list(protein_comparison_to_proceed["Measured"])],alpha=0.4)
                axs[fig_row, fig_col].plot(numpy.log10(x_reg), numpy.log10(predictions), color='red')
                #axs[fig_row, fig_col].plot(x_reg, predictions, color='red')
                axs[fig_row, fig_col].legend(['Identity', "Correlation: {}".format(str(round(regressor.coef_[0][0],2))), 'Data'])
                axs[fig_row, fig_col].plot([numpy.log10(total_min), numpy.log10(total_max)], [numpy.log10(total_min)+2, numpy.log10(total_max)+2], color='black', linewidth=1,alpha=0.6)
                axs[fig_row, fig_col].plot([numpy.log10(total_min), numpy.log10(total_max)], [numpy.log10(total_min)-2, numpy.log10(total_max)-2], color='black', linewidth=1,alpha=0.6)
                #axs[fig_row, fig_col].set_title("Predicted protein numbers per gram dry weight - {}".format(condition))
                axs[fig_row, fig_col].set_title(condition)
                axs[fig_row, fig_col].set_xlabel('Predicted copies per $g_{DW}$ ($Log_{10}$)')
                axs[fig_row, fig_col].set_ylabel('Measured copies per $g_{DW}$ ($Log_{10}$)')
                axs[fig_row, fig_col].set_ylim(math.floor(numpy.log10(min(list(protein_comparison_to_proceed["Measured"])))),math.ceil(numpy.log10(max(list(protein_comparison_to_proceed["Measured"])))))

    plt.show()

def sample_copy_numbers_from_proteome_replicates(Input_data,cols_to_draw_from,target_size=1):
    sample_set=set()
    dimension_too_draw=Input_data.shape[0]
    while len(sample_set)<target_size:
        sample_set.add(tuple(numpy.random.choice(a=cols_to_draw_from,size=dimension_too_draw,replace=True)))
    out=pandas.DataFrame(index=list(Input_data.index))
    out["Gene"]=Input_data["Gene"]
    count=0
    for sample in list(sample_set):
        count+=1
        out["run_{}".format(count)]=[Input_data.loc[list(Input_data.index)[i],sample[i]] for i in list(range(len(sample)))]
    return(out)

def sample_copy_numbers_from_residuals_old(Input_data,replicate_cols,mean_col,replicate_threshold=1,filter_list=[],target_size=1,reps_to_sample=3):
    data_to_use=pandas.DataFrame(columns=Input_data.columns)
    for i in list(Input_data.index):
        finite_count=0
        for j in replicate_cols:
            if not pandas.isna(Input_data.loc[i,j]):
                finite_count+=1
        if finite_count>=replicate_threshold:
            if len(filter_list)>0:
                if i in filter_list:
                    data_to_use.loc[i,:]=Input_data.loc[i,:]
            else:
                data_to_use.loc[i,:]=Input_data.loc[i,:]

    all_residuals=[]
    for i in replicate_cols:
        data_to_use["Log__{}".format(i)]=[numpy.log10(j) for j in list(data_to_use[i])]
    for i in data_to_use.index:
        vals=[]
        for j in replicate_cols:
            if not pandas.isna(data_to_use.loc[i,j]):
                vals.append(numpy.log10(data_to_use.loc[i,j]))
        data_to_use.loc[i,"Log__mean"]=numpy.mean(vals)
    for i in replicate_cols:
        data_to_use["Residual__{}".format(i)]=data_to_use["Log__{}".format(i)]-data_to_use["Log__mean"]
        all_residuals+=list([j for j in list(data_to_use["Residual__{}".format(i)]) if not pandas.isna(j)])

    for i in Input_data.index:
        vals=[]
        for j in replicate_cols:
            if not pandas.isna(Input_data.loc[i,j]):
                vals.append(numpy.log10(Input_data.loc[i,j]))
        Input_data.loc[i,"Log__mean"]=numpy.mean(vals)
    dimension_too_draw=Input_data.shape[0]
    out=pandas.DataFrame(index=list(Input_data.index))
    out["Gene"]=Input_data["Gene"]
    count=0
    df_intermediate=pandas.DataFrame(index=list(Input_data.index))
    df_intermediate["LogMean"]=Input_data["Log__mean"]
    out["mean_noNoise"]=[10**j for j in list(df_intermediate["LogMean"])]
    for run in list(range(target_size)):
        count+=1
        for rep in range(reps_to_sample):
            df_intermediate["LogRes_{}".format(rep+1)]=list(numpy.random.choice(a=all_residuals,size=dimension_too_draw,replace=True))
            df_intermediate["Sampled_{}".format(rep+1)]=df_intermediate["LogMean"]+df_intermediate["LogRes_{}".format(rep+1)]
        for i in list(df_intermediate.index):
            sampled_reps=[]
            for rep in range(reps_to_sample):
                sampled_reps.append(10**df_intermediate.loc[i,"Sampled_{}".format(rep+1)])
            df_intermediate.loc[i,"Sampled_Mean_abs"]=numpy.nanmean(sampled_reps)
        out["run_{}".format(count)]=list(df_intermediate["Sampled_Mean_abs"])
    return(out)

def check_quantile(val,quantiles):
    if not pandas.isna(val):
        for i in range(len(quantiles)):
            if i!=0:
                if (val>quantiles[i]) and (val<=quantiles[i+1]):
                    return(i+1)
            else:
                if (val>=quantiles[i]) and (val<=quantiles[i+1]):
                    return(i+1)
    else:
        return(numpy.nan)

def sample_copy_numbers_from_residuals_quantiles(Input_data,replicate_cols,mean_col,replicate_threshold=1,filter_list=[],target_size=1,reps_to_sample=3,number_quantiles=1,transform_residuals=False,regression_type="lin"):
    dimension_too_draw=Input_data.shape[0]
    out=pandas.DataFrame(index=list(Input_data.index))
    df_intermediate=pandas.DataFrame(index=list(Input_data.index))
    out["Gene"]=Input_data["Gene"]
    out[mean_col]=Input_data[mean_col]

    for i in Input_data.index:
        vals=[]
        for j in replicate_cols:
            Input_data.loc[i,"Log__{}".format(j)]=numpy.log10(Input_data.loc[i,j])
            if not pandas.isna(Input_data.loc[i,j]):
                vals.append(numpy.log10(Input_data.loc[i,j]))
        Input_data.loc[i,"Log__mean"]=numpy.nanmean(vals)
        Input_data.loc[i,"Log__sdev"]=numpy.nanstd(vals)

    if number_quantiles>1:
        quantiles=list(numpy.quantile(a=list(Input_data.loc[pandas.isna(Input_data["Log__mean"])==False,"Log__mean"]),q=[i/number_quantiles for i in list(range(number_quantiles+1))]))
        Input_data["Quantile"]=[check_quantile(val=i,quantiles=quantiles) for i in list(Input_data["Log__mean"])]
    else:
        Input_data["Quantile"]=[1]*Input_data.shape[0]
    out["Quantile"]=Input_data["Quantile"]

    df_intermediate["LogMean"]=Input_data["Log__mean"]
    df_intermediate["LogSdev"]=Input_data["Log__sdev"]
    out["mean_Log_noNoise"]=df_intermediate["LogMean"]
    out["sdev_Log_noNoise"]=df_intermediate["LogSdev"]
    #out["mean_noNoise"]=[10**j for j in list(df_intermediate["LogMean"])]

    data_to_use=pandas.DataFrame(columns=Input_data.columns)
    for i in list(Input_data.index):
        finite_count=0
        for j in replicate_cols:
            if not pandas.isna(Input_data.loc[i,j]):
                finite_count+=1
        Input_data.loc[i,"Number_quantified_replicates"]=finite_count
        out.loc[i,"Number_quantified_replicates"]=finite_count
        if finite_count>=replicate_threshold:
            if not pandas.isna(Input_data.loc[i,"Quantile"]):
                if len(filter_list)>0:
                    if i in filter_list:
                        data_to_use.loc[i,:]=Input_data.loc[i,:]
                else:
                    data_to_use.loc[i,:]=Input_data.loc[i,:]

    if transform_residuals:
        if regression_type=="lin":
        # do linear regression of standard deviation over replicates of protein
            x_reg = numpy.reshape(numpy.array(list(data_to_use["Log__mean"])), (len(list(data_to_use["Log__mean"])), 1))
            y_reg = numpy.reshape(numpy.array(list(data_to_use["Log__sdev"])), (len(list(data_to_use["Log__sdev"])), 1))
            regressor = LinearRegression(fit_intercept=True)
            regressor.fit(x_reg, y_reg)
            slope_sdev=regressor.coef_[0][0]
            offset_sdev=regressor.intercept_[0]
            data_to_use["Fitted_Stdev"]=[offset_sdev+slope_sdev*data_to_use.loc[i,"Log__mean"] for i in data_to_use.index]
            df_intermediate["Fitted_Stdev"]=[offset_sdev+slope_sdev*df_intermediate.loc[i,"LogMean"] for i in df_intermediate.index]
        elif regression_type=="inverse_lin":
            x_reg = numpy.reshape(numpy.array(list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__mean"])), (len(list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__mean"])), 1))
            y_reg = numpy.reshape(numpy.array([1/i for i in list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__sdev"])]), (len(list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__sdev"])), 1))
            regressor = LinearRegression(fit_intercept=True)
            regressor.fit(x_reg, y_reg)
            slope_sdev=regressor.coef_[0][0]
            offset_sdev=regressor.intercept_[0]
            data_to_use["Fitted_Stdev"]=[1/(offset_sdev+slope_sdev*data_to_use.loc[i,"Log__mean"]) for i in data_to_use.index]
            df_intermediate["Fitted_Stdev"]=[1/(offset_sdev+slope_sdev*df_intermediate.loc[i,"LogMean"]) for i in df_intermediate.index]
        out["Fitted_Stdev"]=df_intermediate["Fitted_Stdev"]

    all_residuals={i:[] for i in list(set(list(data_to_use["Quantile"]))) if not pandas.isna(i)}
    #all_residuals[numpy.nan]=[]
    for quantile in list(set(list(data_to_use["Quantile"]))):
        if not pandas.isna(quantile):
            for i in replicate_cols:
                data_to_use.loc[data_to_use["Quantile"]==quantile,"Residual_empirical__{}".format(i)]=data_to_use.loc[data_to_use["Quantile"]==quantile,"Log__{}".format(i)]-data_to_use.loc[data_to_use["Quantile"]==quantile,"Log__mean"]
                if transform_residuals:
                    data_to_use.loc[data_to_use["Quantile"]==quantile,"Residual_empirical__{}".format(i)]/=data_to_use.loc[data_to_use["Quantile"]==quantile,"Fitted_Stdev"]
                all_residuals[quantile]+=list([j for j in list(data_to_use.loc[data_to_use["Quantile"]==quantile,"Residual_empirical__{}".format(i)]) if not pandas.isna(j)])

    count=0
    for i in replicate_cols:
        for  j in data_to_use.index:
            out.loc[j,"Residual_empirical__{}".format(i)]=data_to_use.loc[j,"Residual_empirical__{}".format(i)]
    out2=pandas.DataFrame(index=list(out.index))
    #out2["mean_noNoise"]=out["mean_noNoise"]
    out2["mean_noNoise"]=out[mean_col]
    for run in list(range(target_size)):
        count+=1
        dummyDF_residual=pandas.DataFrame(index=list(Input_data.index))
        dummyDF_sample=pandas.DataFrame(index=list(Input_data.index))
        sample_count=0
        for rep in range(reps_to_sample):
            sample_count+=1
            df_intermediate.loc[(pandas.isna(Input_data["Quantile"])==False)&(Input_data["Number_quantified_replicates"]>=sample_count),"SampledResidual"]=[list(numpy.random.choice(a=all_residuals[Input_data.loc[protein,"Quantile"]],size=1,replace=True))[0] for protein in df_intermediate.loc[(pandas.isna(Input_data["Quantile"])==False)&(Input_data["Number_quantified_replicates"]>=sample_count),:].index]
            if transform_residuals:
                df_intermediate.loc[(pandas.isna(Input_data["Quantile"])==False)&(Input_data["Number_quantified_replicates"]>=sample_count),"SampledResidual"]*=df_intermediate.loc[(pandas.isna(Input_data["Quantile"])==False)&(Input_data["Number_quantified_replicates"]>=sample_count),"Fitted_Stdev"]
            df_intermediate["SampleLogRep_{}__run_{}".format(rep+1,count)]=df_intermediate["SampledResidual"]+df_intermediate["LogMean"]
            out["SampledResidual_{}__run_{}".format(rep+1,count)]=df_intermediate["SampledResidual"]
            out["SampleLogRep_{}__run_{}".format(rep+1,count)]=df_intermediate["SampleLogRep_{}__run_{}".format(rep+1,count)]
            dummyDF_residual["SampledResidual_{}__run_{}".format(rep+1,count)]=df_intermediate["SampledResidual"]
            dummyDF_sample["SampleLogRep_{}__run_{}".format(rep+1,count)]=df_intermediate["SampleLogRep_{}__run_{}".format(rep+1,count)]
        out["MeanSampledResidual__run_{}".format(count)]=dummyDF_residual.mean(axis=1,skipna=True)
        out["MeanSampleLog__run_{}".format(count)]=dummyDF_sample.mean(axis=1,skipna=True)
        out["run_{}".format(count)]=[10**i for i in out["MeanSampleLog__run_{}".format(count)]]
        out2["run_{}".format(count)]=out["run_{}".format(count)]
    out["Mean_of_log_samples"]=out.loc[:,[col for col in out.columns if col.startswith("MeanSampleLog__run_")]].mean(axis=1)
    return(out2)

def estimate_specific_enzyme_efficiencies(rba_session, proteomicsData, flux_bounds, mu, biomass_function=None, target_biomass_function=True, parsimonious_fba=True, only_non_ambigous_proteins=False, chose_isoreaction=False, equalize_identical_enzymes=True, only_identical_reactions_with_twin_enzyme=False, equalize_identical_reactions=True, metabolites_to_ignore=[], impose_on_all_isoreactions=True, zero_on_all_isoreactions=True, condition=None, store_output=True,rxns_to_ignore_when_parsimonious=[]):
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

    old_model = copy.deepcopy(rba_session.model)
    for i in rba_session.model.targets.target_groups._elements_by_id['translation_targets'].concentrations._elements:
        if i.species == 'average_protein_c':
            new_agg = rba.xml.parameters.Aggregate(id_='total_protein', type_='multiplication')
            new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='Total_protein_massfraction_CarbonLimitation'))
            #new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='Total_protein_massfraction_NitrogenLimitation'))
            #new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='Total_protein_massfraction_PhosphorousLimitation'))
            new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='Average_peptide_monomer_inversemillimolarweight'))
            new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='inverse_average_protein_length'))
            rba_session.model.parameters.aggregates._elements.append(new_agg)
            i.value = 'total_protein'
        else:
            rba_session.model.targets.target_groups._elements_by_id['translation_targets'].concentrations._elements.remove(i)
    for i in rba_session.model.targets.target_groups._elements_by_id['transcription_targets'].concentrations._elements:
        if i.species == 'average_mrna_n':
            new_agg = rba.xml.parameters.Aggregate(id_='total_rna', type_='multiplication')
            new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='RNA_massfraction_CarbonLimitation'))
            #new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='RNA_massfraction_NitrogenLimitation'))
            #new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='RNA_massfraction_PhosphorousLimitation'))
            new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='mRNA_inversemillimolarweight'))
            rba_session.model.parameters.aggregates._elements.append(new_agg)
            i.value = 'total_rna'
        else:
            rba_session.model.targets.target_groups._elements_by_id['transcription_targets'].concentrations._elements.remove(i)
    for i in rba_session.model.targets.target_groups._elements_by_id['replication_targets'].concentrations._elements:
        if i.species == 'dna_n':
            new_agg = rba.xml.parameters.Aggregate(id_='total_dna', type_='multiplication')
            new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='DNA_massfraction_CarbonLimitation'))
            #new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='DNA_massfraction_NitrogenLimitation'))
            #new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='DNA_massfraction_PhosphorousLimitation'))
            new_agg.function_references.append(rba.xml.parameters.FunctionReference(function='DNA_inversemillimolarweight'))
            rba_session.model.parameters.aggregates._elements.append(new_agg)
            i.value = 'total_dna'
        else:
            rba_session.model.targets.target_groups._elements_by_id['replication_targets'].concentrations._elements.remove(i)

    rba_session.rebuild_from_model()
    rba_session.set_medium(rba_session.Medium)
    rba_session.add_exchange_reactions()
    rba_session.set_growth_rate(mu)

    if target_biomass_function:
        #rba_session.build_fba_model(objective='targets', maintenanceToBM=True)
        rba_session.build_fba_model(objective='targets', maintenanceToBM=False)
        BMfunction = 'R_BIOMASS_targetsRBA'
    else:
        rba_session.build_fba_model(objective='classic', maintenanceToBM=False)
        BMfunction = biomass_function

    for j in [i for i in rba_session.Medium.keys() if rba_session.Medium[i] == 0]:
        Exrxn = 'R_EX_'+j.split('M_')[-1]+'_e'
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
    rba_session.FBA.set_lb(rxn_LBs)
    rba_session.FBA.set_ub(rxn_UBs)

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
        #ObjectiveValue

    if parsimonious_fba:
        rba_session.FBA.parsimonise(rxns_to_ignore_in_objective=rxns_to_ignore_when_parsimonious)
        rba_session.FBA.set_lb(rxn_LBs)
        rba_session.FBA.set_ub(rxn_UBs)
        rba_session.FBA.set_lb({BMfunction: BMfluxOld})
        rba_session.FBA.set_ub({BMfunction: BMfluxOld})
        rba_session.FBA.solve_lp()

    FluxDistribution = pandas.DataFrame(index=list(rba_session.FBA.SolutionValues.keys()), columns=['FluxValues'])
    FluxDistribution['FluxValues'] = list(rba_session.FBA.SolutionValues.values())

    ProtoIDmap = {}
    for i in rba_session.ModelStructure.ProteinInfo.Elements.keys():
        ProtoID = rba_session.ModelStructure.ProteinInfo.Elements[i]['ProtoID']
        if ProtoID in list(proteomicsData['ID']):
            if not pandas.isna(proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0]):
                if proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0] != 0:
                    if ProtoID in ProtoIDmap.keys():
                        ProtoIDmap[ProtoID]['ModelProteins'].append(i)
                    else:
                        ProtoIDmap.update({ProtoID: {'ModelProteins': [i], 'CopyNumber': proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0]}})

    # NEW METHOD
    # identify all model protein-isoforms, associated with the measured proteins
    # return dict with (proto-)protein IDs as keys and list of associated compartment-isoforms as values
    ProtoProteinMap = rba_session.ModelStructure.ProteinInfo.return_protein_iso_form_map()
    measured_Proteins_Isoform_Map = {p_ID: ProtoProteinMap[p_ID] for p_ID in list(proteomicsData['ID']) if p_ID in list(ProtoProteinMap.keys())}

    # identify all model reactions, associated with the measured proteins
    measured_Proteins_Reaction_Map = {}
    for p_ID in measured_Proteins_Isoform_Map.keys():
        reactions_associated_with_proto_Protein = []
        for isoform in measured_Proteins_Isoform_Map[p_ID]:
            for reaction in rba_session.ModelStructure.ProteinInfo.Elements[isoform]['associatedReactions']:
                #reactions_associated_with_proto_Protein.append(reaction.split('_duplicate')[0])
                reactions_associated_with_proto_Protein.append(reaction)
        if only_non_ambigous_proteins:
            if len(list(reactions_associated_with_proto_Protein)) == 1:
                measured_Proteins_Reaction_Map[p_ID] = list(reactions_associated_with_proto_Protein)
        else:
            measured_Proteins_Reaction_Map[p_ID] = list(set(reactions_associated_with_proto_Protein))
    # choose most likely iso-reaction for each measured-protein associated reaction
    protoRxnDict = {}
    chosen_Isoreactions = {}
    for p_ID in measured_Proteins_Reaction_Map.keys():
        for rxn in measured_Proteins_Reaction_Map[p_ID]:
            rxn_to_split=str(rxn)
            protoRxn = rxn_to_split.split('_duplicate')[0]
            if protoRxn in list(protoRxnDict.keys()):
                if rxn in list(protoRxnDict[protoRxn].keys()):
                    protoRxnDict[protoRxn][rxn] += 1
                else:
                    protoRxnDict[protoRxn].update({rxn: 1})
            else:
                protoRxnDict[protoRxn] = {rxn: 1}
    for prx in protoRxnDict.keys():
        unique_SU_dict={}
        for irx in protoRxnDict[prx].keys():
            enzyme = rba_session.ModelStructure.ReactionInfo.Elements[irx]['Enzyme']
            unique_SU_dict[irx] = len(list(rba_session.ModelStructure.EnzymeInfo.Elements[enzyme]['Subunits'].keys()))
        if chose_isoreaction:
            max_val = max([protoRxnDict[prx][i]/unique_SU_dict[i] for i in protoRxnDict[prx].keys()])
            l = [i for i in protoRxnDict[prx].keys() if protoRxnDict[prx][i]/unique_SU_dict[i] == max_val]
            if len(l)>1:
                max_SU_number=max([unique_SU_dict[i] for i in l])
                l_new=[i for i in l if unique_SU_dict[i]==max_SU_number]
                l=l_new
            l.sort()
            #chosen_Isoreactions[prx] = ';,;'.join([l[0]])
            chosen_Isoreactions[prx] = ';,;'.join(l)
        else:
            l = [i for i in protoRxnDict[prx].keys() if protoRxnDict[prx][i] != 0]
            l.sort()
            chosen_Isoreactions[prx] = ';,;'.join(l)
    chosen_isoreaction_DF = pandas.DataFrame()
    protoreactionsforDF = list(chosen_Isoreactions.keys())
    chosen_isoreaction_DF['ProtoReaction'] = protoreactionsforDF
    chosen_isoreaction_DF['ChosenIsoReaction'] = [chosen_Isoreactions[pr] for pr in protoreactionsforDF]

    # determine model reactions with non-zero flux in FBA#
    overview_out = pandas.DataFrame()
    for i in list(chosen_Isoreactions.keys()):
        all_iso_rxns = chosen_Isoreactions[i].split(';,;')
        if i in list(FluxDistribution.index):
            rxn_flux = FluxDistribution.loc[i, 'FluxValues']
        else:
            rxn_flux = 0
        all_enzyme_concentrations = {}
        for iso_rxn in all_iso_rxns:
            corresponding_enzyme = rba_session.ModelStructure.EnzymeInfo.Elements[rba_session.ModelStructure.ReactionInfo.Elements[iso_rxn]['Enzyme']]
            CompositionDict = {rba_session.ModelStructure.ProteinInfo.Elements[j]['ProtoID']: corresponding_enzyme['Subunits'][j] for j in corresponding_enzyme['Subunits'].keys()}
            CopyNumbers = []
            Stoichiometries = []
            EnzymeNumbers = {}
            for j in CompositionDict.keys():
                if j in ProtoIDmap.keys():
                    CopyNumbers.append(ProtoIDmap[j]['CopyNumber'])
                    Stoichiometries.append(CompositionDict[j])
                    EnzymeNumbers.update({j: ProtoIDmap[j]['CopyNumber']/CompositionDict[j]})

            GM_enzymenumber = gmean(numpy.array(list(EnzymeNumbers.values())))
            if (numpy.isfinite(GM_enzymenumber)) and (GM_enzymenumber != 0):
                all_enzyme_concentrations[iso_rxn] = GM_enzymenumber
        overall_enzyme_concentration = gmean(list(all_enzyme_concentrations.values()))
        if (overall_enzyme_concentration!=0)&(pandas.isna(overall_enzyme_concentration)==False):
            if (rxn_flux !=0)&(pandas.isna(rxn_flux)==False):
                for iso_rxn in all_iso_rxns:
                    overview_out.loc[iso_rxn, 'Comment'] = '1'
                    overview_out.loc[iso_rxn, 'Enzyme_ID'] = rba_session.ModelStructure.ReactionInfo.Elements[iso_rxn]['Enzyme']
                    overview_out.loc[iso_rxn, 'ChosenIsoReaction'] = chosen_Isoreactions[i]
                    overview_out.loc[iso_rxn, 'Proto_Reaction'] = i
                    overview_out.loc[iso_rxn, 'Individual_Isozyme_Concentrations'] = json.dumps(all_enzyme_concentrations)
                    overview_out.loc[iso_rxn, 'CopyNumber'] = overall_enzyme_concentration
                    overview_out.loc[iso_rxn, 'Concentration'] = overall_enzyme_concentration
                    overview_out.loc[iso_rxn, 'Apparent_Concentration'] = overall_enzyme_concentration
                    overview_out.loc[iso_rxn, 'Apparent_Flux'] = rxn_flux
                    overview_out.loc[iso_rxn, 'Flux_all_promisc_rxns'] = rxn_flux
                    overview_out.loc[iso_rxn, 'Flux_FBA'] = rxn_flux

    FluxDistribution.to_csv('Calib_FluxDist_'+condition+'_.csv', sep=';')
    chosen_isoreaction_DF.to_csv('Chosen_IsoRxns_'+condition+'_.csv', sep=';')
    if equalize_identical_enzymes:
        already_handled = []
        identical_enzymes = []
        for i in list(overview_out['Enzyme_ID']):
            identical_set = []
            if i not in already_handled:
                if len(rba_session.ModelStructure.EnzymeInfo.Elements[i]['EnzymesWithIdenticalSubunitComposition']) > 0:
                    identical_set.append(i)
                    already_handled.append(i)
                    for j in rba_session.ModelStructure.EnzymeInfo.Elements[i]['EnzymesWithIdenticalSubunitComposition']:
                        if j not in already_handled:
                            if only_identical_reactions_with_twin_enzyme:
                                R1 = rba_session.ModelStructure.ReactionInfo.Elements[rba_session.ModelStructure.EnzymeInfo.Elements[i]['Reaction']]
                                R2 = rba_session.ModelStructure.ReactionInfo.Elements[rba_session.ModelStructure.EnzymeInfo.Elements[j]['Reaction']]
                                if (len(R1['Compartment_Species']) == 1) & (len(R2['Compartment_Species']) == 1):
                                    reactants1 = {i.rsplit('_{}'.format(R1['Compartment_Species'][0]), 1)[0]: R1['Reactants'][i] for i in list(R1['Reactants'].keys()) if i.rsplit('_{}'.format(R1['Compartment_Species'][0]), 1)[0] not in metabolites_to_ignore}
                                    reactants2 = {i.rsplit('_{}'.format(R2['Compartment_Species'][0]), 1)[0]: R2['Reactants'][i] for i in list(R2['Reactants'].keys()) if i.rsplit('_{}'.format(R2['Compartment_Species'][0]), 1)[0] not in metabolites_to_ignore}
                                    products1 = {i.rsplit('_{}'.format(R1['Compartment_Species'][0]), 1)[0]: R1['Products'][i] for i in list(R1['Products'].keys()) if i.rsplit('_{}'.format(R1['Compartment_Species'][0]), 1)[0] not in metabolites_to_ignore}
                                    products2 = {i.rsplit('_{}'.format(R2['Compartment_Species'][0]), 1)[0]: R2['Products'][i] for i in list(R2['Products'].keys()) if i.rsplit('_{}'.format(R2['Compartment_Species'][0]), 1)[0] not in metabolites_to_ignore}
                                    if (reactants1 == reactants2) & (products1 == products2):
                                        identical_set.append(j)
                                        already_handled.append(j)
                            else:
                                identical_set.append(j)
                                already_handled.append(j)
                    identical_enzymes.append(identical_set)

        ident_set_overview = {}
        for ident_set in identical_enzymes:
            #enzyme_concentration = list(set([overview_out.loc[overview_out['Enzyme_ID'] == enz, 'CopyNumber'].values[0] for enz in ident_set if enz in list(overview_out['Enzyme_ID'])]))[0]
            enzyme_concentration = sum(list(set([overview_out.loc[overview_out['Enzyme_ID'] == enz, 'CopyNumber'].values[0] for enz in ident_set if enz in list(overview_out['Enzyme_ID'])])))
            tot_flux = 0
            for enz in ident_set:
                rxn_id = rba_session.ModelStructure.EnzymeInfo.Elements[enz]['Reaction']
                if rxn_id in list(FluxDistribution.index):
                    tot_flux += abs(FluxDistribution.loc[rxn_id, 'FluxValues'])
                if enz not in list(ident_set_overview.keys()):
                    ident_set_overview[enz] = {'Set_of_Identical_Enzymes': ident_set}

            if (enzyme_concentration!=0)&(pandas.isna(enzyme_concentration)==False):
                if (tot_flux !=0)&(pandas.isna(tot_flux)==False):
                    for enz in ident_set:
                        rxn_id = rba_session.ModelStructure.EnzymeInfo.Elements[enz]['Reaction']
                        if rxn_id in list(FluxDistribution.index):
                            ident_set_overview[enz]['Respective_Enzyme_Level'] = enzyme_concentration * FluxDistribution.loc[rxn_id, 'FluxValues']/tot_flux
                        overview_out.loc[rxn_id, 'Comment'] = '2'
                        overview_out.loc[rxn_id,'Enzyme_ID'] = rba_session.ModelStructure.ReactionInfo.Elements[rxn_id]['Enzyme']
                        overview_out.loc[rxn_id, 'Proto_Reaction'] = rxn_id.split('_duplicate_')[0]
                        if rxn_id in list(FluxDistribution.index):
                            overview_out.loc[rxn_id, 'Concentration'] = enzyme_concentration* FluxDistribution.loc[rxn_id, 'FluxValues']/tot_flux
                        else:
                            overview_out.loc[rxn_id, 'Concentration'] = 0
                        overview_out.loc[rxn_id, 'Apparent_Concentration'] = enzyme_concentration
                        overview_out.loc[rxn_id, 'Apparent_Flux'] = tot_flux
                        overview_out.loc[rxn_id, 'Flux_all_promisc_rxns'] = tot_flux
                        overview_out.loc[rxn_id, 'CopyNumber'] = enzyme_concentration
                        if rxn_id in list(FluxDistribution.index):
                            overview_out.loc[rxn_id, 'Flux_FBA'] = FluxDistribution.loc[rxn_id, 'FluxValues']


    if equalize_identical_reactions:
        RBAmets = [i for i in rba_session.ModelStructure.MetaboliteInfo.Elements.keys() if i.rsplit('_{}'.format(rba_session.ModelStructure.MetaboliteInfo.Elements[i]['Compartment'][0]), 1)[0] not in metabolites_to_ignore]
        RBArxns = list(rba_session.ModelStructure.ReactionInfo.Elements.keys())
        Metabolite_rows = [i for i in list(rba_session.Problem.LP.row_names) if i in RBAmets]
        Rxn_cols = [i for i in list(rba_session.Problem.LP.col_names) if (i in RBArxns) and ('_duplicate' not in i)]
        A = pandas.DataFrame(index=list(rba_session.Problem.LP.row_names), columns=list(rba_session.Problem.LP.col_names), data=rba_session.Problem.LP.A.toarray())
        S = A.loc[Metabolite_rows, Rxn_cols].abs()
        S_t = S.T
        S_t_dup = S_t.loc[S_t.duplicated(keep=False), :]

        Du_index = list(S_t_dup.index)
        ident_sets = {}
        for i in Du_index:
            if rba_session.ModelStructure.ReactionInfo.Elements[i]['Compartment_Species'] != ['e']:
                l = []
                for j in Du_index:
                    if j != i:
                        if (rba_session.ModelStructure.ReactionInfo.Elements[i]['Reactants'] == rba_session.ModelStructure.ReactionInfo.Elements[j]['Reactants']) and (rba_session.ModelStructure.ReactionInfo.Elements[i]['Products'] == rba_session.ModelStructure.ReactionInfo.Elements[j]['Products']):
                            l.append(j)
                        elif (rba_session.ModelStructure.ReactionInfo.Elements[i]['Reactants'] == rba_session.ModelStructure.ReactionInfo.Elements[j]['Products']) and (rba_session.ModelStructure.ReactionInfo.Elements[i]['Products'] == rba_session.ModelStructure.ReactionInfo.Elements[j]['Reactants']):
                            if (rba_session.ModelStructure.ReactionInfo.Elements[i]['Reversible']) or (rba_session.ModelStructure.ReactionInfo.Elements[j]['Reversible']):
                                l.append(j)
                if len(l) >= 1:
                    ident_sets[i] = l
        rxns_already_handled = []
        for i in list(ident_sets.keys()):
            if i not in rxns_already_handled:
                tot_flux = 0
                tot_enzyme_level = 0
                for j in ident_sets[i]:
                    if j not in rxns_already_handled:
                        rxns_already_handled.append(j)
                        if j in list(overview_out.index):
                            tot_enzyme_level += abs(overview_out.loc[j, 'Concentration'])
                        if j in list(FluxDistribution.index):
                            tot_flux += abs(FluxDistribution.loc[j, 'FluxValues'])
                for j in ident_sets[i]:
                    overview_out.loc[i, 'Comment'] = '3'
                    overview_out.loc[j,'Enzyme_ID'] = rba_session.ModelStructure.ReactionInfo.Elements[j]['Enzyme']
                    overview_out.loc[j, 'Proto_Reaction'] = j
                    overview_out.loc[j, 'Apparent_Concentration'] = tot_enzyme_level
                    overview_out.loc[j, 'Apparent_Flux'] = tot_flux
                    if j in list(FluxDistribution.index):
                        overview_out.loc[j, 'Flux_FBA'] = FluxDistribution.loc[j, 'FluxValues']
                    else:
                        overview_out.loc[j, 'Flux_FBA'] = 0

    overview_out['Kapp']=[abs(overview_out.loc[i,'Apparent_Flux']/overview_out.loc[i,'Apparent_Concentration']) for i in overview_out.index]
    overview_out['Flux']=[overview_out.loc[i,'Flux_FBA']/abs(overview_out.loc[i,'Flux_FBA']) for i in overview_out.index]

    if impose_on_all_isoreactions:
        for rx in list(overview_out.index):
            if rx in list(rba_session.ModelStructure.ReactionInfo.Elements.keys()):
                enz = rba_session.ModelStructure.ReactionInfo.Elements[rx]['Enzyme']
            elif rx in list(rba_session.ModelStructure.EnzymeInfo.Elements.keys()):
                enz = rx
            if enz in list(rba_session.ModelStructure.EnzymeInfo.Elements.keys()):
                for isoenz in rba_session.ModelStructure.EnzymeInfo.Elements[enz]['Isozymes']:
                    isorx = rba_session.ModelStructure.EnzymeInfo.Elements[isoenz]['Reaction']
                    if isorx not in overview_out.index:
                        overview_out.loc[isorx, 'Comment'] = '4'
                        overview_out.loc[isorx,'Enzyme_ID'] = isoenz
                        overview_out.loc[isorx, 'Proto_Reaction'] = overview_out.loc[rx, 'Proto_Reaction']
                        overview_out.loc[isorx, 'Kapp'] = overview_out.loc[rx, 'Kapp']
                        overview_out.loc[isorx, 'Flux'] = overview_out.loc[rx, 'Flux']

    if zero_on_all_isoreactions:
        for rx in list(overview_out.index):
            if rx in list(rba_session.ModelStructure.ReactionInfo.Elements.keys()):
                enz = rba_session.ModelStructure.ReactionInfo.Elements[rx]['Enzyme']
            elif rx in list(rba_session.ModelStructure.EnzymeInfo.Elements.keys()):
                enz = rx
            if enz in list(rba_session.ModelStructure.EnzymeInfo.Elements.keys()):
                for isoenz in rba_session.ModelStructure.EnzymeInfo.Elements[enz]['Isozymes']:
                    isorx = rba_session.ModelStructure.EnzymeInfo.Elements[isoenz]['Reaction']
                    if isorx not in overview_out.index:
                        overview_out.loc[isorx, 'Comment'] = '4'
                        overview_out.loc[isorx,'Enzyme_ID'] =isoenz
                        overview_out.loc[isorx, 'Proto_Reaction'] = overview_out.loc[rx, 'Proto_Reaction']
                        overview_out.loc[isorx, 'Kapp'] = 0.0

    rba_session.model = old_model
    rba_session.rebuild_from_model()
    rba_session.set_medium(rba_session.Medium)

    overview_out.sort_index(inplace=True)

    if store_output:
        if condition is not None:
            overview_out.to_csv('SpecKapp_Old_overview_'+condition+'_.csv', sep=';')
        else:
            overview_out.to_csv('SpecKapp_Old_overview_.csv', sep=';')
    #overview_out.to_csv('SpecKapp_5_method_.csv'.format(), sep=';')
    return({"Overview":overview_out,"Flux_Distribution":FluxDistribution})

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
