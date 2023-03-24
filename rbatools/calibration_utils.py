import rba
import copy
import pandas
import time
import numpy
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


def determine_compartment_occupation(Data_input,
                                       Condition,
                                       mass_col='AA_residues',
                                       only_in_model=False,
                                       compartments_to_replace={'DEF':"c"},
                                       compartments_no_original_PG=[],
                                       ribosomal_proteins_as_extra_compartment=True):
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


def build_proteome_overview(input, condition, compartments_to_replace={'DEF':"c", 'DEFA':"c", 'Def':"c"}, compartments_no_original_PG=['n', 'Secreted'], ribosomal_proteins_as_extra_compartment=True):
    out = determine_compartment_occupation(Data_input=input, Condition=condition, compartments_to_replace=compartments_to_replace,
                                           compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=ribosomal_proteins_as_extra_compartment, only_in_model=False)
    out_in_model = determine_compartment_occupation(Data_input=input, Condition=condition, compartments_to_replace=compartments_to_replace,
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


def flux_bounds_from_input(input,rba_session, condition, specific_exchanges=None, specific_directions=None,also_consider_iso_enzmes=True):

    out = pandas.DataFrame(columns=['Reaction_ID', 'LB', 'UB'])

    flux_mean_df = input.loc[input['Type'] == 'ExchangeFlux_Mean', :]
    flux_mean_SE = input.loc[input['Type'] == 'ExchangeFlux_StandardError', :]
    if specific_exchanges is None:
        exchanges_to_set = list(flux_mean_df['ID'])
    else:
        exchanges_to_set = specific_exchanges
    for rx in exchanges_to_set:
        reactions_to_consider=[rx]
        if also_consider_iso_enzmes:
            if rx in rba_session.get_reactions():
                reactions_to_consider+=list(rba_session.get_reaction_information(rx)['Twins'])
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
                        feasible_stati=["optimal","feasible"],
                        try_unscaling_if_sol_status_is_feasible_only_before_unscaling=True,
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
        sol_status=rba_session.Problem.SolutionStatus
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
        sol_status=rba_session.Problem.SolutionStatus
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
        sol_status=rba_session.Problem.SolutionStatus
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
    return({"SolutionStatus":sol_status,"FeasibleRange_def":def_Feasible_Ranges,"FeasibleRange_prok":prok_Feasible_Ranges,"FeasibleRange_euk":euk_Feasible_Ranges,"Mu_def":mumax_def,"Mu_prok":mumax_prok,"Mu_euk":mumax_euk,'Simulation_Results': prok_results, 'Simulation_Results_Euk': euk_results, 'Simulation_Results_DefKapp': def_results, 'Euk_CompSizes': compartment_fractions,"Condition":condition})


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
        axs[ax_indices[count][0],ax_indices[count][1]].scatter([cal_input["RSS_trajectory"].index(min(cal_input["RSS_trajectory"]))+1] , min(cal_input["RSS_trajectory"]) , color="red")
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
    subunit_misprediction_factors_default_efficiency_enzymes={}
    subunit_misprediction_factors_processes={}
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
                        #squared_residuals.append((predicted_protein-measured_protein)**2)
                        #squared_residuals.append(((predicted_protein-measured_protein)*rba_session.get_protein_information(protein=proto_protein_isoform_map[proto_protein_ID][0])["AAnumber"])**2)
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
                                    else:
                                        subunit_misprediction_factors_enzymes[enzyme].append(misprediction_coeff)
                                else:
                                    if enzyme not in subunit_misprediction_factors_default_efficiency_enzymes.keys():
                                        subunit_misprediction_factors_default_efficiency_enzymes[enzyme]=[misprediction_coeff]
                                    else:
                                        subunit_misprediction_factors_default_efficiency_enzymes[enzyme].append(misprediction_coeff)
                            for process in rba_session.get_protein_information(protein=protein)["SupportsProcess"]:
                                if process in process_efficiencies.index:
                                    if process not in subunit_misprediction_factors_processes.keys():
                                        subunit_misprediction_factors_processes[process]=[misprediction_coeff]
                                    else:
                                        subunit_misprediction_factors_processes[process].append(misprediction_coeff)

    enzyme_correction_coefficients={}
    for enzyme in subunit_misprediction_factors_enzymes.keys():
        enzyme_correction_coefficients[enzyme]=numpy.power(numpy.median(subunit_misprediction_factors_enzymes[enzyme]),1/n_th_root_mispred)
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
                    if tolerance is None:
                        enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"]=new_kapp
                    else:
                        if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                            enzyme_efficiencies_out.loc[enzyme_efficiencies_out["Enzyme_ID"]==iso_enzyme_to_consider,"Kapp"]=new_kapp

    if correct_default_kapp_enzymes:
        default_efficiency_enzyme_correction_coefficients={}
        for enzyme in subunit_misprediction_factors_default_efficiency_enzymes.keys():
            if enzyme not in list(enzyme_efficiencies["Enzyme_ID"]):
                default_efficiency_enzyme_correction_coefficients[enzyme]=numpy.power(numpy.median(subunit_misprediction_factors_default_efficiency_enzymes[enzyme]),1/n_th_root_mispred)
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
            process_correction_coefficients[process]=numpy.power(numpy.median(subunit_misprediction_factors_processes[process]),1/n_th_root_mispred)
    for process in process_correction_coefficients.keys():
        correction_coeff=process_correction_coefficients[process]
        old_efficiency=process_efficiencies.loc[process,"Value"]
        new_efficiency=old_efficiency*correction_coeff
        if tolerance is None:
            process_efficiencies_out.loc[process,"Value"]=new_efficiency
        else:
            if abs(numpy.log(tolerance)) <= abs(numpy.log(correction_coeff)):
                process_efficiencies_out.loc[process,"Value"]=new_efficiency

#    return({"sum_of_squared_residuals":sum(squared_residuals),"enzyme_efficiencies":enzyme_efficiencies_out,"process_efficiencies":process_efficiencies_out,"process_efficiency_correction_factors":process_correction_coefficients,"enzyme_efficiency_correction_factors":enzyme_correction_coefficients})
    return({"Sum_of_squared_residuals":sum(squared_residuals),"Kapps":enzyme_efficiencies_out,"ProcessEfficiencies":process_efficiencies_out,"Process_MispredictionFactors":process_correction_coefficients,"Enzyme_MispredictionFactors":enzyme_correction_coefficients})


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


def generate_mean_enzyme_composition_model(rba_session,condition):
    enzymes_already_handled=[]
    for i in rba_session.get_enzymes():
        if i in enzymes_already_handled:
            continue
        all_iso_enzymes=list([i]+rba_session.get_enzyme_information(i)['Isozymes'])
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
                print(j)
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
        enzymes_already_handled+=all_iso_enzymes    
        rba_session.build_model_structure(file_name='/ModelStructure_meancompo_{}.json'.format(condition),print_warnings=False)
    


def calibration_workflow(proteome,
                         condition,
                         reference_condition,
                         gene_ID_column,
                         definition_file,
                         rba_session,
                         process_efficiency_estimation_input=None,
                         spec_kapps=None,
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
                         use_mean_enzyme_composition_for_calibration=False):
    
    correction_settings=machinery_efficiency_correction_settings_from_input(input=definition_file, condition=condition)
    enzyme_efficiency_estimation_settings=enzyme_efficiency_estimation_settings_from_input(input=definition_file, condition=condition)
    t0 = time.time()
    correction_results_compartement_sizes = correction_pipeline(input=proteome,
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
        compartment_densities_and_PGs = extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=correction_results_compartement_sizes,rows_to_exclude=["Ribosomes","Total"]+[i for i in correction_results_compartement_sizes.index if i.startswith("pg_")])
        correction_results_compartement_sizes.to_csv(str('Correction_overview_HackettNielsen_'+condition+'.csv'))
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

    correction_results_compartement_sizes.to_csv(str('Correction_overview_HackettNielsen_corrected_'+condition+'.csv'))
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
    ### define coeff as input ###
    proteome[condition]*=global_protein_scaling_coeff
    ###

    process_efficiencies.to_csv("ProcEffsOrig_{}.csv".format(condition))

    if use_mean_enzyme_composition_for_calibration:
        generate_mean_enzyme_composition_model(rba_session,condition)

    if spec_kapps is None:
        flux_bounds_fba=flux_bounds_from_input(input=definition_file,rba_session=rba_session, condition=condition, specific_exchanges=None, specific_directions=None,also_consider_iso_enzmes=False)
        Specific_Kapps_Results = estimate_specific_enzyme_efficiencies(rba_session=rba_session, 
                                                                               proteomicsData=build_input_proteome_for_specific_kapp_estimation(proteome, condition), 
                                                                               flux_bounds=flux_bounds_fba, 
                                                                               mu=growth_rate_from_input(input=definition_file, condition=condition), 
                                                                               biomass_function=enzyme_efficiency_estimation_settings['biomass_function_in_model'], 
                                                                               target_biomass_function=enzyme_efficiency_estimation_settings['use_target_biomass_function'], 
                                                                               parsimonious_fba=enzyme_efficiency_estimation_settings['parsimonious_fba'], 
                                                                               chose_most_likely_isoreactions=enzyme_efficiency_estimation_settings['chose_most_likely_isoreactions'],
                                                                               impose_on_all_isoreactions=enzyme_efficiency_estimation_settings['impose_on_all_isoreactions'], 
                                                                               zero_on_all_isoreactions=enzyme_efficiency_estimation_settings['zero_on_all_isoreactions'],
                                                                               node_degree_identical_enzyme_network=enzyme_efficiency_estimation_settings['node_degree_identical_enzyme_network'],
                                                                               impose_on_identical_enzymes=enzyme_efficiency_estimation_settings['impose_on_identical_enzymes'],
                                                                               condition=condition, 
                                                                               store_output=True,
                                                                               rxns_to_ignore_when_parsimonious=[])

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
                        for iso_rxn in rba_session.get_reaction_information(fba_rxn)['Twins']:
                            if flux_value<0:
                                fba_flux_directions.update({iso_rxn:{"LB":numpy.nan,"UB":0}})
                            elif flux_value>0:
                                fba_flux_directions.update({iso_rxn:{"LB":0,"UB":numpy.nan}})

        Specific_Kapps.to_csv("Specific_Kapps_Hackett__{}.csv".format(condition), sep=";", decimal=",")
    else:
        print('importing spec kapps')
        Specific_Kapps = pandas.DataFrame()
        Specific_Kapps['Enzyme_ID'] = spec_kapps['ID']
        Specific_Kapps['Kapp'] = spec_kapps[condition]
        Specific_Kapps['Flux'] = spec_kapps[str(condition+'_Flux')]

    spec_kapp_median=Specific_Kapps.loc[(Specific_Kapps['Kapp']!=0)&(pandas.isna(Specific_Kapps['Kapp'])==False),'Kapp'].median()
    Default_Kapps={"default_efficiency":spec_kapp_median,"default_transporter_efficiency":transporter_multiplier*spec_kapp_median}

    flux_bounds_data=flux_bounds_from_input(input=definition_file,rba_session=rba_session, condition=condition, specific_exchanges=None, specific_directions=None,also_consider_iso_enzmes=True)
    Exchanges_to_impose={i:{"LB":flux_bounds_data.loc[i,"LB"],"UB":flux_bounds_data.loc[i,"UB"]} for i in list(flux_bounds_data["Reaction_ID"])}
    if correction_settings['impose_directions_from_fba_during_correction']:
        for rxn in fba_flux_directions.keys():
            if rxn not in Exchanges_to_impose.keys():
                Exchanges_to_impose.update({rxn:fba_flux_directions[rxn]})
        Exchanges_to_impose.update(fba_flux_directions)
    compartment_densities_and_PGs = extract_compsizes_and_pgfractions_from_correction_summary(corrsummary=correction_results_compartement_sizes,rows_to_exclude=["Ribosomes","Total"]+[i for i in correction_results_compartement_sizes.index if i.startswith("pg_")])


    Specific_Kapps_original=Specific_Kapps.copy()
    Default_Kapps_original=Default_Kapps.copy()
    process_efficiencies_original=process_efficiencies.copy()

    if correction_settings['correct_efficiencies']:
        steady_count=0
        increasing_RSS_count=0
        iteration_count=0
        continuation_criterion=True
        previous_RSS=numpy.nan
        efficiencies_over_correction_iterations=[]
        rss_trajectory=[]
        minimum_iteration_number=correction_settings['minimum_iteration_number']
        iteration_limit=correction_settings['iteration_limit']
        steady_limit=correction_settings['steady_rss_limit']
        increasing_RSS_limit=correction_settings['increasing_rss_limit']
        rss_tolerance=correction_settings['rss_tolerance']
        increasing_RSS_factor=correction_settings['increasing_rss_factor']
        while continuation_criterion:
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
                                                             condition_to_look_up="Prokaryotic",
                                                             growth_rate_to_look_up="Mu_prok",
                                                             results_to_look_up="Simulation_Results",
                                                             fixed_mu_when_above_target_mu_in_correction=correction_settings['fixed_growth_rate_global_scaling'],
                                                             n_th_root_mispred=1,
                                                             print_outputs=False,
                                                             adjust_root=correction_settings['abjust_root_of_correction_coeffs_global_scaling'])

            
            Simulation_results=results_global_scaling["simulation_results"]
            Specific_Kapps=results_global_scaling["specific_kapps"]
            Default_Kapps=results_global_scaling["default_kapps"]
            process_efficiencies=results_global_scaling["process_efficiencies"]

            if len(list(Simulation_results["Simulation_Results"].keys()))!=0:
                efficiencies_over_correction_iterations.append({"Specific_Kapps":Specific_Kapps.copy(),"Default_Kapps":Default_Kapps.copy(),"Process_Efficiencies":process_efficiencies.copy()})

                KappCorrectionResults=efficiency_correction(enzyme_efficiencies=Specific_Kapps,
                                                            simulation_results=Simulation_results["Simulation_Results"],
                                                            protein_data=build_input_proteome_for_specific_kapp_estimation(proteome, condition),
                                                            rba_session=rba_session,
                                                            condition_to_look_up="Prokaryotic",
                                                            default_enzyme_efficiencies=Default_Kapps,
                                                            tolerance=None,
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
                previous_RSS=current_RSS

            if current_RSS>rss_trajectory[0]:
                continuation_criterion=False
            if steady_count>=steady_limit:
                continuation_criterion=False
            elif iteration_count>=iteration_limit:
                continuation_criterion=False
            elif increasing_RSS_count>=increasing_RSS_limit:
                continuation_criterion=False

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
                                                             condition_to_look_up="Prokaryotic",
                                                             growth_rate_to_look_up="Mu_prok",
                                                             results_to_look_up="Simulation_Results",
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
                                                         condition_to_look_up="Prokaryotic",
                                                         growth_rate_to_look_up="Mu_prok",
                                                         results_to_look_up="Simulation_Results",
                                                         fixed_mu_when_above_target_mu_in_correction=False)
            
        Simulation_results=results_global_scaling["simulation_results"]
        Specific_Kapps_to_return=results_global_scaling["specific_kapps"]
        Default_Kapps_to_return=results_global_scaling["default_kapps"]
        process_efficiencies_to_return=results_global_scaling["process_efficiencies"]
        rss_trajectory=None

    #if use_mean_enzyme_composition_for_calibration:
        #rba_session.reload_model()
        #if enzyme_efficiency_estimation_settings['impose_on_all_isoreactions']:
        #elif enzyme_efficiency_estimation_settings['zero_on_all_isoreactions']:
        #os.remove('{}/ModelStructure_meancompo_{}.json'.format(rba_session.xml_dir,condition))

    if print_outputs:
        print("")
        print("Runtime - {} : {}".format(condition,time.time() - t0))
        print("")

    return({"RSS_trajectory":rss_trajectory,
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
                              adjust_root=True):
    
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
                                             print_output=False,
                                             apply_model=False,
                                             transporter_multiplier=transporter_multiplier,
                                             start_val=0,
                                             Mu_approx_precision=mu_approx_precision,
                                             max_mu_in_dichotomy=2*mu_measured)

    mumax_predicted=simulation_results[growth_rate_to_look_up]
    print("Global start - {}:{}".format(condition,mumax_predicted))
    predicted_growth_rates=[mumax_predicted]

    if mumax_predicted == 0:
        mu_misprediction_factor=10
    else:
        mu_misprediction_factor=mu_measured/mumax_predicted
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
                                                    print_output=False,
                                                    apply_model=False,
                                                    transporter_multiplier=transporter_multiplier,
                                                    start_val=mu_measured,
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
                print("Measured: {} - Predicted: {} - mispred coeff: {} - root: {} - runs_of_sign: {} -- status: {}".format(mu_measured,mumax_predicted,mu_misprediction_factor,n_th_root_mispred,runs_of_sign,simulation_results["SolutionStatus"]))

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
                                                        print_output=False,
                                                        apply_model=False,
                                                        transporter_multiplier=transporter_multiplier)
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
                                                    print_output=False,
                                                    apply_model=False,
                                                    transporter_multiplier=transporter_multiplier,
                                                    start_val=mu_measured,
                                                    Mu_approx_precision=mu_approx_precision,
                                                    max_mu_in_dichotomy=2*mu_measured)
            mumax_predicted=simulation_results[growth_rate_to_look_up]
            if print_outputs:
                print("Measured: {} - Predicted: {} - mispred coeff: {} - root: {} - runs_of_sign: {}-- status: {}".format(mu_measured,mumax_predicted,mu_misprediction_factor,n_th_root_mispred,runs_of_sign,simulation_results["SolutionStatus"]))

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
                                                print_output=False,
                                                apply_model=False,
                                                transporter_multiplier=transporter_multiplier,
                                                start_val=mu_measured,
                                                Mu_approx_precision=mu_approx_precision,
                                                max_mu_in_dichotomy=2*mu_measured)

    return({"specific_kapps":specific_kapps_out,"default_kapps":default_kapps_out,"process_efficiencies":process_efficiencies_out,"correction_factor":best_cumulative_correction_factor,"simulation_results":simulation_results})


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


def determine_calibration_flux_distribution(rba_session,
                                            mu,
                                            flux_bounds,
                                            biomass_function,
                                            target_biomass_function,
                                            parsimonious_fba,
                                            rxns_to_ignore_when_parsimonious
                                            ):
    rba_session.rebuild_from_model()
    rba_session.set_medium(rba_session.Medium)
    rba_session.add_exchange_reactions()
    rba_session.set_growth_rate(mu)
    if target_biomass_function:
        #original_density_constraint_signs=rba_session.Problem.get_constraint_types(constraints=[i for i in rba_session.get_density_constraints() if i in rba_session.Problem.LP.row_names])
        original_medium = copy.deepcopy(rba_session.Medium)
        rba_session.set_medium({i:100.0 for i in original_medium.keys()})
        #rba_session.Problem.set_constraint_types({i:"E" for i in rba_session.get_density_constraints() if i in rba_session.Problem.LP.row_names})
        derive_bm_from_rbasolution=False
        solved=rba_session.solve()
        if solved:
            derive_bm_from_rbasolution=True
        else:
            #rba_session.Problem.set_constraint_types(original_density_constraint_signs)
            solved2=rba_session.solve()
            if solved2:
                derive_bm_from_rbasolution=True
        rba_session.set_medium(original_medium)
        rba_session.build_fba_model(rba_derived_biomass_function=True,from_rba_solution=derive_bm_from_rbasolution)
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

    if parsimonious_fba:
        rba_session.FBA.parsimonise(rxns_to_ignore_in_objective=rxns_to_ignore_when_parsimonious)
        rba_session.FBA.set_lb(rxn_LBs)
        rba_session.FBA.set_ub(rxn_UBs)
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


def determine_machinery_concentration(rba_session,machinery_composition,proteomicsData,proto_proteins=False):
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


def pre_select_iso_reactions(measured_proteins_reaction_map,rba_session,chose_most_quantified,keep_isorxns_specific_to_quantified_proteins=False):
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
            unique_subunit_dict[irx] = len(list(rba_session.ModelStructure.EnzymeInfo.Elements[enzyme]['Subunits'].keys()))
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


def determine_reactions_associated_with_measured_proto_protein(measured_proteins_isoform_map,rba_session):
    # identify all model reactions, associated with the measured proteins
    out = {}
    for p_ID in measured_proteins_isoform_map.keys():
        reactions_associated_with_proto_protein = []
        for isoform in measured_proteins_isoform_map[p_ID]:
            for reaction in rba_session.ModelStructure.ProteinInfo.Elements[isoform]['associatedReactions']:
                reactions_associated_with_proto_protein.append(reaction)
        out[p_ID] = list(set(reactions_associated_with_proto_protein))
    return(out)


def estimate_specific_enzyme_efficiencies(rba_session, 
                                          proteomicsData, 
                                          flux_bounds, 
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
                                          rxns_to_ignore_when_parsimonious=[]):
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
                                                             biomass_function=biomass_function,
                                                             target_biomass_function=target_biomass_function,
                                                             parsimonious_fba=parsimonious_fba,
                                                             rxns_to_ignore_when_parsimonious=rxns_to_ignore_when_parsimonious
                                                             )
    FluxDistribution.to_csv('Calib_FluxDist_'+condition+'_.csv', sep=';')

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
                                                     keep_isorxns_specific_to_quantified_proteins=False)
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
                model_enzyme_concentration=determine_machinery_concentration(rba_session=rba_session,
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
            respective_concentration=determine_machinery_concentration(rba_session=rba_session,
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
            overview_out.to_csv('SpecKapp_Network_overview_'+condition+'_.csv', sep=';')
        else:
            overview_out.to_csv('SpecKapp_Network_overview_.csv', sep=';')
    # 9: ...#
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

