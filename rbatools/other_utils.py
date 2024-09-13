import pandas

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
                    out.loc[i, 'Absolute_Reference_restored'] = CopyNumber_match/(2**FoldChange_match)
    for gene in list(out['ID']):
        Abs_Ref = out.loc[gene, 'Absolute_Reference_restored']
        for condition in conditions_in_fold_change_data_to_restore:
            if condition in list(fold_changes.columns):
                out.loc[gene, condition] = Abs_Ref * (2**fold_changes.loc[fold_changes['Gene'] == gene, condition].values[0])
    return(out)


def add_annotations_to_proteome(input, ID_column, annotations):
    for i in input.index:
        if input.loc[i, ID_column] in annotations.index:
            input.loc[i, 'AA_residues'] = annotations.loc[input.loc[i, ID_column], 'AA_residues']
            input.loc[i, 'Location'] = annotations.loc[input.loc[i, ID_column], 'Location']
            input.loc[i, 'InModel'] = annotations.loc[input.loc[i, ID_column], 'InModel']
            input.loc[i, 'IsRibosomal'] = annotations.loc[input.loc[i, ID_column], 'IsRibosomal']
    return(input)


def extract_proteomes_from_calibration_results(calib_results):
    proteomes=pandas.DataFrame()
    for calib_result in calib_results:
        condition=calib_result["Condition"]
        proteome=calib_result["Proteome"]
        for i in list(proteome.index):
            proteomes.loc[i,condition]=proteome.loc[i,"copy_number"]
    return(proteomes)


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
