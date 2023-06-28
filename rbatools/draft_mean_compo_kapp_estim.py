
def pre_select_iso_reactions(measured_Proteins_Reaction_Map,rba_session,chose_most_quantified):
    # choose most likely iso-reaction for each measured-protein associated reaction
    protoRxnDict = {}
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
    out = {}
    for prx in protoRxnDict.keys():
        unique_SU_dict={}
        for irx in protoRxnDict[prx].keys():
            enzyme = rba_session.ModelStructure.ReactionInfo.Elements[irx]['Enzyme']
            unique_SU_dict[irx] = len(list(rba_session.ModelStructure.EnzymeInfo.Elements[enzyme]['Subunits'].keys()))
        if chose_most_quantified:
            max_val = max([protoRxnDict[prx][i]/unique_SU_dict[i] for i in protoRxnDict[prx].keys()])
            list_isorxns = [i for i in protoRxnDict[prx].keys() if protoRxnDict[prx][i]/unique_SU_dict[i] == max_val]
            if len(list_isorxns)>1:
                max_SU_number=max([unique_SU_dict[i] for i in list_isorxns])
                selected=[i for i in list_isorxns if unique_SU_dict[i]==max_SU_number]
        else:
            selected = [i for i in protoRxnDict[prx].keys() if protoRxnDict[prx][i] != 0]
        selected.sort()
        out[prx] = selected
    return(out)

def determine_reactions_associated_with_measured_proto_protein(measured_Proteins_Isoform_Map,rba_session):
    # identify all model reactions, associated with the measured proteins
    out = {}
    for p_ID in measured_Proteins_Isoform_Map.keys():
        reactions_associated_with_proto_Protein = []
        for isoform in measured_Proteins_Isoform_Map[p_ID]:
            for reaction in rba_session.ModelStructure.ProteinInfo.Elements[isoform]['associatedReactions']:
                reactions_associated_with_proto_Protein.append(reaction)
        out[p_ID] = list(set(reactions_associated_with_proto_Protein))
    return(out)

def estimate_specific_enzyme_efficiencies_refactored(rba_session, 
                                                     proteomicsData, 
                                                     flux_bounds, 
                                                     mu, 
                                                     biomass_function=None, 
                                                     target_biomass_function=True, 
                                                     parsimonious_fba=True, 
                                                     chose_most_quantified=True, 
                                                     equalize_identical_enzymes=True, 
                                                     impose_on_all_isoreactions=False, 
                                                     zero_on_all_isoreactions=True, 
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

    FluxDistribution=determine_calibration_flux_distribution(rba_session=rba_session,
                                                             mu=mu,
                                                             flux_bounds=flux_bounds,
                                                             biomass_function=biomass_function,
                                                             target_biomass_function=target_biomass_function,
                                                             parsimonious_fba=parsimonious_fba,
                                                             rxns_to_ignore_when_parsimonious=rxns_to_ignore_when_parsimonious
                                                             )
    FluxDistribution.to_csv('Calib_FluxDist_'+condition+'_.csv', sep=';')

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
    measured_Proteins_Reaction_Map = determine_reactions_associated_with_measured_proto_protein(measured_Proteins_Isoform_Map=measured_Proteins_Isoform_Map,
                                                                                                rba_session=rba_session)

    chosen_Isoreactions=pre_select_iso_reactions(measured_Proteins_Reaction_Map=measured_Proteins_Reaction_Map,
                                                 rba_session=rba_session,
                                                 chose_most_quantified=chose_most_quantified)

    # determine model reactions with non-zero flux in FBA#
    overview_out = pandas.DataFrame()
    for i in list(chosen_Isoreactions.keys()):
        all_iso_rxns = chosen_Isoreactions[i]
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
        ### ??? ###
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
                            identical_set.append(j)
                            already_handled.append(j)
                    identical_enzymes.append(identical_set)

        ident_set_overview = {}
        for ident_set in identical_enzymes:
            ## ??? ###
            enzyme_concentration = list(set([overview_out.loc[overview_out['Enzyme_ID'] == enz, 'CopyNumber'].values[0] for enz in ident_set if enz in list(overview_out['Enzyme_ID'])]))[0]
            #enzyme_concentration = sum(list(set([overview_out.loc[overview_out['Enzyme_ID'] == enz, 'CopyNumber'].values[0] for enz in ident_set if enz in list(overview_out['Enzyme_ID'])])))
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

