import pandas
import numpy
import json
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.linear_model import LinearRegression

from rbatools.regression_utils import  do_regression , lin_predictions , quad_predictions


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



def extract_predicted_exchange_fluxes(inputs=[],result_object='Simulation_Results', run='Prokaryotic', metabolite=''):
    out = []
    for i in inputs:
        try:
            out.append(abs(i[result_object]['ExchangeFluxes'].loc[metabolite, run]))
        except:
            out.append(numpy.nan)
    return(out)



def extract_predicted_growth_rates(inputs=[],result_object='Simulation_Results', run='Prokaryotic'):
    out = []
    for i in inputs:
        try:
            out.append(i[result_object]['Mu'].loc['Mu', run])
        except:
            out.append(numpy.nan)
    return(out)



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



def generate_proto_protein_map(rba_session):
    protomap={}
    for i in rba_session.get_proteins():
        protoID=rba_session.get_protein_information(protein=i)["ProtoID"]
        if protoID in protomap.keys():
            protomap[protoID].append(i)
        else:
            protomap[protoID]=[i]
    return(protomap)



def plot_predicted_fluxes(simulation_outputs,types=['Fixed_PG_Eukaryotic_fixed_sizes','Fixed_PG_Eukaryotic','Eukaryotic_fixed_sizes',"DefaultKapp","Prokaryotic","Eukaryotic"],input_definition=None):

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
                Glc_Hackett_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_Mean')&(input_definition['ID'] == 'R_EX_glc__D_e'), condition].values[0])
            except:
                Glc_Hackett_list.append(numpy.nan)
            try:
                EtOH_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'Flux_Range_Mean')&(input_definition['ID'] == 'R_EX_etoh_e'), condition].values[0])
            except:
                EtOH_Hackett_list.append(numpy.nan)
            try:
                Ac_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'Flux_Range_Mean')&(input_definition['ID'] == 'R_EX_ac_e'), condition].values[0])
            except:
                Ac_Hackett_list.append(numpy.nan)
            try:
                Acald_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'Flux_Range_Mean')&(input_definition['ID'] == 'R_EX_acald_e'), condition].values[0])
            except:
                Acald_Hackett_list.append(numpy.nan)
            try:
                Glyc_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'Flux_Range_Mean')&(input_definition['ID'] == 'R_EX_glyc_e'), condition].values[0])
            except:
                Glyc_Hackett_list.append(numpy.nan)
            try:
                Lac_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'Flux_Range_Mean')&(input_definition['ID'] == 'R_EX_lac__D_e'), condition].values[0])
            except:
                Lac_Hackett_list.append(numpy.nan)
            try:
                Succ_Hackett_list.append(-input_definition.loc[(input_definition['Type'] == 'Flux_Range_Mean')&(input_definition['ID'] == 'R_EX_succ_e'), condition].values[0])
            except:
                Succ_Hackett_list.append(numpy.nan)

            try:
                Glc_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_StandardError')&(input_definition['ID'] == 'R_EX_glc__D_e_SE'), condition].values[0])
            except:
                Glc_Hackett_SE_list.append(numpy.nan)
            try:
                EtOH_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_StandardError')&(input_definition['ID'] == 'R_EX_etoh_e_SE'), condition].values[0])
            except:
                EtOH_Hackett_SE_list.append(numpy.nan)
            try:
                Ac_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_StandardError')&(input_definition['ID'] == 'R_EX_ac_e_SE'), condition].values[0])
            except:
                Ac_Hackett_SE_list.append(numpy.nan)
            try:
                Acald_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_StandardError')&(input_definition['ID'] == 'R_EX_acald_e_SE'), condition].values[0])
            except:
                Acald_Hackett_SE_list.append(numpy.nan)
            try:
                Glyc_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_StandardError')&(input_definition['ID'] == 'R_EX_glyc_e_SE'), condition].values[0])
            except:
                Glyc_Hackett_SE_list.append(numpy.nan)
            try:
                Lac_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_StandardError')&(input_definition['ID'] == 'R_EX_lac__D_e_SE'), condition].values[0])
            except:
                Lac_Hackett_SE_list.append(numpy.nan)
            try:
                Succ_Hackett_SE_list.append(input_definition.loc[(input_definition['Type'] == 'Flux_Range_StandardError')&(input_definition['ID'] == 'R_EX_succ_e_SE'), condition].values[0])
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
    Mus_predicted_euk_fixed = extract_predicted_growth_rates(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes')
    Mus_predicted_fixed_pg_euk = extract_predicted_growth_rates(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic')
    Mus_predicted_fixed_pg_euk_fixed = extract_predicted_growth_rates(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes')

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

    Glc_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_glc__D')
    EtOH_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_etoh')
    O2_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_o2')
    CO2_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_co2')
    Ac_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_ac')
    Glycerol_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_glyc')
    Acald_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_acald')
    Lac_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_lac__D')
    Succ_Exchange_predicted_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_Euk_fixed', run='Eukaryotic_fixed_sizes', metabolite='M_succ')

    Glc_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_glc__D')
    EtOH_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_etoh')
    O2_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_o2')
    CO2_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_co2')
    Ac_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_ac')
    Glycerol_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_glyc')
    Acald_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_acald')
    Lac_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_lac__D')
    Succ_Exchange_predicted_fixed_pg_euk = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk', run='Fixed_PG_Eukaryotic', metabolite='M_succ')

    Glc_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_glc__D')
    EtOH_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_etoh')
    O2_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_o2')
    CO2_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_co2')
    Ac_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_ac')
    Glycerol_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_glyc')
    Acald_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_acald')
    Lac_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_lac__D')
    Succ_Exchange_predicted_fixed_pg_euk_fixed = extract_predicted_exchange_fluxes(inputs=simulation_outputs,result_object='Simulation_Results_fixed_pg_Euk_fixed', run='Fixed_PG_Eukaryotic_fixed_sizes', metabolite='M_succ')


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

    Glc_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_glc__D_e',bound_type="Min")
    EtOH_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_etoh_e',bound_type="Min")
    O2_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_o2_e',bound_type="Min")
    Ac_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_ac_e',bound_type="Min")
    Glycerol_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_glyc_e',bound_type="Min")
    Acald_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_acald_e',bound_type="Min")
    Lac_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_lac__D_e',bound_type="Min")
    Succ_VAmin_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_succ_e',bound_type="Min")
    Glc_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_glc__D_e',bound_type="Max")
    EtOH_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_etoh_e',bound_type="Max")
    O2_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_o2_e',bound_type="Max")
    Ac_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_ac_e',bound_type="Max")
    Glycerol_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_glyc_e',bound_type="Max")
    Acald_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_acald_e',bound_type="Max")
    Lac_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_lac__D_e',bound_type="Max")
    Succ_VAmax_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_euk_fixed', variable='R_EX_succ_e',bound_type="Max")

    Glc_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_glc__D_e',bound_type="Min")
    EtOH_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_etoh_e',bound_type="Min")
    O2_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_o2_e',bound_type="Min")
    Ac_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_ac_e',bound_type="Min")
    Glycerol_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_glyc_e',bound_type="Min")
    Acald_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_acald_e',bound_type="Min")
    Lac_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_lac__D_e',bound_type="Min")
    Succ_VAmin_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_succ_e',bound_type="Min")
    Glc_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_glc__D_e',bound_type="Max")
    EtOH_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_etoh_e',bound_type="Max")
    O2_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_o2_e',bound_type="Max")
    Ac_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_ac_e',bound_type="Max")
    Glycerol_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_glyc_e',bound_type="Max")
    Acald_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_acald_e',bound_type="Max")
    Lac_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_lac__D_e',bound_type="Max")
    Succ_VAmax_fixed_pg_euk=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk', variable='R_EX_succ_e',bound_type="Max")

    Glc_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_glc__D_e',bound_type="Min")
    EtOH_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_etoh_e',bound_type="Min")
    O2_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_o2_e',bound_type="Min")
    Ac_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_ac_e',bound_type="Min")
    Glycerol_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_glyc_e',bound_type="Min")
    Acald_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_acald_e',bound_type="Min")
    Lac_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_lac__D_e',bound_type="Min")
    Succ_VAmin_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_succ_e',bound_type="Min")
    Glc_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_glc__D_e',bound_type="Max")
    EtOH_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_etoh_e',bound_type="Max")
    O2_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_o2_e',bound_type="Max")
    Ac_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_ac_e',bound_type="Max")
    Glycerol_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_glyc_e',bound_type="Max")
    Acald_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_acald_e',bound_type="Max")
    Lac_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_lac__D_e',bound_type="Max")
    Succ_VAmax_fixed_pg_euk_fixed=extract_feasible_bounds(inputs=simulation_outputs,feasible_range_object='FeasibleRange_fixed_pg_euk_fixed', variable='R_EX_succ_e',bound_type="Max")

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
    if "Eukaryotic_fixed_sizes" in types:
        axs[0, 0].scatter(Mu_Hackett, Mus_predicted_euk_fixed, color='steelblue')
        legendlist.append("Euk. (fixed)")
    if "Fixed_PG_Eukaryotic" in types:
        axs[0, 0].scatter(Mu_Hackett, Mus_predicted_fixed_pg_euk, color='green')
        legendlist.append("Euk_PG.")
    if "Fixed_PG_Eukaryotic_fixed_sizes" in types:
        axs[0, 0].scatter(Mu_Hackett, Mus_predicted_fixed_pg_euk_fixed, color='pink')
        legendlist.append("Euk_Pg. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[0, 1].scatter(Mus_predicted_euk_fixed, Glc_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
    if "Fixed_PG_Eukaryotic" in types:
        axs[0, 1].scatter(Mus_predicted_fixed_pg_euk, Glc_Exchange_predicted_fixed_pg_euk, color='green', alpha=0.8)
        legendlist.append("Euk_PG.")
    if "Fixed_PG_Eukaryotic_fixed_sizes" in types:
        axs[0, 1].scatter(Mus_predicted_fixed_pg_euk_fixed, Glc_Exchange_predicted_fixed_pg_euk_fixed, color='pink', alpha=0.8)
        legendlist.append("Euk_PG. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[0, 1].plot([Mus_predicted_euk_fixed[i]]*2,[Glc_VAmin_euk_fixed[i],Glc_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[0, 1].scatter(Mus_predicted_euk_fixed, Glc_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[0, 1].scatter(Mus_predicted_euk_fixed, Glc_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[0, 2].scatter(Mus_predicted_euk_fixed, O2_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[0, 2].plot([Mus_predicted_euk_fixed[i]]*2,[O2_VAmin_euk_fixed[i],O2_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[0, 2].scatter(Mus_predicted_euk_fixed, O2_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[0, 2].scatter(Mus_predicted_euk_fixed, O2_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[1, 0].scatter(Mus_predicted_euk_fixed, EtOH_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
    if "Fixed_PG_Eukaryotic" in types:
        axs[1, 0].scatter(Mus_predicted_fixed_pg_euk, EtOH_Exchange_predicted_fixed_pg_euk, color='green', alpha=0.8)
        legendlist.append("Euk_PG.")
    if "Fixed_PG_Eukaryotic_fixed_sizes" in types:
        axs[1, 0].scatter(Mus_predicted_fixed_pg_euk_fixed, EtOH_Exchange_predicted_fixed_pg_euk_fixed, color='pink', alpha=0.8)
        legendlist.append("Euk_PG. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[1, 0].plot([Mus_predicted_euk_fixed[i]]*2,[EtOH_VAmin_euk_fixed[i],EtOH_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[1, 0].scatter(Mus_predicted_euk_fixed, EtOH_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[1, 0].scatter(Mus_predicted_euk_fixed, EtOH_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[1, 1].scatter(Mus_predicted_euk_fixed, Ac_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[1, 1].plot([Mus_predicted_euk_fixed[i]]*2,[Ac_VAmin_euk_fixed[i],Ac_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[1, 1].scatter(Mus_predicted_euk_fixed, Ac_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[1, 1].scatter(Mus_predicted_euk_fixed, Ac_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[1, 2].scatter(Mus_predicted_euk_fixed, Glycerol_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[1, 2].plot([Mus_predicted_euk_fixed[i]]*2,[Glycerol_VAmin_euk_fixed[i],Glycerol_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[1, 2].scatter(Mus_predicted_euk_fixed, Glycerol_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[1, 2].scatter(Mus_predicted_euk_fixed, Glycerol_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[2, 0].scatter(Mus_predicted_euk_fixed, Acald_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[2, 0].plot([Mus_predicted_euk_fixed[i]]*2,[Acald_VAmin_euk_fixed[i],Acald_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[2, 0].scatter(Mus_predicted_euk_fixed, Acald_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[2, 0].scatter(Mus_predicted_euk_fixed, Acald_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[2, 1].scatter(Mus_predicted_euk_fixed, Lac_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[2, 1].plot([Mus_predicted_euk_fixed[i]]*2,[Lac_VAmin_euk_fixed[i],Lac_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[2, 1].scatter(Mus_predicted_euk_fixed, Lac_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[2, 1].scatter(Mus_predicted_euk_fixed, Lac_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
    if "Eukaryotic_fixed_sizes" in types:
        axs[2, 2].scatter(Mus_predicted_euk_fixed, Succ_Exchange_predicted_euk_fixed, color='steelblue', alpha=0.8)
        legendlist.append("Euk. (fixed)")
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
    if "Eukaryotic_fixed_sizes" in types:
        for i in range(len(Mus_predicted_euk_fixed)):
            axs[2, 2].plot([Mus_predicted_euk_fixed[i]]*2,[Succ_VAmin_euk_fixed[i],Succ_VAmax_euk_fixed[i]],color="steelblue",alpha=0.1)
        axs[2, 2].scatter(Mus_predicted_euk_fixed, Succ_VAmin_euk_fixed, color='steelblue',marker=6, alpha=0.5)
        axs[2, 2].scatter(Mus_predicted_euk_fixed, Succ_VAmax_euk_fixed, color='steelblue',marker=7, alpha=0.5)
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
        n_rows=2
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
        axs[ax_indices[count][0],ax_indices[count][1]].scatter([cal_input["RSS_trajectory"].index(min(cal_input["RSS_trajectory"]))+1] , min(cal_input["RSS_trajectory"]) , color="green",alpha=0.7)
        axs[ax_indices[count][0],ax_indices[count][1]].scatter([cal_input["Chosen_RSS_index"]+1] , cal_input["RSS_trajectory"][cal_input["Chosen_RSS_index"]] , color="red",alpha=0.7)
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
    elif type=="Eukaryotic_fixed_sizes":
        results_object="Simulation_Results_Euk_fixed"
    elif type=="DefaultKapp":
        results_object="Simulation_Results_DefKapp"


def plot_protein_protein_comparison(predicted_proteomes,measured_proteomes,conditions):
    plot_dimensions={1:(1,1),2:(1,2),3:(1,3),4:(2,3),5:(2,3),6:(2,3),7:(3,3),8:(3,3),9:(3,3)}
    plot_indices={1:(1,1),2:(1,2),3:(1,3),4:(2,1),5:(2,2),6:(2,3),7:(3,1),8:(3,2),9:(3,3)}
    number_conditions=len(conditions)
    fig, axs = plt.subplots(plot_dimensions[number_conditions][0], plot_dimensions[number_conditions][1], figsize=(28, 7), sharex=True)
    condition_count=0

    for condition in conditions:
        condition_count+=1
        fig_row=plot_indices[condition_count][0]-1
        fig_col=plot_indices[condition_count][1]-1
        #protein_comparison=pandas.DataFrame()
        if condition in list(predicted_proteomes.columns):
            if condition in list(measured_proteomes.columns):
                respective_measured_proteome=measured_proteomes.loc[pandas.isna(measured_proteomes[condition])==False]
                respective_predicted_proteome=predicted_proteomes.loc[pandas.isna(predicted_proteomes[condition])==False]
                x=[6.023e20*respective_predicted_proteome.loc[i,condition] for i in respective_predicted_proteome.index if i in respective_measured_proteome.index]
                y=[6.023e20*respective_measured_proteome.loc[i,condition] for i in respective_predicted_proteome.index if i in respective_measured_proteome.index]

                regression_results=do_linear_regression_on_proteome_prediction(x=x,
                                                                               y=y,
                                                                               fit_intercept=False)
                x_reg=regression_results['X_regression']
                #y_reg=regression_results['Y_regression']
                regressor=regression_results['Regressor']
                predictions=regression_results['Prediction']

                total_max=max([max(x),max(y)])
                total_min=min([min(x),min(y)])
                #[numpy.log10(i) for i in list(protein_comparison_to_prodeed["Predicted"])]
                axs[fig_row, fig_col].plot([numpy.log10(total_min), numpy.log10(total_max)], [numpy.log10(total_min), numpy.log10(total_max)], color='green', linewidth=3,alpha=0.6)
                axs[fig_row, fig_col].scatter([numpy.log10(i) for i in x],[numpy.log10(i) for i in y],alpha=0.4)
                axs[fig_row, fig_col].plot(numpy.log10(x_reg), numpy.log10(predictions), color='red')
                #axs[fig_row, fig_col].plot(x_reg, predictions, color='red')
                axs[fig_row, fig_col].legend(['Identity', "Correlation: {}".format(str(round(regressor.coef_[0][0],2))), 'Data'])
                axs[fig_row, fig_col].plot([numpy.log10(total_min), numpy.log10(total_max)], [numpy.log10(total_min)+2, numpy.log10(total_max)+2], color='black', linewidth=1,alpha=0.6)
                axs[fig_row, fig_col].plot([numpy.log10(total_min), numpy.log10(total_max)], [numpy.log10(total_min)-2, numpy.log10(total_max)-2], color='black', linewidth=1,alpha=0.6)
                #axs[fig_row, fig_col].set_title("Predicted protein numbers per gram dry weight - {}".format(condition))
                axs[fig_row, fig_col].set_title(condition)
                axs[fig_row, fig_col].set_xlabel('Predicted copies per $g_{DW}$ ($Log_{10}$)')
                axs[fig_row, fig_col].set_ylabel('Measured copies per $g_{DW}$ ($Log_{10}$)')
                axs[fig_row, fig_col].set_ylim(math.floor(numpy.log10(min(y))),math.ceil(numpy.log10(max(y))))

    plt.show()


def do_linear_regression_on_proteome_prediction(x,y,fit_intercept):
    x_reg = numpy.reshape(numpy.array(x), (len(x), 1))
    y_reg = numpy.reshape(numpy.array(y), (len(y), 1))
    regressor = LinearRegression(fit_intercept=fit_intercept)
    regressor.fit(x_reg, y_reg)
    predictions = regressor.predict(x_reg)
    return({'Regressor':regressor,
            'Prediction':predictions,
            'X_regression':x_reg,
            'Y_regression':y_reg})

def regression_on_specific_enzyme_efficiencies(Spec_Kapps,min_kapp,max_kapp,conditions,growth_rates,only_lin=False,impose_on_isoenzymes=True,monotonous_quadratic=False,fill_in_missing_conditions=False,permit_quadratic_model=True):
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
            regression_results=do_regression(x_to_fit=x_s,y_to_fit=y_s,x_to_plot=x_to_plot,max_val=max_kapp,min_val=min_kapp,monotonous_quadratic=monotonous_quadratic,total_x_range=(min(x),max(x)),permit_quadratic_model=permit_quadratic_model)
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


def regression_on_compartment_sizes(Comp_sizes,conditions,growth_rates,monotonous_quadratic=False,permit_quadratic_model=True,plotting=True):
    if plotting:
        pp = PdfPages('CompSize_Plots_refined3.pdf')
    out=pandas.DataFrame()
    for comp in Comp_sizes.index:
        x=[growth_rates[i] for i in conditions]
        y=list(Comp_sizes.loc[comp,conditions])
        x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
        regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=1,min_val=0,monotonous_quadratic=monotonous_quadratic,total_x_range=(0,1),permit_quadratic_model=permit_quadratic_model)
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
        if plotting:
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
    if plotting:
        pp.close()
    return(out)


def regression_on_pg_fractions(PG_sizes,conditions,growth_rates,monotonous_quadratic=False,permit_quadratic_model=True,plotting=True):
    if plotting:
        pp = PdfPages('PGfraction_Plots_refined3.pdf')
    out=pandas.DataFrame()
    for comp in PG_sizes.index:
        x=[growth_rates[i]for i in conditions]
        y=list(PG_sizes.loc[comp,conditions])

        x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
        regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=1,min_val=0,monotonous_quadratic=monotonous_quadratic,total_x_range=(0,1),permit_quadratic_model=permit_quadratic_model)
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
        if plotting:
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
    if plotting:
        pp.close()
    return(out)


def regression_on_process_efficiencies(Process_efficiencies,min_efficiency,max_efficiency,conditions,growth_rates,monotonous_quadratic=False,permit_quadratic_model=True):
    pp = PdfPages('ProcessEfficiencies_Plots_refined3.pdf')
    out=pandas.DataFrame()
    for process in Process_efficiencies.index:
        x=[growth_rates[i]for i in conditions]
        y=list(Process_efficiencies.loc[process,conditions])

        x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
        regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=max_efficiency,min_val=min_efficiency,monotonous_quadratic=monotonous_quadratic,total_x_range=(min(x),max(x)),permit_quadratic_model=permit_quadratic_model)
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


def regression_on_default_enzyme_efficiencies(default_kapps,min_kapp,max_kapp,conditions,growth_rates,monotonous_quadratic=False,permit_quadratic_model=True):
    pp = PdfPages('DefKapp_Plots_refined3.pdf')
    out=pandas.DataFrame()
    x=[growth_rates[i]for i in conditions]
    y=list(default_kapps.loc[conditions,"Default Kapp"])

    x_to_plot=(numpy.linspace(x[0], x[-1], 1000))
    regression_results=do_regression(x_to_fit=x,y_to_fit=y,x_to_plot=x_to_plot,max_val=max_kapp,min_val=min_kapp,monotonous_quadratic=monotonous_quadratic,total_x_range=(min(x),max(x)),permit_quadratic_model=permit_quadratic_model)
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

