import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *
import time
from multiprocessing import Pool , cpu_count
from joblib import Parallel, delayed
import math

def apply_single_knock_outs(input_dict):
    genes_to_KO=input_dict["genes_to_KO"]
    Mumax_WT=input_dict["Mumax_WT"]
    condition=input_dict["condition"]
    Input_Data=input_dict["input_data"]
    rba_session = SessionRBA(input_dict["xml_dir"])
    rba_session.add_exchange_reactions()

    compartment_sizes_from_calibration=input_dict["compartment_sizes_from_calibration"]
    pg_fractions_from_calibration=input_dict["pg_fractions_from_calibration"]
    process_efficiencies_from_calibration=input_dict["process_efficiencies_from_calibration"]
    default_kapps_from_calibration=input_dict["default_kapps_from_calibration"]
    specific_kapps_from_calibration=input_dict["specific_kapps_from_calibration"]
    transporter_multiplier=input_dict["transporter_multiplier"]

    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes_from_calibration,pg_fractions=pg_fractions_from_calibration,condition=condition)
    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
    # Process efficiencies & Def/Spec Kapps
    process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies_from_calibration,condition=condition,parameter_name_suffix='_apparent_efficiency')
    Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=default_kapps_from_calibration,condition=condition,transporter_multiplier=transporter_multiplier)
    Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=specific_kapps_from_calibration,condition=condition)

    inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
    rba_session.rebuild_from_model()
    rba_session.set_medium(medium_concentrations_from_input(input=Input_Data, condition=condition))

    single_ko_mus={}
    for gene in genes_to_KO:
        try:
            rba_session.apply_gene_knock_out(gene=gene)
            rba_session.set_growth_rate(Mumax_WT)
            wt_feasible=rba_session.solve()
            if wt_feasible:
                single_ko_mus.update({gene:Mumax_WT})
            else:
                rba_session.set_growth_rate(0.0)
                zero_feasible=rba_session.solve()
                if zero_feasible:
                    single_ko_mus.update({gene:rba_session.find_max_growth_rate(max_value=Mumax_WT)})
                else:
                    single_ko_mus.update({gene:0.0})
            rba_session.undo_gene_knock_out(gene=gene)            
        except:
            single_ko_mus.update({gene:numpy.nan})
        print("{}: {} / {}".format(gene,single_ko_mus[gene],Mumax_WT))
    return(single_ko_mus)

def apply_double_knock_outs(input_dict):
    genes_to_KO=input_dict["genes_to_KO"]
    Mumax_WT=input_dict["Mumax_WT"]
    condition=input_dict["condition"]
    Input_Data=input_dict["input_data"]
    Singele_KO_results=input_dict["Single_KO_results"]
    rba_session = SessionRBA(input_dict["xml_dir"])
    rba_session.add_exchange_reactions()

    compartment_sizes_from_calibration=input_dict["compartment_sizes_from_calibration"]
    pg_fractions_from_calibration=input_dict["pg_fractions_from_calibration"]
    process_efficiencies_from_calibration=input_dict["process_efficiencies_from_calibration"]
    default_kapps_from_calibration=input_dict["default_kapps_from_calibration"]
    specific_kapps_from_calibration=input_dict["specific_kapps_from_calibration"]
    transporter_multiplier=input_dict["transporter_multiplier"]

    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes_from_calibration,pg_fractions=pg_fractions_from_calibration,condition=condition)
    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
    # Process efficiencies & Def/Spec Kapps
    process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies_from_calibration,condition=condition,parameter_name_suffix='_apparent_efficiency')
    Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=default_kapps_from_calibration,condition=condition,transporter_multiplier=transporter_multiplier)
    Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=specific_kapps_from_calibration,condition=condition)

    inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
    rba_session.rebuild_from_model()
    rba_session.set_medium(medium_concentrations_from_input(input=Input_Data, condition=condition))

    double_ko_mus={}
    for gene_pair in genes_to_KO:
        gene="__+__".join(gene_pair)
        try:
            #print({gene_pair[0]:float(Singele_KO_results.loc[gene_pair[0],"Mu_max_KO"]),gene_pair[1]:float(Singele_KO_results.loc[gene_pair[1],"Mu_max_KO"])})
            min_mu_of_single_KO=min([float(Singele_KO_results.loc[gene_pair[0],"Mu_max_KO"]),float(Singele_KO_results.loc[gene_pair[1],"Mu_max_KO"])])
            rba_session.apply_gene_knock_out(gene=gene_pair[0])
            rba_session.apply_gene_knock_out(gene=gene_pair[1])
            #rba_session.set_growth_rate(Mumax_WT)
            #wt_feasible=rba_session.solve()
            #if wt_feasible:
            #    double_ko_mus.update({gene:Mumax_WT})
            rba_session.set_growth_rate(min_mu_of_single_KO)
            min_mu_of_single_KO_feasible=rba_session.solve()
            if min_mu_of_single_KO_feasible:
                double_ko_mus.update({gene:min_mu_of_single_KO})
            else:
                rba_session.set_growth_rate(0.0)
                zero_feasible=rba_session.solve()
                if zero_feasible:
                    #double_ko_mus.update({gene:rba_session.find_max_growth_rate(max_value=Mumax_WT)})
                    double_ko_mus.update({gene:rba_session.find_max_growth_rate(max_value=min_mu_of_single_KO)})
                else:
                    double_ko_mus.update({gene:0.0})
            rba_session.undo_gene_knock_out(gene=gene_pair[0])            
            rba_session.undo_gene_knock_out(gene=gene_pair[1])            
        except:
            double_ko_mus.update({gene:numpy.nan})
        print("{}: {} / {}".format(gene,double_ko_mus[gene],Mumax_WT))
    return(double_ko_mus)

def main(number_chunks=2,transporter_multiplier=1,condition='Hackett_C01',xml_dir='../Yeast_iMM904_RBA_model',n_parallel_processes=None,compartments_to_consider=None,double_KOs=False,single_KOs_to_import=None,file_name="",threshold_for_zero_mu=None):
    Input_Data = pandas.read_csv('../DataSetsYeastRBACalibration/Calibration_InputDefinition_plus_Nlim_Frick_fluxes.csv', sep=';', decimal=',', index_col=0)

    default_kapps_from_calibration=pandas.read_csv("../default_kapps_refactored_WF.csv",index_col=0)
    specific_kapps_from_calibration=pandas.read_csv("../specific_kapps_refactored_WF.csv",index_col=0)
    process_efficiencies_from_calibration=pandas.read_csv("../process_efficiencies_refactored_WF.csv",index_col=0)
    pg_fractions_from_calibration=pandas.read_csv("../pg_fractions_refactored_WF.csv",index_col=0)
    compartment_sizes_from_calibration = pandas.read_csv("../compartment_sizes_refactored_WF.csv",index_col=0)

    ##########################################################################################

    rba_session = SessionRBA(xml_dir)
    rba_session.add_exchange_reactions()

    compartment_densities_and_PGs=generate_compartment_size_and_pg_input(compartment_sizes=compartment_sizes_from_calibration,pg_fractions=pg_fractions_from_calibration,condition=condition)
    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
        rba_session.model.parameters.functions._elements_by_id[str('fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density'].values[0]
        rba_session.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction'].values[0]
    # Process efficiencies & Def/Spec Kapps
    process_efficiencies_to_inject=generate_process_efficiency_input(process_efficiencies=process_efficiencies_from_calibration,condition=condition,parameter_name_suffix='_apparent_efficiency')
    Default_Kapps_to_inject=generate_default_kapp_input(default_kapps=default_kapps_from_calibration,condition=condition,transporter_multiplier=transporter_multiplier)
    Specific_Kapps_to_inject=generate_specific_kapp_input(specific_kapps=specific_kapps_from_calibration,condition=condition)

    inject_estimated_efficiencies_into_model(rba_session, specific_kapps=Specific_Kapps_to_inject, default_kapps=Default_Kapps_to_inject, process_efficiencies=process_efficiencies_to_inject)
    rba_session.rebuild_from_model()
    print(medium_concentrations_from_input(input=Input_Data, condition=condition))
    rba_session.set_medium(medium_concentrations_from_input(input=Input_Data, condition=condition))

    Mumax_WT = rba_session.find_max_growth_rate()

    if threshold_for_zero_mu is not None:
        zero_mu_multiplier=threshold_for_zero_mu
    else:
        zero_mu_multiplier=0


    if single_KOs_to_import is None:
        enzymes_with_zero_kapp=list(Specific_Kapps_to_inject.loc[Specific_Kapps_to_inject["Kapp"]==0,"Enzyme_ID"])
        all_metabolic_genes=[]
        for enzyme in rba_session.get_enzymes():
            if enzyme not in enzymes_with_zero_kapp:
                SUs=list(rba_session.get_enzyme_information(enzyme=enzyme)["Subunits"].keys())
                for su in SUs:
                    comp=rba_session.get_protein_information(protein=su)["Compartment"]
                    if su not in all_metabolic_genes:
                        if not su.startswith("average_protein"):
                            if compartments_to_consider is None:
                                all_metabolic_genes.append(su)
                            else:
                                for compartment in compartments_to_consider:
                                    if comp==compartment:
                                        all_metabolic_genes.append(su)

        #metabolic_genes=all_metabolic_genes[:100]
        metabolic_genes=all_metabolic_genes
        if number_chunks is not None:
            if number_chunks==1:
                gene_lists_to_KO=[metabolic_genes]
            else:
                gene_lists_to_KO=[metabolic_genes[i:i+math.ceil(len(metabolic_genes)/number_chunks)] for i in range(0,len(metabolic_genes),math.ceil(len(metabolic_genes)/number_chunks))]
        else:
            num_cores=cpu_count()
            if num_cores==1:
                gene_lists_to_KO=[metabolic_genes]
            else:
                gene_lists_to_KO=[metabolic_genes[i:i+math.ceil(len(metabolic_genes)/num_cores)] for i in range(0,len(metabolic_genes),math.ceil(len(metabolic_genes)/num_cores))]

        print("{} single KOs to perform".format(len(metabolic_genes)))

        input_dicts=[]
        for gene_list_to_KO in gene_lists_to_KO:
            dict_to_add={"xml_dir":xml_dir,
                        "genes_to_KO":gene_list_to_KO,
                        "condition":condition,
                        "input_data":Input_Data,
                        "compartment_sizes_from_calibration":compartment_sizes_from_calibration,
                        "pg_fractions_from_calibration":pg_fractions_from_calibration,
                        "process_efficiencies_from_calibration":process_efficiencies_from_calibration,
                        "default_kapps_from_calibration":default_kapps_from_calibration,
                        "specific_kapps_from_calibration":specific_kapps_from_calibration,
                        "transporter_multiplier":1,
                        "Mumax_WT":Mumax_WT}
            input_dicts.append(dict_to_add)

        if n_parallel_processes!=1:
            if n_parallel_processes is None:
                num_cores=cpu_count()
                n_jobs=min(num_cores,len(input_dicts))
            else:
                n_jobs=n_parallel_processes
            print("Jobs:{}".format(n_jobs))

            pool=Pool(n_jobs)

            #simulation_results_pool=pool.map_async(simulate_run,input_dicts)
            #simulation_results=simulation_results_pool.get()

            single_KO_results=pool.imap_unordered(apply_single_knock_outs,input_dicts)

            #simulation_results=Parallel(n_jobs=n_jobs)(delayed(simulate_run)(input_dict) for input_dict in input_dicts)
        else:
            single_KO_results=[]
            for input_dict in input_dicts:
                single_KO_results.append(apply_single_knock_outs(input_dict))
        
        output_single=pandas.DataFrame()
        for i in single_KO_results:
            for j in i.keys():
                output_single.loc[j,"ID"]=j
                output_single.loc[j,"Mu_max_KO"]=i[j]
                output_single.loc[j,"Mu_max_WT"]=Mumax_WT
                output_single.loc[j,"Relative_Mu"]=i[j]/Mumax_WT
    
        output_single.to_csv("Single_Gene_KOs_{}.csv".format(file_name),sep=";",decimal=",")
    else:
        print("Importing single KO data from: {}".format(single_KOs_to_import))
        output_single=pandas.read_csv(single_KOs_to_import,sep=";",decimal=",",index_col=0)

    if double_KOs:
        genes_for_double_KOs=list(output_single.loc[output_single["Relative_Mu"]>zero_mu_multiplier,"ID"])
        gene_pairs={}
        #count=0
        for i in genes_for_double_KOs:
            for j in genes_for_double_KOs:
                if i!=j:
                    if "{}__+__{}".format(j,i) not in list(gene_pairs.keys()):
                        gene_pairs["{}__+__{}".format(i,j)]=[i,j]
                        #count+=1
                #if count>=80:
                #    break
            #if count>=80:
            #    break
        gene_pair_list=list(gene_pairs.values())
        if number_chunks is not None:
            if number_chunks==1:
                doublegene_lists_to_KO=[gene_pair_list]
            else:
                doublegene_lists_to_KO=[gene_pair_list[i:i+math.ceil(len(gene_pair_list)/number_chunks)] for i in range(0,len(gene_pair_list),math.ceil(len(gene_pair_list)/number_chunks))]
        else:
            num_cores=cpu_count()
            if num_cores==1:
                doublegene_lists_to_KO=[gene_pair_list]
            else:
                doublegene_lists_to_KO=[gene_pair_list[i:i+math.ceil(len(gene_pair_list)/num_cores)] for i in range(0,len(gene_pair_list),math.ceil(len(gene_pair_list)/num_cores))]
        print("{} double KOs to perform".format(len(gene_pair_list)))

        input_dicts=[]
        for doublegene_list_to_KO in doublegene_lists_to_KO:
            dict_to_add={"xml_dir":xml_dir,
                        "genes_to_KO":doublegene_list_to_KO,
                        "condition":condition,
                        "input_data":Input_Data,
                        "compartment_sizes_from_calibration":compartment_sizes_from_calibration,
                        "pg_fractions_from_calibration":pg_fractions_from_calibration,
                        "process_efficiencies_from_calibration":process_efficiencies_from_calibration,
                        "default_kapps_from_calibration":default_kapps_from_calibration,
                        "specific_kapps_from_calibration":specific_kapps_from_calibration,
                        "transporter_multiplier":1,
                        "Mumax_WT":Mumax_WT,
                        "Single_KO_results":output_single}
            input_dicts.append(dict_to_add)

        if n_parallel_processes!=1:
            if n_parallel_processes is None:
                num_cores=cpu_count()
                n_jobs=min(num_cores,len(input_dicts))
            else:
                n_jobs=n_parallel_processes
            print("Jobs:{}".format(n_jobs))

            pool=Pool(n_jobs)

            #simulation_results_pool=pool.map_async(simulate_run,input_dicts)
            #simulation_results=simulation_results_pool.get()

            double_KO_results=pool.imap_unordered(apply_double_knock_outs,input_dicts)

            #simulation_results=Parallel(n_jobs=n_jobs)(delayed(simulate_run)(input_dict) for input_dict in input_dicts)
        else:
            double_KO_results=[]
            for input_dict in input_dicts:
                double_KO_results.append(apply_double_knock_outs(input_dict))
        
        output_double=pandas.DataFrame()
        for i in double_KO_results:
            for j in i.keys():
                output_double.loc[j,"ID"]=j
                output_double.loc[j,"Mu_max_KO"]=i[j]
                output_double.loc[j,"Mu_max_WT"]=Mumax_WT
                output_double.loc[j,"Relative_Mu"]=i[j]/Mumax_WT
    
        output_double.to_csv("Double_Gene_KOs_{}.csv".format(file_name),sep=";",decimal=",")


if __name__ == "__main__":
    t0=time.time()
    main(number_chunks=None,
             transporter_multiplier=1,
             condition='Mean_01',
             xml_dir='../Yeast_iMM904_RBA_model',
             #compartments_to_consider=None,
             compartments_to_consider=["c","mIM","m"],
             #compartments_to_consider=["m"],
             n_parallel_processes=None,
             double_KOs=True,
             #single_KOs_to_import=None,
             single_KOs_to_import="Single_Gene_KOs_c_m_mIM.csv",
             file_name="c_m_mIM",
             threshold_for_zero_mu=0.005)
    print("Total time: {}".format(time.time()-t0))
