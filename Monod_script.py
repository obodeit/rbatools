import pandas
from rbatools.rba_session import SessionRBA
from rbatools.calibration_utils import *
import time
from multiprocessing import Pool , cpu_count
from joblib import Parallel, delayed
import math

def main(transporter_multiplier=1,
         condition='Hackett_C01',
         xml_dir='../Yeast_iMM904_RBA_model'):
    
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

    glcs=[(i+1)*0.001 for i in range(35)]
    #glcs=[(i+1)*0.5 for i in range(3)]
    Results_glc_screen = {"Mu":[],"M_glc__D":[],"M_o2":[],"M_etoh":[],"M_ac":[],"M_acald":[]} # generate empty structure to store results of interest

    for c in glcs:
        print({'M_glc__D': c})
        rba_session.set_medium({'M_glc__D': c})
        mumax = rba_session.find_max_growth_rate()
        print(mumax)
        J_ex = rba_session.return_exchange_fluxes() # obtain corresponding exchange-fluxes
        for i in Results_glc_screen.keys():
            if i == "Mu":
                Results_glc_screen["Mu"].append(mumax)
            else:
                if i in J_ex.keys():
                    Results_glc_screen[i].append(J_ex[i])
                else:
                    Results_glc_screen[i].append(0)


    # plot results of glucose screen
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    ax.plot(glcs, Results_glc_screen["Mu"], linewidth=5, alpha=0.6, color="dimgray")
    ax.legend(["Monod curve"], loc="upper left", fontsize=12, frameon=False)
    ax.plot(glcs, Results_glc_screen["Mu"], linewidth=13, alpha=0.6, color="gainsboro")
    ax.plot(glcs, Results_glc_screen["Mu"], linewidth=10, alpha=0.6, color="lightgray")
    ax.plot(glcs, Results_glc_screen["Mu"], linewidth=7, alpha=0.6, color="silver")
    ax.plot(glcs, Results_glc_screen["Mu"], linewidth=5, alpha=0.6, color="gray")
    ax.plot(glcs, Results_glc_screen["Mu"], linewidth=2, alpha=0.6, color="dimgray")
    ax.plot(glcs, Results_glc_screen["Mu"], linewidth=1, alpha=0.6, color="black")
    ax2=ax.twinx()
    ax2.plot(glcs, [i for i in Results_glc_screen["M_etoh"]], linewidth=4,
            color="mediumseagreen", alpha=0.9)
    ax2.plot(glcs, [-i for i in Results_glc_screen["M_o2"]], linewidth=4, 
            color="indianred", alpha=0.9, linestyle=(0, (1, 1)))
    ax2.plot(glcs, [-i for i in Results_glc_screen["M_glc__D"]], linewidth=4, 
            color="steelblue", alpha=0.9, linestyle=(0, (3, 1, 1, 1)))
    #ax2.plot(glcs, [i for i in Results_glc_screen["M_acald"]], linewidth=4, 
    #        color="steelblue", alpha=0.9, linestyle=(0, (3, 1, 1, 1)))
    ax2.legend(["$EtOH$ excretion","$O_2$ uptake","$Glucose$ uptake"], 
            loc="lower right", fontsize=12, frameon=False)
    #ax2.legend(["$EtOH$ excretion","$O_2$ uptake","$Acetaldehyde$ excretion"], 
    #        loc="lower right", fontsize=12, frameon=False)

    ax.set_title("Glucose dependent growth-rate and exchange fluxes", fontsize=15)
    ax.set_xlabel("Glucose concentration (mM)", fontsize=12)
    #ax.set_ylim(-0.01,0.69)
    #ax.set_xlim(0.004,0.075)
    ax.tick_params(axis='y', labelcolor="dimgrey")
    ax2.tick_params(axis='y', labelcolor="black")
    ax.set_ylabel("Growth rate ($h^{-1}$)", color="dimgrey", fontsize=12)
    ax2.set_ylabel('Exchange fluxes ($\\frac{mmol}{h \\times g_{DW}}$)', 
                color="black", fontsize=12)
    #ax2.set_ylim(-0.2,15.6)
    plt.tight_layout()
    plt.savefig("Monod_curve_{}.png".format(condition),dpi=400)
    plt.show()

if __name__ == "__main__":
    t0=time.time()
    main(transporter_multiplier=1,
             condition='Hackett_C03',
             xml_dir='../Yeast_models/Yeast_iMM904_RBA_model')
    print("Total time: {}".format(time.time()-t0))
