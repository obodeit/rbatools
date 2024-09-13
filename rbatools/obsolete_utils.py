import pandas

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

