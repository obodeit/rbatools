
import pandas
import numpy
#from sklearn.linear_model import LinearRegression
from rbatools.regression_utils import  do_lin_regression

def check_quantile(val,quantiles):
    """
    _summary_

    Parameters
    ----------
    val : _type_
        _description_
    quantiles : _type_
        _description_
    """
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


def sample_copy_numbers_from_residuals_quantiles(Input_data,
                                                 replicate_cols,
                                                 mean_col,
                                                 replicate_threshold=1,
                                                 target_size=1,
                                                 start_sample=0,
                                                 number_quantiles=1,
                                                 transform_residuals=False,
                                                 regression_type="lin",
                                                 start_run_id=0):
    """
    _summary_

    Parameters
    ----------
    Input_data : _type_
        _description_
    replicate_cols : _type_
        _description_
    mean_col : _type_
        _description_
    replicate_threshold : int, optional
        _description_, by default 1
    target_size : int, optional
        _description_, by default 1
    number_quantiles : int, optional
        _description_, by default 1
    transform_residuals : bool, optional
        _description_, by default False
    regression_type : str, optional
        _description_, by default "lin"
    start_run_id : int, optional
        _description_, by default 0
    """

    samples_df=pandas.DataFrame(index=list(Input_data.index))
    samples_df["Gene"]=Input_data["Gene"]

    sampled_replicates_df=pandas.DataFrame(index=list(Input_data.index))
    sampled_replicates_df["Gene"]=Input_data["Gene"]

    empirical_data_df=pandas.DataFrame(index=list(Input_data.index))

    for i in Input_data.index:
        vals=[]
        for j in replicate_cols:
            empirical_data_df.loc[i,"Log__{}".format(j)]=numpy.log10(Input_data.loc[i,j])
            vals.append(numpy.log10(Input_data.loc[i,j]))
        empirical_data_df.loc[i,"Log__mean"]=numpy.nanmean(vals)
        empirical_data_df.loc[i,"Log__sdev"]=numpy.nanstd(vals)

    if number_quantiles>1:
        quantiles=list(numpy.quantile(a=list(empirical_data_df.loc[pandas.isna(empirical_data_df["Log__mean"])==False,"Log__mean"]),q=[i/number_quantiles for i in list(range(number_quantiles+1))]))
        empirical_data_df["Quantile"]=[check_quantile(val=i,quantiles=quantiles) for i in list(empirical_data_df["Log__mean"])]
    else:
        empirical_data_df["Quantile"]=[1]*empirical_data_df.shape[0]
    
    for i in list(Input_data.index):
        finite_count=0
        for j in replicate_cols:
            if not pandas.isna(Input_data.loc[i,j]):
                finite_count+=1
        empirical_data_df.loc[i,"Number_quantified_replicates"]=finite_count

    if transform_residuals:
        if regression_type=="lin":
            lin_regression_results=do_lin_regression(x_to_fit=list(empirical_data_df.loc[(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),"Log__mean"]),
                                                     y_to_fit=list(empirical_data_df.loc[(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),"Log__sdev"]),
                                                     fit_intercept=True)
            slope_sdev=lin_regression_results["A"]
            offset_sdev=lin_regression_results["B"]
            empirical_data_df["Fitted_Stdev"]=[offset_sdev+slope_sdev*empirical_data_df.loc[i,"Log__mean"] for i in empirical_data_df.index]
        elif regression_type=="inverse_lin":
            lin_regression_results=do_lin_regression(x_to_fit=list(empirical_data_df.loc[(empirical_data_df["Log__sdev"]!=0)&(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),"Log__mean"]),
                                                     y_to_fit=[1/i for i in list(empirical_data_df.loc[(empirical_data_df["Log__sdev"]!=0)&(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),"Log__sdev"])],
                                                     fit_intercept=True)
            slope_sdev=lin_regression_results["A"]
            offset_sdev=lin_regression_results["B"]
            empirical_data_df["Fitted_Stdev"]=[1/(offset_sdev+slope_sdev*empirical_data_df.loc[i,"Log__mean"]) for i in empirical_data_df.index]

    all_residuals={i:[] for i in list(set(list(empirical_data_df["Quantile"]))) if not pandas.isna(i)}
    for quantile in list(set(list(empirical_data_df["Quantile"]))):
        if not pandas.isna(quantile):
            for i in replicate_cols:
                empirical_data_df.loc[empirical_data_df["Quantile"]==quantile,"Residual_empirical__{}".format(i)]=empirical_data_df.loc[empirical_data_df["Quantile"]==quantile,"Log__{}".format(i)]-empirical_data_df.loc[empirical_data_df["Quantile"]==quantile,"Log__mean"]
                residual_col_prefix="Residual_empirical__"
                if transform_residuals:
                    residual_col_prefix="Residual_empirical_scaled__"
                    empirical_data_df.loc[empirical_data_df["Quantile"]==quantile,"Residual_empirical_scaled__{}".format(i)]=empirical_data_df.loc[empirical_data_df["Quantile"]==quantile,"{}{}".format(residual_col_prefix,i)]/empirical_data_df.loc[empirical_data_df["Quantile"]==quantile,"Fitted_Stdev"]
                
                
                all_residuals[quantile]+=list([j for j in list(empirical_data_df.loc[empirical_data_df["Quantile"]==quantile,"{}{}".format(residual_col_prefix,i)]) if not pandas.isna(j)])

    if mean col in list(Input_data.columns):
        samples_df["mean_noNoise"]=Input_data[mean_col]
        samples_df["log_mean_noNoise"]=[10**ifor i in list(empirical_data_df.loc["Log__mean"])]

    count=start_run_id
    for run in list(range(start_sample,start_sample+target_size)):
        count+=1
        intermediate_sampling_DF=pandas.DataFrame(index=list(empirical_data_df.index))

        for rep in range(len(replicate_cols)):
            intermediate_sampling_DF.loc[(pandas.isna(empirical_data_df["Quantile"])==False)&(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),"sampled_residual"]=[list(numpy.random.choice(a=all_residuals[empirical_data_df.loc[protein,"Quantile"]],size=1,replace=True))[0] for protein in intermediate_sampling_DF.loc[(pandas.isna(empirical_data_df["Quantile"])==False)&(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),:].index]
            if transform_residuals:
                intermediate_sampling_DF.loc[(pandas.isna(empirical_data_df["Quantile"])==False)&(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),"sampled_residual"]*=intermediate_sampling_DF.loc[(pandas.isna(empirical_data_df["Quantile"])==False)&(empirical_data_df["Number_quantified_replicates"]>=replicate_threshold),"Fitted_Stdev"]
            
            intermediate_sampling_DF["sample_{}_log_rep_{}".format(rep+1,count)]=intermediate_sampling_DF["sampled_residual"]+empirical_data_df["Log__mean"]
            sampled_replicates_df["sample_{}_log_rep_{}".format(rep+1,count)]=intermediate_sampling_DF["sample_{}_log_rep_{}".format(rep+1,count)]

        sampled_replicates_df["sample_{}_log_mean".format(count)]=list(intermediate_sampling_DF.mean(axis=1,skipna=True))
        samples_df["sample_{}".format(count)]=[10**i for i in list(intermediate_sampling_DF.mean(axis=1,skipna=True))]
        
    return({"empirical_data":empirical_data_df,"sampled_replicates":sampled_replicates_df,"samples":samples_df})

