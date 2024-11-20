
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


def sample_copy_numbers_from_residuals_quantiles(Input_data,replicate_cols,mean_col,replicate_threshold=1,filter_list=[],target_size=1,reps_to_sample=3,number_quantiles=1,transform_residuals=False,regression_type="lin",start_run_id=0,mean_no_noise=True,sample_mean=True):
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
    filter_list : list, optional
        _description_, by default []
    target_size : int, optional
        _description_, by default 1
    reps_to_sample : int, optional
        _description_, by default 3
    number_quantiles : int, optional
        _description_, by default 1
    transform_residuals : bool, optional
        _description_, by default False
    regression_type : str, optional
        _description_, by default "lin"
    start_run_id : int, optional
        _description_, by default 0
    mean_no_noise : bool, optional
        _description_, by default True
    sample_mean : bool, optional
        _description_, by default True
    """
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
            #x_reg = numpy.reshape(numpy.array(list(data_to_use["Log__mean"])), (len(list(data_to_use["Log__mean"])), 1))
            #y_reg = numpy.reshape(numpy.array(list(data_to_use["Log__sdev"])), (len(list(data_to_use["Log__sdev"])), 1))
            #regressor = LinearRegression(fit_intercept=True)
            #regressor.fit(x_reg, y_reg)
            #slope_sdev=regressor.coef_[0][0]
            #offset_sdev=regressor.intercept_[0]
            lin_regression_results=do_lin_regression(x_to_fit=list(data_to_use["Log__mean"]),y_to_fit=list(data_to_use["Log__sdev"]),fit_intercept=True)
            slope_sdev=lin_regression_results["A"]
            offset_sdev=lin_regression_results["B"]
            data_to_use["Fitted_Stdev"]=[offset_sdev+slope_sdev*data_to_use.loc[i,"Log__mean"] for i in data_to_use.index]
            df_intermediate["Fitted_Stdev"]=[offset_sdev+slope_sdev*df_intermediate.loc[i,"LogMean"] for i in df_intermediate.index]
        elif regression_type=="inverse_lin":
            #x_reg = numpy.reshape(numpy.array(list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__mean"])), (len(list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__mean"])), 1))
            #y_reg = numpy.reshape(numpy.array([1/i for i in list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__sdev"])]), (len(list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__sdev"])), 1))
            #regressor = LinearRegression(fit_intercept=True)
            #regressor.fit(x_reg, y_reg)
            #slope_sdev=regressor.coef_[0][0]
            #offset_sdev=regressor.intercept_[0]
            lin_regression_results=do_lin_regression(x_to_fit=list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__mean"]),y_to_fit=[1/i for i in list(data_to_use.loc[data_to_use["Log__sdev"]!=0,"Log__sdev"])],fit_intercept=True)
            slope_sdev=lin_regression_results["A"]
            offset_sdev=lin_regression_results["B"]
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

    count=start_run_id
    if mean_no_noise:
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
        out["sample_{}".format(count)]=[10**i for i in out["MeanSampleLog__run_{}".format(count)]]
        out2["sample_{}".format(count)]=out["sample_{}".format(count)]
    if sample_mean:
        out["Mean_of_log_samples"]=out.loc[:,[col for col in out.columns if col.startswith("MeanSampleLog__run_")]].mean(axis=1)
    out2["Gene"]=list(out2.index)
    return(out2)

