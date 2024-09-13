import pandas
import numpy
import json
import scipy.signal
from scipy.optimize import curve_fit

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


def do_regression(x_to_fit,y_to_fit,x_to_plot,max_val,min_val,monotonous_quadratic=False,total_x_range=(0,1),permit_quadratic_model=True):
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
    if permit_quadratic_model:
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
    else:
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

