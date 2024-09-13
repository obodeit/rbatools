from func_timeout import func_timeout, FunctionTimedOut
from multiprocessing import Pool , cpu_count

def long_running_function(duration):
    import time
    time.sleep(duration)
    return("Completed")

def run_with_timeout(input):
    run_id=list(input.keys())[0]
    print("Running: {}".format(run_id))
    #duration=input[run_id]["Sleep"]
    try:
        result = func_timeout(input[run_id]["Limit"], long_running_function, args=(input[run_id]["Sleep"],))
        return({run_id:result})
    except FunctionTimedOut:
        #print(f"Function did not complete within {timeout} seconds.")
        return({run_id:None})

def main(timeout_duration,n_parallel_processes):
    inputs=[{1:{"Sleep":1,"Limit":timeout_duration}},
               {2:{"Sleep":12,"Limit":timeout_duration}},
               {3:{"Sleep":10,"Limit":timeout_duration}},
               {4:{"Sleep":11,"Limit":timeout_duration}},
               {5:{"Sleep":3,"Limit":timeout_duration}},
               {6:{"Sleep":7,"Limit":timeout_duration}},
               {7:{"Sleep":19,"Limit":timeout_duration}},
               {8:{"Sleep":1,"Limit":timeout_duration}},
               {11:{"Sleep":3,"Limit":timeout_duration}},
               {12:{"Sleep":12,"Limit":timeout_duration}},
               {13:{"Sleep":10,"Limit":timeout_duration}},
               {14:{"Sleep":11,"Limit":timeout_duration}},
               {15:{"Sleep":3,"Limit":timeout_duration}},
               {16:{"Sleep":7,"Limit":timeout_duration}},
               {17:{"Sleep":19,"Limit":timeout_duration}},
               {18:{"Sleep":1,"Limit":timeout_duration}}]

    if n_parallel_processes!=1:
        if n_parallel_processes is None:
            num_cores=cpu_count()
            n_jobs=min(num_cores,len(inputs))
        else:
            n_jobs=min(n_parallel_processes,len(inputs))            
        pool=Pool(n_jobs)
        results=pool.imap_unordered(run_with_timeout,inputs)
    else:
        results=[run_with_timeout(input_dict) for input_dict in inputs]

    print([i for i in results])


if __name__ == "__main__":
    main(timeout_duration = 11,n_parallel_processes=8)
