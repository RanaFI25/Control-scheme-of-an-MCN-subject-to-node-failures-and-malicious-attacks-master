function [ctrlLaw,tf] = Gr_Tf(ctrlLaw)
    s = [1,1,1,2,3,4];
    d = [2,3,4,5,5,5];
    sched_func = {[1,2],[1,3],[1,4],[2,5],[3,5],[4,5]};
    [paths,G] = createGraph(s,d,50);
    [tf_empty,transfer_function] = transferFunctionCalculate(paths, sched_func,G);
    
    if tf_empty
        ctrlLaw = 0;
        tf = 0; 
    else
        tf = transfer_function;
    end
end



