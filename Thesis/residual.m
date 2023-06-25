function r = residual(measured_Value, estimated_Value) %residuals
    diff = measured_Value - estimated_Value;
    
    tresholdUpperBound = 1; % upperbound
    tresholdLowerBound = -1; %lowerbound
    if all(diff > tresholdUpperBound) & all(diff < tresholdLowerBound)
        display("Out of treshold");
        r = 1;
    else r = 0;
    end
end
