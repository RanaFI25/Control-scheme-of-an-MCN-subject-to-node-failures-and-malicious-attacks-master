function [controllable,observable] = contAndObsv(all_models,B_matrix,C_matrix)
len_models = length(all_models);
% Controllabiltiy
for curr_model = 1:len_models
    if rank(ctrb(all_models{curr_model},B_matrix{curr_model})) == size(all_models{curr_model},1)
        %disp(['Model ' , num2str(curr_model) ,' is controllable']);
        controllable = 1;
    else
        %disp(['Model ' , num2str(curr_model) ,' not controllable']);
        controllable = 0;
        break;
    end
end

%Observability


for curr_model = 1:len_models
    if rank(obsv(all_models{curr_model},C_matrix{curr_model})) == size(all_models{curr_model},1)
        %disp(['Model ' , num2str(curr_model) ,' is observable']);
        observable = 1; 
    else
        %disp(['Model ' , num2str(curr_model) ,' not observable']);
        observable = 0; 
        break;
    end
end

end

