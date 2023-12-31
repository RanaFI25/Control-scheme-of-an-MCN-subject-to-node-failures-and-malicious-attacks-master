function [tf_empty, transfer_function] = transferFunctionCalculate(paths,sched_func,G)

if isempty(paths) % if there is no path then no need of computing the transfer function
    tf_empty = 1;
    return;
else
    max_num_hops = max(cellfun(@(path) numel(path)-1, paths)); 
    maxDR = max_num_hops;
    pathsCopy = paths{1};
    sched_length = length(sched_func);
    delay = 1;
    iteration = false;
    delay_all_paths = {};
    

  
        % Loop through each path of the graph. 
        for i = 1: length(paths)
            pathsCopy = route_path(paths{i}); % Formating paths from 1 4 5 to cells {[1 4]} {[4 5]} corresponding to edges instead. 
            delay_all_paths{i} = [delaychecker(pathsCopy,sched_func), i]; % Computing delays/Hops for current paths by using defined scheduling function
        end

    %Sorting paths from least delay to maximum delay. 
    first_elems = cellfun(@(x) x(1), delay_all_paths);
    [~, idx] = sort(first_elems);
    sorted_A = delay_all_paths(idx);
    least_delay = sorted_A{1}(1);
    max_delay = sorted_A{length(sorted_A)}(1);
    
    %Creating dictionary which includes delays and weights. 
    dic_delays = dictionary;
    dic_weight = dictionary;
    for curr_delay = least_delay:max_delay
        first_elements = cellfun(@(x) x(1), sorted_A);
        output_array = sorted_A(first_elements == curr_delay); 

        %Now output array correspond to the path we could take in
        %curr_delay/curr_hops. we want to check what paths has been taken to compute product of weights.
        second_element = cellfun(@(x) x(2), output_array);
        nr_paths = length(second_element);
        curr_paths = {};
        for k = 1: nr_paths
            curr_paths{k} = route_path(paths{second_element(k)});

        end

        % As mentioned in report that each path can be taken only once
        % thats why this section is created to check we don't have taken
        % any duplicate paths beacuse it will otherwise impact the weights.
        
        seen_before = false(length(curr_paths),length(curr_paths{1}));
        sz3 = size(seen_before);
        for i = 1: sz3(1)
            len = sz3(2);
            for j = 1: len
                if i <= length(curr_paths) && j <= length(curr_paths{i})
                    curr = curr_paths{i}(1,j);
                    curr_value = 1;
                    ny_length = length(curr_paths);
                    for kk = 1: ny_length

                        if (any(cellfun(@(c) isequal(c, curr{1}), curr_paths{curr_value})) && curr_value ~= i)
                            seen_before(i,j) = true;
                            curr_paths(curr_value) = [];
                            curr_value= curr_value-1;
                        else
                            %not exist
                        end

                        curr_value= curr_value+1;
                    end

                end

            end
        end
        dic_delays(curr_delay) = length(curr_paths);
        % Computes the product of weight for current path.
        total_weight = 0;
        for len_paths = 1:length(curr_paths)
            curr_weight = 1;
            for len_p = 1:length(curr_paths{len_paths})
                curr_route = curr_paths{len_paths}(len_p);
                edge_idx = findedge(G, curr_route{1}(1), curr_route{1}(2));
                edge_weight = G.Edges.Weight(edge_idx);
                curr_weight = curr_weight * edge_weight;
            end
            total_weight = total_weight+curr_weight; %Sum the weight of all paths.  
        end
        dic_weight(curr_delay) = total_weight;

    end

    z = tf('z');
    transfer_function = 0;
    possible_path = length(inedges(G,paths{1}(length(paths{1}))));  
    keys = dic_delays.keys;
    for i = 1:numel(keys)
        key = keys(i);
        value = dic_delays(key);
        weight_value = dic_weight(key);
        transfer_function = transfer_function + (weight_value/z^key);%*(value/possible_path); % Defining the transfer function where numerator is the total weight and denominator is the z^delay. 
    end
    transfer_function = minreal(transfer_function); % Minimal realization of tranfer function.  
    tf_empty = 0;
end

end

