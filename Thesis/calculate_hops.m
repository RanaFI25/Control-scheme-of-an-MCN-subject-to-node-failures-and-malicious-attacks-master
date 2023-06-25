function [total_hops, path] = calculate_hops(periods, src, dest)
% Inputs:
%   periods: cell array of scheduling periods, each period represented as a
%            matrix with rows as source-destination pairs
%   src: source node
%   dest: destination node
% Outputs:
%   total_hops: total number of hops between src and dest
%   path: a vector containing the intermediate nodes in the path from src to dest

path = [];
total_hops = 0;

% Check if src and dest are in the same period
same_period = false;
for p = 1:length(periods)
    period = periods{p};
    if any(period(:,1)==src) && any(period(:,2)==dest)
        same_period = true;
        period_idx = p;
        break;
    end
end

if same_period
    period = periods{period_idx};
    [total_hops, path] = calculate_hops_in_period(period, src, dest);
else
    % Find the period that contains the source node
    src_period_idx = 0;
    for p = 1:length(periods)
        period = periods{p};
        if any(period(:,1)==src)
            src_period_idx = p;
            break;
        end
    end
    if src_period_idx == 0
        error('Error: could not find source node in any period');
    end
    
    % Find the period that contains the destination node
    dest_period_idx = 0;
    for p = 1:length(periods)
        period = periods{p};
        if any(period(:,2)==dest)
            dest_period_idx = p;
            break;
        end
    end
    if dest_period_idx == 0
        error('Error: could not find destination node in any period');
    end
    
    % Calculate hops within the source period
    src_period = periods{src_period_idx};
    [hops_src, path_src] = calculate_hops_in_period(src_period, src, src_period(end,2));
    
    % Calculate hops within the destination period
    dest_period = periods{dest_period_idx};
    [hops_dest, path_dest] = calculate_hops_in_period(dest_period, dest_period(1,1), dest);
    
    % Calculate hops between the source and destination periods
    intermediate_hops = 0;
    intermediate_path = [];
    for p = src_period_idx+1:dest_period_idx-1
        period = periods{p};
        intermediate_hops = intermediate_hops + size(period,1);
        intermediate_path = [intermediate_path; period(:,1)];
    end
    
    % Combine the hops and paths
    total_hops = hops_src + intermediate_hops + hops_dest;
    path = [path_src; intermediate_path; path_dest(2:end)];
end

end
