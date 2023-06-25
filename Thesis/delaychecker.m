function delay = delaychecker(pathsCopy,sched_func)
   % Looping through from source to destination and check if we can
   % make hops towards destination that corresponds to delay.
    pos = 1; 
    sep_delay = {};
        
     % Removing from pathsCopy in every iteration. when its empty we have transmitted one packet. 
    while ~isempty(pathsCopy)    
        
        sz1 = size(sched_func{pos}); % size of each sublist in cell array 
        
        for ll = 1:sz1(1) % iteration from 1 to the length of each sublist
            if ~isempty(pathsCopy) && all(pathsCopy{1} == sched_func{pos}(ll,:))   %check if current route exist in any sub list 
                pathsCopy(1) = []; % removing from original path. 
                %sched_funcb{pos}(ll,:) = [];
                if(numel(sep_delay) >= pos && ~isempty(sep_delay))
                  
                     sep_delay{pos} = sep_delay{pos} + 1;
                else
                     sep_delay{pos} = 1;
                    
                end
            end
        end
            if pos == length(sched_func) % Now first iteration is done and we may have a hop. 
                pos = 1;  % Resetting current position so that we can go through again and made more hops towards destination  
                max_delay = max([sep_delay{:}]); % Find maximum number of hop during this iteration 
                [sep_delay{:}] = deal(max_delay);
            else
                pos = pos + 1; 
            end
    end
    delay = max([sep_delay{:}]);
end