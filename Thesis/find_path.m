function path = find_path(period, src, des)
    % Initialize the path to be an empty array
    path = [];
    
    % Iterate over each sender in this period
    for i = 1:size(period, 1)
        % Check if the sender is the source
        if period(i, 1) == src
            % Find the receiver in this period
            receiver = period(i, 2);
            % Check if the receiver is the destination
            if receiver == des
                % Add the source and destination to the path
                path = [src, des];
                % Exit the loop
                break;
            else
                % Find the nodes in the path from the receiver to the destination
                subpath = find_path(period, receiver, des);
                % Check if a path was found
                if ~isempty(subpath)
                    % Add the source and the nodes in the subpath to the path
                    path = [src, subpath];
                    % Exit the loop
                    break;
                end
            end
        end
    end
end