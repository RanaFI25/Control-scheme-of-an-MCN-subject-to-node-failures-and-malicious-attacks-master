function clear_path = route_path(p)
routing = {};
for i = 1:length(p)-1
    routing{i} = p(i:i+1);
end
clear_path = routing;
end



