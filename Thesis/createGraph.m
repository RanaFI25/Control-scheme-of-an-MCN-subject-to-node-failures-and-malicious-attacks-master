function [paths,G] = createGraph(s,d,brokenNode,weight)
%s = [1,1,1,2,3,4];
%d = [2,3,4,5,5,5];


src = s(1);
dst = d(length(d));

names = {'V1', 'V2', 'V3','V4', 'V5'};%,'V6', 'V7', 'V8','V9', 'V10'};


%weight = 2*rand(1,length(s)) - 1;

G= digraph(s,d,weight,names);
paths = allpaths(G, src, dst);

paths_cell = num2cell(paths);

p = plot(G, 'Layout','auto','EdgeLabel', G.Edges.Weight);
p.MarkerSize = 10;
p.EdgeFontSize = 10;
p.ArrowSize = 10;
p.NodeFontSize =12;

p.EdgeColor = 'Black';
p.NodeColor = 'blue';
highlight(p, src, 'NodeColor', 'r');
highlight(p, dst, 'NodeColor', 'g');
camroll(90);






newCellArray = {}; % initialize an empty cell array
for i = 1:length(paths)
    if ~ismember(brokenNode, paths{i})
        newCellArray{end+1} = paths{i}; % add the cell to the new cell array
    end
end
disp(newCellArray);
paths = newCellArray;



% Count the number of hops in the shortest path
max_num_hops = max(cellfun(@(path) numel(path)-1, paths));

% Display the result
fprintf('The path between nodes %d and %d has %d hops.\n', src, dst, max_num_hops);

end

