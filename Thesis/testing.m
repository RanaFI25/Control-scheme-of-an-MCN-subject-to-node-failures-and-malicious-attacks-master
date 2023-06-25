clear all 
clc 
close all

s = [1,1,2,2,3,3,4,4,5,5];
d = [2,3,4,5,2,5,5,6,4,6];

src = 1; 
dst = 6;
names = {'Vc', 'V2', 'V3', 'V4', 'V5', 'Vu'};
weight = ones(1,length(s));
G= digraph(s,d,weight,names);
paths = allpaths(G, 1, 6);

paths_cell = num2cell(paths);
 
% Display the paths
disp('All paths between nodes 1 and 6:');
disp(paths_cell);
p = plot(G, 'Layout','force');
p.MarkerSize = 8;
p.EdgeColor = 'Black';
p.NodeColor = 'yellow';
%highlight(p, {'Vc', 'V2', 'V4', 'V5'}, {'V2', 'V4', 'V5', 'Vu'}, "EdgeColor","red");
% highlight(p, {'Vc', 'V2', 'V3', 'V4'}, {'V2', 'V3', 'V4', 'Vu'}, "EdgeColor","magenta");
% highlight(p, {'Vc', 'V2', 'V3', 'V5', 'V4'}, {'V2', 'V3', 'V5', 'V4', 'Vu'}, "EdgeColor","cyan");
% highlight(p, {'Vc', 'V2', 'V3', 'V5'}, {'V2', 'V3', 'V5', 'Vu'}, "EdgeColor","red");
% highlight(p, {'Vc', 'V2', 'V4', 'V5'}, {'V2', 'V4', 'V5', 'Vu'}, "EdgeColor","green");
% highlight(p, {'Vc', 'V2', 'V4'}, {'V2', 'V4', 'Vu'}, "EdgeColor","blue");
highlight(p, [src, dst], 'NodeColor', 'r');
camroll(-90);

brokenNode = 10;

newCellArray = {}; % initialize an empty cell array
for i = 1:length(paths)
    if ~ismember(brokenNode, paths{i})
        newCellArray{end+1} = paths{i}; % add the cell to the new cell array
    end
end
disp(newCellArray);
%paths = newCellArray;