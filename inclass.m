clear all; close all ;clc; 

adjacency_matrix = [
    0,0,1,0,0,0;
    0,0,2,0,0,5;
    1,2,0,3,4,5;
    0,0,3,0,5,0;
    0,0,4,5,0,0;
    0,5,5,0,0,0;];

communities = [3, 4, 1, 5, 4, 4];

%% plot graph
G = graph(double(adjacency_matrix));

% Plot the graph with nodes colored based on labels
figure;
h = plot(G, 'NodeLabel', {}, 'Layout', 'force', 'MarkerSize', 12);
custom_colors = [
    1 0 0; % Red for unique number 1
    0 1 0; % Green for unique number 2
    0 0 1; % Blue for unique number 3
    1 1 0; % Yellow for unique number 4
    1 0 1; % Magenta for unique number 5
    0 1 1  % Cyan for unique number 6
];
 % Red (1) and blue (2) colors for nodes
node_colors = communities + 1; % Adjust labels to match colormap indices (1-based)
h.NodeCData = node_colors;

% Add node labels as text
node_labels = arrayfun(@num2str, 1:numel(communities), 'UniformOutput', false);
text(h.XData, h.YData, node_labels, 'FontSize', 10, 'HorizontalAlignment', 'center');

%% calculate modularity
n = size(adjacency_matrix, 1);
m = sum(adjacency_matrix(:)) / 2;
degrees = sum(adjacency_matrix, 2);

modularity_value = 0;

unweighted_adjacency_matrix = adjacency_matrix>0;
unweighted_degrees = sum(unweighted_adjacency_matrix, 2);
unweighted_modularity_value = 0;

for i = 1:n
    for j = 1:n
        if (abs(communities(i) - communities(j))<=1) && (abs(adjacency_matrix(i, j))>0)
            modularity_value = modularity_value + (adjacency_matrix(i, j) - degrees(i)*degrees(j)/(2*m));
            % what will happen when we consider unweighted adjacency ? - can
            % we even consider it ?
            unweighted_modularity_value = unweighted_modularity_value + (unweighted_adjacency_matrix(i, j) - unweighted_degrees(i)*unweighted_degrees(j)/(2*m));
        end
    end
end

modularity_value = modularity_value / (2*m);
unweighted_modularity_value = unweighted_modularity_value / (2*m);

% Display the modularity value
disp(['Modularity: ' num2str(modularity_value)]);
disp(['Modularity Unweighted: ' num2str(unweighted_modularity_value)]);

