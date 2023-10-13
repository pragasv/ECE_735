close all; clear all; clc; 

karate_data = load('zachary_dataset.mat');
adjacency_matrix = karate_data.adj;
communities = karate_data.label; 

%% plot graph
G = graph(double(adjacency_matrix));

% Plot the graph with nodes colored based on labels
figure;
h = plot(G, 'NodeLabel', {}, 'Layout', 'force', 'MarkerSize', 12);
colormap([1 0 0; 0 0 1]); % Red (1) and blue (2) colors for nodes
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
        if (communities(i) == communities(j)) && (abs(adjacency_matrix(i, j))>0)
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

%% community_louvain :
%       gamma,
%           resolution parameter (optional)
%               gamma>1,        detects smaller modules
%               0<=gamma<1,     detects larger modules
%               gamma=1,        classic modularity (default)
%       M0,     
%           initial community affiliation vector (optional)
%       B,
%           objective-function type or custom objective matrix (optional)
%           'modularity',       modularity (default)
%           'potts',            Potts-model Hamiltonian (for binary networks)
%           'negative_sym',     symmetric treatment of negative weights
%           'negative_asym',    asymmetric treatment of negative weights
%           B,                  custom objective-function matrix
%
%           Note: see Rubinov and Sporns (2011) for a discussion of
%           symmetric vs. asymmetric treatment of negative weights.


% play around with gamma and see if we can retrive our original groupings
gamma=1;
[communities, modularity_value] = community_louvain(adjacency_matrix, gamma);
unique_values = unique(communities);
custom_colors = [
    1 0 0; % Red for unique number 1
    0 1 0; % Green for unique number 2
    0 0 1; % Blue for unique number 3
    1 1 0; % Yellow for unique number 4
    1 0 1; % Magenta for unique number 5
    0 1 1  % Cyan for unique number 6
];


% Plot the graph with nodes colored based on labels
figure;
h = plot(G, 'NodeLabel', {}, 'Layout', 'force', 'MarkerSize', 12);
node_colors = communities + 1; % Adjust labels to match colormap indices (1-based)
h.NodeCData = node_colors;

% Add node labels as text
node_labels = arrayfun(@num2str, 1:numel(communities), 'UniformOutput', false);
text(h.XData, h.YData, node_labels, 'FontSize', 10, 'HorizontalAlignment', 'center');

disp(['louvain Modularity: ' num2str(modularity_value)]);

