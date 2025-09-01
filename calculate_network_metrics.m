
function [clustering, path_length, global_efficiency, modularity, degree] = calculate_network_metrics(adj_matrix)
    % Calculate basic network metrics
    n = size(adj_matrix, 1);
    
    % Clustering coefficient
    clustering = mean(clustering_coef_bu(adj_matrix));
    
    % Path length
    D = distance_bin(adj_matrix);
    path_length = mean(D(D~=Inf & ~eye(n)));
    
    % Global efficiency
    global_efficiency = efficiency_bin(adj_matrix);
    
    % Modularity
    [~, modularity] = community_louvain(adj_matrix);
    
    % Node degree
    degree = sum(adj_matrix, 2);
end

