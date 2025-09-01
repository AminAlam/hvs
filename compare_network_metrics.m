
function compare_network_metrics(plv_values_active, plv_values_control, stim_types, days, threshold_method)
    % Input:
    % plv_values_active, plv_values_control: 5D matrices (stim × day × subject × ch × ch)
    % stim_types: cell array of stimulus type labels
    % days: cell array of day labels
    % threshold_method: 'proportion' or 'absolute' (default: 'proportion')
    
    if nargin < 5
        threshold_method = 'proportion';
    end
    
    [n_stim, n_days, n_subj_active, n_ch, ~] = size(plv_values_active);
    n_subj_control = size(plv_values_control, 3);
    
    % Initialize arrays for network metrics
    metrics_active = struct();
    metrics_control = struct();
    
    % Define network metrics to calculate
    metrics_active.clustering = zeros(n_stim, n_days, n_subj_active);
    metrics_active.path_length = zeros(n_stim, n_days, n_subj_active);
    metrics_active.global_efficiency = zeros(n_stim, n_days, n_subj_active);
    metrics_active.modularity = zeros(n_stim, n_days, n_subj_active);
    metrics_active.degree = zeros(n_stim, n_days, n_subj_active, n_ch);
    
    metrics_control.clustering = zeros(n_stim, n_days, n_subj_control);
    metrics_control.path_length = zeros(n_stim, n_days, n_subj_control);
    metrics_control.global_efficiency = zeros(n_stim, n_days, n_subj_control);
    metrics_control.modularity = zeros(n_stim, n_days, n_subj_control);
    metrics_control.degree = zeros(n_stim, n_days, n_subj_control, n_ch);
    
    % Calculate network metrics for each condition and subject
    for stim = 1:n_stim
        for day = 1:n_days
            % Process active group
            for subj = 1:n_subj_active
                adj_matrix = squeeze(plv_values_active(stim, day, subj, :, :));
                
                % Threshold the adjacency matrix
                if strcmp(threshold_method, 'proportion')
                    thresh = prctile(adj_matrix(:), 70);  % Keep top 30% of connections
                else
                    thresh = 0.5;  % Absolute threshold
                end
                adj_binary = adj_matrix > thresh;
                
                % Calculate network metrics
                [metrics_active.clustering(stim,day,subj), ...
                 metrics_active.path_length(stim,day,subj), ...
                 metrics_active.global_efficiency(stim,day,subj), ...
                 metrics_active.modularity(stim,day,subj), ...
                 metrics_active.degree(stim,day,subj,:)] = calculate_network_metrics(adj_binary);
            end
            
            % Process control group
            for subj = 1:n_subj_control
                adj_matrix = squeeze(plv_values_control(stim, day, subj, :, :));
                
                % Threshold the adjacency matrix
                if strcmp(threshold_method, 'proportion')
                    thresh = prctile(adj_matrix(:), 70);
                else
                    thresh = 0.5;
                end
                adj_binary = adj_matrix > thresh;
                
                % Calculate network metrics
                [metrics_control.clustering(stim,day,subj), ...
                 metrics_control.path_length(stim,day,subj), ...
                 metrics_control.global_efficiency(stim,day,subj), ...
                 metrics_control.modularity(stim,day,subj), ...
                 metrics_control.degree(stim,day,subj,:)] = calculate_network_metrics(adj_binary);
            end
        end
    end

    metrics_control.path_length(isnan(metrics_control.path_length)) = 0;
    metrics_active.path_length(isnan(metrics_active.path_length)) = 0;
    metrics_active.global_efficiency(isnan(metrics_active.global_efficiency)) = 0;
    metrics_control.global_efficiency(isnan(metrics_control.global_efficiency)) = 0;
    metrics_active.modularity(isnan(metrics_active.modularity)) = 0;
    metrics_control.modularity(isnan(metrics_control.modularity)) = 0;
    
    % Perform statistical comparisons
    % metric_names = {'Clustering Coefficient', 'Path Length', 'Global Efficiency', 'Modularity'};
    metric_names = {'Clustering Coefficient', 'Path Length', 'Global Efficiency', 'Modularity'};
    p_values = zeros(n_stim, n_days, length(metric_names));

    color_set1 = [230, 159, 0]/255;  % Orange
    color_set2 = [86, 180, 233]/255; % Sky Blue
    
    % Statistical testing for each metric
    for stim = 1:n_stim
        for day = 1:n_days
            % Compare clustering coefficient
            [p_values(stim,day,1)] = ranksum(squeeze(metrics_active.clustering(stim,day,:)), ...
                                           squeeze(metrics_control.clustering(stim,day,:)));
            % Compare path length
            [p_values(stim,day,2)] = ranksum(squeeze(metrics_active.path_length(stim,day,:)), ...
                                           squeeze(metrics_control.path_length(stim,day,:)));
            % Compare global efficiency
            [p_values(stim,day,3)] = ranksum(squeeze(metrics_active.global_efficiency(stim,day,:)), ...
                                           squeeze(metrics_control.global_efficiency(stim,day,:)));
            % Compare modularity
            [p_values(stim,day,4)] = ranksum(squeeze(metrics_active.modularity(stim,day,:)), ...
                                           squeeze(metrics_control.modularity(stim,day,:)));
        end
    end
    
    % Visualize results
    for metric = 1:length(metric_names)
        subplot(2,2,metric);
        
        % Extract data for current metric
        switch metric
            case 1
                active_data = metrics_active.clustering;
                control_data = metrics_control.clustering;
            case 2
                active_data = metrics_active.path_length;
                control_data = metrics_control.path_length;
            case 3
                active_data = metrics_active.global_efficiency;
                control_data = metrics_control.global_efficiency;
            case 4
                active_data = metrics_active.modularity;
                control_data = metrics_control.modularity;
        end
        
        % Plot data
        for stim = 3:n_stim
            subplot(2,2,metric);
            
            % Calculate means and standard errors
            active_mean = squeeze(mean(active_data(stim,:,:), 3));
            control_mean = squeeze(mean(control_data(stim,:,:), 3));
            active_se = squeeze(std(active_data(stim,:,:), [], 3)) / sqrt(n_subj_active);
            control_se = squeeze(std(control_data(stim,:,:), [], 3)) / sqrt(n_subj_control);
            
            % Plot with error bars
            errorbar((1:n_days), active_mean, active_se, '-o', 'LineWidth', 1.5, 'Color', color_set1)
            hold on;
            errorbar((1:n_days), control_mean, control_se, '--s', 'LineWidth', 1.5, 'Color', color_set2)
        end
        
        title(metric_names{metric}, 'Interpreter', 'tex', 'FontSize', 16);
        xlabel('Recording Day', 'Interpreter', 'tex', 'FontSize', 14);
        ylabel('Value', 'Interpreter', 'tex', 'FontSize', 14);
        stim_types_control = [];
        stim_types_active = [];
        for stim_type = stim_types(3:end)
            stim_types_control = [stim_types_control, strcat(strrep(string(stim_type(1)), '_', ' '), " - control")];
            stim_types_active = [stim_types_active, strcat(strrep(string(stim_type(1)), '_', ' '), " - active")];
        end
        legend([stim_types_control, stim_types_active], 'Location', 'best', 'Interpreter', 'tex', 'FontSize', 12)

        recording_days = {"Day 1", "Day 5", "Day 19"};
        
        grid on;
        ax = gca;
        ax.XTick = 1:length(recording_days);
        ax.XTickLabel = recording_days;
        ax.FontSize = 12;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
    end
    sgtitle('Network Metrics Comparison Between Groups', 'Interpreter', 'tex', 'FontSize', 18);
end
