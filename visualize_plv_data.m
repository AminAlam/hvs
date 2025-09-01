
function visualize_plv_data(plv_values, stim_types, days, subjects)
    % Input:
    % plv_values: 5D matrix (stim_type × day_no × subject_no × no_ch × no_ch)
    % stim_types: cell array of stimulus type labels
    % days: array of day numbers or cell array of day labels
    % subjects: array of subject numbers or cell array of subject labels
    
    % 1. Function to plot PLV matrix for a specific condition
    function plot_plv_matrix(plv_matrix, title_str, min_value, max_value)
        imagesc(plv_matrix);
        colorbar;
        colormap('parula');
        % clim([0 1]);  % Set color limits between 0 and 1
        clim([min_value, max_value])
        axis square;
        title(title_str, 'Interpreter', 'tex');
        xlabel('Channel','interpreter', 'tex');
        ylabel('Channel','interpreter', 'tex');
        
        ch_list = {'FP1', 'FP2', 'C3', 'C4', 'TP7', 'TP8', 'O1', 'O2'};
        ax = gca;
        ax.XTick = 1:length(ch_list);
        ax.YTick = 1:length(ch_list);
        ax.XTickLabel = ch_list;
        ax.YTickLabel = ch_list;
        ax.FontSize = 12;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
    end
    
    % 2. Function to plot mean PLV across subjects for each condition and day
    function plot_mean_plv_conditions(plv_values)
        [n_stim, n_days, ~, n_ch, ~] = size(plv_values);
        mean_plv = squeeze(mean(plv_values, 3));  % Average across subjects
        min_value = min(mean_plv, [], 'all');
        max_value = max(mean_plv, [], 'all');
        
        % figure('Position', [100 100 1200 300]);
        for stim = 1:n_stim
            for day = 1:n_days
                subplot(n_stim, n_days, (day-1)*n_stim + stim);
                plot_plv_matrix(squeeze(mean_plv(stim, day, :, :)), ...
                    sprintf('Stim: %s, Day: %s', strrep(stim_types{stim}, '_', ' '), strrep(days{day}, '_', ' ')), min_value, max_value);
            end
        end
        sgtitle('Mean PLV Across Subjects','interpreter', 'tex');
    end
    
    % Call all visualization functions
    plot_mean_plv_conditions(plv_values);
    % plot_plv_distributions(plv_values);
    % plot_temporal_evolution(plv_values);
end
