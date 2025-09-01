
function compare_plv_groups(plv_values_active, plv_values_control, stim_types, days, alpha)
    % Input:
    % plv_values_active: 5D matrix for active group
    % plv_values_control: 5D matrix for control group
    % stim_types: cell array of stimulus type labels
    % days: cell array of day labels
    % alpha: significance level (default 0.05)
    
    if nargin < 5
        alpha = 0.05;
    end
    
    [n_stim, n_days, ~, n_ch, ~] = size(plv_values_active);
    
    % Initialize matrices to store statistics
    p_values = zeros(n_stim, n_days, n_ch, n_ch);
    stat_diff = zeros(n_stim, n_days, n_ch, n_ch);  % Store effect sizes
    sig_mask = zeros(n_stim, n_days, n_ch, n_ch);   % Binary mask for significant differences
    
    % Perform statistical tests for each condition and channel pair
    for stim = 1:n_stim
        for day = 1:n_days
            for ch1 = 1:n_ch
                for ch2 = 1:n_ch
                    if ch1 ~= ch2  % Skip diagonal elements
                        % Extract PLV values for this condition and channel pair
                        active_vals = squeeze(plv_values_active(stim, day, :, ch1, ch2));
                        control_vals = squeeze(plv_values_control(stim, day, :, ch1, ch2));
                        
                        % Perform Mann-Whitney U test (non-parametric)
                        [p, ~, stats] = ranksum(active_vals, control_vals);
                        
                        % Store results
                        p_values(stim, day, ch1, ch2) = p;
                        stat_diff(stim, day, ch1, ch2) = mean(active_vals) - mean(control_vals);
                        sig_mask(stim, day, ch1, ch2) = p < alpha;
                    end
                end
            end
        end
    end
    
    % Apply FDR correction for multiple comparisons
    for stim = 1:n_stim
        for day = 1:n_days
            p_temp = squeeze(p_values(stim, day, :, :));
            [sig_mask_fdr, ~] = fdr_correction(p_temp, alpha);
            sig_mask(stim, day, :, :) = sig_mask_fdr;
        end
    end

    min_value = min(stat_diff, [], 'all');
    max_value = max(stat_diff, [], 'all');
    
    % Visualization functions
    % 1. Plot difference matrices with significance markers
    for day = 1:n_days
        for stim = 1:n_stim
            subplot(n_stim, n_days, (day-1)*n_stim + stim);
            
            % Plot difference matrix
            diff_mat = squeeze(stat_diff(stim, day, :, :));
            sig_mat = squeeze(sig_mask(stim, day, :, :));
            
            imagesc(diff_mat);
            colormap('parula');
            colorbar;
            
            % Mark significant differences
            hold on;
            [row, col] = find(sig_mat);
            row_without_diag = [];
            col_without_diag = [];
            for i = 1:length(row)
                if row(i) ~= col(i)
                    row_without_diag = [row_without_diag, row(i)];
                    col_without_diag = [col_without_diag, col(i)];
                end
            end
            plot(col_without_diag, row_without_diag, 'k*', 'MarkerSize', 3);
            
            title(sprintf('Stim: %s, Day: %s', strrep(stim_types{stim}, '_', '-'), strrep(days{day}, '_', '-')), 'Interpreter', 'tex', 'FontSize', 16);
            xlabel('Channel', 'interpreter', 'tex', 'FontSize', 14);
            ylabel('Channel', 'interpreter', 'tex', 'FontSize', 14);
            axis square;
            clim([min_value, max_value])
            
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
    end
    sgtitle('PLV Differences (Active - Control) with Significant Connections (*)',  'interpreter', 'tex', 'FontSize', 18);
end
