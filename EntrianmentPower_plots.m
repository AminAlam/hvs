%% Power of 60Hz
clc
close all

path2save = fullfile("results", "Power_60");

no_channels = 8;
% target_freq_range = [59.36, 59.38];
target_freq_range = [59.28, 59.42];
% target_freq_range = [59.36, 59.40];
% target_freq_range = [59.95, 60.05];
whole_freq_range = [2, 80];

avg_power_over_days.active = zeros(numel(active_subjects), numel(stim_types), numel(recording_days), no_channels);
avg_power_over_days.control = zeros(numel(control_subjects), numel(stim_types), numel(recording_days), no_channels);

for subject_type = ["control", "active"]
    subject_counter = 0;
    for subject = eval(strcat(subject_type, "_subjects"))
        subject_counter = subject_counter+1;
        recording_day_counter = 0;
        for recording_day = recording_days
            recording_day_counter = recording_day_counter+1;
            stim_type_counter = 0;
            for stim_type = stim_types
                stim_type_counter = stim_type_counter+1;

                subject  = string(subject(1));
                recording_day = string(recording_day(1));
                stim_type = string(stim_type(1));

                data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
                recording_id = char(strcat(subject, "_", recording_day, "_", stim_type));
                target_power = [];
                whole_power = [];
                for ch = 1:size(data, 1)
                    target_power = [target_power, calcPower(data(ch, :)', fs, target_freq_range)];
                    whole_power = [whole_power, calcPower(data(ch, :)', fs, whole_freq_range)];
                end
                target_power(:, 5) = mean(target_power(:, [1:4, 6:end]), 2);
                whole_power(:, 5) = mean(whole_power(:, [1:4, 6:end]), 2);
                normalized_target_power = target_power./whole_power;

                data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).target_power = target_power;
                data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).whole_power = whole_power;
                data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).normalized_target_power = normalized_target_power;

                avg_power_over_days.(sprintf("%s",subject_type))(subject_counter, stim_type_counter, recording_day_counter, :) = normalized_target_power;
            end
        end
    end
end

disp('avg_power_over_days matrix made')

colorMatrix = [
    0, 0, 0;
    1, 0, 0;   % Red
    0, 1, 0;   % Green
    0, 0, 1;   % Blue   
    1, 1, 0;
    1, 0, 1;
    0, 1, 1;
    1, 1, 1;
    0.5, 1, 0;
    1, 0.5, 0;
    1, 0, 0.5;
    0.5, 1, 1;
    1, 0.5, 1;
    1, 1, 0.5;
];

shapesMatrix = {"-", "--", "-."};
legends_ = {};

% Individual subjects plot
figure
hold on
subject_color_counter = 0;
for subject_type = ["control", "active"]
    subject_counter = 0;
    for subject = eval(strcat(subject_type, "_subjects"))
        subject_counter = subject_counter+1;
        recording_day_counter = 0;
        subject_color_counter = subject_color_counter+1;
        for recording_day = recording_days
            recording_day_counter = recording_day_counter+1;
            stim_type_counter = 0;

            data_points = squeeze(avg_power_over_days.(sprintf("%s",subject_type))(subject_counter, :, recording_day_counter, :));
            data_points = data_points./data_points(1, :);
            data_points = mean(data_points, 2);
            plot(data_points, 'color', colorMatrix(subject_color_counter, :), 'LineWidth', 2, 'LineStyle', shapesMatrix{recording_day_counter})
            legends_{end+1} = [char(subject_type), '_', char(recording_day{1})];
            set(gca, 'YScale', 'log')
        end
    end
end

xticks([1 2 3]);
xticklabels({'No Light', 'Constant Light', 'Stimulus Light'});
legend(legends_)

% % Mean of the subjects plot
% figure
% hold on
% subject_color_counter = 0;
% for subject_type = ["control", "active"]
%     recording_day_counter = 0;
%     subject_color_counter = subject_color_counter+1;
%     for recording_day = recording_days
%         recording_day_counter = recording_day_counter+1;
%         stim_type_counter = 0;
% 
%         data_points = squeeze(avg_power_over_days.(sprintf("%s",subject_type))(:, :, recording_day_counter, :));
%         % data_points = data_points./data_points(1, :);
%         % data_points = mean(data_points, [1,2]);
%         % plot(data_points, 'color', colorMatrix(subject_color_counter, :), 'LineWidth', 2, 'LineStyle', shapesMatrix{recording_day_counter})
% 
%     end
% end
% 
% path2save = fullfile("results", "topoplot");
% 
% 
% % Topograph plot
% for subject_type = ["control", "active"]
%     subject_counter = 0;
%     for subject = eval(strcat(subject_type, "_subjects"))
%         subject_counter = subject_counter+1;
%         recording_day_counter = 0;
%         subject_color_counter = subject_color_counter+1;
%         for recording_day = recording_days
%             fig = figure;
%             recording_day_counter = recording_day_counter+1;
%             stim_type_counter = 0;
% 
%             data_points = squeeze(avg_power_over_days.(sprintf("%s",subject_type))(subject_counter, :, recording_day_counter, :));
%             % data_points = data_points./data_points(1, :);
% 
%             max_value = max(data_points, [], 'all');
%             min_value = 0;
% 
%             subplot(1,3,1)
%             plot_topography(ch_list, data_points(1,:))
%             clim([min_value, max_value])
%             colormap parula
% 
%             subplot(1,3,2)
%             plot_topography(ch_list, data_points(2,:))
%             clim([min_value, max_value])
%             colormap parula
% 
%             subplot(1,3,3)
%             plot_topography(ch_list, data_points(3,:))
%             clim([min_value, max_value])
%             colormap parula
% 
%             recording_day = string(recording_day(1));
%             subject_type = string(subject_type(1));
%             subject = string(subject(1));
% 
%             sgtitle(strcat(subject_type, "_", recording_day, "_", subject), 'interpreter','tex');
%             recording_id = strcat(subject_type, recording_day);
%             image_path_dir = fullfile(path2save, subject);
%             image_path = fullfile(image_path_dir, strcat(recording_id, ".pdf"));
%             if ~exist(image_path_dir, 'dir')
%                 mkdir(image_path_dir);
%             end
%             exportgraphics(fig, image_path, 'Resolution', 300);
%         end
%     end
% end


%% Descriptive analysis for power of 60Hz (Mean and Std of different conditions)
clc
close all

path2save = fullfile("results", "Power_60");

group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};

results = {};

for subject_type = ["control", "active"]
    recording_day_counter = 0;
    for recording_day = recording_days
        recording_day_counter = recording_day_counter+1;
        subject_counter = 0;
        data_group = zeros(numel(eval(strcat(subject_type, "_subjects"))), 3, no_channels);
        for subject = eval(strcat(subject_type, "_subjects"))
            subject_counter = subject_counter+1;
            data_points = squeeze(avg_power_over_days.(sprintf("%s",subject_type))(subject_counter, :, recording_day_counter, :));
            data_group(subject_counter, :, :) = data_points./data_points(1,:);
        end
        % data_group = squeeze(mean(data_group, 1));
        light_condition_counter = 0;
        for light_condition = light_conditions
            light_condition_counter = light_condition_counter+1;
            group_specific_light = squeeze(data_group(:, light_condition_counter, :));
            mean_group = mean(group_specific_light, "all");
            std_group = std(group_specific_light, [], "all");
            min_group = min(group_specific_light, [], "all");
            max_group = max(group_specific_light, [], "all");
            results = [results; {char(subject_type), day_labels{recording_day_counter}, ...
                light_conditions{light_condition_counter}, mean_group, std_group, min_group, max_group}];
        end
    end
end

results_table = cell2table(results, ...
    'VariableNames', {'Group', 'Day', 'LightCondition', 'MeanGroup', 'StdGroup', 'MinGroup', 'MaxGroup'});
writetable(results_table, fullfile(path2save, "descriptive_analysis.csv"));
disp('Results saved to group_statistics.csv');

%% Group AVG Topo Plots
clc
close all
path2save = fullfile("results", "topoplot");
group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};
no_channels = 8;
max_value = 0;
min_value = 0;
for subject_type = ["control", "active"]
    recording_day_counter = 0;
    for recording_day = recording_days
        recording_day_counter = recording_day_counter+1;
        subject_counter = 0;
        data_group = zeros(numel(eval(strcat(subject_type, "_subjects"))), 3, no_channels);
        for subject = eval(strcat(subject_type, "_subjects"))
            subject_counter = subject_counter+1;
            stim_type_counter = 0;
            data_points = squeeze(avg_power_over_days.(sprintf("%s",subject_type))(subject_counter, :, recording_day_counter, :));
            data_group(subject_counter, :, :) = data_points./data_points(1,:);
        end
        data_group = squeeze(mean(data_group, 1));

        if max(data_group, [], 'all') > max_value
            max_value = max(mag2db(data_group), [], 'all');
        end
        
        if min(data_group, [], 'all') < min_value
            min_value = min(mag2db(data_group), [], 'all');
        end

    end
end

% min_value = 0;
% max_value = 4;

for subject_type = ["control", "active"]
    fig = figure;
    i = 1;
    set(gcf, 'Position', get(0, 'Screensize'));
    recording_day_counter = 0;
    for recording_day = recording_days
        recording_day_counter = recording_day_counter+1;
        subject_counter = 0;
        data_group = zeros(numel(eval(strcat(subject_type, "_subjects"))), 3, no_channels);
        for subject = eval(strcat(subject_type, "_subjects"))
            subject_counter = subject_counter+1;
            stim_type_counter = 0;
            data_points = squeeze(avg_power_over_days.(sprintf("%s",subject_type))(subject_counter, :, recording_day_counter, :));
            data_group(subject_counter, :, :) = data_points./data_points(1,:);
        end
        data_group = squeeze(mean(data_group, 1));
        data_group = mag2db(data_group);
        recording_day = string(recording_day(1));

        subplot(3,3,i)
        plot_topography(ch_list, data_group(1,:))
        title(strcat("No light", " | ", day_labels{recording_day_counter}), 'interpreter', 'tex', 'FontSize', fontsize_labels);
        clim([min_value, max_value])
        colormap parula
        i = i+1;

        subplot(3,3,i)
        plot_topography(ch_list, data_group(2,:))
        title(strcat("Constant light", " | ", day_labels{recording_day_counter}), 'interpreter', 'tex', 'FontSize', fontsize_labels);
        clim([min_value, max_value])
        colormap parula
        i = i+1;

        subplot(3,3,i)
        plot_topography(ch_list, data_group(3,:))
        title(strcat("Stimulus light", " | ", day_labels{recording_day_counter}), 'interpreter', 'tex', 'FontSize', fontsize_labels);
        clim([min_value, max_value])
        colormap parula
        colorbar_handle = colorbar;
        ylabel(colorbar_handle, 'Normalized Power (dB)', 'interpreter', 'latex', 'FontSize', fontsize_axis); % Customize the label and font size
        i = i+1;
    end
    stim_type = string(stim_type(1));
    image_path_dir = fullfile(path2save);
    if ~exist(image_path_dir, 'dir')
        mkdir(image_path_dir);
    end
    ax = gca;
    ax.FontSize = fontsize_axis;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    sgtitle(strcat(subject_type), 'interpreter', 'tex', 'FontSize', 18);
    image_path = fullfile(image_path_dir, strcat("Topoplots", "_", subject_type, "_", stim_type, ".pdf"));
    exportgraphics(fig, image_path, 'Resolution', 300);
    image_path = fullfile(image_path_dir, strcat("Topoplots", "_", subject_type, "_", stim_type, ".png"));
    exportgraphics(fig, image_path, 'Resolution', 300);
end
%% QQ plots
clc
close all
data = cat(5, avg_power_over_days.control, avg_power_over_days.active); % subjects * light_condition * days * channels * subject_type
data = data./data(:,1, :, :, :);
path2save = fullfile("results", "qq_plot");

group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};
num_days = size(data, 3);
num_light_conditions = size(data, 2);
num_groups = size(data, 5);

% Loop over the days and create box plots | control vs active each and stim
% type
min_value = 0;
max_value = 0;

for day_idx = 1:num_days
    % Extract data for the current day
    for light_idx = 2:num_light_conditions
        data_control = squeeze(data(:,light_idx,day_idx,:,1));  % Control condition, averaging over subjects and channels
        data_active = squeeze(data(:,light_idx,day_idx,:,2));    % Active condition, averaging over subjects and channels        
        data_control_vec = data_control(:);
        data_active_vec = data_active(:);
        if max(data_control_vec) > max_value
            max_value = max(data_control_vec);
        end
        if max(data_active_vec) > max_value
            max_value = max(data_active_vec);
        end
    end
end

for day_idx = 1:num_days
    % Extract data for the current day
    for light_idx = 2:num_light_conditions
        data_control = squeeze(data(:,light_idx,day_idx,:,1));  % Control condition, averaging over subjects and channels
        data_active = squeeze(data(:,light_idx,day_idx,:,2));  % Active condition, averaging over subjects and channels        
        data_control_vec = data_control(:);
        data_active_vec = data_active(:);

        % Test for normality
        [h_control, p_control] = kstest((data_control_vec - mean(data_control_vec)) / std(data_control_vec)); % K-S test for control
        [h_active, p_active] = kstest((data_active_vec - mean(data_active_vec)) / std(data_active_vec));     % K-S test for active
        
        % Determine if the data is normally distributed
        if h_control == 0
            control_title = sprintf('Control (Normal, p = %.3f)', p_control);
        else
            control_title = sprintf('Control (Not Normal, p = %.3f)', p_control);
        end
        
        if h_active == 0
            active_title = sprintf('Active (Normal, p = %.3f)', p_active);
        else
            active_title = sprintf('Active (Not Normal, p = %.3f)', p_active);
        end

        % Create Q-Q plots
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        subplot(2,1,1);
        qqplot(data_control_vec);
        title(control_title,'interpreter', 'tex', 'FontSize', fontsize_title);

        ax = gca;
        ax.FontSize = fontsize_axis;
        ax.LineWidth = 1.5;
        ax.Box = 'on';

        subplot(2,1,2);
        qqplot(data_active_vec);
        title(active_title,'interpreter', 'tex', 'FontSize', fontsize_title);
                
        sgtitle([light_conditions{light_idx} ', ' day_labels{day_idx}],'interpreter', 'tex', 'FontSize', fontsize_title);
        ax = gca;
        ax.FontSize = fontsize_axis;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
        
        image_path_dir = fullfile(path2save);
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end

        image_path = fullfile(image_path_dir, strcat("qq_plot", "_day", num2str(day_idx), "_", num2str(light_idx), ".pdf"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, strcat("qq_plot", "_day", num2str(day_idx), "_", num2str(light_idx), ".png"));
        exportgraphics(fig, image_path, 'Resolution', 300);
    end
end
%% Box plots (all channels of each subjects)
% clc
% close all
% data = cat(5, avg_power_over_days.control, avg_power_over_days.active); % subjects * light_condition * days * channels * subject_type
% 
% group_labels = {'Control', 'Active'};
% day_labels = {'Day 1', 'Day 5', 'Day 19'};
% light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};
% num_days = size(data, 3);
% num_light_conditions = size(data, 2);
% num_groups = size(data, 5);
% 
% % Find global min and max values
% min_value = 0;
% max_value = 0;
% for day_idx = 1:num_days
%     for light_idx = 1:num_light_conditions
%         data_control = squeeze(data(:,light_idx,day_idx,:,1));
%         data_active = squeeze(data(:,light_idx,day_idx,:,2));
%         max_value = max([max_value; data_control(:); data_active(:)]);
%     end
% end
% 
% % Control vs Active comparisons
% for day_idx = 1:num_days
%     for light_idx = 1:num_light_conditions
%         data_control = squeeze(data(:,light_idx,day_idx,:,1));
%         data_active = squeeze(data(:,light_idx,day_idx,:,2));
% 
%         figure;
%         hold on;
% 
%         % Create box plots
%         scatter(repmat(1, length(data_control(:)), 1), data_control(:), 'b', 'jitter', 'on', 'jitterAmount', 0.05);
%         scatter(repmat(2, length(data_active(:)), 1), data_active(:), 'r', 'jitter', 'on', 'jitterAmount', 0.05);
%         boxplot([data_control(:), data_active(:)], 'Labels', {'Control', 'Active'});
% 
%         % Perform Mann-Whitney U test (non-parametric alternative to t-test)
%         [p_value, ~] = ranksum(data_control(:), data_active(:));
% 
%         % Add significance bracket
%         bracket_height = max_value * 1.1;  % Place bracket above the maximum value
%         add_significance_bracket(1, 2, bracket_height, p_value);
% 
%         title([light_conditions{light_idx} ', ' day_labels{day_idx}]);
%         ylabel('Normalized Entrainment Power');
%         ylim([min_value, max_value*2]);  % Extend y-axis to accommodate brackets
%     end
% end
% 
% % Time effect comparisons
% for group_idx = 1:num_groups
%     for light_idx = 1:num_light_conditions
%         figure;
%         hold on;
% 
%         day_1 = data(:,light_idx,1,:,group_idx);
%         day_5 = data(:,light_idx,2,:,group_idx);
%         day_19 = data(:,light_idx,3,:,group_idx);
% 
%         % Create box plots
%         scatter(repmat(1, length(day_1(:)), 1), day_1(:), 'b', 'jitter', 'on', 'jitterAmount', 0.05);
%         scatter(repmat(2, length(day_5(:)), 1), day_5(:), 'r', 'jitter', 'on', 'jitterAmount', 0.05);
%         scatter(repmat(3, length(day_19(:)), 1), day_19(:), 'g', 'jitter', 'on', 'jitterAmount', 0.05);
%         boxplot([day_1(:), day_5(:), day_19(:)], 'Labels', {'Day 1', 'Day 5', 'Day 19'});
% 
%         % Perform Kruskal-Wallis test between all pairs
%         % Day 1 vs Day 5
%         [p_value_1_5, ~] = ranksum(day_1(:), day_5(:));
%         % Day 5 vs Day 19
%         [p_value_5_19, ~] = ranksum(day_5(:), day_19(:));
%         % Day 1 vs Day 19
%         [p_value_1_19, ~] = ranksum(day_1(:), day_19(:));
% 
%         % Add significance brackets at different heights
%         bracket_height_1 = max_value * 1.1;
%         bracket_height_2 = max_value * 1.2;
%         bracket_height_3 = max_value * 1.3;
% 
%         add_significance_bracket(1, 2, bracket_height_1, p_value_1_5);    % Day 1 vs 5
%         add_significance_bracket(2, 3, bracket_height_2, p_value_5_19);   % Day 5 vs 19
%         add_significance_bracket(1, 3, bracket_height_3, p_value_1_19);   % Day 1 vs 19
% 
%         title([group_labels{group_idx} ', ' light_conditions{light_idx}]);
%         ylabel('Normalized Entrainment Power');
%         ylim([min_value, max_value*2]);  % Extend y-axis to accommodate all brackets
%     end
% end

%% Complete analysis with box plots and normalized time course plots
% Box plots with subject-specific markers
clc
close all

path2save = fullfile("results", "entrainment_power_box_plots");

% Set default font sizes
fontsize_title = 16;
fontsize_labels = 14;
fontsize_axis = 12;

% Define marker styles for subjects
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
num_marker_styles = length(marker_styles);

% Get data dimensions
num_subjects = size(avg_power_over_days.control, 1);
num_channels = size(avg_power_over_days.control, 4);
num_light_conditions = size(avg_power_over_days.control, 2);
num_days = size(avg_power_over_days.control, 3);

% Reshape data for easier subject-level analysis
control_data = reshape(avg_power_over_days.control, [num_subjects, num_light_conditions, num_days, num_channels]);
active_data = reshape(avg_power_over_days.active, [num_subjects, num_light_conditions, num_days, num_channels]);
control_data = control_data./control_data(:, 1, :, :);
active_data = active_data./active_data(:, 1, :, :);
% Labels
group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};

% Find global min and max for consistent y-axis limits
all_data = cat(1, control_data(:), active_data(:));
global_min = min(all_data);
global_max = max(all_data);
y_range = [global_min, global_max*1.4]; % Extend upper limit for significance brackets
y_range_avg_channels = [global_min, global_max*0.8]; % Extend upper limit for significance brackets
%%
table_row_idx = 1;
p_values_table = {};
% Create subject-level plots for each light condition and day
for light_idx = 1:num_light_conditions
    for day_idx = 1:num_days
        fig = figure('Position', [100, 100, 800, 600]);
        set(gcf, 'Position', get(0, 'Screensize'));
        hold on;
        
        % Get data for current condition
        control_current = squeeze(control_data(:, light_idx, day_idx, :));
        active_current = squeeze(active_data(:, light_idx, day_idx, :));
        
        % Create box plots first (to be in background)
        boxplot([control_current(:), active_current(:)], ...
            [ones(size(control_current(:))); 2*ones(size(active_current(:)))], ...
            'Labels', {'Control', 'Active'}, 'Width', 0.7);
        h = findobj(gca, 'Tag', 'Box');
        patch(get(h(1), 'XData'), get(h(1), 'YData'), color_set2, 'FaceAlpha', 0.2);
        patch(get(h(2), 'XData'), get(h(2), 'YData'), color_set1, 'FaceAlpha', 0.2);
        
        % Plot individual subject data points with specific markers
        for subj = 1:num_subjects
            % Control group
            marker_idx = mod(subj-1, num_marker_styles) + 1;
            scatter(ones(size(control_current(subj,:))) + randn(size(control_current(subj,:)))*0.05, ...
                control_current(subj,:), 100, ...
                marker_styles{marker_idx}, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', color_set1, ...
                'MarkerFaceAlpha', 0.6);
            
            % Active group
            scatter(2*ones(size(active_current(subj,:))) + randn(size(active_current(subj,:)))*0.05, ...
                active_current(subj,:), 100, ...
                marker_styles{marker_idx}, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', color_set2, ...
                'MarkerFaceAlpha', 0.6);
        end
        
        % Perform statistical test
        [p_value, ~] = ranksum(control_current(:), active_current(:));
        add_significance_bracket(1, 2, global_max*1.05, p_value);

        % Store p-values in table
        p_values_table{table_row_idx, 1} = light_conditions{light_idx};
        p_values_table{table_row_idx, 2} = day_labels{day_idx};
        p_values_table{table_row_idx, 3} = p_value;
        table_row_idx = table_row_idx + 1;
        
        % Formatting
        title([light_conditions{light_idx} ', ' day_labels{day_idx}], ...
            'FontSize', fontsize_title, 'FontWeight', 'bold','interpreter', 'tex');
        ylabel('Normalized Entrainment Power', 'FontSize', fontsize_labels,'interpreter', 'tex');
        xlabel('Groups', 'FontSize', fontsize_labels,'interpreter', 'tex');
        ylim(y_range);
        
        ax = gca;
        ax.FontSize = fontsize_axis;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
            
        image_path_dir = fullfile(path2save, "subject_level_each_light_day_pair");
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        image_path = fullfile(image_path_dir, strcat("BoxPlot", "_", light_conditions{light_idx}, "_", day_labels{day_idx}, ".pdf"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, strcat("BoxPlot", "_", light_conditions{light_idx}, "_", day_labels{day_idx}, ".png"));
        exportgraphics(fig, image_path, 'Resolution', 300);
            
        % Adjust figure position to accommodate legend
        % set(gcf, 'Position', [100, 100, 1000, 600]);
    end
end

% Save as CSV in the same directory
csv_path = fullfile(path2save, "subject_level_each_light_day_pair", "statistical_analysis_results.csv");
p_values_table_headers = {'Light Condition', 'Day', 'p-value (control vs active)'};
p_values_table_final = cell2table(p_values_table, 'VariableNames', p_values_table_headers);
writetable(p_values_table_final, csv_path);

%%
% Time course plots
for light_idx = 1:num_light_conditions
    % figure('Position', [100, 100, 1200, 600]);
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    
    % Plot control group
    % subplot(1, 2, 1);
    hold on;
    
    % Individual trajectories
    for subj = 1:num_subjects
        control_trajectory = squeeze(mean(control_data(subj, light_idx, :, :), 4));
        plot(1:num_days, control_trajectory, '-o', ...
            'Color', [color_set_control, 0.3], 'LineWidth', 1.5, ...
            'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'none');
    end
    
    % Mean trajectory
    mean_control = squeeze(mean(mean(control_data(:, light_idx, :, :), 4), 1));
    plot(1:num_days, mean_control, '-o', ...
        'Color', color_set_control, 'LineWidth', 3, ...
        'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'none', ...
        'MarkerSize', 12);
    % 
    % title(['Control Group - ' light_conditions{light_idx}], ...
    %     'FontSize', fontsize_title, 'FontWeight', 'bold');
    xlabel('Time Point', 'FontSize', fontsize_labels);
    ylabel('Mean of Normalized Entrainment Power', 'FontSize', fontsize_labels);
    set(gca, 'XTick', 1:num_days, 'XTickLabel', day_labels);
    ylim(y_range_avg_channels);
     
    % 
    % % Plot active group
    % subplot(1, 2, 2);
    % hold on;
    
    % Individual trajectories
    for subj = 1:num_subjects
        active_trajectory = squeeze(mean(active_data(subj, light_idx, :, :), 4));
        plot(1:num_days, active_trajectory, '-o', ...
            'Color', [color_set_active, 0.3], 'LineWidth', 1.5, ...
            'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'none');
    end
    
    % Mean trajectory
    mean_active = squeeze(mean(mean(active_data(:, light_idx, :, :), 4), 1));
    plot(1:num_days, mean_active, '-o', ...
        'Color', color_set_active, 'LineWidth', 3, ...
        'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'none', ...
        'MarkerSize', 12);
    
    % title(['Active Group - ' light_conditions{light_idx}], ...
    %     'FontSize', fontsize_title, 'FontWeight', 'bold');
    xlabel('Time Point', 'FontSize', fontsize_labels, 'interpreter', 'tex');
    ylabel('Mean of Normalized Entrainment Power', 'FontSize', fontsize_labels, 'interpreter', 'tex');
    set(gca, 'XTick', 1:num_days, 'XTickLabel', day_labels);
    ylim(y_range_avg_channels);

    ax = gca;
    ax.FontSize = fontsize_axis;
    ax.LineWidth = 1.5;
    ax.Box = 'on';

    % Add a dummy plot to represent single channels (for legend)
    dummySingle1 = plot(nan, nan, 'Color', color_set_control, 'LineWidth', 1.5);
    dummySingle2 = plot(nan, nan, 'Color', color_set_active, 'LineWidth', 1.5);

    legend([dummySingle1, dummySingle2], {'Control', 'Active'}, 'FontSize', 14, 'Location', 'best');
    
     
    image_path_dir = fullfile(path2save, "subject_level_over_time");
    
    if ~exist(image_path_dir, 'dir')
        mkdir(image_path_dir);
    end
    sgtitle([light_conditions{light_idx}], ...
    'FontSize', fontsize_title+2, 'FontWeight', 'bold','interpreter', 'tex');
    image_path = fullfile(image_path_dir, strcat("BoxPlot", "_", light_conditions{light_idx}, ".pdf"));
    exportgraphics(fig, image_path, 'Resolution', 300);
    image_path = fullfile(image_path_dir, strcat("BoxPlot", "_", light_conditions{light_idx}, ".png"));
    exportgraphics(fig, image_path, 'Resolution', 300);
    
end
%%
% Day comparisons for each group and light condition
table_row_idx = 1;
p_values_table = {};
for group_idx = 1:2
    for light_idx = 1:num_light_conditions
        % figure('Position', [100, 100, 800, 600]);
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        hold on;
        
        if group_idx == 1
            current_data = control_data;
            group_color = color_set1;
        else
            current_data = active_data;
            group_color = color_set2;
        end
        
        % Get data for each day
        day_1_data = squeeze(current_data(:, light_idx, 1, :));
        day_5_data = squeeze(current_data(:, light_idx, 2, :));
        day_19_data = squeeze(current_data(:, light_idx, 3, :));
        
        % Create scatter plots with jitter
        
        % Plot individual subject data points with specific markers
        for subj = 1:num_subjects
            % Control group
            marker_idx = mod(subj-1, num_marker_styles) + 1;
        
            scatter(ones(size(day_1_data(subj, :))), day_1_data(subj, :), 50, marker_styles{marker_idx}, ...
                'MarkerFaceColor', color_set1, 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.05);
            scatter(2*ones(size(day_5_data(subj, :))), day_5_data(subj, :), 50, marker_styles{marker_idx}, ...
                'MarkerFaceColor', color_set2, 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.05);
            scatter(3*ones(size(day_19_data(subj, :))), day_19_data(subj, :), 50, marker_styles{marker_idx}, ...
                'MarkerFaceColor', color_set3, 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.05);
        end

        % Create box plots
        boxplot([day_1_data(:), day_5_data(:), day_19_data(:)], ...
            'Labels', {'Day 1', 'Day 5', 'Day 19'});
        h = findobj(gca, 'Tag', 'Box');
        patch(get(h(1), 'XData'), get(h(1), 'YData'), color_set3, 'FaceAlpha', 0.2);
        patch(get(h(2), 'XData'), get(h(2), 'YData'), color_set2, 'FaceAlpha', 0.2);
        patch(get(h(3), 'XData'), get(h(2), 'YData'), color_set1, 'FaceAlpha', 0.2);
        
        % Statistical tests
        [p_value_1_5, ~] = ranksum(day_1_data(:), day_5_data(:));
        [p_value_5_19, ~] = ranksum(day_5_data(:), day_19_data(:));
        [p_value_1_19, ~] = ranksum(day_1_data(:), day_19_data(:));

        % Store p-values in table
        p_values_table{table_row_idx, 1} = group_labels{group_idx};
        p_values_table{table_row_idx, 2} = light_conditions{light_idx};
        p_values_table{table_row_idx, 3} = "Day 1 vs Day 5";
        p_values_table{table_row_idx, 4} = p_value_1_5;
        table_row_idx = table_row_idx + 1;

        p_values_table{table_row_idx, 1} = group_labels{group_idx};
        p_values_table{table_row_idx, 2} = light_conditions{light_idx};
        p_values_table{table_row_idx, 3} = "Day 5 vs Day 19";
        p_values_table{table_row_idx, 4} = p_value_5_19;
        table_row_idx = table_row_idx + 1;

        p_values_table{table_row_idx, 1} = group_labels{group_idx};
        p_values_table{table_row_idx, 2} = light_conditions{light_idx};
        p_values_table{table_row_idx, 3} = "Day 1 vs Day 19";
        p_values_table{table_row_idx, 4} = p_value_1_19;
        table_row_idx = table_row_idx + 1;
        
        % Add significance brackets
        add_significance_bracket(1, 2, global_max*1.05, p_value_1_5);
        add_significance_bracket(2, 3, global_max*1.1, p_value_5_19);
        add_significance_bracket(1, 3, global_max*1.15, p_value_1_19);
        
        % Formatting
        title([group_labels{group_idx} ' - ' light_conditions{light_idx}], ...
            'FontSize', fontsize_title, 'FontWeight', 'bold', 'interpreter', 'tex');
        ylabel('Normalized Entrainment Power', 'FontSize', fontsize_labels, 'interpreter', 'tex');
        xlabel('Time Points', 'FontSize', fontsize_labels, 'interpreter', 'tex');
        ylim(y_range);
        
        ax = gca;
        ax.FontSize = fontsize_axis;
        ax.LineWidth = 1.5;
        ax.Box = 'on';

        image_path_dir = fullfile(path2save, "day_comparisons_each_group_light_pair");
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        image_path = fullfile(image_path_dir, strcat("BoxPlot", "_", light_conditions{light_idx}, "_", group_labels{group_idx}, ".pdf"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, strcat("BoxPlot", "_", light_conditions{light_idx}, "_", group_labels{group_idx}, ".png"));
        exportgraphics(fig, image_path, 'Resolution', 300);

    end
end

% Save as CSV in the same directory
csv_path = fullfile(path2save, "light_comparisons", "statistical_analysis_results.csv");
p_values_table_headers = {'Group', 'Light condition', 'Day X vs Day y', 'p-value (Day x vs Day y)'};
p_values_table_final = cell2table(p_values_table, 'VariableNames', p_values_table_headers);
writetable(p_values_table_final, csv_path);
%% Active vs control over days in the same plot
clc;
close all;
% Initialize table for statistical results
table_row_idx = 1;
p_values_table = {};

% Light condition labels
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};

for light_idx = 1:3  % Loop over light conditions
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    
    % Collect data for all days and groups
    data_all_days = [];
    group_labels_all = [];
    day_labels_all = [];
    
    for day_idx = 1:num_days  % Loop over days
        control_data_flat = squeeze(control_data(:, light_idx, day_idx, :));
        active_data_flat = squeeze(active_data(:, light_idx, day_idx, :));
        control_data_flat = control_data_flat(:);
        active_data_flat = active_data_flat(:);

        % Normality check (Shapiro-Wilk test)
        [h_control, p_control] = swtest(control_data_flat, 0.05);  % Test for Control group
        [h_active, p_active] = swtest(active_data_flat, 0.05);    % Test for Active group

        p_values_table{table_row_idx, 1} = light_conditions{light_idx};
        p_values_table{table_row_idx, 2} = day_labels{day_idx};
        p_values_table{table_row_idx, 3} = 'Control Group Normality';
        p_values_table{table_row_idx, 4} = 'Shapiro-Wilk Test';
        p_values_table{table_row_idx, 5} = 'N/A'; % No test statistic available for Shapiro-Wilk
        p_values_table{table_row_idx, 6} = length(control_data_flat); % Sample size
        p_values_table{table_row_idx, 7} = p_control;
        table_row_idx = table_row_idx + 1;

        p_values_table{table_row_idx, 1} = light_conditions{light_idx};
        p_values_table{table_row_idx, 2} = day_labels{day_idx};
        p_values_table{table_row_idx, 3} = 'Active Group Normality';
        p_values_table{table_row_idx, 4} = 'Shapiro-Wilk Test';
        p_values_table{table_row_idx, 5} = 'N/A'; % No test statistic available for Shapiro-Wilk
        p_values_table{table_row_idx, 6} = length(active_data_flat); % Sample size
        p_values_table{table_row_idx, 7} = p_active;
        table_row_idx = table_row_idx + 1;
        
        % Append to all data
        data_all_days = [data_all_days; control_data_flat; active_data_flat];
        group_labels_all = [group_labels_all; repmat({'Control'}, length(control_data_flat), 1); ...
                                           repmat({'Active'}, length(active_data_flat), 1)];
        day_labels_all = [day_labels_all; repmat({day_labels{day_idx}}, length(control_data_flat) + length(active_data_flat), 1)];
        
        % Inter-group statistical analysis (Control vs. Active for this day)
        [p_value_inter_group, h_stat_inter_group, stats_inter_group] = ranksum(control_data_flat, active_data_flat);
        
        % Record inter-group p-value and other details
        p_values_table{table_row_idx, 1} = light_conditions{light_idx};
        p_values_table{table_row_idx, 2} = day_labels{day_idx};
        p_values_table{table_row_idx, 3} = 'Control vs Active';
        p_values_table{table_row_idx, 4} = 'Ranksum Test';
        p_values_table{table_row_idx, 5} = stats_inter_group.ranksum; % Test statistic
        p_values_table{table_row_idx, 6} = length(control_data_flat) + length(active_data_flat); % Sample size
        p_values_table{table_row_idx, 7} = p_value_inter_group;
        table_row_idx = table_row_idx + 1;
        
        % Add significance annotation for inter-group comparison
        if p_value_inter_group < 0.05
            add_significance_bracket((day_idx-1)*2+1, (day_idx-1)*2+2, global_max*1.05, p_value_inter_group);
        end
    end
    
    % Create boxplot for current light condition
    boxplot(data_all_days, {day_labels_all, group_labels_all}, ...
        'Labels', strcat(day_labels_all, '-', group_labels_all), ...
        'LabelOrientation', 'inline');

    hold on;

    h = findobj(gca, 'Tag', 'Box');

    for day_idx = 1:num_days  % Loop over days
        control_data_flat = squeeze(control_data(:, light_idx, day_idx, :));
        active_data_flat = squeeze(active_data(:, light_idx, day_idx, :));

        set(h((day_idx-1)*2+1), 'Color', color_set_active, 'LineWidth', 1.5);
        set(h((day_idx-1)*2+2), 'Color', color_set_control, 'LineWidth', 1.5);
        
        for subj = 1:size(control_data, 1)
            marker_idx = mod(subj-1, num_marker_styles) + 1;
            scatter((day_idx-1)*2+ones(size(control_data_flat(subj, :))), control_data_flat(subj, :), 50, marker_styles{marker_idx}, ...
            'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.1);
        end
        for subj = 1:size(active_data, 1)
            marker_idx = mod(subj-1, num_marker_styles) + 1;
            scatter((day_idx-1)*2+1+ones(size(control_data_flat(subj, :))), active_data_flat(subj, :), 50, marker_styles{marker_idx}, ...
            'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.1);
        end
    end
    
    % Add title and formatting
    title([light_conditions{light_idx} ' Comparison Across Days'], ...
        'FontSize', fontsize_title, 'FontWeight', 'bold');
    ylabel('Normalized Entrainment Power', 'FontSize', fontsize_labels, 'Interpreter', 'tex');
    xlabel('Day - Group', 'FontSize', fontsize_labels, 'Interpreter', 'tex');
    
    ax = gca;
    ax.FontSize = fontsize_axis;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    
    % Intra-group comparisons across days
    for group_idx = 1:2
        if group_idx == 1
            current_data = control_data;
            group_name = 'Control';
        else
            current_data = active_data;
            group_name = 'Active';
        end
        
        % Extract data across days for the current light condition
        group_data = [];
        day_labels_in_group = [];
        for day_idx = 1:num_days
            day_data = squeeze(current_data(:, light_idx, day_idx, :));
            group_data = [group_data; day_data(:)];
            day_labels_in_group = [day_labels_in_group; repmat({day_labels{day_idx}}, length(day_data(:)), 1)];
        end
        
        % Perform Kruskal-Wallis test
        [p_kw, tbl_kw, stats_kw] = kruskalwallis(group_data, day_labels_in_group, 'off');
        chi2_stat = tbl_kw{2, 5}; % Extract the chi-squared statistic from the Kruskal-Wallis table
        
        % Record Kruskal-Wallis results
        p_values_table{table_row_idx, 1} = light_conditions{light_idx};
        p_values_table{table_row_idx, 2} = group_name;
        p_values_table{table_row_idx, 3} = 'Across Days';
        p_values_table{table_row_idx, 4} = 'Kruskal-Wallis Test';
        p_values_table{table_row_idx, 5} = chi2_stat; % Chi-squared statistic
        p_values_table{table_row_idx, 6} = length(group_data); % Sample size
        p_values_table{table_row_idx, 7} = p_kw;
        table_row_idx = table_row_idx + 1;
        
        % Post-hoc Dunn's test if significant
        if p_kw < 0.05
            c = multcompare(stats_kw, 'Display', 'off');
            for i = 1:size(c, 1)
                p_value_pair = c(i, 6);  % p-value
                if p_value_pair < 0.05
                    day1 = day_labels{c(i, 1)};
                    day2 = day_labels{c(i, 2)};
                    % Record p-value
                    p_values_table{table_row_idx, 1} = light_conditions{light_idx};
                    p_values_table{table_row_idx, 2} = group_name;
                    p_values_table{table_row_idx, 3} = [day1 ' vs ' day2];
                    p_values_table{table_row_idx, 4} = "Dunn's Test";
                    p_values_table{table_row_idx, 5} = 'N/A'; % No specific test statistic for Dunn's
                    p_values_table{table_row_idx, 6} = 'N/A'; % Dunn's test uses ranks from KW
                    p_values_table{table_row_idx, 7} = p_value_pair;
                    table_row_idx = table_row_idx + 1;

                    % Add significance annotation
                    x1 = find(strcmp(day_labels, day1));
                    x2 = find(strcmp(day_labels, day2));
                    if strcmp(group_name, 'Control')
                        add_significance_bracket((x1-1)*2+1, (x2-1)*2+1, global_max*(1.15+randi(100)/500), p_value_pair);
                    else
                        add_significance_bracket((x1-1)*2+2, (x2-1)*2+2, global_max*(1.15+randi(100)/500), p_value_pair);
                    end
                    
                end
            end
        end
    end
    ylim(y_range);
    % Save figure
    image_path_dir = fullfile(path2save, "light_comparisons");
    if ~exist(image_path_dir, 'dir')
        mkdir(image_path_dir);
    end
    image_path = fullfile(image_path_dir, strcat("BoxPlot_", light_conditions{light_idx}, ".png"));
    exportgraphics(fig, image_path, 'Resolution', 300);
    image_path = fullfile(image_path_dir, strcat("BoxPlot_", light_conditions{light_idx}, ".pdf"));
    exportgraphics(fig, image_path, 'Resolution', 300);
end

% Save p-values with additional statistical details
csv_path = fullfile(path2save, "light_comparisons", "statistical_analysis_results.csv");
p_values_table_headers = {'Light Condition', 'Day/Group', 'Comparison', 'Test Type', 'Test Statistic', 'Sample Size', 'p-value'};
p_values_table_final = cell2table(p_values_table, 'VariableNames', p_values_table_headers);
writetable(p_values_table_final, csv_path);

%% Active vs control over days seperated by brain region
clc;
close all;

% Define brain regions based on channel locations
brain_regions = struct();
brain_regions.Frontal = [1, 2];    % FP1, FP2
brain_regions.Central = [3, 4];    % C3, C4  
brain_regions.Temporal = [5, 6];   % TP7, TP8
brain_regions.Occipital = [7, 8];  % O1, O2

region_names = fieldnames(brain_regions);
num_regions = length(region_names);

% Create path for saving results
path2save = fullfile("results", "entrainment_brain_region_analysis");
if ~exist(path2save, 'dir')
    mkdir(path2save);
end

% Set up colors and font sizes
color_set_control = [81, 90, 95]/255;
color_set_active = [77,139,145]/255;
color_set1 = [230, 159, 0]/255;  % Orange
color_set2 = [86, 180, 233]/255; % Sky Blue
color_set3 = [0, 158, 115]/255;  % Bluish Green

fontsize_title = 16;
fontsize_labels = 14;
fontsize_axis = 12;

% Extract and organize data by brain regions - keep all individual channels instead of averaging
control_data_regions = zeros(numel(control_subjects), numel(stim_types), numel(recording_days), sum(cellfun(@length, struct2cell(brain_regions))));
active_data_regions = zeros(numel(active_subjects), numel(stim_types), numel(recording_days), sum(cellfun(@length, struct2cell(brain_regions))));

% Store individual channels from each brain region instead of averaging
channel_to_region_map = [];
region_channel_labels = {};
global_channel_idx = 1;

for region_idx = 1:num_regions
    region_channels = brain_regions.(region_names{region_idx});
    
    for ch_idx = 1:length(region_channels)
        % Control group - store individual channel data
        control_data_regions(:, :, :, global_channel_idx) = avg_power_over_days.control(:, :, :, region_channels(ch_idx));
        
        % Active group - store individual channel data  
        active_data_regions(:, :, :, global_channel_idx) = avg_power_over_days.active(:, :, :, region_channels(ch_idx));
        
        % Keep track of which region each channel belongs to
        channel_to_region_map = [channel_to_region_map, region_idx];
        region_channel_labels{global_channel_idx} = sprintf('%s_Ch%d', region_names{region_idx}, region_channels(ch_idx));
        
        global_channel_idx = global_channel_idx + 1;
    end
end

% Normalize by no-light condition
control_data_regions = control_data_regions ./ control_data_regions(:, 1, :, :);
active_data_regions = active_data_regions ./ active_data_regions(:, 1, :, :);

% Find global min and max for consistent scaling
all_data_regions = cat(5, control_data_regions, active_data_regions);
global_min = min(all_data_regions(:));
global_max = max(all_data_regions(:));
y_range = [global_min, global_max*1.4];

% Define marker styles for subjects
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
num_marker_styles = length(marker_styles);

% Regional comparison plots for each light condition and day
p_values_table = {};
table_row_idx = 1;

for light_idx = 2:numel(stim_types)  % Skip no_light condition
    for day_idx = 1:numel(recording_days)
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        
        for region_idx = 1:num_regions
            subplot(2, 2, region_idx);
            hold on;
            
            % Extract data for current condition and region - all channels from this region
            region_channel_indices = find(channel_to_region_map == region_idx);
            control_region_data = squeeze(control_data_regions(:, light_idx, day_idx, region_channel_indices));
            active_region_data = squeeze(active_data_regions(:, light_idx, day_idx, region_channel_indices));
            
            % Flatten data while keeping track of subject identity
            control_flat = [];
            control_subject_labels = [];
            active_flat = [];
            active_subject_labels = [];
            
            % Control group - keep track of which subject each data point belongs to
            for subj = 1:size(control_region_data, 1)
                subj_data = control_region_data(subj, :);
                control_flat = [control_flat, subj_data];
                control_subject_labels = [control_subject_labels, subj * ones(1, length(subj_data))];
            end
            
            % Active group - keep track of which subject each data point belongs to
            for subj = 1:size(active_region_data, 1)
                subj_data = active_region_data(subj, :);
                active_flat = [active_flat, subj_data];
                active_subject_labels = [active_subject_labels, subj * ones(1, length(subj_data))];
            end
            
            % Create box plots
            boxplot([control_flat(:), active_flat(:)], ...
                [ones(size(control_flat(:))); 2*ones(size(active_flat(:)))], ...
                'Labels', {'Control', 'Active'}, 'Width', 0.7);
            
            h = findobj(gca, 'Tag', 'Box');
            patch(get(h(1), 'XData'), get(h(1), 'YData'), color_set2, 'FaceAlpha', 0.2);
            patch(get(h(2), 'XData'), get(h(2), 'YData'), color_set1, 'FaceAlpha', 0.2);
            
            % Plot individual channel data points with subject-specific markers
            for i = 1:length(control_flat)
                subj = control_subject_labels(i);
                marker_idx = mod(subj-1, num_marker_styles) + 1;
                
                % Control group
                scatter(1 + randn(1)*0.05, control_flat(i), 60, ...
                    marker_styles{marker_idx}, ...
                    'MarkerEdgeColor', 'none', ...
                    'MarkerFaceColor', color_set_control, ...
                    'MarkerFaceAlpha', 0.7);
            end
            
            for i = 1:length(active_flat)
                subj = active_subject_labels(i);
                marker_idx = mod(subj-1, num_marker_styles) + 1;
                
                % Active group
                scatter(2 + randn(1)*0.05, active_flat(i), 60, ...
                    marker_styles{marker_idx}, ...
                    'MarkerEdgeColor', 'none', ...
                    'MarkerFaceColor', color_set_active, ...
                    'MarkerFaceAlpha', 0.7);
            end
            
            % Statistical test
            [p_value, ~] = ranksum(control_flat, active_flat);
            add_significance_bracket(1, 2, global_max*1.05, p_value);
            
            % Store p-values in table
            p_values_table{table_row_idx, 1} = light_conditions{light_idx};
            p_values_table{table_row_idx, 2} = day_labels{day_idx};
            p_values_table{table_row_idx, 3} = region_names{region_idx};
            p_values_table{table_row_idx, 4} = p_value;
            table_row_idx = table_row_idx + 1;
            
            % Formatting
            title(sprintf('%s Region', region_names{region_idx}), ...
                'FontSize', fontsize_title-2, 'FontWeight', 'bold', 'Interpreter', 'tex');
            ylabel('Normalized Entrainment Power', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
            xlabel('Groups', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
            ylim(y_range);
            
            ax = gca;
            ax.FontSize = fontsize_axis-1;
            ax.LineWidth = 1.5;
            ax.Box = 'on';
        end
        
        sgtitle(sprintf('Brain Regional Analysis: %s, %s', ...
            light_conditions{light_idx}, day_labels{day_idx}), ...
            'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');
        
        % Save figure
        image_path_dir = fullfile(path2save, "regional_comparisons");
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        image_path = fullfile(image_path_dir, sprintf("Regional_Analysis_%s_%s.pdf", ...
            light_conditions{light_idx}, day_labels{day_idx}));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, sprintf("Regional_Analysis_%s_%s.png", ...
            light_conditions{light_idx}, day_labels{day_idx}));
        exportgraphics(fig, image_path, 'Resolution', 300);
    end
end

% Time course analysis by brain region
for light_idx = 2:numel(stim_types)  % Skip no_light condition
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    
    for region_idx = 1:num_regions
        subplot(2, 2, region_idx);
        hold on;
        
        % Extract data for this region across all days
        region_channel_indices = find(channel_to_region_map == region_idx);
        
        % Calculate means and standard errors across days using all channels
        control_means = [];
        active_means = [];
        control_se = [];
        active_se = [];
        
        for day_idx = 1:numel(recording_days)
            control_day_data = squeeze(control_data_regions(:, light_idx, day_idx, region_channel_indices));
            active_day_data = squeeze(active_data_regions(:, light_idx, day_idx, region_channel_indices));
            
            % Flatten all channel data for this day
            control_flat = control_day_data(:);
            active_flat = active_day_data(:);
            
            control_means = [control_means, mean(control_flat)];
            active_means = [active_means, mean(active_flat)];
            control_se = [control_se, std(control_flat) / sqrt(length(control_flat))];
            active_se = [active_se, std(active_flat) / sqrt(length(active_flat))];
        end
        
        % Plot individual subject trajectories (faded) - average across channels per subject
        for subj = 1:size(control_data_regions, 1)
            control_trajectory = squeeze(mean(control_data_regions(subj, light_idx, :, region_channel_indices), 4));
            plot(1:numel(recording_days), control_trajectory, '-', ...
                'Color', [color_set_control, 0.3], 'LineWidth', 1);
        end
        
        for subj = 1:size(active_data_regions, 1)
            active_trajectory = squeeze(mean(active_data_regions(subj, light_idx, :, region_channel_indices), 4));
            plot(1:numel(recording_days), active_trajectory, '-', ...
                'Color', [color_set_active, 0.3], 'LineWidth', 1);
        end
        
        % Plot mean trajectories with error bars
        errorbar(1:numel(recording_days), control_means, control_se, '-o', ...
            'Color', color_set_control, 'LineWidth', 3, 'MarkerSize', 8, ...
            'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'none');
        
        errorbar(1:numel(recording_days), active_means, active_se, '-s', ...
            'Color', color_set_active, 'LineWidth', 3, 'MarkerSize', 8, ...
            'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'none');
        
        % Formatting
        title(sprintf('%s Region', region_names{region_idx}), ...
            'FontSize', fontsize_title-2, 'FontWeight', 'bold', 'Interpreter', 'tex');
        xlabel('Time Point', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
        ylabel('Normalized Entrainment Power', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
        
        set(gca, 'XTick', 1:numel(recording_days), 'XTickLabel', day_labels);
        ylim([global_min*0.8, global_max*0.4]);
        grid on;
        
        if region_idx == 1
            legend('Control Individual', 'Active Individual', 'Control Mean', 'Active Mean', ...
                'Location', 'best', 'FontSize', fontsize_axis-1, 'Interpreter', 'tex');
        end
        
        ax = gca;
        ax.FontSize = fontsize_axis-1;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
    end
    
    sgtitle(sprintf('Temporal Evolution by Brain Region: %s', light_conditions{light_idx}), ...
        'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Save figure
    image_path_dir = fullfile(path2save, "temporal_evolution");
    if ~exist(image_path_dir, 'dir')
        mkdir(image_path_dir);
    end
    image_path = fullfile(image_path_dir, sprintf("Temporal_Evolution_%s.pdf", light_conditions{light_idx}));
    exportgraphics(fig, image_path, 'Resolution', 300);
    image_path = fullfile(image_path_dir, sprintf("Temporal_Evolution_%s.png", light_conditions{light_idx}));
    exportgraphics(fig, image_path, 'Resolution', 300);
end

% Summary heatmap showing effect sizes by region
effect_sizes = zeros(numel(stim_types)-1, numel(recording_days), num_regions);
p_values_matrix = zeros(numel(stim_types)-1, numel(recording_days), num_regions);

for light_idx = 2:numel(stim_types)
    for day_idx = 1:numel(recording_days)
        for region_idx = 1:num_regions
            % Get all channels for this region
            region_channel_indices = find(channel_to_region_map == region_idx);
            control_data = squeeze(control_data_regions(:, light_idx, day_idx, region_channel_indices));
            active_data = squeeze(active_data_regions(:, light_idx, day_idx, region_channel_indices));
            
            % Flatten all channel data
            control_flat = control_data(:);
            active_flat = active_data(:);
            
            % Calculate Cohen's d (effect size)
            pooled_std = sqrt(((length(control_flat)-1)*var(control_flat) + ...
                              (length(active_flat)-1)*var(active_flat)) / ...
                              (length(control_flat) + length(active_flat) - 2));
            cohens_d = (mean(active_flat) - mean(control_flat)) / pooled_std;
            
            effect_sizes(light_idx-1, day_idx, region_idx) = cohens_d;
            
            % Get p-value
            [p_val, ~] = ranksum(control_flat, active_flat);
            p_values_matrix(light_idx-1, day_idx, region_idx) = p_val;
        end
    end
end

% Create heatmap figure
fig = figure;
set(gcf, 'Position', get(0, 'Screensize'));

for region_idx = 1:num_regions
    subplot(2, 2, region_idx);
    
    region_effects = squeeze(effect_sizes(:, :, region_idx));
    region_pvals = squeeze(p_values_matrix(:, :, region_idx));
    
    imagesc(region_effects);
    colormap('abyss');
    colorbar;
    clim([-1, 1]);
    
    % Mark significant effects
    hold on;
    [row, col] = find(region_pvals < 0.05);
    if ~isempty(row)
        plot(col, row, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    title(sprintf('%s Region', region_names{region_idx}), ...
        'FontSize', fontsize_title-2, 'FontWeight', 'bold', 'Interpreter', 'tex');
    xlabel('Day', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
    ylabel('Light Condition', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
    
    set(gca, 'XTick', 1:numel(recording_days), 'XTickLabel', day_labels);
    set(gca, 'YTick', 1:numel(stim_types)-1, 'YTickLabel', light_conditions(2:end));
    
    ax = gca;
    ax.FontSize = fontsize_axis-1;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
end

sgtitle('Effect Sizes (Cohen''s d) by Brain Region (* p < 0.05)', ...
    'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');

% Save heatmap figure
image_path_dir = fullfile(path2save, "effect_size_heatmaps");
if ~exist(image_path_dir, 'dir')
    mkdir(image_path_dir);
end
image_path = fullfile(image_path_dir, "Regional_Effect_Sizes.pdf");
exportgraphics(fig, image_path, 'Resolution', 300);
image_path = fullfile(image_path_dir, "Regional_Effect_Sizes.png");
exportgraphics(fig, image_path, 'Resolution', 300);

% Save statistical results
csv_path = fullfile(path2save, "regional_statistical_analysis.csv");
p_values_table_headers = {'Light Condition', 'Day', 'Brain Region', 'p-value (control vs active)'};
p_values_table_final = cell2table(p_values_table, 'VariableNames', p_values_table_headers);
writetable(p_values_table_final, csv_path);

% Save effect size data
effect_size_results = {};
result_idx = 1;
for light_idx = 2:numel(stim_types)
    for day_idx = 1:numel(recording_days)
        for region_idx = 1:num_regions
            effect_size_results{result_idx, 1} = light_conditions{light_idx};
            effect_size_results{result_idx, 2} = day_labels{day_idx};
            effect_size_results{result_idx, 3} = region_names{region_idx};
            effect_size_results{result_idx, 4} = effect_sizes(light_idx-1, day_idx, region_idx);
            effect_size_results{result_idx, 5} = p_values_matrix(light_idx-1, day_idx, region_idx);
            result_idx = result_idx + 1;
        end
    end
end

effect_size_headers = {'Light Condition', 'Day', 'Brain Region', 'Effect Size (Cohens d)', 'p-value'};
effect_size_table = cell2table(effect_size_results, 'VariableNames', effect_size_headers);
csv_path_effects = fullfile(path2save, "regional_effect_sizes.csv");
writetable(effect_size_table, csv_path_effects);

fprintf('Brain region analysis completed. Results saved to: %s\n', path2save);

% Summary statistics by region
fprintf('\n=== BRAIN REGION ANALYSIS SUMMARY ===\n');
for region_idx = 1:num_regions
    fprintf('\n%s Region:\n', region_names{region_idx});
    region_channels = brain_regions.(region_names{region_idx});
    ch_list = {'FP1', 'FP2', 'C3', 'C4', 'TP7', 'TP8', 'O1', 'O2'};
    fprintf('Channels: %s\n', strjoin(ch_list(region_channels), ', '));
    
    % Find most significant effects
    region_pvals = squeeze(p_values_matrix(:, :, region_idx));
    [min_p, min_idx] = min(region_pvals(:));
    [light_idx, day_idx] = ind2sub(size(region_pvals), min_idx);
    
    if min_p < 0.05
        fprintf('Most significant effect: %s, %s (p = %.4f)\n', ...
            light_conditions{light_idx+1}, day_labels{day_idx}, min_p);
        fprintf('Effect size: %.3f\n', effect_sizes(light_idx, day_idx, region_idx));
    else
        fprintf('No significant effects found\n');
    end
end

% Time-resolved 60Hz Power Analysis (Normalized by no_light baseline)
clc
close all

path2save = fullfile("results", "Time_resolved_60Hz_power");
if ~exist(path2save, 'dir')
    mkdir(path2save);
end

% Analysis parameters
no_channels = 8;
num_stimulus_parts = 5;  % Divide stimulus into n equal parts

% Initialize storage for normalized power data
% Structure: [constant_light + 10 stimulus parts] x [subjects] x [channels] x [days] x [groups]
time_points = num_stimulus_parts + 1;  % constant_light + 10 stimulus parts
normalized_power_data = struct();
normalized_power_data.control = zeros(time_points, numel(control_subjects), no_channels, numel(recording_days));
normalized_power_data.active = zeros(time_points, numel(active_subjects), no_channels, numel(recording_days));

% Process each subject and condition
for subject_type = ["control", "active"]
    fprintf('Processing %s subjects...\n', subject_type);
    subject_counter = 0;
    subjects = eval(strcat(subject_type, "_subjects"));
    
    for subject = subjects
        subject_counter = subject_counter + 1;
        subject = string(subject(1));
        fprintf('  Subject: %s\n', subject);
        
        recording_day_counter = 0;
        for recording_day = recording_days
            recording_day_counter = recording_day_counter + 1;
            recording_day = string(recording_day(1));
            
            % Get no_light data for normalization
            if isfield(data_struct.(sprintf(subject)), sprintf(recording_day)) && ...
               isfield(data_struct.(sprintf(subject)).(sprintf(recording_day)), 'no_light')
                
                no_light_data = data_struct.(sprintf(subject)).(sprintf(recording_day)).no_light.data_clean;
                
                % Calculate total 60Hz power in no_light condition for normalization
                baseline_power = zeros(1, no_channels);
                for ch = 1:no_channels
                    baseline_power(ch) = calcPower(no_light_data(ch, :)', fs, target_freq_range);
                end
                
                % Process constant_light condition (time point 1)
                if isfield(data_struct.(sprintf(subject)).(sprintf(recording_day)), 'constant_light')
                    constant_light_data = data_struct.(sprintf(subject)).(sprintf(recording_day)).constant_light.data_clean;
                    
                    % Calculate total 60Hz power in constant_light
                    constant_power = zeros(1, no_channels);
                    for ch = 1:no_channels
                        constant_power(ch) = calcPower(constant_light_data(ch, :)', fs, target_freq_range);
                    end
                    
                    % Normalize by baseline and store
                    normalized_power_data.(sprintf(subject_type))(1, subject_counter, :, recording_day_counter) = ...
                        constant_power ./ baseline_power;
                end
                
                % Process stimulus_light condition (time points 2-11)
                if isfield(data_struct.(sprintf(subject)).(sprintf(recording_day)), 'stimulus_light')
                    stimulus_data = data_struct.(sprintf(subject)).(sprintf(recording_day)).stimulus_light.data_clean;
                    [~, total_samples] = size(stimulus_data);
                    
                    % Divide stimulus into 10 equal parts
                    samples_per_part = floor(total_samples / num_stimulus_parts);
                    
                    for part = 1:num_stimulus_parts
                        % Calculate start and end indices for this part
                        start_idx = (part - 1) * samples_per_part + 1;
                        if part == num_stimulus_parts
                            end_idx = total_samples;  % Include remaining samples in last part
                        else
                            end_idx = part * samples_per_part;
                        end
                        
                        % Extract data for this part
                        part_data = stimulus_data(:, start_idx:end_idx);
                        
                        % Calculate 60Hz power for this part
                        part_power = zeros(1, no_channels);
                        for ch = 1:no_channels
                            part_power(ch) = calcPower(part_data(ch, :)', fs, target_freq_range);
                        end
                        
                        % Normalize by baseline and store
                        normalized_power_data.(sprintf(subject_type))(part + 1, subject_counter, :, recording_day_counter) = ...
                            part_power ./ baseline_power;
                    end
                end
                
                fprintf('    Completed day %s\n', recording_day);
            else
                fprintf('    Warning: No baseline data found for %s %s\n', subject, recording_day);
            end
        end
    end
end

% Create plots for each recording day
fprintf('Creating plots...\n');

% Define colors
color_set_control = [81, 90, 95]/255;
color_set_active = [77,139,145]/255;

% Define marker styles for subjects
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
num_marker_styles = length(marker_styles);

% Find global Y-axis limits across all data
global_min = inf;
global_max = -inf;
global_max_with_error = -inf;

for day_idx = 1:numel(recording_days)
    control_data_day = squeeze(normalized_power_data.control(:, :, :, day_idx));
    active_data_day = squeeze(normalized_power_data.active(:, :, :, day_idx));
    
    % Find min and max for this day
    control_valid = control_data_day(~isnan(control_data_day) & control_data_day > 0);
    active_valid = active_data_day(~isnan(active_data_day) & active_data_day > 0);
    
    if ~isempty(control_valid)
        global_min = min(global_min, min(control_valid));
        global_max = max(global_max, max(control_valid));
    end
    
    if ~isempty(active_valid)
        global_min = min(global_min, min(active_valid));
        global_max = max(global_max, max(active_valid));
    end
    
    % Calculate max with error bars for each time point
    for time_point = 1:time_points
        control_time_data = squeeze(control_data_day(time_point, :, :));
        control_time_data = control_time_data(:);
        control_time_data = control_time_data(~isnan(control_time_data) & control_time_data > 0);
        
        active_time_data = squeeze(active_data_day(time_point, :, :));
        active_time_data = active_time_data(:);
        active_time_data = active_time_data(~isnan(active_time_data) & active_time_data > 0);
        
        if ~isempty(control_time_data)
            control_upper = mean(control_time_data) + std(control_time_data) + 1;
            global_max_with_error = max(global_max_with_error, control_upper);
        end
        
        if ~isempty(active_time_data)
            active_upper = mean(active_time_data) + std(active_time_data) + 1;
            global_max_with_error = max(global_max_with_error, active_upper);
        end
    end
end

% Add some padding to the limits
y_range_padding = (global_max - global_min) * 0.05;
y_limits = [global_min - y_range_padding, global_max + y_range_padding];
y_limits_with_error = [global_min - y_range_padding, global_max_with_error];

fprintf('Global Y-axis limits: [%.3f, %.3f]\n', y_limits(1), y_limits(2));
fprintf('Global Y-axis limits with error bars: [%.3f, %.3f]\n', y_limits_with_error(1), y_limits_with_error(2));

for day_idx = 1:numel(recording_days)
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    
    % Extract data for this day
    control_data_day = squeeze(normalized_power_data.control(:, :, :, day_idx));  % time_points x subjects x channels
    active_data_day = squeeze(normalized_power_data.active(:, :, :, day_idx));    % time_points x subjects x channels
    
    % Reshape data for plotting: each time point has num_channels * num_subjects data points
    x_positions = [];
    control_values = [];
    active_values = [];
    
    for time_point = 1:time_points
        % Control group data for this time point
        control_time_data = squeeze(control_data_day(time_point, :, :));  % subjects x channels
        control_time_data = control_time_data(:);  % Flatten to vector
        control_time_data = control_time_data(~isnan(control_time_data) & control_time_data > 0);  % Remove invalid values
        
        % Active group data for this time point  
        active_time_data = squeeze(active_data_day(time_point, :, :));    % subjects x channels
        active_time_data = active_time_data(:);    % Flatten to vector
        active_time_data = active_time_data(~isnan(active_time_data) & active_time_data > 0);    % Remove invalid values
        
        % Create x positions with slight jitter for visibility
        control_x = time_point - 0.1 + randn(length(control_time_data), 1) * 0.02;
        active_x = time_point + 0.1 + randn(length(active_time_data), 1) * 0.02;
        
        % Plot data points
        scatter(control_x, control_time_data, 50, 'o', ...
            'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.6);
        scatter(active_x, active_time_data, 50, 's', ...
            'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.6);
        
        % Store for group statistics
        x_positions = [x_positions; control_x; active_x];
        control_values = [control_values; control_time_data];
        active_values = [active_values; active_time_data];
        
        % Calculate and plot group means
        if ~isempty(control_time_data)
            control_mean = median(control_time_data);
            plot(time_point - 0.1, control_mean, 'o', 'MarkerSize', 12, ...
                'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
        end
        
        if ~isempty(active_time_data)
            active_mean = median(active_time_data);
            plot(time_point + 0.1, active_mean, 's', 'MarkerSize', 12, ...
                'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
        end
    end
    
    % Connect means with lines
    control_means = [];
    active_means = [];
    valid_time_points = [];
    
    for time_point = 1:time_points
        control_time_data = squeeze(control_data_day(time_point, :, :));
        control_time_data = control_time_data(:);
        control_time_data = control_time_data(~isnan(control_time_data) & control_time_data > 0);
        
        active_time_data = squeeze(active_data_day(time_point, :, :));
        active_time_data = active_time_data(:);
        active_time_data = active_time_data(~isnan(active_time_data) & active_time_data > 0);
        
        if ~isempty(control_time_data) && ~isempty(active_time_data)
            control_means = [control_means, median(control_time_data)];
            active_means = [active_means, median(active_time_data)];
            valid_time_points = [valid_time_points, time_point];
        end
    end
    
    if length(valid_time_points) > 1
        plot(valid_time_points - 0.1, control_means, '-', 'Color', color_set_control, 'LineWidth', 2);
        plot(valid_time_points + 0.1, active_means, '-', 'Color', color_set_active, 'LineWidth', 2);
    end
    
    % Formatting
    xlabel('Time Points', 'FontSize', 16, 'Interpreter', 'tex');
    ylabel('Normalized 60Hz Power', 'FontSize', 16, 'Interpreter', 'tex');
    title(sprintf('Normalized 60Hz Power Over Time - %s', strrep(char(recording_days{day_idx}), '_', ' ')), ...
        'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Set x-axis labels
    x_labels = {'Constant Light'};
    for i = 1:num_stimulus_parts
        x_labels{end+1} = sprintf('Stimulus %d - %d min', ceil((i-1) * 20/num_stimulus_parts), ceil(i * 20/num_stimulus_parts));
    end
    
    xlim([0.5, time_points + 0.5]);
    set(gca, 'XTick', 1:time_points, 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
    
    % Set consistent Y-axis limits
    ylim(y_limits);
    
    % Add legend
    legend('Control Individual', 'Active Individual', 'Control Median', 'Active Median', ...
        'Control Trajectory', 'Active Trajectory', ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'tex');
    
    grid on;
    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    
    % Save figure
    filename = sprintf('Normalized_60Hz_TimeCourse_%s', char(recording_days{day_idx}));
    exportgraphics(fig, fullfile(path2save, [filename '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(path2save, [filename '.pdf']), 'Resolution', 300);
    
    close(fig);
end

% Create error bar plots for each recording day
fprintf('Creating error bar plots...\n');

for day_idx = 1:numel(recording_days)
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    
    % Extract data for this day
    control_data_day = squeeze(normalized_power_data.control(:, :, :, day_idx));
    active_data_day = squeeze(normalized_power_data.active(:, :, :, day_idx));
    
    % Calculate means and standard deviations for each time point
    control_means = [];
    control_stds = [];
    active_means = [];
    active_stds = [];
    valid_time_points = [];
    
    for time_point = 1:time_points
        % Control group data for this time point
        control_time_data = squeeze(control_data_day(time_point, :, :));
        control_time_data = control_time_data(:);
        control_time_data = control_time_data(~isnan(control_time_data) & control_time_data > 0);
        
        % Active group data for this time point  
        active_time_data = squeeze(active_data_day(time_point, :, :));
        active_time_data = active_time_data(:);
        active_time_data = active_time_data(~isnan(active_time_data) & active_time_data > 0);
        
        if ~isempty(control_time_data) && ~isempty(active_time_data)
            control_means = [control_means, mean(control_time_data)];
            control_stds = [control_stds, std(control_time_data)];
            active_means = [active_means, mean(active_time_data)];
            active_stds = [active_stds, std(active_time_data)];
            valid_time_points = [valid_time_points, time_point];
        end
    end
    
    % Plot error bars
    if ~isempty(valid_time_points)
        errorbar(valid_time_points - 0.1, control_means, control_stds, 'o-', ...
            'Color', color_set_control, 'LineWidth', 3, 'MarkerSize', 8, ...
            'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'k', 'CapSize', 8);
        errorbar(valid_time_points + 0.1, active_means, active_stds, 's-', ...
            'Color', color_set_active, 'LineWidth', 3, 'MarkerSize', 8, ...
            'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'k', 'CapSize', 8);
    end
    
    % Formatting
    xlabel('Time Points', 'FontSize', 16, 'Interpreter', 'tex');
    ylabel('Normalized 60Hz Power', 'FontSize', 16, 'Interpreter', 'tex');
    title(sprintf('Normalized 60Hz Power Over Time (Mean ± SD) - %s', strrep(char(recording_days{day_idx}), '_', ' ')), ...
        'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Set x-axis labels
    x_labels = {'Constant Light'};
    for i = 1:num_stimulus_parts
        x_labels{end+1} = sprintf('Stimulus %d - %d min', ceil((i-1) * 20/num_stimulus_parts), ceil(i * 20/num_stimulus_parts));
    end
    
    xlim([0.5, time_points + 0.5]);
    set(gca, 'XTick', 1:time_points, 'XTickLabel', x_labels, 'XTickLabelRotation', 45);
    
    % Set consistent Y-axis limits for error bar plots
    ylim(y_limits_with_error);
    
    % Add legend
    legend('Control (Mean ± SD)', 'Active (Mean ± SD)', ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'tex');
    
    grid on;
    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    
    % Save figure
    filename = sprintf('Normalized_60Hz_TimeCourse_ErrorBars_%s', char(recording_days{day_idx}));
    exportgraphics(fig, fullfile(path2save, [filename '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(path2save, [filename '.pdf']), 'Resolution', 300);
    
    close(fig);
end

% Statistical analysis
fprintf('Performing statistical analysis...\n');
stats_results = {};
result_idx = 1;

for day_idx = 1:numel(recording_days)
    control_data_day = squeeze(normalized_power_data.control(:, :, :, day_idx));
    active_data_day = squeeze(normalized_power_data.active(:, :, :, day_idx));
    
    for time_point = 1:time_points
        % Extract data for this time point
        control_time_data = squeeze(control_data_day(time_point, :, :));
        control_time_data = control_time_data(:);
        control_time_data = control_time_data(~isnan(control_time_data) & control_time_data > 0);
        
        active_time_data = squeeze(active_data_day(time_point, :, :));
        active_time_data = active_time_data(:);
        active_time_data = active_time_data(~isnan(active_time_data) & active_time_data > 0);
        
        % Perform statistical test
        if length(control_time_data) >= 3 && length(active_time_data) >= 3
            [p_value, ~] = ranksum(control_time_data, active_time_data);
            
            % Determine time point label
            if time_point == 1
                time_label = 'Constant Light';
            else
                time_label = sprintf('%d%% to %d%% of Stimulus', (time_point - 1) *100/num_stimulus_parts, time_point *100/num_stimulus_parts);
            end
            
            % Store results
            stats_results{result_idx, 1} = strrep(char(recording_days{day_idx}), '_', ' ');
            stats_results{result_idx, 2} = time_label;
            stats_results{result_idx, 3} = median(control_time_data);
            stats_results{result_idx, 4} = std(control_time_data);
            stats_results{result_idx, 5} = median(active_time_data);
            stats_results{result_idx, 6} = std(active_time_data);
            stats_results{result_idx, 7} = p_value;
            stats_results{result_idx, 8} = length(control_time_data);
            stats_results{result_idx, 9} = length(active_time_data);
            
            result_idx = result_idx + 1;
        end
    end
end

% Save statistical results
if ~isempty(stats_results)
    headers = {'Day', 'Time_Point', 'Control_Median', 'Control_Std', 'Active_Median', 'Active_Std', ...
              'p_value', 'N_Control', 'N_Active'};
    stats_table = cell2table(stats_results, 'VariableNames', headers);
    csv_path = fullfile(path2save, 'time_course_statistics.csv');
    writetable(stats_table, csv_path);
    fprintf('Statistical results saved to: %s\n', csv_path);
end

fprintf('Time-resolved 60Hz power analysis completed!\n');