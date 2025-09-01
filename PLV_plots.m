%% PLV - inter channels
clc
close all
no_channels = 8;
active_subjects = {"SU1_24_JEz",...
    "SU1_24_4V9",...
    "SU1_24_SOH",...
    "SU1_24_jWX",...
    "SU1_24_Q0R",...
    "SU1_24_jp6"};
control_subjects = {"SU1_24_UzK",...
    "SU1_24_n4X",...
    "SU1_24_RDj",...
    "SU1_24_gTo",...
    "SU1_24_nub",...
    "SU1_24_u1W"};
recording_days = {"day_1", ...
    "day_5", ...
    "day_19"};
stim_types = {"no_light",...
    "constant_light",...
    "stimulus_light"};


plv_values_active =  zeros(numel(stim_types), numel(recording_days), numel(active_subjects), no_channels, no_channels); % 3*3*6*8*8
plv_values_control =  zeros(numel(stim_types), numel(recording_days), numel(control_subjects), no_channels, no_channels); % 3*3*6*8*8

for subject_type = ["control", "active"]
    subject_counter = 0;
    for subject = eval(strcat(subject_type, "_subjects"))
        subject_counter = subject_counter + 1;
        recording_day_counter = 0;
        for recording_day = recording_days
            recording_day_counter = recording_day_counter + 1;
            stim_type_counter = 0;
            for stim_type = stim_types
                stim_type_counter = stim_type_counter + 1;
                
                subject = string(subject(1));
                recording_day = string(recording_day(1));
                stim_type = string(stim_type(1));
                
                eegData = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
                [plv_timeresolved, time_vector] = calculate_timeresolved_plv(eegData, fs, [59.3, 59.4], 10, 9);
                plv = squeeze(mean(plv_timeresolved, 3));
                eval(strcat("plv_values_", subject_type, "(stim_type_counter, recording_day_counter, subject_counter, :, :) = plv;"));
                % visualize_timeresolved_plv(plv_timeresolved, time_vector, [7, 8])
                % title(strcat(recording_day, "_", stim_type,"_",  subject_type), 'Interpreter', 'tex')

            end
        end
    end
end

stim_type_counter = 0;
for stim_type = stim_types
    stim_type_counter = stim_type_counter + 1;
    plv_values_active(stim_type_counter, :, :, :, :) = plv_values_active(stim_type_counter, :, :, :, :)./plv_values_active(1, :, :, :, :);
    plv_values_control(stim_type_counter, :, :, :, :) = plv_values_control(stim_type_counter, :, :, :, :)./plv_values_control(1, :, :, :, :);
end

% for i = 2:size(plv_values_active, 4)
%     for j = 1:i-1
%         plv_values_control(:, :, :, i, j) = 0;
%         plv_values_active(:, :, :, i, j) = 0;
%     end
% end

data = cat(6, plv_values_control, plv_values_active); % light_condition * days * subjects * channels * subject_type


data_disc = "light_condition * days * subjects * channels * subject_type";
path2save = fullfile("results", "PLV");
save(fullfile("results", "PLV", 'PLV_Values'), 'data', 'data_disc' , '-V7.3');


%% Descriptive analysis for power of 60Hz (Mean and Std of different conditions)
clc
close all

path2save = fullfile("results", "PLV");

group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};

results = {};

subject_type_no = 0;
for subject_type = ["control", "active"]
    subject_type_no = subject_type_no+1;
    recording_day_counter = 0;
    for recording_day = recording_days
        recording_day_counter = recording_day_counter+1;
        subject_counter = 0;
        data_group = zeros(size(data, 1), size(data, 3), size(data, 4), size(data, 5)); % light_condition * subjects * channels * channels
        for subject = eval(strcat(subject_type, "_subjects"))
            subject_counter = subject_counter+1;
            data_points = data(:, recording_day_counter, subject_counter, :, :, subject_type_no);
            data_group(:, subject_counter, :, :) = data_points; % light_condition * subjects * channels * channels
        end
        % data_group = squeeze(mean(data_group, 1));
        light_condition_counter = 0;
        for light_condition = light_conditions
            light_condition_counter = light_condition_counter+1;
            group_specific_light = squeeze(data_group(light_condition_counter, :, :, :));
            under_diag_values = [];
            for i = 2:size(group_specific_light, 2)
                for j = 1:i-1
                    under_diag_values = [group_specific_light(:, i, j), under_diag_values];
                end
            end
            mean_group = mean(under_diag_values, "all");
            std_group = std(under_diag_values, [], "all");
            min_group = min(under_diag_values, [], "all");
            max_group = max(under_diag_values, [], "all");
            results = [results; {char(subject_type), day_labels{recording_day_counter}, ...
                light_conditions{light_condition_counter}, mean_group, std_group, min_group, max_group}];
        end
    end
end

results_table = cell2table(results, ...
    'VariableNames', {'Group', 'Day', 'LightCondition', 'MeanGroup', 'StdGroup', 'MinGroup', 'MaxGroup'});
writetable(results_table, fullfile(path2save, "descriptive_analysis.csv"));
disp('Results saved to group_statistics.csv');

%% PLV plots
clc
close all

% Set default font sizes
fontsize_title = 16;
fontsize_labels = 14;
fontsize_axis = 12;

% Define colorblind-friendly colors
color_set1 = [230, 159, 0]/255;  % Orange
color_set2 = [86, 180, 233]/255; % Sky Blue
color_set3 = [0, 158, 115]/255;  % Bluish Green

% Define marker styles for subjects
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
num_marker_styles = length(marker_styles);
num_days = 3;
num_light_conditions = 3;
path2save = fullfile("results", "PLV");

%%
fig = figure;
set(gcf, 'Position', get(0, 'Screensize'));
visualize_plv_data(plv_values_active, stim_types, recording_days, active_subjects)
ax = gca;
ax.FontSize = fontsize_axis;
ax.LineWidth = 1.5;
ax.Box = 'on';
image_path_dir = fullfile(path2save);
if ~exist(image_path_dir, 'dir')
    mkdir(image_path_dir);
end
image_path = fullfile(image_path_dir, strcat("PLV_between_channels_active", ".pdf"));
exportgraphics(fig, image_path, 'Resolution', 300);
image_path = fullfile(image_path_dir, strcat("PLV_between_channels_active", ".png"));
exportgraphics(fig, image_path, 'Resolution', 300);

fig = figure;
set(gcf, 'Position', get(0, 'Screensize'));
visualize_plv_data(plv_values_control, stim_types, recording_days, control_subjects)
ax = gca;
ax.FontSize = fontsize_axis;
ax.LineWidth = 1.5;
ax.Box = 'on';
image_path_dir = fullfile(path2save);
if ~exist(image_path_dir, 'dir')
    mkdir(image_path_dir);
end
image_path = fullfile(image_path_dir, strcat("PLV_between_channels_control", ".pdf"));
exportgraphics(fig, image_path, 'Resolution', 300);
image_path = fullfile(image_path_dir, strcat("PLV_between_channels_control", ".png"));
exportgraphics(fig, image_path, 'Resolution', 300);

fig = figure;
set(gcf, 'Position', get(0, 'Screensize'));
compare_plv_groups(plv_values_active, plv_values_control, stim_types, recording_days, 0.05)
ax = gca;
ax.FontSize = fontsize_axis;
ax.LineWidth = 1.5;
ax.Box = 'on';
image_path_dir = fullfile(path2save);
if ~exist(image_path_dir, 'dir')
    mkdir(image_path_dir);
end
image_path = fullfile(image_path_dir, strcat("PLV_difference_control_and_active", ".pdf"));
exportgraphics(fig, image_path, 'Resolution', 300);
image_path = fullfile(image_path_dir, strcat("PLV_difference_control_and_active", ".png"));
exportgraphics(fig, image_path, 'Resolution', 300);

fig = figure;
set(gcf, 'Position', get(0, 'Screensize'));
compare_network_metrics(plv_values_active, plv_values_control, stim_types, recording_days, 'proportion')
image_path_dir = fullfile(path2save);
if ~exist(image_path_dir, 'dir')
    mkdir(image_path_dir);
end
image_path = fullfile(image_path_dir, strcat("PLV_matrix_connectivity", ".pdf"));
exportgraphics(fig, image_path, 'Resolution', 300);
image_path = fullfile(image_path_dir, strcat("PLV_matrix_connectivity", ".png"));
exportgraphics(fig, image_path, 'Resolution', 300);

% Control vs Active comparisons (mean of PLV of channels)
data = cat(6, plv_values_control, plv_values_active); % light_condition * days * subjects * channels * channels * subject_type

min_value = 0.4;
max_value = 0;

% Marker styles for subjects
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
num_marker_styles = length(marker_styles);

for day_idx = 1:num_days
    for light_idx = 1:num_light_conditions
        data_control = mean(squeeze(data(light_idx, day_idx, :, :, :, 1)), [2, 3]);
        data_active = mean(squeeze(data(light_idx, day_idx, :, :, :, 2)), [2, 3]);
        max_value = max([max_value; data_control(:); data_active(:)]);
        % min_value = min([min_value; data_control(:); data_active(:)]);
    end
end

p_values_table = {};
table_row_idx = 1;

for day_idx = 1:num_days
    for light_idx = 1:num_light_conditions
        data_control = squeeze(data(light_idx, day_idx, :, :, :, 1));
        data_active = squeeze(data(light_idx, day_idx, :, :, :, 2));

        under_diag_values_active = [];
        under_diag_values_control = [];

        % Extract values below the diagonal
        for i = 2:size(data_control, 2)
            for j = 1:i-1
                under_diag_values_active = [data_active(:, i, j), under_diag_values_active];
                under_diag_values_control = [data_control(:, i, j), under_diag_values_control];
            end
        end
        data_control = mean(under_diag_values_control, 1);
        data_active = mean(under_diag_values_active, 1);

        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        hold on;

        % Shape-code and scatter plots with jitter for each subject
        for subj = 1:length(data_control)
            marker_idx = mod(subj - 1, num_marker_styles) + 1;
            
            % Control group
            scatter(1 + randn(size(data_control(subj))) * 0.05, ...
                data_control(subj), 100, ...
                marker_styles{marker_idx}, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', color_set1, ...
                'MarkerFaceAlpha', 0.6);
            
            % Active group
            scatter(2 + randn(size(data_active(subj))) * 0.05, ...
                data_active(subj), 100, ...
                marker_styles{marker_idx}, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', color_set2, ...
                'MarkerFaceAlpha', 0.6);
        end

        % Box plots
        boxplot([data_control(:), data_active(:)], 'Labels', {'Control', 'Active'}, 'Width', 0.7);
        h = findobj(gca, 'Tag', 'Box');
        patch(get(h(1), 'XData'), get(h(1), 'YData'), color_set2, 'FaceAlpha', 0.2);
        patch(get(h(2), 'XData'), get(h(2), 'YData'), color_set1, 'FaceAlpha', 0.2);

        % Perform Mann-Whitney U test
        [p_value, ~] = ranksum(data_control(:), data_active(:));

        % Add significance bracket
        add_significance_bracket(1, 2, max_value * 1.1, p_value);

        % Formatting
        title([light_conditions{light_idx} ', ' day_labels{day_idx}], ...
            'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');
        ylabel('PLV', 'FontSize', fontsize_labels, 'Interpreter', 'tex');
        xlabel('Groups', 'FontSize', fontsize_labels, 'Interpreter', 'tex');
        ylim([min_value, max_value * 1.4]);

        ax = gca;
        ax.FontSize = fontsize_axis;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
        
        image_path_dir = fullfile(path2save);
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        image_path = fullfile(image_path_dir, strcat("Box_plot_", light_conditions{light_idx}, day_labels{day_idx}, ".pdf"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, strcat("Box_plot_", light_conditions{light_idx}, day_labels{day_idx}, ".png"));
        exportgraphics(fig, image_path, 'Resolution', 300);

        p_values_table{table_row_idx, 1} = light_conditions{light_idx};
        p_values_table{table_row_idx, 2} = day_labels{day_idx};
        p_values_table{table_row_idx, 3} = p_value;
        table_row_idx = table_row_idx + 1;
    end
end

csv_path = fullfile(path2save, "subject_level_each_light_day_pair", "statistical_analysis_results.csv");
p_values_table_headers = {'Light Condition', 'Day', 'p-value (control vs active)'};
p_values_table_final = cell2table(p_values_table, 'VariableNames', p_values_table_headers);
writetable(p_values_table_final, csv_path);


%% Active vs control over days in the same plot
clc
close all
% Initialize table for statistical results
table_row_idx = 1;
p_values_table = {};
path2save = fullfile("results", "PLV");
num_days = 3;
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
num_marker_styles = length(marker_styles);

control_data = plv_values_control;
active_data = plv_values_active;
y_range = [0.4, 1.4];
global_max = 1;


for light_idx = 1:3  % Loop over light conditions
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    
    % Collect data for all days and groups
    data_all_days = [];
    group_labels_all = [];
    day_labels_all = [];
    
    for day_idx = 1:num_days  % Loop over days
        control_data_flat = squeeze(control_data(light_idx, day_idx, :, :, :));
        active_data_flat = squeeze(active_data(light_idx, day_idx, :, :, :));

        under_diag_values_active = [];
        under_diag_values_control = [];
        % Extract values below the diagonal
        for i = 2:size(control_data_flat, 2)
            for j = 1:i-1
                under_diag_values_active = [active_data_flat(:, i, j), under_diag_values_active];
                under_diag_values_control = [control_data_flat(:, i, j), under_diag_values_control];
            end
        end

        control_data_flat = under_diag_values_control(:);
        active_data_flat = under_diag_values_active(:);

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
        [p_value_inter_group, h_stat_inter_group] = ranksum(control_data_flat, active_data_flat);

        % Record inter-group p-value
        p_values_table{table_row_idx, 1} = light_conditions{light_idx};
        p_values_table{table_row_idx, 2} = day_labels{day_idx};
        p_values_table{table_row_idx, 3} = 'Control vs Active';
        p_values_table{table_row_idx, 4} = 'Ranksum Test';
        p_values_table{table_row_idx, 5} = h_stat_inter_group; % Test statistic (Ranksum)
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

    hold on

    h = findobj(gca, 'Tag', 'Box');

    for day_idx = 1:num_days  % Loop over days
        control_data_flat = squeeze(control_data(light_idx, day_idx, :, :, :));
        active_data_flat = squeeze(active_data(light_idx, day_idx, :, :, :));

        under_diag_values_active = [];
        under_diag_values_control = [];
        % Extract values below the diagonal
        for i = 2:size(control_data_flat, 2)
            for j = 1:i-1
                under_diag_values_active = [active_data_flat(:, i, j), under_diag_values_active];
                under_diag_values_control = [control_data_flat(:, i, j), under_diag_values_control];
            end
        end
        control_data_flat = under_diag_values_control;
        active_data_flat = under_diag_values_active;

        set(h((day_idx-1)*2+1), 'Color', color_set_active);
        set(h((day_idx-1)*2+1), 'LineWidth', 1.5);

        set(h((day_idx-1)*2+2), 'Color', color_set_control);
        set(h((day_idx-1)*2+2), 'LineWidth', 1.5);
        % patch(get(h((day_idx-1)*2+1), 'XData'), get(h((day_idx-1)*2+1), 'YData'), color_set_active, 'FaceAlpha', 0.5);
        % patch(get(h((day_idx-1)*2+2), 'XData'), get(h((day_idx-1)*2+2), 'YData'), color_set_control, 'FaceAlpha', 0.5);
        
        for subj = 1:size(control_data, 3)
            marker_idx = mod(subj-1, num_marker_styles) + 1;
            scatter((day_idx-1)*2+ones(size(control_data_flat(subj, :))), control_data_flat(subj, :), 50, marker_styles{marker_idx}, ...
            'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.1);
        end
        for subj = 1:size(active_data, 3)
            marker_idx = mod(subj-1, num_marker_styles) + 1;
            scatter((day_idx-1)*2+1+ones(size(active_data_flat(subj, :))), active_data_flat(subj, :), 50, marker_styles{marker_idx}, ...
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
            day_data = squeeze(current_data(light_idx, day_idx, :, :, :));
            under_diag_values = [];
            % Extract values below the diagonal
            for i = 2:size(day_data, 2)
                for j = 1:i-1
                    under_diag_values = [day_data(:, i, j), under_diag_values];
                end
            end
            day_data = under_diag_values;
            group_data = [group_data; day_data(:)];
            day_labels_in_group = [day_labels_in_group; repmat({day_labels{day_idx}}, length(day_data(:)), 1)];
        end
        
        % Perform Kruskal-Wallis test
        [p_kw, tbl_kw, stats_kw] = kruskalwallis(group_data, day_labels_in_group, 'off');
        chi2_stat = tbl_kw{2, 5}; % Extract the chi-squared statistic
        
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
            c = multcompare(stats, 'Display', 'off');
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
                    p_values_table{table_row_idx, 6} = 'N/A'; % Sample size not applicable for pairwise comparisons
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

% Save p-values as CSV
csv_path = fullfile(path2save, "light_comparisons", "statistical_analysis_results.csv");
p_values_table_headers = {'Light Condition', 'Day/Group', 'Comparison', 'Test Type', 'Test Statistic', 'Sample Size', 'p-value'};
p_values_table_final = cell2table(p_values_table, 'VariableNames', p_values_table_headers);
writetable(p_values_table_final, csv_path);


%% PLV - Active vs control over days seperated by brain region
clc;
close all;

% Define brain regions based on channel locations
brain_regions = struct();
brain_regions.Frontal = [1, 2];     % FP1, FP2
brain_regions.Central = [3, 4];     % C3, C4  
brain_regions.Temporal = [5, 6];    % TP7, TP8
brain_regions.Occipital = [7, 8];   % O1, O2

region_names = fieldnames(brain_regions);
num_regions = length(region_names);

% Create path for saving results
path2save = fullfile("results", "PLV_brain_region_analysis");
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

% Extract and organize PLV data by brain regions
% PLV data structure: light_condition * days * subjects * channels * channels
num_region_pairs = num_regions * (num_regions - 1) / 2 + num_regions; % Including within-region and between-region pairs

% Initialize storage for regional PLV data - keep all channel pairs instead of averaging
% Structure: light_condition * days * subjects * region_pairs * individual_channel_pairs
plv_data_regions = struct();

% Calculate maximum number of channel pairs per region pair
max_pairs_per_region = 0;
for r1 = 1:num_regions
    region_channels = brain_regions.(region_names{r1});
    if length(region_channels) > 1
        num_within_pairs = length(region_channels) * (length(region_channels) - 1) / 2;
        max_pairs_per_region = max(max_pairs_per_region, num_within_pairs);
    end
end

for r1 = 1:num_regions
    for r2 = r1+1:num_regions
        region1_channels = brain_regions.(region_names{r1});
        region2_channels = brain_regions.(region_names{r2});
        num_between_pairs = length(region1_channels) * length(region2_channels);
        max_pairs_per_region = max(max_pairs_per_region, num_between_pairs);
    end
end

% Initialize storage with maximum possible channel pairs
plv_data_regions.control = NaN(numel(stim_types), numel(recording_days), numel(control_subjects), num_region_pairs, max_pairs_per_region);
plv_data_regions.active = NaN(numel(stim_types), numel(recording_days), numel(active_subjects), num_region_pairs, max_pairs_per_region);

% Create region pair labels
region_pair_labels = {};
region_pair_idx = 1;

% Within-region pairs
for r1 = 1:num_regions
    region_pair_labels{region_pair_idx} = sprintf('%s-Within', region_names{r1});
    region_pair_idx = region_pair_idx + 1;
end

% Between-region pairs
for r1 = 1:num_regions
    for r2 = r1+1:num_regions
        region_pair_labels{region_pair_idx} = sprintf('%s-%s', region_names{r1}, region_names{r2});
        region_pair_idx = region_pair_idx + 1;
    end
end

% Extract PLV values for each region pair - keep all individual channel pairs
for light_idx = 1:numel(stim_types)
    for day_idx = 1:numel(recording_days)
        % Process control group
        for subj = 1:numel(control_subjects)
            plv_matrix = squeeze(plv_values_control(light_idx, day_idx, subj, :, :));
            region_pair_idx = 1;
            
            % Within-region PLV (keep all channel pairs within each region)
            for r1 = 1:num_regions
                region_channels = brain_regions.(region_names{r1});
                if length(region_channels) > 1
                    pair_counter = 1;
                    for ch1_idx = 1:length(region_channels)
                        for ch2_idx = ch1_idx+1:length(region_channels)
                            ch1 = region_channels(ch1_idx);
                            ch2 = region_channels(ch2_idx);
                            plv_data_regions.control(light_idx, day_idx, subj, region_pair_idx, pair_counter) = plv_matrix(ch1, ch2);
                            pair_counter = pair_counter + 1;
                        end
                    end
                else
                    plv_data_regions.control(light_idx, day_idx, subj, region_pair_idx, 1) = 1; % Single channel case
                end
                region_pair_idx = region_pair_idx + 1;
            end
            
            % Between-region PLV (keep all channel pairs between regions)
            for r1 = 1:num_regions
                for r2 = r1+1:num_regions
                    region1_channels = brain_regions.(region_names{r1});
                    region2_channels = brain_regions.(region_names{r2});
                    
                    pair_counter = 1;
                    for ch1 = region1_channels
                        for ch2 = region2_channels
                            plv_data_regions.control(light_idx, day_idx, subj, region_pair_idx, pair_counter) = plv_matrix(ch1, ch2);
                            pair_counter = pair_counter + 1;
                        end
                    end
                    region_pair_idx = region_pair_idx + 1;
                end
            end
        end
        
        % Process active group
        for subj = 1:numel(active_subjects)
            plv_matrix = squeeze(plv_values_active(light_idx, day_idx, subj, :, :));
            region_pair_idx = 1;
            
            % Within-region PLV
            for r1 = 1:num_regions
                region_channels = brain_regions.(region_names{r1});
                if length(region_channels) > 1
                    pair_counter = 1;
                    for ch1_idx = 1:length(region_channels)
                        for ch2_idx = ch1_idx+1:length(region_channels)
                            ch1 = region_channels(ch1_idx);
                            ch2 = region_channels(ch2_idx);
                            plv_data_regions.active(light_idx, day_idx, subj, region_pair_idx, pair_counter) = plv_matrix(ch1, ch2);
                            pair_counter = pair_counter + 1;
                        end
                    end
                else
                    plv_data_regions.active(light_idx, day_idx, subj, region_pair_idx, 1) = 1; % Single channel case
                end
                region_pair_idx = region_pair_idx + 1;
            end
            
            % Between-region PLV
            for r1 = 1:num_regions
                for r2 = r1+1:num_regions
                    region1_channels = brain_regions.(region_names{r1});
                    region2_channels = brain_regions.(region_names{r2});
                    
                    pair_counter = 1;
                    for ch1 = region1_channels
                        for ch2 = region2_channels
                            plv_data_regions.active(light_idx, day_idx, subj, region_pair_idx, pair_counter) = plv_matrix(ch1, ch2);
                            pair_counter = pair_counter + 1;
                        end
                    end
                    region_pair_idx = region_pair_idx + 1;
                end
            end
        end
    end
end

% Normalize by no-light condition
plv_data_regions.control = plv_data_regions.control ./ plv_data_regions.control(1, :, :, :, :);
plv_data_regions.active = plv_data_regions.active ./ plv_data_regions.active(1, :, :, :, :);

% Find global min and max for consistent scaling
all_data_regions = cat(6, plv_data_regions.control, plv_data_regions.active);
all_data_regions = all_data_regions(~isnan(all_data_regions));
global_min = min(all_data_regions(:));
global_max = max(all_data_regions(:));
y_range = [global_min*0.95, global_max*1.15];

% Define marker styles for subjects
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
num_marker_styles = length(marker_styles);

% % Regional comparison plots for each light condition and day
% p_values_table = {};
% table_row_idx = 1;
% 
% for light_idx = 2:numel(stim_types)  % Skip no_light condition
%     for day_idx = 1:numel(recording_days)
%         fig = figure;
%         set(gcf, 'Position', get(0, 'Screensize'));
% 
%         % Create subplots for each region pair (3x3 grid)
%         num_plots = min(num_region_pairs, 9);  % Limit to 9 subplots
%         subplot_rows = ceil(sqrt(num_plots));
%         subplot_cols = ceil(num_plots / subplot_rows);
% 
%         for pair_idx = 1:num_plots
%             subplot(subplot_rows, subplot_cols, pair_idx);
%             hold on;
% 
%             % Extract data for current condition and region pair - flatten all channel pairs per subject
%             control_region_data = squeeze(plv_data_regions.control(light_idx, day_idx, :, pair_idx, :));
%             active_region_data = squeeze(plv_data_regions.active(light_idx, day_idx, :, pair_idx, :));
% 
%             % Remove NaN values and reshape for plotting
%             control_flat = [];
%             control_subject_labels = [];
%             active_flat = [];
%             active_subject_labels = [];
% 
%             % Control group - keep track of which subject each data point belongs to
%             for subj = 1:size(control_region_data, 1)
%                 subj_data = control_region_data(subj, :);
%                 subj_data = subj_data(~isnan(subj_data));
%                 control_flat = [control_flat, subj_data];
%                 control_subject_labels = [control_subject_labels, subj * ones(1, length(subj_data))];
%             end
% 
%             % Active group - keep track of which subject each data point belongs to
%             for subj = 1:size(active_region_data, 1)
%                 subj_data = active_region_data(subj, :);
%                 subj_data = subj_data(~isnan(subj_data));
%                 active_flat = [active_flat, subj_data];
%                 active_subject_labels = [active_subject_labels, subj * ones(1, length(subj_data))];
%             end
% 
%             % Create box plots
%             if ~isempty(control_flat) && ~isempty(active_flat)
%                 boxplot([control_flat(:), active_flat(:)], ...
%                     [ones(size(control_flat(:))); 2*ones(size(active_flat(:)))], ...
%                     'Labels', {'Control', 'Active'}, 'Width', 0.7);
% 
%                 h = findobj(gca, 'Tag', 'Box');
%                 if length(h) >= 2
%                     patch(get(h(1), 'XData'), get(h(1), 'YData'), color_set2, 'FaceAlpha', 0.2);
%                     patch(get(h(2), 'XData'), get(h(2), 'YData'), color_set1, 'FaceAlpha', 0.2);
%                 end
% 
%                 % Plot individual channel data points with subject-specific markers
%                 for i = 1:length(control_flat)
%                     subj = control_subject_labels(i);
%                     marker_idx = mod(subj-1, num_marker_styles) + 1;
% 
%                     % Control group
%                     scatter(1 + randn(1)*0.05, control_flat(i), 40, ...
%                         marker_styles{marker_idx}, ...
%                         'MarkerEdgeColor', 'none', ...
%                         'MarkerFaceColor', color_set_control, ...
%                         'MarkerFaceAlpha', 0.7);
%                 end
% 
%                 for i = 1:length(active_flat)
%                     subj = active_subject_labels(i);
%                     marker_idx = mod(subj-1, num_marker_styles) + 1;
% 
%                     % Active group
%                     scatter(2 + randn(1)*0.05, active_flat(i), 40, ...
%                         marker_styles{marker_idx}, ...
%                         'MarkerEdgeColor', 'none', ...
%                         'MarkerFaceColor', color_set_active, ...
%                         'MarkerFaceAlpha', 0.7);
%                 end
% 
%                 % Statistical test
%                 [p_value, ~] = ranksum(control_flat, active_flat);
%                 add_significance_bracket(1, 2, global_max*1.05, p_value);
% 
%                 % Store p-values in table
%                 p_values_table{table_row_idx, 1} = light_conditions{light_idx};
%                 p_values_table{table_row_idx, 2} = day_labels{day_idx};
%                 p_values_table{table_row_idx, 3} = region_pair_labels{pair_idx};
%                 p_values_table{table_row_idx, 4} = p_value;
%                 table_row_idx = table_row_idx + 1;
%             end
% 
%             % Formatting
%             title(sprintf('%s', region_pair_labels{pair_idx}), ...
%                 'FontSize', fontsize_title-4, 'FontWeight', 'bold', 'Interpreter', 'tex');
%             ylabel('Normalized PLV', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
%             xlabel('Groups', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
%             ylim(y_range);
% 
%             ax = gca;
%             ax.FontSize = fontsize_axis-2;
%             ax.LineWidth = 1.5;
%             ax.Box = 'on';
%         end
% 
%         sgtitle(sprintf('PLV Brain Regional Analysis: %s, %s', ...
%             light_conditions{light_idx}, day_labels{day_idx}), ...
%             'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');
% 
%         % Save figure
%         image_path_dir = fullfile(path2save, "regional_comparisons");
%         if ~exist(image_path_dir, 'dir')
%             mkdir(image_path_dir);
%         end
%         image_path = fullfile(image_path_dir, sprintf("PLV_Regional_Analysis_%s_%s.pdf", ...
%             light_conditions{light_idx}, day_labels{day_idx}));
%         exportgraphics(fig, image_path, 'Resolution', 300);
%         image_path = fullfile(image_path_dir, sprintf("PLV_Regional_Analysis_%s_%s.png", ...
%             light_conditions{light_idx}, day_labels{day_idx}));
%         exportgraphics(fig, image_path, 'Resolution', 300);
%     end
% end

% Regional comparison plots for each light condition and day
p_values_table = {};
table_row_idx = 1;

for light_idx = 2:numel(stim_types)  % Skip no_light condition
    for day_idx = 1:numel(recording_days)
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        
        num_plots = min(num_region_pairs, 12);  % Limit to 9 subplots
        hold on;
        x_positions = [];
        x_ticks = {};
        for pair_idx = 5:num_plots
            pair_idx_adjusted = pair_idx-4;
            
            % Extract data for current condition and region pair - flatten all channel pairs per subject
            control_region_data = squeeze(plv_data_regions.control(light_idx, day_idx, :, pair_idx, :));
            active_region_data = squeeze(plv_data_regions.active(light_idx, day_idx, :, pair_idx, :));
            
            % Remove NaN values and reshape for plotting
            control_flat = [];
            control_subject_labels = [];
            active_flat = [];
            active_subject_labels = [];
            
            % Control group - keep track of which subject each data point belongs to
            for subj = 1:size(control_region_data, 1)
                subj_data = control_region_data(subj, :);
                subj_data = subj_data(~isnan(subj_data));
                control_flat = [control_flat, subj_data];
                control_subject_labels = [control_subject_labels, subj * ones(1, length(subj_data))];
            end
            
            % Active group - keep track of which subject each data point belongs to
            for subj = 1:size(active_region_data, 1)
                subj_data = active_region_data(subj, :);
                subj_data = subj_data(~isnan(subj_data));
                active_flat = [active_flat, subj_data];
                active_subject_labels = [active_subject_labels, subj * ones(1, length(subj_data))];
            end
            
            % Create box plots
            if ~isempty(control_flat) && ~isempty(active_flat)
                base_x = (pair_idx_adjusted - 1) * 3;  % Space between region pairs
                control_x = base_x + 1;
                active_x = base_x + 2;
                x_positions(end+1) = base_x+1.5;
                x_ticks{end+1} = region_pair_labels{pair_idx};

                boxplot([control_flat(:), active_flat(:)], ...
                    [ones(size(control_flat(:))); 2*ones(size(active_flat(:)))], ...
                    'Labels', {'Control', 'Active'}, 'Width', 0.7, 'Positions', [control_x, active_x]);
                
                h = findobj(gca, 'Tag', 'Box');
                if length(h) >= 2
                    patch(get(h(1), 'XData'), get(h(1), 'YData'), color_set_active, 'FaceAlpha', 0.2);
                    patch(get(h(2), 'XData'), get(h(2), 'YData'), color_set_control, 'FaceAlpha', 0.2);
                end
                
                % Plot individual channel data points with subject-specific markers
                for i = 1:length(control_flat)
                    subj = control_subject_labels(i);
                    marker_idx = mod(subj-1, num_marker_styles) + 1;
                    
                    % Control group
                    scatter(control_x + randn(1)*0.05, control_flat(i), 40, ...
                        marker_styles{marker_idx}, ...
                        'MarkerEdgeColor', 'none', ...
                        'MarkerFaceColor', color_set_control, ...
                        'MarkerFaceAlpha', 0.7);
                end
                
                for i = 1:length(active_flat)
                    subj = active_subject_labels(i);
                    marker_idx = mod(subj-1, num_marker_styles) + 1;
                    
                    % Active group
                    scatter(active_x + randn(1)*0.05, active_flat(i), 40, ...
                        marker_styles{marker_idx}, ...
                        'MarkerEdgeColor', 'none', ...
                        'MarkerFaceColor', color_set_active, ...
                        'MarkerFaceAlpha', 0.7);
                end
                
                % Statistical test
                [p_value, ~] = ranksum(control_flat, active_flat);
                add_significance_bracket(control_x, active_x, global_max*1.05, p_value);
                
                % Store p-values in table
                p_values_table{table_row_idx, 1} = light_conditions{light_idx};
                p_values_table{table_row_idx, 2} = day_labels{day_idx};
                p_values_table{table_row_idx, 3} = region_pair_labels{pair_idx};
                p_values_table{table_row_idx, 4} = p_value;
                table_row_idx = table_row_idx + 1;
            end
            
            % Formatting
            % title(sprintf('%s', region_pair_labels{pair_idx}), ...
            %     'FontSize', fontsize_title-4, 'FontWeight', 'bold', 'Interpreter', 'tex');
            ylabel('Normalized PLV', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
            xlabel('Groups', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
            ylim(y_range);
            
            ax = gca;
            ax.FontSize = fontsize_axis-2;
            ax.LineWidth = 1.5;
            ax.Box = 'on';
        end

         set(gca, 'XTick', x_positions, 'XTickLabel', x_ticks, ...
            'XTickLabelRotation', 45);
        
        sgtitle(sprintf('PLV Brain Regional Analysis: %s, %s', ...
            light_conditions{light_idx}, day_labels{day_idx}), ...
            'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');
        
        % Save figure
        image_path_dir = fullfile(path2save, "regional_comparisons");
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        image_path = fullfile(image_path_dir, sprintf("PLV_Regional_Analysis_%s_%s.pdf", ...
            light_conditions{light_idx}, day_labels{day_idx}));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, sprintf("PLV_Regional_Analysis_%s_%s.png", ...
            light_conditions{light_idx}, day_labels{day_idx}));
        exportgraphics(fig, image_path, 'Resolution', 300);
    end
end

%% Time course analysis by brain region pairs for PLV
for light_idx = 2:numel(stim_types)  % Skip no_light condition
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    
    % Create subplots for region pairs (3x3 grid)
    num_plots = min(num_region_pairs, 9);
    subplot_rows = ceil(sqrt(num_plots));
    subplot_cols = ceil(num_plots / subplot_rows);
    
    for pair_idx = 1:num_plots
        subplot(subplot_rows, subplot_cols, pair_idx);
        hold on;
        
        % Calculate means and standard errors across days using all channel pairs
        control_means = [];
        active_means = [];
        control_se = [];
        active_se = [];
        
        for day_idx = 1:numel(recording_days)
            % Get all channel pairs for this day
            control_day_data = squeeze(plv_data_regions.control(light_idx, day_idx, :, pair_idx, :));
            active_day_data = squeeze(plv_data_regions.active(light_idx, day_idx, :, pair_idx, :));
            
            % Flatten and remove NaN
            control_flat = control_day_data(~isnan(control_day_data));
            active_flat = active_day_data(~isnan(active_day_data));
            
            control_means = [control_means, mean(control_flat)];
            active_means = [active_means, mean(active_flat)];
            control_se = [control_se, std(control_flat) / sqrt(length(control_flat))];
            active_se = [active_se, std(active_flat) / sqrt(length(active_flat))];
        end
        
        % Plot individual subject trajectories (faded) - average across channel pairs per subject
        for subj = 1:size(plv_data_regions.control, 3)
            control_trajectory = [];
            for day_idx = 1:numel(recording_days)
                subj_data = squeeze(plv_data_regions.control(light_idx, day_idx, subj, pair_idx, :));
                subj_data = subj_data(~isnan(subj_data));
                control_trajectory = [control_trajectory, mean(subj_data)];
            end
            plot(1:numel(recording_days), control_trajectory, '-', ...
                'Color', [color_set_control, 0.3], 'LineWidth', 1);
        end
        
        for subj = 1:size(plv_data_regions.active, 3)
            active_trajectory = [];
            for day_idx = 1:numel(recording_days)
                subj_data = squeeze(plv_data_regions.active(light_idx, day_idx, subj, pair_idx, :));
                subj_data = subj_data(~isnan(subj_data));
                active_trajectory = [active_trajectory, mean(subj_data)];
            end
            plot(1:numel(recording_days), active_trajectory, '-', ...
                'Color', [color_set_active, 0.3], 'LineWidth', 1);
        end
        
        % Plot mean trajectories with error bars
        errorbar(1:numel(recording_days), control_means, control_se, '-o', ...
            'Color', color_set_control, 'LineWidth', 2, 'MarkerSize', 6, ...
            'MarkerFaceColor', color_set_control, 'MarkerEdgeColor', 'none');
        
        errorbar(1:numel(recording_days), active_means, active_se, '-s', ...
            'Color', color_set_active, 'LineWidth', 2, 'MarkerSize', 6, ...
            'MarkerFaceColor', color_set_active, 'MarkerEdgeColor', 'none');
        
        % Formatting
        title(sprintf('%s', region_pair_labels{pair_idx}), ...
            'FontSize', fontsize_title-4, 'FontWeight', 'bold', 'Interpreter', 'tex');
        xlabel('Time Point', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
        ylabel('Normalized PLV', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
        
        set(gca, 'XTick', 1:numel(recording_days), 'XTickLabel', day_labels);
        ylim([global_min*0.9, global_max*1.1]);
        grid on;
        
        if pair_idx == 1
            legend('Control Individual', 'Active Individual', 'Control Mean', 'Active Mean', ...
                'Location', 'best', 'FontSize', fontsize_axis-2, 'Interpreter', 'tex');
        end
        
        ax = gca;
        ax.FontSize = fontsize_axis-2;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
    end
    
    sgtitle(sprintf('PLV Temporal Evolution by Brain Region Pairs: %s', light_conditions{light_idx}), ...
        'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Save figure
    image_path_dir = fullfile(path2save, "temporal_evolution");
    if ~exist(image_path_dir, 'dir')
        mkdir(image_path_dir);
    end
    image_path = fullfile(image_path_dir, sprintf("PLV_Temporal_Evolution_%s.pdf", light_conditions{light_idx}));
    exportgraphics(fig, image_path, 'Resolution', 300);
    image_path = fullfile(image_path_dir, sprintf("PLV_Temporal_Evolution_%s.png", light_conditions{light_idx}));
    exportgraphics(fig, image_path, 'Resolution', 300);
end

% Summary heatmap showing PLV effect sizes by region pairs
effect_sizes = zeros(numel(stim_types)-1, numel(recording_days), num_region_pairs);
p_values_matrix = zeros(numel(stim_types)-1, numel(recording_days), num_region_pairs);

for light_idx = 2:numel(stim_types)
    for day_idx = 1:numel(recording_days)
        for pair_idx = 1:num_region_pairs
            % Get all channel pairs for this condition
            control_data = squeeze(plv_data_regions.control(light_idx, day_idx, :, pair_idx, :));
            active_data = squeeze(plv_data_regions.active(light_idx, day_idx, :, pair_idx, :));
            
            % Flatten and remove NaN
            control_flat = control_data(~isnan(control_data));
            active_flat = active_data(~isnan(active_data));
            
            if ~isempty(control_flat) && ~isempty(active_flat)
                % Calculate Cohen's d (effect size)
                pooled_std = sqrt(((length(control_flat)-1)*var(control_flat) + ...
                                  (length(active_flat)-1)*var(active_flat)) / ...
                                  (length(control_flat) + length(active_flat) - 2));
                cohens_d = (mean(active_flat) - mean(control_flat)) / pooled_std;
                
                effect_sizes(light_idx-1, day_idx, pair_idx) = cohens_d;
                
                % Get p-value
                [p_val, ~] = ranksum(control_flat, active_flat);
                p_values_matrix(light_idx-1, day_idx, pair_idx) = p_val;
            end
        end
    end
end

% Create heatmap figure
fig = figure;
set(gcf, 'Position', get(0, 'Screensize'));

% Plot each region pair as a separate subplot
num_plots = min(num_region_pairs, 9);
subplot_rows = ceil(sqrt(num_plots));
subplot_cols = ceil(num_plots / subplot_rows);

for pair_idx = 1:num_plots
    subplot(subplot_rows, subplot_cols, pair_idx);
    
    pair_effects = squeeze(effect_sizes(:, :, pair_idx));
    pair_pvals = squeeze(p_values_matrix(:, :, pair_idx));
    
    imagesc(pair_effects);
    colormap('abyss');
    colorbar;
    clim([-1, 1]);
    
    % Mark significant effects
    hold on;
    [row, col] = find(pair_pvals < 0.05);
    if ~isempty(row)
        plot(col, row, 'k*', 'MarkerSize', 8, 'LineWidth', 2);
    end
    
    title(sprintf('%s', region_pair_labels{pair_idx}), ...
        'FontSize', fontsize_title-4, 'FontWeight', 'bold', 'Interpreter', 'tex');
    xlabel('Day', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
    ylabel('Light Condition', 'FontSize', fontsize_labels-2, 'Interpreter', 'tex');
    
    set(gca, 'XTick', 1:numel(recording_days), 'XTickLabel', day_labels);
    set(gca, 'YTick', 1:numel(stim_types)-1, 'YTickLabel', light_conditions(2:end));
    
    ax = gca;
    ax.FontSize = fontsize_axis-2;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
end

sgtitle('PLV Effect Sizes (Cohen''s d) by Brain Region Pairs (* p < 0.05)', ...
    'FontSize', fontsize_title, 'FontWeight', 'bold', 'Interpreter', 'tex');

% Save heatmap figure
image_path_dir = fullfile(path2save, "effect_size_heatmaps");
if ~exist(image_path_dir, 'dir')
    mkdir(image_path_dir);
end
image_path = fullfile(image_path_dir, "PLV_Regional_Effect_Sizes.pdf");
exportgraphics(fig, image_path, 'Resolution', 300);
image_path = fullfile(image_path_dir, "PLV_Regional_Effect_Sizes.png");
exportgraphics(fig, image_path, 'Resolution', 300);

% Save PLV regional statistical results
csv_path = fullfile(path2save, "plv_regional_statistical_analysis.csv");
p_values_table_headers = {'Light Condition', 'Day', 'Brain Region Pair', 'p-value (control vs active)'};
p_values_table_final = cell2table(p_values_table, 'VariableNames', p_values_table_headers);
writetable(p_values_table_final, csv_path);

% Save effect size data
effect_size_results = {};
result_idx = 1;
for light_idx = 2:numel(stim_types)
    for day_idx = 1:numel(recording_days)
        for pair_idx = 1:num_region_pairs
            effect_size_results{result_idx, 1} = light_conditions{light_idx};
            effect_size_results{result_idx, 2} = day_labels{day_idx};
            effect_size_results{result_idx, 3} = region_pair_labels{pair_idx};
            effect_size_results{result_idx, 4} = effect_sizes(light_idx-1, day_idx, pair_idx);
            effect_size_results{result_idx, 5} = p_values_matrix(light_idx-1, day_idx, pair_idx);
            result_idx = result_idx + 1;
        end
    end
end

effect_size_headers = {'Light Condition', 'Day', 'Brain Region Pair', 'Effect Size (Cohens d)', 'p-value'};
effect_size_table = cell2table(effect_size_results, 'VariableNames', effect_size_headers);
csv_path_effects = fullfile(path2save, "plv_regional_effect_sizes.csv");
writetable(effect_size_table, csv_path_effects);

fprintf('PLV brain region analysis completed. Results saved to: %s\n', path2save);

% PLV Summary statistics by region pairs
fprintf('\n=== PLV BRAIN REGION ANALYSIS SUMMARY ===\n');
for pair_idx = 1:num_region_pairs
    fprintf('\n%s:\n', region_pair_labels{pair_idx});
    
    % Find most significant effects
    pair_pvals = squeeze(p_values_matrix(:, :, pair_idx));
    [min_p, min_idx] = min(pair_pvals(:));
    [light_idx, day_idx] = ind2sub(size(pair_pvals), min_idx);
    
    if min_p < 0.05
        fprintf('Most significant effect: %s, %s (p = %.4f)\n', ...
            light_conditions{light_idx+1}, day_labels{day_idx}, min_p);
        fprintf('Effect size: %.3f\n', effect_sizes(light_idx, day_idx, pair_idx));
    else
        fprintf('No significant effects found\n');
    end
end