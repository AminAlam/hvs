%% Load Clean files
clc
clear
close all

data_path = "/Users/amin/Documents/Syntropic/HVS/Recordings";
list_of_subjects = {"SU1_24_JEz",...
                      "SU1_24_4V9",...
                      "SU1_24_SOH",...
                      "SU1_24_UzK", ...
                      "SU1_24_jWX",...
                      "SU1_24_n4X",...
                      "SU1_24_RDj",...
                      "SU1_24_gTo",...
                      "SU1_24_Q0R",...
                      "SU1_24_nub",...
                      "SU1_24_u1W",...
                      "SU1_24_jp6"};

recording_days = {"day_1", ...
                "day_5", ...
                "day_19"};

stim_types = {"no_light", ...
            "constant_light", ... 
            "stimulus_light"};

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


fft_y0 = -30;
fft_y1 = -3;

% Define colorblind-friendly colors
color_set1 = [230, 159, 0]/255;  % Orange
color_set2 = [86, 180, 233]/255; % Sky Blue
color_set3 = [0, 158, 115]/255;  % Bluish Green

color_set_active = [77,139,145]/255;
color_set_control = [81, 90, 95]/255;

group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};

ch_list = {'FP1', 'FP2', 'C3', 'C4', 'TP7', 'TP8', 'O1', 'O2'};

micro_volt_constant = 1e-6;

fs = 250;

fontsize_title = 20;
fontsize_labels = 16;
fontsize_axis = 14;

for subject = list_of_subjects
    for recording_day = recording_days
        for stim_type = stim_types
            subject = string(subject(1));
            recording_day = string(recording_day(1));
            stim_type = string(stim_type(1));

            filepath = char(fullfile(data_path, subject));
            filename = char(strcat(recording_day, "_", stim_type, "_ICA", ".set"));

            EEG = pop_loadset('filename',filename,'filepath',filepath);
            data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean = EEG.data;
        end
    end
end




%% BrainNet Viewer
clc
close all
data = cat(6, plv_values_control, plv_values_active); % light_condition * days * subjects * channels * subject_type

for group_indx = 1:2
    fig_dummy = figure;
    fig_main = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot_rows = num_days;
    subplot_cols = num_light_conditions;
    
    tiledlayout(subplot_rows,subplot_cols, 'Padding', 'none', 'TileSpacing', 'compact'); 
    
    for day_idx = 1:num_days
        for light_idx = 1:num_light_conditions
            % Calculate subplot position
            
            % Process data
            data_ = mean(squeeze(data(light_idx, day_idx, :, :, :, group_indx)), 1);
            % Prepare edge file
            edge_data = squeeze(data_);
            edge_file = 'electrodes.edge';
            writematrix(edge_data, edge_file, 'Delimiter', ' ', 'FileType', 'text');
            % Prepare node file (you might need to create this similarly)
            node_file = 'electrodes.node';
            % BrainNet Viewer path
            brainnet_path = '../BrainNetViewer_20191031';
            surface_file = fullfile(brainnet_path, 'Data', 'SurfTemplate', 'BrainMesh_ICBM152.nv');
    
            figure(fig_dummy);
            BrainNet_MapCfg(surface_file, node_file, edge_file, 'BrainNet_conf.mat', 'electrode_network.png');
            
            figure(fig_main);
            nexttile
            imshow('electrode_network.png')
            title(sprintf('Day %d, Light Condition %d', day_idx, light_idx));
        end
    end

    image_path_dir = fullfile(path2save, "PLV_on_brain");
    if ~exist(image_path_dir, 'dir')
        mkdir(image_path_dir);
    end
    image_path = fullfile(image_path_dir, strcat(group_labels{group_indx}, ".png"));
    exportgraphics(fig_main, image_path, 'Resolution', 300);
    image_path = fullfile(image_path_dir, strcat(group_labels{group_indx}, ".pdf"));
    exportgraphics(fig_main, image_path, 'Resolution', 300);
end

% Adjust figure properties
sgtitle('Brain Network Visualization');
%%
edge_data = [
    0, 0.8, 0.5, 0.2, 0.1, 0, 0, 0;
    0.8, 0, 0.7, 0.3, 0.1, 0, 0, 0;
    0.5, 0.7, 0, 0.9, 0.6, 0, 0, 0;
    0.2, 0.3, 0.9, 0, 0.8, 0, 0, 0;
    0.1, 0.1, 0.6, 0.8, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0
];
edge_file = 'electrodes.edge';
writematrix(edge_data, edge_file, 'Delimiter', ' ', 'FileType', 'text');
node_file = 'electrodes.node';
% Step 4: Visualize with BrainNet Viewer
brainnet_path = '../BrainNetViewer_20191031'; % Replace with your BrainNet Viewer path
surface_file = fullfile(brainnet_path, 'Data', 'SurfTemplate', 'BrainMesh_ICBM152.nv');

BrainNet_MapCfg(surface_file, node_file, edge_file, 'electrode_network.png', 'Viewtype', "3");
