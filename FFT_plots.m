%% FFTs
clc
close all

path2save = fullfile("results", "FFTs");

for subject = list_of_subjects
    for recording_day = recording_days
        for stim_type = stim_types
            subject = string(subject(1));
            recording_day = string(recording_day(1));
            stim_type = string(stim_type(1));

            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            recording_id = char(strcat(subject, "_", recording_day, "_", stim_type));

            fig = figure;
            set(gcf, 'Position', get(0, 'Screensize'));
            plotFFT(data', fs, recording_id)
            image_path_dir = fullfile(path2save, subject);
            image_path = fullfile(image_path_dir, strcat(recording_id, ".pdf"));
            image_path = fullfile(image_path_dir, strcat(recording_id, ".png"));
            ylim([-40, 15])
            if ~exist(image_path_dir, 'dir')
                mkdir(image_path_dir);
            end
            exportgraphics(fig, image_path, 'Resolution', 300);
            close all;
        end
    end
end

%% Group AVG FFT Plots
clc
close all

path2save = fullfile("results", "FFTs");

fs = 250;

group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};

% active
day_counter = 0;
for recording_day = recording_days
    recording_day = string(recording_day(1));
    data_active = [];
    day_counter = day_counter+1;
    stim_couter = 0;
    for stim_type = stim_types
        stim_couter = stim_couter+1;
        stim_type = string(stim_type(1));
        for subject = active_subjects
            subject = string(subject(1));
            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            data_active = [data_active, data];
        end
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        plotFFT(data_active', fs, strcat("Avg FFT active", " | ", light_conditions{stim_couter}, " | ", day_labels{day_counter}))
        image_path_dir = fullfile(path2save, strcat("AVG_results", "_", stim_type, "_", recording_day));
        image_path = fullfile(image_path_dir, strcat("Avg FFT active", ".pdf"));
        image_path = fullfile(image_path_dir, strcat("Avg FFT active", ".png"));
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end

        ylim([-70, 15])
        exportgraphics(fig, image_path, 'Resolution', 300);
    
    end
end

% control
day_counter = 0;
for recording_day = recording_days
    recording_day = string(recording_day(1));
    data_active = [];
    day_counter = day_counter+1;
    stim_couter = 0;
    for stim_type = stim_types
        stim_couter = stim_couter+1;
        stim_type = string(stim_type(1));
        for subject = control_subjects
            subject = string(subject(1));
            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            data_active = [data_active, data];
        end
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        plotFFT(data_active', fs, strcat("Avg FFT control", " | ", light_conditions{stim_couter}, " | ", day_labels{day_counter}))
        image_path_dir = fullfile(path2save, strcat("AVG_results", "_", stim_type, "_", recording_day));
        image_path = fullfile(image_path_dir, strcat("Avg FFT control", ".pdf"));
        image_path = fullfile(image_path_dir, strcat("Avg FFT control", ".png"));
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end

        ylim([fft_y0, fft_y1])
        exportgraphics(fig, image_path, 'Resolution', 300);
    
    end
end

%% FFT active on top of control
clc
close all
path2save = fullfile("results", "FFTs");
fft_y0 = -30;
fft_y1 = 0;

day_counter = 0;
for recording_day = recording_days
    recording_day = string(recording_day(1));
    day_counter = day_counter+1;
    stim_couter = 0;
    for stim_type = stim_types
        data_active = [];
        data_control = [];
        stim_couter = stim_couter+1;
        stim_type = string(stim_type(1));
        for subject = active_subjects
            subject = string(subject(1));
            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            notch_freq = 50;
            quality_factor = 30;
            notchFilter = designfilt('bandstopiir', ...
                                     'FilterOrder', 2, ...
                                     'HalfPowerFrequency1', notch_freq - notch_freq/quality_factor, ...
                                     'HalfPowerFrequency2', notch_freq + notch_freq/quality_factor, ...
                                     'SampleRate', fs);
            
            [b, a] = tf(notchFilter);
            for ch = 1:size(data, 1)  % Loop over channels (rows)
                data(ch, :) = filtfilt(b, a, data(ch, :)); % Zero-phase filtering
            end
            data_active = [data_active, data];
        end

        for subject = control_subjects
            subject = string(subject(1));
            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            notch_freq = 50;
            quality_factor = 30;
            notchFilter = designfilt('bandstopiir', ...
                                     'FilterOrder', 2, ...
                                     'HalfPowerFrequency1', notch_freq - notch_freq/quality_factor, ...
                                     'HalfPowerFrequency2', notch_freq + notch_freq/quality_factor, ...
                                     'SampleRate', fs);
            
            [b, a] = tf(notchFilter);
            for ch = 1:size(data, 1)  % Loop over channels (rows)
                data(ch, :) = filtfilt(b, a, data(ch, :)); % Zero-phase filtering
            end
            data_control = [data_control, data];
        end
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        plotFFT_plain(data_active', fs, color_set_active)
        plotFFT_plain(data_control', fs, color_set_control)

        % patch([48, 52, 52, 48],[fft_y0, fft_y0, fft_y1, fft_y1] , [0.2, 0.2, 0.2]);

        title(strcat("Active vs Control FFT | ",light_conditions{stim_couter}, " | ", day_labels{day_counter}), 'interpreter', 'tex', 'FontSize', 20);
        ax = gca;
        ax.FontSize = 14;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
        legend("Active", "Control", 'interpreter', 'tex', 'FontSize', 14)

        image_path_dir = fullfile(path2save, strcat("AVG_results", "_", stim_type, "_", recording_day));
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        ylim([fft_y0, fft_y1])

        image_path = fullfile(image_path_dir, strcat("Avg FFT", ".pdf"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, strcat("Avg FFT", ".png"));
        exportgraphics(fig, image_path, 'Resolution', 300);
    end
end