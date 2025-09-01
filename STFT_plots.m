%% STFTs
clc
close all
group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};
path2save = fullfile("results", "STFTs");
for subject = list_of_subjects
    subject = string(subject(1));
    for recording_day = recording_days
        min_value = 100;
        max_value = 0;
        recording_day = string(recording_day(1));
        for stim_type = stim_types
            stim_type = string(stim_type(1));
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

            % data = data*1e-6;
            % Calculate STFT for all channels and take the mean
            [num_channels, num_samples] = size(data);
            stft_tmp = stft(data(1,:), fs, 'Window', hamming(256,'periodic'), 'OverlapLength', 128);
            all_stfts = zeros(size(stft_tmp,1), size(stft_tmp,2), num_channels); % Adjust these dimensions based on your STFT parameters
            
            for chan = 1:num_channels
                stft_tmp = stft(data(chan,:), fs, 'Window', hamming(256,'periodic'), 'OverlapLength', 128);
                all_stfts(:,:,chan) = abs(stft_tmp);
            end

            all_stfts = mag2db(all_stfts);
            
            % Calculate mean STFT across channels
            stft_data = mean(all_stfts, 3);
            
            min_tmp = min(stft_data, [], 'all');
            max_tmp = max(stft_data, [], 'all');
            if min_tmp < min_value
                min_value = min_tmp;
            end
            if max_tmp > max_value
                max_value = max_tmp;
            end
        end
        
        for stim_type = stim_types
            stim_type = string(stim_type(1));
            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            % data = data*1e-6;
            recording_id = char(strcat(strrep(subject, '_', ' '), " | ", strrep(recording_day, '_', ' '), " | ", strrep(stim_type, '_', ' ')));
            
            % Calculate mean STFT for plotting
            [num_channels, ~] = size(data);
            stft_tmp = stft(data(1,:), fs, 'Window', hamming(256,'periodic'), 'OverlapLength', 128);
            all_stfts = zeros(size(stft_tmp,1), size(stft_tmp,2), num_channels); % Adjust these dimensions based on your STFT parameters
            
            
            for chan = 1:num_channels
                [stft_tmp, f, t] = stft(data(chan,:), fs, 'Window', hamming(256,'periodic'), 'OverlapLength', 128);
                all_stfts(:,:,chan) = mag2db(abs(stft_tmp));
            end
            
            mean_stft = mean(all_stfts, 3);
            std_stft = std(all_stfts, 0, 3);
            
            % Create figure and plot
            fig = figure;
            set(gcf, 'Position', get(0, 'Screensize'));
            
            % %%%%%% Plot the mean STFT %%%%%%
            imagesc(t, f, mean_stft);
            axis xy;
            colormap parula;
            clim([min_value, max_value]);
            ylim([0, 70]);
            
            % Add labels and title
            xlabel('Time (s)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
            ylabel('Frequency (Hz)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
            title(sprintf('Mean STFT across channels | %s', recording_id), 'interpreter','tex', 'FontSize', fontsize_title);
            colorbar_handle = colorbar; % Create colorbar
            % Add label to the colorbar
            ylabel(colorbar_handle, 'Amplitude (dB)', 'interpreter', 'tex', 'FontSize', fontsize_axis); % Customize the label and font size
    
            ax = gca;
            ax.FontSize = fontsize_axis;
            ax.LineWidth = 1.5;
            ax.Box = 'on';
            
            % Save the figure
            image_path_dir = fullfile(path2save, subject);
            image_path = fullfile(image_path_dir, strcat(recording_id, "_mean.pdf"));
            image_path = fullfile(image_path_dir, strcat(recording_id, "_mean.png"));
            if ~exist(image_path_dir, 'dir')
                mkdir(image_path_dir);
            end
            exportgraphics(fig, image_path, 'Resolution', 300);

            % %%%%%% Plot the std STFT %%%%%%
            fig = figure;
            set(gcf, 'Position', get(0, 'Screensize'));
            imagesc(t, f, std_stft);
            axis xy;
            colormap parula;
            clim([0, 10]);
            ylim([0, 70]);
            

            % Add labels and title
            xlabel('Time (s)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
            ylabel('Frequency (Hz)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
            title(sprintf('Std STFT across channels | %s', recording_id), 'interpreter','tex', 'FontSize', fontsize_title);
            colorbar_handle = colorbar; % Create colorbar
            % Add label to the colorbar
            ylabel(colorbar_handle, 'Amplitude (dB)', 'interpreter', 'tex', 'FontSize', fontsize_axis); % Customize the label and font size

            ax = gca;
            ax.FontSize = fontsize_axis;
            ax.LineWidth = 1.5;
            ax.Box = 'on';

            % Save the figure
            image_path_dir = fullfile(path2save, subject);
            image_path = fullfile(image_path_dir, strcat(recording_id, "_std.pdf"));
            exportgraphics(fig, image_path, 'Resolution', 300);
            image_path = fullfile(image_path_dir, strcat(recording_id, "_std.png"));
            exportgraphics(fig, image_path, 'Resolution', 300);

            close all;
        end
    end
end
%% Full recording STFT
group_labels = {'Control', 'Active'};
day_labels = {'Day 1', 'Day 5', 'Day 19'};
light_conditions = {'No Light', 'Constant Light', 'Stimulus Light'};
path2save = fullfile("results", "STFTs");
min_value = 100;
max_value = 0;
for subject = list_of_subjects
    subject = string(subject(1));
    for recording_day = recording_days
        recording_day = string(recording_day(1));
        data_all = [];
        for stim_type = stim_types
            stim_type = string(stim_type(1));
            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            data = data*micro_volt_constant;
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
            data_all = [data_all, data];
        end

        % Calculate STFT for all channels and take the mean
        [num_channels, num_samples] = size(data_all);
        stft_tmp = stft(data_all(1,:), fs, 'Window', hamming(256,'periodic'), 'OverlapLength', 128);
        all_stfts = zeros(size(stft_tmp,1), size(stft_tmp,2), num_channels); % Adjust these dimensions based on your STFT parameters
        
        for chan = 1:num_channels
            stft_tmp = stft(data_all(chan,:), fs, 'Window', hamming(256,'periodic'), 'OverlapLength', 128);
            all_stfts(:,:,chan) = abs(stft_tmp);
        end
        all_stfts = mag2db(all_stfts);
        
        % Calculate mean STFT across channels
        stft_data = mean(all_stfts, 3);
        
        min_tmp = min(stft_data, [], 'all');
        max_tmp = max(stft_data, [], 'all');
        if min_tmp < min_value
            min_value = min_tmp;
        end
        if max_tmp > max_value
            max_value = max_tmp;
        end
        
    end
end
min_value = -100;
max_value = -55;

for subject = list_of_subjects
    subject = string(subject(1));
    for recording_day = recording_days
        recording_day = string(recording_day(1));
        data_all = [];
        t_change = [];
        for stim_type = stim_types
            stim_type = string(stim_type(1));
            data = data_struct.(sprintf(subject)).(sprintf(recording_day)).(sprintf(stim_type)).data_clean;
            data = data*micro_volt_constant;
            data_all = [data_all, data];
            t_change = [t_change, length(data_all)];
        end

        recording_id = char(strcat(strrep(subject, '_', ' '), " | ", strrep(recording_day, '_', ' '), " | ", strrep("all", '_', ' ')));
        
        % Calculate mean STFT for plotting
        [num_channels, ~] = size(data_all);
        stft_tmp = stft(data_all(1,:), fs, 'Window', hamming(512,'periodic'), 'OverlapLength', 32);
        all_stfts = zeros(size(stft_tmp,1), size(stft_tmp,2), num_channels); % Adjust these dimensions based on your STFT parameters
        
        
        for chan = 1:num_channels
            [stft_tmp, f, t] = stft(data_all(chan,:), fs, 'Window', hamming(512,'periodic'), 'OverlapLength', 32);
            all_stfts(:,:,chan) = mag2db(abs(stft_tmp).^2)/2;
        end
        
        mean_stft = mean(all_stfts, 3);
        std_stft = std(all_stfts, 0, 3);
        
        % Create figure and plot
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        
        % Plot the mean STFT
        imagesc(t, f, mean_stft);
        axis xy;
        colormap parula;
        clim([min_value, max_value]);
        ylim([0, 70]);
        hold on
        for t_ = t_change
            plot([t_/fs, t_/fs], [0, 70], 'LineWidth', 2, 'Color', 'k')
        end
        
        % Add labels and title
        xlabel('Time (s)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
        ylabel('Frequency (Hz)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
        title(sprintf('Mean STFT across channels | %s', recording_id), 'interpreter','tex', 'FontSize', fontsize_title);
        colorbar_handle = colorbar;
        ylabel(colorbar_handle, 'Power (dB)', 'interpreter', 'tex', 'FontSize', fontsize_axis); % Customize the label and font size
        
        ax = gca;
        ax.FontSize = fontsize_axis;
        ax.LineWidth = 1.5;
        ax.Box = 'on';
        
        % Save the figure
        image_path_dir = fullfile(path2save, 'full_recording');
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        image_path = fullfile(image_path_dir, strcat(recording_id, "_mean.pdf"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, strcat(recording_id, "_mean.png"));
        exportgraphics(fig, image_path, 'Resolution', 300);



        % %%%%%% Plot the std STFT %%%%%%
        fig = figure;
        set(gcf, 'Position', get(0, 'Screensize'));
        imagesc(t, f, std_stft);
        axis xy;
        colormap parula;
        clim([0, 15]);
        ylim([0, 70]);

        % Add labels and title
        xlabel('Time (s)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
        ylabel('Frequency (Hz)', 'interpreter', 'tex', 'FontSize', fontsize_labels);
        title(sprintf('Std STFT across channels | %s', recording_id), 'interpreter','tex', 'FontSize', fontsize_title);
        colorbar_handle = colorbar; % Create colorbar
        % Add label to the colorbar
        ylabel(colorbar_handle, 'Amplitude (dB)', 'interpreter', 'tex', 'FontSize', fontsize_axis); % Customize the label and font size

        ax = gca;
        ax.FontSize = fontsize_axis;
        ax.LineWidth = 1.5;
        ax.Box = 'on';

        % Save the figure
        image_path_dir = fullfile(path2save, 'full_recording');
        if ~exist(image_path_dir, 'dir')
            mkdir(image_path_dir);
        end
        image_path = fullfile(image_path_dir, strcat(recording_id, "_std.pdf"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        image_path = fullfile(image_path_dir, strcat(recording_id, "_std.png"));
        exportgraphics(fig, image_path, 'Resolution', 300);
        close all;
    end
end
