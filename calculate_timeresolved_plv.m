
function [plv_timeresolved, time_vector] = calculate_timeresolved_plv(data, fs, freq_range, window_size, window_step)
    % Input:
    % data: channels x timepoints matrix
    % fs: sampling frequency
    % freq_range: [low_freq high_freq] frequency range to analyze
    % window_size: size of sliding window in seconds
    % window_step: step size for sliding window in seconds
    
    % Convert window parameters from seconds to samples
    window_samples = round(window_size * fs);
    step_samples = round(window_step * fs);
    
    % Get dimensions
    [no_channels, time_points] = size(data);
    
    % Calculate number of windows
    no_windows = floor((time_points - window_samples) / step_samples) + 1;
    
    % Generate time vector for plotting
    time_vector = ((window_samples/2) : step_samples : (window_samples/2 + (no_windows-1)*step_samples)) / fs;
    
    % Filter the data in the frequency band of interest
    [b, a] = butter(2, freq_range/(fs/2), 'bandpass');
    filtered_data = filtfilt(b, a, data.').';
    
    % Initialize time-resolved PLV matrix
    plv_timeresolved = zeros(no_channels, no_channels, no_windows);
    
    % Calculate phase using Hilbert transform for all channels
    phase_data = angle(hilbert(filtered_data.'));  % time_points x channels
    
    % Calculate PLV for each time window
    for win = 1:no_windows
        % Define window indices
        start_idx = (win-1)*step_samples + 1;
        end_idx = start_idx + window_samples - 1;
        
        % Extract phases for current window
        window_phases = phase_data(start_idx:end_idx, :);
        
        % Calculate PLV between all channel pairs for this window
        for chan1 = 1:no_channels
            for chan2 = chan1+1:no_channels
                % Calculate phase difference between channels
                phase_diff = window_phases(:,chan1) - window_phases(:,chan2);
                
                % Calculate PLV for this window
                plv_value = abs(mean(exp(1i * phase_diff)));
                
                % Store the PLV value (symmetric matrix)
                plv_timeresolved(chan1, chan2, win) = plv_value;
                plv_timeresolved(chan2, chan1, win) = plv_value;
            end
        end
        
        % Set diagonal to 1 (PLV of a channel with itself)
        plv_timeresolved(:,:,win) = plv_timeresolved(:,:,win) + eye(no_channels);
    end
end
