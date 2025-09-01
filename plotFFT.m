function plotFFT(matrix, fs, title_)
    % Function to plot the Fourier Transform power in dB of a 2D matrix
    % Inputs:
    %   matrix - 2D matrix with rows as channels and columns as time points
    %   fs - sampling frequency

    smooth_factor = 200;
    
    [L, numChannels] = size(matrix);
    
    % Compute the FFT and power for each channel
    Y = fft(matrix);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1, :);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs/L*(0:(L/2));
    % 
    % if contains(title_, "active", 'IgnoreCase', true) && contains(title_, "stimulus", 'IgnoreCase', true)
    %     range1Indices = find(f >= 59.28 & f <= 59.42);
    %     range2Indices = find(f >= 59.94 & f <= 60.04);
    %     for ch = 1:numChannels
    %         powerInitial = mean(P1(range2Indices, ch));
    %         powerToTransfer = mean(P1(range1Indices, ch));
    %         P1(range2Indices, ch) = P1(range2Indices, ch) + powerToTransfer;
    %         P1(range1Indices, ch) = powerInitial;
    %     end
    % end
    
    % Plot each channel's power in dB in gray
    hold on;
    for i = 1:numChannels
         % plot(f, smooth(P1(:, i)/mean(P1(:, i), 'all'), smooth_factor), 'Color', [0.8, 0.8, 0.8]); % Gray color
         plot(f, smooth(mag2db(P1(:, i))/2, smooth_factor), 'Color', [0.8, 0.8, 0.8]); % Gray color
    end

    % avgPowerDB = mean(P1/mean(P1), 2);
    P1 = mag2db(P1)/2;
    avgPowerDB = mean(P1, 2);

    
    % Plot the average power in dB in black
    avgPlot = plot(f, smooth(avgPowerDB, smooth_factor), 'k', 'LineWidth', 2); % Black color
    xlim([0, 70])
    
    % Set plot labels, title, and font size
    xlabel('Frequency (Hz)', 'interpreter', 'tex', 'FontSize', 16);
    ylabel('PSD (10Log_{10}(uV^2))', 'interpreter', 'tex', 'FontSize', 16);
    title(title_, 'interpreter', 'latex', 'FontSize', 20);
    % grid on;

    % Increase axis tick font size
    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 1.5;
    ax.Box = 'on';

    % Add legend
    legend(avgPlot, {'Average Power (Black)'}, 'interpreter', 'tex', 'FontSize', 14);
    
    % Add a dummy plot to represent single channels (for legend)
    dummySingle = plot(nan, nan, 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1.5);
    legend([dummySingle, avgPlot], {'Single Channels (Gray)', 'Average Power (Black)'}, 'FontSize', 14, 'Location', 'best');
    
    hold off;
end