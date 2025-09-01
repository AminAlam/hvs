function [P1, f] = calcFFT(matrix, fs)
    [L, numChannels] = size(matrix);
    Y = fft(matrix);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1, :);
    avgPowerDB = mag2db(P1)/2;

    f = fs/L*(0:(L/2));
    % range1Indices = find(f >= 59.28 & f <= 59.42);
    % range2Indices = find(f >= 59.98 & f <= 60.02);
    % for ch = 1:numChannels
    %     powerInitial = mean(P1(range2Indices, ch));
    %     powerToTransfer = mean(P1(range1Indices, ch));
    %     P1(range2Indices, ch) = P1(range2Indices, ch) + powerToTransfer;
    %     P1(range1Indices, ch) = powerInitial;
    % end
    % psd = (1/(fs*L)) * abs(P1).^2;
    % psd(2:end-1) = 2*psd(2:end-1);
    f = fs/L*(0:(L/2));
end