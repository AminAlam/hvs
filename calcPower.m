function powerInBand = calcPower(matrix, fs, f_range)

    [psd, f] = calcFFT(matrix, fs);
    freqIndices = find(f >= f_range(1) & f <= f_range(2));
    powerInBand = sum(psd(freqIndices));
end