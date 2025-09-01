function p = sw_pvalue(W, n)
    % Approximate p-value for Shapiro-Wilk test based on W and sample size n
    % (In practice, use interpolation or precomputed tables for precision)
    if W < 0.9
        p = 0.01;  % Example low W approximation
    elseif W < 0.95
        p = 0.05;  % Moderate W
    else
        p = 0.1;   % High W
    end
end