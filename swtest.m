
function [h, p] = swtest(x, alpha)
    % Shapiro-Wilk Test for Normality
    % Input:
    %   x - data vector
    %   alpha - significance level (default: 0.05)
    % Output:
    %   h - test result (0: normal, 1: not normal)
    %   p - p-value of the test
    
    if nargin < 2
        alpha = 0.05;
    end
    
    x = sort(x(:));  % Sort data
    n = length(x);   % Sample size
    
    if n < 3 || n > 5000
        error('Sample size must be between 3 and 5000 for SW test.');
    end
    
    % Shapiro-Wilk coefficients (approximated for simplicity)
    coeffs = swcoeffs(n);  % Function to compute W coefficients
    a = coeffs.a;          % Vector of coefficients
    
    % Compute W statistic
    x_mean = mean(x);
    S = sum((x - x_mean).^2);               % Total sum of squares
    W = (sum(a .* x).^2) / S;               % W statistic
    
    % Approximate p-value (look-up or interpolation)
    p = sw_pvalue(W, n);                    % Function to compute p-value
    
    % Hypothesis result
    h = (p < alpha);
end
