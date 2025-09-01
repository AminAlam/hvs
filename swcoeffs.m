
function coeffs = swcoeffs(n)
    % Compute Shapiro-Wilk coefficients for sample size n
    % (Use approximations for simplicity; full implementation requires look-up tables)
    coeffs.a = norminv((1:n) / (n + 1))';  % Normal quantiles
    coeffs.a = coeffs.a / sqrt(sum(coeffs.a.^2));  % Normalize
end

