function Z = zscore_columns(X)
    % Z-score each column of a 2D matrix
    % X is the input matrix
    % Z is the matrix where each column is z-scored
    
    % Get the mean and standard deviation of each column
    col_mean = mean(X);
    col_std = std(X);
    
    % Z-score calculation
    Z = (X - col_mean) ./ col_std;
end

