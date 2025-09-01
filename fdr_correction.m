
% Helper function for FDR correction
function [sig_mask, p_adj] = fdr_correction(p_values, alpha)
    % Perform FDR correction on p-values
    p_vec = p_values(:);
    [p_sorted, idx] = sort(p_vec);
    m = length(p_vec);
    
    % Calculate critical values
    i = (1:m)';
    crit_p = (i/m) * alpha;
    
    % Find largest p-value meeting critical value
    k = find(p_sorted <= crit_p, 1, 'last');
    
    if isempty(k)
        thresh = 0;
    else
        thresh = p_sorted(k);
    end
    
    % Create mask and adjusted p-values
    sig_mask = p_values <= thresh;
    size(p_values)
    size(repmat(1:m, [m 1]))
    p_adj = p_values * m ./ repmat(1:sqrt(m), [sqrt(m) 1]);
end
