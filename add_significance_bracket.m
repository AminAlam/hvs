
% Function to add significance brackets
function add_significance_bracket(x1, x2, y, p_value)
    hold on;
    % Plot the bracket
    plot([x1 x1 x2 x2], [y y y y], 'k-');
    % Add significance annotation
    if p_value < 0.001
        text((x1+x2)/2, y*1.03, '***', 'HorizontalAlignment', 'center', 'FontSize', 16);
    elseif p_value < 0.01
        text((x1+x2)/2, y*1.03, '**', 'HorizontalAlignment', 'center', 'FontSize', 16);
    elseif p_value < 0.05
        text((x1+x2)/2, y*1.03, '*', 'HorizontalAlignment', 'center', 'FontSize', 16);
    else
        text((x1+x2)/2, y*1.03, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 16);
    end
end

