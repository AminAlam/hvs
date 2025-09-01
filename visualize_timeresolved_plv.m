
% Example visualization function
function visualize_timeresolved_plv(plv_timeresolved, time_vector, channel_pair)
    % Input:
    % plv_timeresolved: channels x channels x time windows matrix
    % time_vector: time points corresponding to each window
    % channel_pair: [chan1 chan2] pair of channels to visualize
    
    % Extract PLV time course for specified channel pair
    plv_timecourse = squeeze(plv_timeresolved(channel_pair(1), channel_pair(2), :));
    
    % Create figure
    figure;
    
    % Plot PLV time course
    plot(time_vector, plv_timecourse, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('PLV');
    grid on;
    ylim([0 1]);
    
    % Add mean PLV line
    hold on;
    plot(time_vector([1 end]), [mean(plv_timecourse) mean(plv_timecourse)], '--r', 'LineWidth', 1);
    legend('PLV', 'Mean PLV');
end
