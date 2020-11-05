function [current_vector] = stimulusPlot(current, start, s_duration, t_duration, step)
    stimulus_current = current;
    stimulus_start = start;
    stimulus_duration = s_duration;
    total_time = t_duration;
    current_vector = zeros(1,total_time / step);
    current_vector (stimulus_start/step:(stimulus_start + stimulus_duration) / step) = stimulus_current;
    % plot(current_vector);
end