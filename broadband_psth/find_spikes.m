function [spikes] = find_spikes(amp_signal, sd_size, time)

[total_chan_num, num_data] = size(amp_signal);
spikes = struct;

disp(['Creating spike table...']);
for current_chan_num = 1:total_chan_num
    
    temp_spikes = [];
    signal_average = mean(amp_signal(current_chan_num, :));     %takes mean of entire signal, including windows of stim
    signal_sd = std(amp_signal(current_chan_num, :)) * sd_size;       %takes SD of entire signal, including windows of stim
    threshold = signal_average + signal_sd;

    spike_time_index = find(diff(amp_signal(current_chan_num,:) > threshold) == 1);
    spike_time = zeros(1,length(spike_time_index));
    
    for i=1:length(spike_time_index)
        spike_time(i)=time(spike_time_index(i));
    end
    
    
%     for current_sample = 1:length(time)
%         if amp_signal(current_chan_num, current_sample) > (signal_average + signal_sd) 
%             temp_spikes = [temp_spikes, time(1,current_sample)]; 
%         end    
%     end
%     
    spikes(current_chan_num).spike_table = spike_time;
    
end

end

