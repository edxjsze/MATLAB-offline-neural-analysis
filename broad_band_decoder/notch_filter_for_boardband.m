function notch_filtered_map = notch_filter_for_boardband(board_band_map, notch_filter_frequency, ...
    notch_filter_bandwidth, sample_rate, use_notch_bandstop)   
    notch_filtered_map = [];
    [tot_channel, ~] = size(board_band_map);
    if use_notch_bandstop
        stopband = [notch_filter_frequency - notch_filter_bandwidth/2 ...
            notch_filter_frequency + notch_filter_bandwidth/2];
        for channel_index = 1:1     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%channel_index = 1:tot_channel
            notch_filtered_data = bandstop(board_band_map{channel_index, 2}, stopband, sample_rate);
            notch_filtered_map = [notch_filtered_map; {board_band_map{channel_index, 1}}, ...
                {notch_filtered_data}];
        end
    else
        for channel_index = 1:1     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%channel_index = 1:tot_channel
            notch_filtered_data = notch_filter(board_band_map{channel_index, 2}, sample_rate, ...
                notch_filter_frequency, notch_filter_bandwidth);
            notch_filtered_map = [notch_filtered_map; {board_band_map{channel_index, 1}}, ...
                {notch_filtered_data}];
        end
    end
end