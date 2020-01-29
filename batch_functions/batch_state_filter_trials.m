function [psth_path] = batch_state_filter_trials(data_path, kalman_path, config)

    if config.use_trajectories
        psth_path = [data_path, '/psth'];
        kalman_list = get_file_list(kalman_path, '.mat', config.ignore_sessions);
        psth_list = get_file_list(psth_path, '.mat', config.ignore_sessions);
        [gpfa_path, ~] = create_dir(data_path, 'gpfa');
    else
        kalman_list = get_file_list(kalman_path, '.mat', config.ignore_sessions);
        psth_list = get_file_list(data_path, '.mat', config.ignore_sessions);
    end

    for file_i = 1:length(kalman_list)
        filename = kalman_list(file_i).name;
        kalman_file = fullfile(kalman_path, filename);
        [~, ~, ~, session_num, session_date, ~] = get_filename_info(erase(filename, '.mat'));
        psth_file_struct = psth_list(contains({psth_list.name}, num2str(session_num)) & ...
            contains({psth_list.name}, num2str(session_date)));
        psth_file = fullfile(psth_file_struct.folder, psth_file_struct.name);

        %% load state table from kalman mat
        load(kalman_file, 'measurements');
        load(psth_file, 'psth_struct', 'labeled_data', 'event_ts');

        %TODO update all_events field in psth_struct and trials in state table with desired trials
        [psth_struct, state_struct] = remove_trials(measurements, psth_struct, labeled_data);

        if config.use_trajectories
            save(kalman_file, 'state_struct', '-append');
            save(fullfile(gpfa_path, psth_file_struct.name), 'psth_struct', 'labeled_data', 'event_ts');
            psth_path = gpfa_path;
        else
            save(kalman_file, 'state_struct', '-append');
            save(psth_file, 'psth_struct', 'labeled_data', 'event_ts', '-append');
        end
    end

end