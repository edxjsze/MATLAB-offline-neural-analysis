function sep_analysis_path = do_sep_analysis(animal_name, slice_path, config)
    sep_analysis_start = tic;
    fprintf('SEP analysis for %s \n', animal_name);
    analysis_vars = {'baseline_window_start', 'baseline_window_end', 'standard_deviation_coefficient', ...
        'early_response_start', 'early_response_end', 'late_response_start', 'late_response_end', 'ignore_sessions'};
    analysis_log = make_struct_log(config, analysis_vars);
    [sep_analysis_path, failed_path] = create_dir(slice_path, 'sep_analysis');
    file_list = get_file_list(slice_path, '.mat', config.ignore_sessions);
    export_params(sep_analysis_path, 'sep_analysis', failed_path, config);
    for file_index = 1:length(file_list)
        try
            %% Load file contents
            file = [slice_path, '/', file_list(file_index).name];
            [~, filename, ~] = fileparts(file);
            load(file, 'sliced_signal', 'sep_window', 'sep_log');
            %extract the file name
            %% Check sliced variables to make sure they are not empty
            empty_vars = check_variables(file, sliced_signal);
            if empty_vars
                continue
            end            
            
            %% Average sliced data into SEP
            sep_data = average_sliced_data(sliced_signal, config.trial_range);
            
            %% Apply sep analysis

            sep_analysis_results = cal_sep_analysis(animal_name, sep_data, sep_window,...
                config.baseline_window_start, config.baseline_window_end, config.standard_deviation_coefficient, ...
                config.early_response_start, config.early_response_end, config.late_response_start, config.late_response_end);
            
            %% Apply sep region analysis
            % (These analyses are updated if changes are made in the GUI)
            sep_analysis_results = region_sep_analysis(sep_analysis_results);

            %% Saving outputs
            matfile = fullfile(sep_analysis_path, ['analysis_', filename, '.mat']);
            %% Check output to make sure there are no issues with the output
            empty_vars = check_variables(matfile, sep_analysis_results);
            if empty_vars
                continue
            end
            %% Save file if all variables are not empty
                    save(matfile, '-v7.3', 'sep_analysis_results', 'analysis_log', 'sep_log');
        catch ME
            handle_ME(ME, failed_path, filename);
        end
    end
    fprintf('Finished SEP analysis for %s. It took %s s\n', ...
        animal_name, num2str(toc(sep_analysis_start)));
end