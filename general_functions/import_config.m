function [config] = import_config(animal_path, analysis_type)

    %% Grabs config file and creates labels
    animal_csv_path = [animal_path, '/*.csv'];
    csv_files = dir(animal_csv_path);
    for csv = 1:length(csv_files)
        csv_file = fullfile(animal_path, csv_files(csv).name);
        if contains(csv_files(csv).name, [analysis_type, '_config.csv'])
            config_table = readtable(csv_file);
        end
    end

    if exist('config_table') == 0
        error('Must have config file to run analysis')
    end

    config = struct;
    for i = 1:height(config_table)
        current_name = config_table.key{i};
        current_value = config_table.value{i};
        value = convert_string(current_value);
        config.(current_name) = value;
    end
end

function [value] = convert_string(string_value)
    if all(ismember(string_value, '0123456789+-.eEdD,; '))
        value = str2num(string_value);
    elseif strcmpi(string_value, 'true')
        value = true;
    elseif strcmpi(string_value, 'false')
        value = false;
    else
        value = string_value;
    end
end