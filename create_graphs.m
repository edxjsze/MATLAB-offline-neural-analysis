function [] = create_graphs()

    %% Extracting data from excel sheet and standard deviation shading on line plots
    %% credit is due towards Nate Bridges.

    %% Colors
    max_color = 255;
    dark_green = [0 (102/max_color) 0];
    light_green = [(102/max_color) (153/max_color) 0];
    black = [0 0 0];
    burgandy = [(102/max_color) 0 0]; % burgandy
    orange = [(255/max_color) (51/max_color) 0]; % orange
    % purple = [(102/max_color) 0 (255/max_color)];
    % lavender = [(153/max_color) (153/max_color) (255/max_color)];    
    close all;
    [filename, file_path] = uigetfile(pwd);
    filename = fullfile(file_path, filename);
    [dir_path, ~, ~] = fileparts(filename);
    [~, sheets, ~] = xlsfinfo(filename);

    graph_path = strcat(dir_path, '/graphs');
    if ~exist(graph_path, 'dir')
       mkdir(dir_path, 'graphs');
    end

    axis_limits = [];
    for i = 1:length(sheets)
        sheet_name = sheets{i};
        [data,txt,raw] = xlsread(filename, sheet_name);
        learning_mean = [];
        non_learning_mean = [];
        control_mean = [];
        learning_std = [];
        non_learning_std = [];
        control_std = [];
        learning_std_error = [];
        non_learning_std_error = [];
        control_std_error = [];
        overall_mean = [];
        % outputTable = table;
        if contains(sheet_name, 'Bar')
            figure('visible','on');
            if contains(sheet_name, 'Performance')
                all_means = [data(3) data(13) data(23); data(8) data(18) data(28)];
                all_std = [data(4) data(14) data(24); data(9) data(19) data(29)];
                upper_bounds = all_means + all_std;
                lower_bounds = all_means - all_std;
                categories = categorical({'Early','Late'});
                
                ax = axes;
                b = bar(all_means, 'BarWidth', 1);
                xticks(ax,[1 2]);
                xticklabels(ax,{ 'Early', 'Late'});
                
                hold on;
                %% Adds error bars
                groups = size(all_means, 1);
                bars = size(all_means, 2);
                groupwidth = min(0.8, bars/(bars + 1.5));
                for k = 1:bars
                    center = (1:groups) - groupwidth/2 + (2*k-1) * groupwidth / (2*bars);
                    errorbar(center, all_means(:,k), all_std(:,k), 'k', 'linestyle', 'none');
                end
                current_graph = gca;
                current_graph.Clipping = 'off';

                %% Set axes limits
                sheet_strs = strsplit(sheet_name, '_');
                current_sheet = sheet_strs{1};
                all_sheets = axis_limits(:,1);
                limit_index = find(contains(all_sheets, current_sheet));
                ylim([axis_limits{limit_index, 2} axis_limits{limit_index, 4}])
                yticks([axis_limits{limit_index, 2} axis_limits{limit_index, 3} axis_limits{limit_index, 4}]);

                %% Sets colors for bars
                for k = 1:size(all_means,2)
                    if mod(k, 3) == 0
                        b(k).FaceColor = [0 0 0];
                        b(k-1).FaceColor = light_green; % green
                        b(k-2).FaceColor = dark_green; % magenta
                    end
                end

                %% Creates Legends
                lg = legend('Learning','Non-learning','Control');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(sheet_name);
                ytickformat('%.2f')
                hold off;
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            elseif contains(sheet_name, 'SynRed')
                all_means = [data(23) data(3); data(28) data(8)];
                all_std = [data(24) data(4); data(29) data(9)];
                upper_bounds = all_means + all_std;
                lower_bounds = all_means - all_std;
                categories = categorical({'Early','Late'});
                
                ax = axes;
                b = bar(all_means, 'BarWidth', 1);
                xticks(ax,[1 2]);
                xticklabels(ax,{ 'Early', 'Late'});
                
                hold on;
                %% Adds error bars
                groups = size(all_means, 1);
                bars = size(all_means, 2);
                groupwidth = min(0.8, bars/(bars + 1.5));
                for k = 1:bars
                    center = (1:groups) - groupwidth/2 + (2*k-1) * groupwidth / (2*bars);
                    errorbar(center, all_means(:,k), all_std(:,k), 'k', 'linestyle', 'none');
                end
                current_graph = gca;
                current_graph.Clipping = 'off';

                %% Set axes limits
                sheet_strs = strsplit(sheet_name, '_');
                current_sheet = sheet_strs{1};
                all_sheets = axis_limits(:,1);
                limit_index = find(contains(all_sheets, current_sheet));
                ylim([axis_limits{limit_index, 2} axis_limits{limit_index, 4}])
                yticks([axis_limits{limit_index, 2} axis_limits{limit_index, 3} axis_limits{limit_index, 4}]);

                %% Sets colors for bars
                for k = 1:size(all_means,2)
                    if mod(k, 2) == 0
                        b(k).FaceColor = 'b';
                        b(k-1).FaceColor = 'r';
                    end
                end

                %% Creates Legends
                lg = legend('Direct','Indirect');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(sheet_name);
                ytickformat('%.2f')
                hold off;
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            elseif contains(sheet_name, 'Timng')
                %% Early
                all_means = [data(3) data(13); data(43) data(53)];
                all_std = [data(4) data(14); data(44) data(54)];
                upper_bounds = all_means + all_std;
                lower_bounds = all_means - all_std;
                categories = categorical({'Indirect','Direct'});
                
                ax = axes;
                b = bar(all_means, 'BarWidth', 1);
                xticks(ax,[1 2]);
                xticklabels(ax,{ 'Indirect', 'Direct'});
                
                hold on;
                %% Adds error bars
                groups = size(all_means, 1);
                bars = size(all_means, 2);
                groupwidth = min(0.8, bars/(bars + 1.5));
                for k = 1:bars
                    center = (1:groups) - groupwidth/2 + (2*k-1) * groupwidth / (2*bars);
                    errorbar(center, all_means(:,k), all_std(:,k), 'k', 'linestyle', 'none');
                end
                current_graph = gca;
                current_graph.Clipping = 'off';

                %% Set axes limits
                % max and min are taken twice since the first call returns a 3x1 matrix with the max and min
                % value across the categories. The second call then returns the max and min for the overall graph.
                upper_limit = max(max(upper_bounds));
                lower_limit = min(min(lower_bounds));
                % Add a buffer to the bounds so that the error bars are not on the exact edge of the graph
                upper_limit = upper_limit + (upper_limit) * 0.05;
                lower_limit = lower_limit + (lower_limit) * 0.05;
                y_midpoint = (upper_limit + 0)/2;
                ylim([min([0 lower_limit]) upper_limit]);
                yticks([0 y_midpoint upper_limit]);

                %% Sets colors for bars
                for k = 1:size(all_means,2)
                    if mod(k, 2) == 0
                        b(k).FaceColor = burgandy;
                        b(k-1).FaceColor = orange;
                    end
                end

                %% Creates Legends
                lg = legend('Timing','Count');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(['EARLY ' sheet_name]);
                hold off;
                graph_name = ['EARLY_', sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));

                %% Late
                figure('visible','on'); 
                all_means = [data(8) data(18); data(48) data(58)];
                all_std = [data(9) data(19); data(49) data(59)];
                upper_bounds = all_means + all_std;
                lower_bounds = all_means - all_std;
                
                ax = axes;
                b = bar(all_means, 'BarWidth', 1);
                xticks(ax,[1 2]);
                xticklabels(ax,{ 'Indirect', 'Direct'});
                
                hold on;
                %% Adds error bars
                groups = size(all_means, 1);
                bars = size(all_means, 2);
                groupwidth = min(0.8, bars/(bars + 1.5));
                for k = 1:bars
                    center = (1:groups) - groupwidth/2 + (2*k-1) * groupwidth / (2*bars);
                    errorbar(center, all_means(:,k), all_std(:,k), 'k', 'linestyle', 'none');
                end
                current_graph = gca;
                current_graph.Clipping = 'off';

                %% Set axes limits
                % max and min are taken twice since the first call returns a 3x1 matrix with the max and min
                % value across the categories. The second call then returns the max and min for the overall graph.
                upper_limit = max(max(upper_bounds));
                lower_limit = min(min(lower_bounds));
                % Add a buffer to the bounds so that the error bars are not on the exact edge of the graph
                upper_limit = upper_limit + (upper_limit) * 0.05;
                lower_limit = lower_limit + (lower_limit) * 0.05;
                y_midpoint = (upper_limit + 0)/2;
                ylim([min([0 lower_limit]) upper_limit]);
                yticks([0 y_midpoint upper_limit]);

                %% Sets colors for bars
                for k = 1:size(all_means,2)
                    if mod(k, 2) == 0
                        b(k).FaceColor = burgandy;
                        b(k-1).FaceColor = orange;
                    end
                end

                %% Creates Legends
                lg = legend('Timing','Count');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(['LATE ' sheet_name]);
                hold off;
                graph_name = ['LATE_', sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            else
                all_means = [data(23) data(3) data(13) ; data(28) data(8) data(18)];
                all_std = [data(24) data(4) data(14); data(29) data(9) data(19)];
                upper_bounds = all_means + all_std;
                lower_bounds = all_means - all_std;
                categories = categorical({'Early','Late'});
                
                ax = axes;
                b = bar(all_means, 'BarWidth', 1);
                xticks(ax,[1 2]);
                xticklabels(ax,{ 'Early', 'Late'});
                
                hold on;
                %% Adds error bars
                groups = size(all_means, 1);
                bars = size(all_means, 2);
                groupwidth = min(0.8, bars/(bars + 1.5));
                for k = 1:bars
                    center = (1:groups) - groupwidth/2 + (2*k-1) * groupwidth / (2*bars);
                    errorbar(center, all_means(:,k), all_std(:,k), 'k', 'linestyle', 'none');
                end
                current_graph = gca;
                current_graph.Clipping = 'off';

                %% Set axes limits
                sheet_strs = strsplit(sheet_name, '_');
                current_sheet = sheet_strs{1};
                all_sheets = axis_limits(:,1);
                limit_index = find(contains(all_sheets, current_sheet));
                ylim([axis_limits{limit_index, 2} axis_limits{limit_index, 4}])
                yticks([axis_limits{limit_index, 2} axis_limits{limit_index, 3} axis_limits{limit_index, 4}]);

                %% Sets colors for bars
                for k = 1:size(all_means,2)
                    if mod(k, 3) == 0
                        b(k).FaceColor = [0 0 0];
                        b(k-1).FaceColor = [0 0 1];
                        b(k-2).FaceColor = [1 0 0];
                    end
                end

                %% Creates Legends
                lg = legend('Direct','Indirect','Control');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(sheet_name);
                ytickformat('%.2f')
                hold off;
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            end
        elseif contains(sheet_name, '_Line')
            figure('visible','on');
            if contains(sheet_name, 'Performance')
                %% Extract relevant variables
                dayCol=1; % column in spreadsheet where "Day" located  
                animalTypeCol=2; % column in spreadsheet where "Animal_Type" located
                meanCol=3; % column in spreadsheet where "Mean" text located
                measureCol=5; % column in spreadsheet where values to be plotted located

                % day
                dayCell=cellfun(@(x) str2double(x) ,raw(:,dayCol), 'UniformOutput',false);
                days=sort(repmat([dayCell{~isnan([dayCell{:}])}]',3,1));

                % mean values
                meanCell=cellfun(@(x) strcmpi(x,'Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                meanLogical=[meanCell{:}];
                means=[raw{[meanCell{:}],measureCol}]';

                % animal type (e.g. learning)
                animalTypeCell=cellfun(@(x) (strcmpi(x,'Learining')|strcmpi(x,'Non-learning')|strcmpi(x,'Control'))...
                ,raw(:,animalTypeCol),...
                'UniformOutput',false);
                animalTypeLogical=[animalTypeCell{:}];
                animalTypes=raw([animalTypeCell{:}],animalTypeCol);

                % standard error
                standardErrorCell=cellfun(@(x) strcmpi(x,'Std. Error of Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                standardErrors=[raw{[standardErrorCell{:}],measureCol}]';

                % compile into table for processing 
                % fprintf('%d %d %d %d\n', length(days), length(animalTypes), length(means), length(standardErrors));
                
                outputTable=table(days,animalTypes,means,standardErrors);
                open('outputTable')
                
                %% Parse data by relevant group (e.g. learing, nonlearning, control)

                % learning animals
                learning_x=outputTable.days(strcmpi(outputTable.animalTypes,'Learining'));  % x-values used in plotting 
                learning_y=outputTable.means(strcmpi(outputTable.animalTypes,'Learining')); % y-values used in plotting 
                learning_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Learining')); % error bar-values used in plotting 

                % non-learning animals
                nonlearning_x=outputTable.days(strcmpi(outputTable.animalTypes,'Non-learning')); % x-values used in plotting 
                nonlearning_y=outputTable.means(strcmpi(outputTable.animalTypes,'Non-learning')); % y-values used in plotting 
                nonlearning_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Non-learning')); % error bar-values used in plotting 

                % control animals 
                control_x=outputTable.days(strcmpi(outputTable.animalTypes,'Control')); % x-values used in plotting 
                control_y=outputTable.means(strcmpi(outputTable.animalTypes,'Control')); % y-values used in plotting 
                control_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Control')); % error bar-values used in plotting 

                y_global_min = min(outputTable.means - outputTable.standardErrors);
                y_global_max = max(outputTable.means + outputTable.standardErrors);
                y_midpoint = (y_global_max + y_global_min)/2;
                %% Plot figure 

                % generate line with shading
                close
                figure
                hold on;
                plot(learning_x, learning_y, 'Color', dark_green, 'LineWidth', 4);
                plot(nonlearning_x, nonlearning_y, 'Color', light_green, 'LineWidth', 4);
                plot(control_x, control_y, 'Color', 'k', 'LineWidth', 4);
                transparency=0.2; % transparency of background
                [l,p] = boundedline(learning_x, learning_y, learning_e, ...
                    nonlearning_x, nonlearning_y, nonlearning_e, 'b',...
                    control_x, control_y, control_e, 'cmap', [dark_green; light_green; 0 0 0],...
                    'transparency', transparency);


                % figure properties 
                innerlineWidth=3;  % width of mean line 

                % apply properties to all lines in figure 
                for lineProperty=1:length(l)
                    l(lineProperty).LineWidth=innerlineWidth;
                end

                % axis labels 
                ylabel('Performance Change');
                xlabel('Day');
                axis tight;
                yticks([y_global_min 0 y_global_max]);
                axis_limits = [axis_limits; sheet_name, {y_global_min}, {0}, {y_global_max}];
                %% Legend
                lg = legend('Learning','Non-learning','Control');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title('Offline Performance');
                ytickformat('%.2f')
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
                
                % for point = 1:length(data)
                %     if mod(point, 15) == 5
                %         %Learning
                %         learning_mean = [learning_mean; data(point-2)];
                %         learning_std = [learning_std; data(point-1)];
                %         learning_std_error = [learning_std_error; data(point)];
                %     elseif mod(point, 15) == 10
                %         %Non-learning
                %         non_learning_mean = [non_learning_mean; data(point-2)];
                %         non_learning_std = [non_learning_std; data(point-1)];
                %         non_learning_std_error = [non_learning_std_error; data(point)];
                %     elseif mod(point, 15) == 0
                %         %Control
                %         control_mean = [control_mean; data(point-2)];
                %         control_std = [control_std; data(point-1)];
                %         control_std_error = [control_std_error; data(point)];
                %     end
                % end

                % % x = [[1:length(learning_mean)]';[length(learning_mean):-1:1]'];
                % % close all
                % % hold on
                % % fill(x,[learning_neg;flipud(learning_pos)],[.9 .9 .9],'linestyle','none');
                % % plot(learning_mean)
                
                % learning_pos = learning_mean + learning_std_error;
                % learning_neg = learning_mean - learning_std_error;
                % non_learning_pos = non_learning_mean + non_learning_std_error;
                % non_learning_neg = non_learning_mean - non_learning_std_error;
                % control_pos = control_mean + control_std_error;
                % control_neg = control_mean - control_std_error;
                % x = [[1:length(learning_mean)]';[length(learning_mean):-1:1]'];
                % fill(x,[learning_neg;flipud(learning_pos)],[.9 .9 .9],'linestyle','none');
                % hold on;
                % % fill(x,[non_learning_neg;flipud(non_learning_pos)],[.5 .5 .5],'linestyle','none');
                % % fill(x,[control_neg;flipud(control_pos)],[.9 .9 .9],'linestyle','none');
                % plot(learning_mean);
                % % plot(non_learning_mean, 'Color', 'm');
                % % plot(control_mean, 'k');
                % hold off;


                % plot(non_learning_mean, 'Color', 'm');
                % title(sheet_name);
            elseif contains(sheet_name, 'Timng')
                %% Extract relevant variables
                meanCol=4; % column in spreadsheet where "Mean" text located
                measureCol=6; % column in spreadsheet where values to be plotted located

                % day
                days = repmat(([1:26]-1), [1,2])';
                % days=sort(dayCell);

                % mean values
                meanCell=cellfun(@(x) strcmpi(x,'Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                meanLogical=[meanCell{:}];
                means=[raw{[meanCell{:}],measureCol}]';

                % animal type (e.g. learning)
                direct_labels = repmat({'Direct'}, [52,1]);
                indirect_labels = repmat({'Indirect'}, [52,1]);

                % standard error
                standardErrorCell=cellfun(@(x) strcmpi(x,'Std. Error of Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                standardErrors=[raw{[standardErrorCell{:}],measureCol}]';

                % compile into table for processing 
                indirect_means = means(1:52);
                indirect_errors = standardErrors(1:52);
                direct_means = means((end-51):end);
                direct_errors = standardErrors((end-51):end);
                indirect_table = table(days, indirect_labels, indirect_means, indirect_errors);
                direct_table = table(days, direct_labels, direct_means, direct_errors);
                %% Create Direct graphs
                open('direct_table')

                %% Parse data by relevant group (e.g. learing, nonlearning, control)

                %% direct
                timing_x=direct_table.days(1:26);  % x-values used in plotting 
                timing_y=direct_table.direct_means(1:26); % y-values used in plotting 
                timing_e=direct_table.direct_errors(1:26); % error bar-values used in plotting 

                count_x=direct_table.days(27:end);  % x-values used in plotting 
                count_y=direct_table.direct_means(27:end); % y-values used in plotting 
                count_e=direct_table.direct_errors(27:end); % error bar-values used in plotting 

                y_global_min = min(direct_table.direct_means - direct_table.direct_errors);
                y_global_max = max(direct_table.direct_means + direct_table.direct_errors);
                y_midpoint = (y_global_max + y_global_min)/2;

                %% Plot figure 

                % generate line with shading
                close
                figure
                hold on;
                transparency=0.2; % transparency of background
                plot(timing_x, timing_y, 'Color', orange, 'LineWidth', 4);
                plot(count_x, count_y, 'Color', burgandy, 'LineWidth', 4);
                [l,p] = boundedline(timing_x, timing_y, timing_e,...
                    count_x, count_y, count_e, 'cmap', [orange; burgandy],'transparency', transparency);

                % figure properties 
                innerlineWidth=4;  % width of mean line 

                % apply properties to all lines in figure 
                for lineProperty=1:length(l)
                    l(lineProperty).LineWidth=innerlineWidth;
                end

                axis tight;
                yticks([y_global_min y_midpoint y_global_max]);
                axis_limits = [axis_limits; sheet_name, {y_global_min}, {y_midpoint}, {y_global_max}];
                %% Legend
                lg = legend('Timing','Count');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(['DIRECT ' sheet_name]);
                ytickformat('%.2f')
                graph_name = ['DIRECT_' sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));

                %% Indirect
                open('indirect_table')
                timing_x=indirect_table.days(1:26);  % x-values used in plotting 
                timing_y=indirect_table.indirect_means(1:26); % y-values used in plotting 
                timing_e=indirect_table.indirect_errors(1:26); % error bar-values used in plotting 

                count_x=indirect_table.days(27:end);  % x-values used in plotting 
                count_y=indirect_table.indirect_means(27:end); % y-values used in plotting 
                count_e=indirect_table.indirect_errors(27:end); % error bar-values used in plotting

                y_global_min = min(indirect_table.indirect_means - indirect_table.indirect_errors);
                y_global_max = max(indirect_table.indirect_means + indirect_table.indirect_errors);
                y_midpoint = (y_global_max + y_global_min)/2;
                %% Plot figure 

                % generate line with shading
                close
                figure
                hold on;
                transparency=0.2; % transparency of background
                plot(timing_x, timing_y, 'Color', orange, 'LineWidth', 4);
                plot(count_x, count_y, 'Color', burgandy, 'LineWidth', 4);
                [l,p] = boundedline(timing_x, timing_y, timing_e,...
                    count_x, count_y, count_e, 'cmap', [orange; burgandy],'transparency', transparency);

                % figure properties 
                innerlineWidth=4;  % width of mean line 

                % apply properties to all lines in figure 
                for lineProperty=1:length(l)
                    l(lineProperty).LineWidth=innerlineWidth;
                end

                axis tight;
                yticks([y_global_min y_midpoint y_global_max]);
                axis_limits = [axis_limits; sheet_name, {y_global_min}, {y_midpoint}, {y_global_max}];
                %% Legend
                lg = legend('Timing','Count');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(['INDIRECT ' sheet_name]);
                ytickformat('%.2f')
                graph_name = ['INDIRECT_' sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            elseif contains(sheet_name, 'SynRed')
                %% Extract relevant variables
                meanCol=3; % column in spreadsheet where "Mean" text located
                measureCol=5; % column in spreadsheet where values to be plotted located

                % day
                days = repmat(([1:26]-1), [1,2])';
                % days=sort(dayCell);

                % mean values
                meanCell=cellfun(@(x) strcmpi(x,'Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                meanLogical=[meanCell{:}];
                means=[raw{[meanCell{:}],measureCol}]';

                % animal type (e.g. learning)
                animalTypes = [repmat({'Indirect'}, [26,1]); repmat({'Direct'}, [26,1])];

                % standard error
                standardErrorCell=cellfun(@(x) strcmpi(x,'Std. Error of Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                standardErrors=[raw{[standardErrorCell{:}],measureCol}]';

                indirect_means = means(1:26);
                indirect_errors = standardErrors(1:26);
                direct_means = means((end-25):end);
                direct_errors = standardErrors((end-25):end);
                means = [indirect_means; direct_means];
                standardErrors = [indirect_errors; direct_errors];

                % compile into table for processing
                outputTable=table(days,animalTypes,means,standardErrors);
                open('outputTable');

                %% Parse data by relevant group (e.g. learing, nonlearning, control)

                % direct animals
                direct_x=outputTable.days(strcmpi(outputTable.animalTypes,'Direct'));  % x-values used in plotting 
                direct_y=outputTable.means(strcmpi(outputTable.animalTypes,'Direct')); % y-values used in plotting 
                direct_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Direct')); % error bar-values used in plotting 

                % non-learning animals
                indirect_x=outputTable.days(strcmpi(outputTable.animalTypes,'Indirect')); % x-values used in plotting 
                indirect_y=outputTable.means(strcmpi(outputTable.animalTypes,'Indirect')); % y-values used in plotting 
                indirect_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Indirect')); % error bar-values used in plotting

                y_global_min = min(outputTable.means - outputTable.standardErrors);
                y_global_max = max(outputTable.means + outputTable.standardErrors);
                y_midpoint = (y_global_max + y_global_min)/2;
                %% Plot figure 
                
                % generate line with shading
                close
                figure
                hold on;
                plot(direct_x, direct_y, 'Color', 'r', 'LineWidth', 4);
                plot(indirect_x, indirect_y, 'Color', 'b', 'LineWidth', 4);
                transparency=0.2; % transparency of background
                [l,p] = boundedline(direct_x, direct_y, direct_e, 'r',...
                indirect_x, indirect_y, indirect_e, 'b', 'transparency', transparency);
                
                % figure properties 
                innerlineWidth=4;  % width of mean line 
                
                % apply properties to all lines in figure 
                for lineProperty=1:length(l)
                    l(lineProperty).LineWidth=innerlineWidth;
                end
                
                %% Legend
                axis tight;
                yticks([y_global_min 0 y_global_max]);
                axis_limits = [axis_limits; sheet_name, {y_global_min}, {0}, {y_global_max}];
                lg = legend('Direct','Indirect');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(sheet_name);
                ytickformat('%.2f')
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            else
                %% Extract relevant variables
                meanCol=3; % column in spreadsheet where "Mean" text located
                measureCol=5; % column in spreadsheet where values to be plotted located

                % day
                days = repmat(([1:26]-1), [1,3])';
                % days=sort(dayCell);

                % mean values
                meanCell=cellfun(@(x) strcmpi(x,'Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                meanLogical=[meanCell{:}];
                means=[raw{[meanCell{:}],measureCol}]';

                % animal type (e.g. learning)
                animalTypes = [repmat({'Indirect'}, [26,1]); repmat({'Control'}, [26,1]); repmat({'Direct'}, [26,1])];

                % standard error
                standardErrorCell=cellfun(@(x) strcmpi(x,'Std. Error of Mean'),raw(:,meanCol),...
                'UniformOutput',false);
                standardErrors=[raw{[standardErrorCell{:}],measureCol}]';

                % compile into table for processing 
                outputTable=table(days,animalTypes,means,standardErrors);
                open('outputTable');

                %% Parse data by relevant group (e.g. learing, nonlearning, control)

                % direct animals
                direct_x=outputTable.days(strcmpi(outputTable.animalTypes,'Direct'));  % x-values used in plotting 
                direct_y=outputTable.means(strcmpi(outputTable.animalTypes,'Direct')); % y-values used in plotting 
                direct_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Direct')); % error bar-values used in plotting 

                % non-learning animals
                indirect_x=outputTable.days(strcmpi(outputTable.animalTypes,'Indirect')); % x-values used in plotting 
                indirect_y=outputTable.means(strcmpi(outputTable.animalTypes,'Indirect')); % y-values used in plotting 
                indirect_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Indirect')); % error bar-values used in plotting 

                % control animals 
                control_x=outputTable.days(strcmpi(outputTable.animalTypes,'Control')); % x-values used in plotting 
                control_y=outputTable.means(strcmpi(outputTable.animalTypes,'Control')); % y-values used in plotting 
                control_e=outputTable.standardErrors(strcmpi(outputTable.animalTypes,'Control')); % error bar-values used in plotting 

                y_global_min = min(outputTable.means - outputTable.standardErrors);
                y_global_max = max(outputTable.means + outputTable.standardErrors);
                y_midpoint = (y_global_max + y_global_min)/2;
                %% Plot figure 

                % generate line with shading
                close
                figure
                hold on;
                plot(direct_x, direct_y, 'Color', 'r', 'LineWidth', 4);
                plot(indirect_x, indirect_y, 'Color', 'b', 'LineWidth', 4);
                plot(control_x, control_y, 'Color', 'k', 'LineWidth', 4);
                transparency=0.2; % transparency of background
                [l,p] = boundedline(direct_x, direct_y, direct_e, 'r',...
                    indirect_x, indirect_y, indirect_e, 'b', control_x, control_y, control_e, 'k', ...
                    'transparency', transparency);

                % figure properties 
                innerlineWidth=4;  % width of mean line 

                % apply properties to all lines in figure 
                for lineProperty=1:length(l)
                    l(lineProperty).LineWidth=innerlineWidth;
                end

                %% Legend
                axis tight;
                yticks([y_global_min 0 y_global_max]);
                axis_limits = [axis_limits; sheet_name, {y_global_min}, {0}, {y_global_max}];
                lg = legend('Direct','Indirect','Control');
                legend('boxoff');
                lg.Location = 'BestOutside';
                lg.Orientation = 'Horizontal';

                title(sheet_name);
                ytickformat('%.2f')
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            end
        elseif contains(sheet_name, 'Correlation')
            if ~contains(sheet_name, 'Change')
                scatter(data(:,2), data(:,1), 'k', 'filled');
                hold on;
                coefficents = polyfit(data(:,2), data(:,1), 1);
                y_fit = polyval(coefficents, xlim);
                plot(xlim, y_fit, 'r');
                ytickformat('%.2f');
                title(sheet_name);
                hold off;
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));

            elseif contains(sheet_name, 'Change')
                days = (1:5:30)-1
                direct = data(1,:)
                control = data(2,:);
                indirect = data(3,:);
                plot(days, direct, 'r', 'LineWidth', 4);
                hold on
                plot(days, indirect, 'b', 'LineWidth', 4);
                plot(days, control, 'k', 'LineWidth', 4);
                y_limits = ylim;
                y_midpoint = (y_limits(1) + y_limits(2)) / 2;
                yticks([y_limits(1) y_midpoint y_limits(2)]);
                ytickformat('%.2f');
                title(sheet_name);
                hold off;
                graph_name = [sheet_name '.png'];
                saveas(gcf, fullfile(graph_path, graph_name));
            end
        else
            disp(sheet_name);
        end
    end
end