%%TODO:
%Compare dat to example to ensure that format is correct
%Fix events for days which don't have 400 total trials
%Remove low firing neurons, ex: FR < 5 spikes / second
%% Animal categories
learning = ['PRAC03', 'TNC16', 'RAVI19', 'RAVI20', 'RAVI019', 'RAVI020'];
non_learning = ['LC02', 'TNC06', 'TNC12', 'TNC25'];
control = ['TNC01', 'TNC03', 'TNC04', 'TNC14'];
right_direct = ['RAVI19', 'PRAC03', 'LC02', 'TNC12'];
left_direct = ['RAVI20', 'TNC16', 'TNC25', 'TNC06''.csv'];
%trial_range = [1,100]; %Work in progress

% Include Animal folder names in ignored_animals if you do not want to run
% them
ignored_animals = [];

% % With this you will need to select the .plx file that I send you, and it
% % will take the spike times and event times
%[file,path]=uigetfile('*.plx');

original_path = uigetdir(pwd);
animal_list = dir(original_path);
%%
if length(animal_list) > 2
    for animal = 3: length(animal_list)
        animal_name = animal_list(animal).name; %file
        animal_path = [animal_list(animal).folder, '\', animal_name, '\']; %path
        if ~isempty(ignored_animals) && contains(ignored_animals, animal_name)
            continue;
        else
            animal_name_path = dir(animal_path);
            for day = 3:length(animal_name_path)
                
                %Clear some variables to ensure they don't affect trials
                %going forward. (Probably not needed, will test later)
                clear dat spk seqTrain result event1_neural_traj event3_neural_traj event4_neural_traj event6_neural_traj estParams Spikes confusion_matrix correct1 correct2 correct3 correct4...
                    ConfusionValues
                path = animal_path;
                file = animal_name_path(day).name;
                pathfile = char(strcat(path,file));
                disp(file);
                [pathstr,name,ext] = fileparts(pathfile);
                plx = '.plx';
                check = ne(ext,plx)';
                if any(check) ==1
                    continue;
                else
                    try
                        starttime = tic;
                        
                        datafile = char(strcat(path,file));
                        [tscounts, wfcounts, evcounts, slowcounts] = plx_info(datafile,1);
                        
                        %%
                        % Fixed: Problems with PRAC 03 day 23
                        % Modified above code to be more robust, should work going forward,
                        % may need to revisit if encountering unexpected problems.
                        %
                        %
                        [nunits1, nchannels1] = size( tscounts );
                        allts = cell(nunits1, nchannels1);
                        for iunit = 0:nunits1-1   % starting with unit 0 (unsorted)
                            for ich = 1:nchannels1-1
                                if ( tscounts( iunit+1 , ich+1 ) > 0 )
                                    % get the timestamps for this channel and unit
                                    [nts, allts{iunit+1,ich}] = plx_ts(datafile, ich , iunit );
                                end
                            end
                        end
                        svStrobed=[];
                        svdummy=[];
                        % and finally the events
                        [u,nevchannels] = size( evcounts );
                        if ( nevchannels > 0 )
                            % need the event chanmap to make any sense of these
                            [u,evchans] = plx_event_chanmap(datafile);
                            for iev = 1:nevchannels
                                if ( evcounts(iev) > 0 )
                                    evch = evchans(iev);
                                    if ( evch == 257 )
                                        [nevs{iev}, tsevs{iev}, svStrobed] = plx_event_ts(datafile, evch);
                                    else
                                        [nevs{iev}, tsevs{iev}, svdummy{iev}] = plx_event_ts(datafile, evch);
                                    end
                                end
                            end
                        end
                        % %
                        events=[];
                        j=0;
                        if length(svStrobed)>1
                            events = tsevs{1,17};
                            events = [svStrobed,events];
                        else
                            for i=1:length(evcounts)
                                if evcounts(i) >= 100
                                    [nevs{i}, tsevs{i}, svdummy] = plx_event_ts(datafile, i);
                                    j=j+1;
                                    eventsingle(1:evcounts(i),1)=j;
                                    events= [events;eventsingle,tsevs{i}];
                                    eventsingle=[];
                                end
                            end
                            for i=1:length(events)
                                if events(i,1)==2
                                    events(i,1)=3;
                                elseif events(i,1)==3
                                    events(i,1)=4;
                                elseif events(i,1)==4
                                    events(i,1)=6;
                                end
                            end
                        end
                        %
                        %% Removes Doubles and Triples from events
                        i=1;
                        while i <= length(events)-1
                            if abs(events(i,2)-events(i+1,2)) < 2
                                events(i+1,:) = [];
                            end
                            i = i+1;
                        end
                        i=1;
                        while i <= length(events)-1
                            if abs(events(i,2)-events(i+1,2)) < 2
                                events(i+1,:) = [];
                            end
                            i = i+1;
                        end
                        i=1;
                        while i <= length(events)-1
                            if abs(events(i,2)-events(i+1,2)) < 2
                                events(i+1,:) = [];
                            end
                            i = i+1;
                        end
                        
                        spk = [];
                        
                        %% Choosing a Hemisphere
                        % Separate Right & Left Hemispheres
                        %Change hemisphere between 1 for Right direct or 2 for left direct
                        if contains(right_direct,animal_list(animal).name)
                            hemisphere = 1;
                        else
                            hemisphere = 2;
                        end
                        
                        if hemisphere == 1
                            % When plotting graphs label figures with specific figure name, NOT
                            % just Figure(#)
                            hemiSide = ' Right Hemisphere'; % Used to label Figure as 'Right Hemisphere'
                            for i = 1:16
                                for j = 2:5
                                    if length(allts{j,i}) >= 1 ;
                                        spk = [spk,allts(j,i)];
                                    end
                                end
                            end
                            
                        else
                            hemiSide = ' Left Hemisphere'; % Used to label Figure as 'Left Hemisphere'
                            for i = 17:32
                                for j = 2:5
                                    if length(allts{j,i}) >= 1 ;
                                        spk = [spk,allts(j,i)];
                                    end
                                end
                            end
                        end
                        
                        spiketimes=spk;
                        
                        %%
                        %Organize Events into individual variables that hold all timestamps of a
                        %single event type
                        %For angle data, right is positive, left is negative
                        
                        event1=[]; %Right Fast
                        event3=[]; %Right Slow
                        event4=[]; %Left  Fast
                        event6=[]; %Left  Slow
                        for i=1:length(events)
                            if events(i,1)==1
                                event1=[event1; events(i,2)];
                            elseif events(i,1)==3
                                event3=[event3; events(i,2)];
                            elseif events(i,1)==4
                                event4= [event4; events(i,2)];
                            else
                                event6= [event6; events(i,2)];
                            end
                        end
                        
                        newspikes=[];
                        
                        for i=1:length(spiketimes)
                            for j=1:length(spiketimes{1,i})
                                newspikes(i,j)=spiketimes{1,i}(j);
                            end
                        end
                        
                        %edge endpoints are extended to solve a problem using discretize.
                        
                        
                        %totalrelspikes is the (400 trials)x(Bins*Neurons) matrix which has each event trial for each
                        %neuron with data put into 100 bins (-0.2 : 0.2) seconds.
                        %Binned every 1 ms(see edge above)
                        %
                        % % [relspikes1]= Eventspiketimes(event1, newspikes, edge);
                        % % [relspikes3]= Eventspiketimes(event3, newspikes, edge);
                        % % [relspikes4]= Eventspiketimes(event4, newspikes, edge);
                        % % [relspikes6]= Eventspiketimes(event6, newspikes, edge);
                        
                        % [relspikes1]= event_spike_times(event1, newspikes, event_window);
                        % [relspikes3]= event_spike_times(event3, newspikes, event_window);
                        % [relspikes4]= event_spike_times(event4, newspikes, event_window);
                        % [relspikes6]= event_spike_times(event6, newspikes, event_window);
                        
                        %
                        %
                        % %This next one is important-has spikes in bins
                        %
                        % reltotalspikes = [relspikes1;relspikes3;relspikes4;relspikes6];
                        %
                        %
                        % %% Uses graphPSTH to create a single PSTH for a neuron given an event (KEEP Code - Commented out for faster running purposes)
                        % % Uses graphtestPSTH to create a PSTH which represents every trial except one, so that the left out trial can be used to test a decoder.
                        %
                        % % % neuronnum= input('Enter Neuron Number(1-#neurons):');
                        % % % graphPSTH(neuronnum,reltotalspikes,edge);
                        %
                        % %% Creates a matrix of spikes to perform pca on
                        % %sums each column (bin) for every neuron on a per event basis.
                        %
                        % count1 = sum(relspikes1);
                        % count2 = sum(relspikes3);
                        % count3 = sum(relspikes4);
                        % count4 = sum(relspikes6);
                        %
                        % %% Normalize counts
                        %
                        % [normcount1] = normcount(count1, edge);
                        % [normcount2] = normcount(count2, edge);
                        % [normcount3] = normcount(count3, edge);
                        % [normcount4] = normcount(count4, edge);
                        %
                        %
                        %
                        % %%  PCA
                        %
                        % pcacount1=[];
                        % pcacount2=[];
                        % pcacount3=[];
                        % pcacount4=[];
                        %
                        %
                        % for i=1:dimensions
                        %     pcacount1(1:(length(edge)-1),i)=normcount1((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        %     pcacount2(1:(length(edge)-1),i)=normcount2((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        %     pcacount3(1:(length(edge)-1),i)=normcount3((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        %     pcacount4(1:(length(edge)-1),i)=normcount4((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        % end
                        %
                        %
                        %
                        % pcatotalcount=[];
                        % coeff=[]; %PC by column
                        % pcascore=[];
                        % latent=[]; % eigenvalues
                        % tsquared=[]; % Hotelling's T-squared statistic for each observation
                        % explained=[]; % eigenvalue percentage of total variance
                        % mu=[]; % estimated mean of each variable
                        %
                        % %pcatotalcount is a matrix with the counts such that the bins have been
                        % %transposed so that pca can be performed
                        %
                        % pcatotalcount=[pcacount1;pcacount2;pcacount3;pcacount4];
                        % [coeff,pcascore,latent,tsquared,explained,mu]=pca(pcatotalcount);
                        %
                        % %% Plot PCA from 1 to N (KEEP Code - Commented out for faster running purposes)
                        % %Explained is percentage of total variance in each PC
                        %
                        % % %  figure;
                        % % %  plot(explained) %Figure 1: This plot will show how much variance is explained by each Principal Component (PC)
                        % % %  figure
                        % % %  for i =1:length(explained)
                        % % %      sumexplained(i) = sum(explained(1:i));
                        % % %  end
                        % % %  plot(sumexplained) %Figure 2: This plot will show what percentage (out of 100) of the variance is explained by by using each PC from 1 to x.
                        %
                        %  %% Graph PCs (KEEP Code - Commented out for faster running purposes)
                        % %Sumexplained shows a cumulative sum of the variance percentage
                        % %Choose N for how many PC to graph
                        %
                        % % % N=input('How Many PCs would you like to graph (1-# of neurons)? ');
                        % % % for i=1:N
                        % % %     figure
                        % % %     hold on
                        % % %     pca1=pcascore(1:400,i);
                        % % %     pca2=pcascore(401:800,i);
                        % % %     pca3=pcascore(801:1200,i);
                        % % %     pca4=pcascore(1201:1600,i);
                        % % %     plot(pca1);
                        % % %     plot(pca2);
                        % % %     plot(pca3);
                        % % %     plot(pca4);
                        % % %     legend('Tilt 1','Tilt 3','Tilt 4','Tilt 6');
                        % % %     title([file,' PC',num2str(i)]);
                        % % %     hold off
                        % % % end
                        %
                        % %pca(1,2,3,4) will graph the pc score relative to an event(tilt 1,3,4,6)
                        % %(1:Right Fast, 3:Right Slow, 4:Left Fast, 6:Left Slow)
                        %
                        % %Compare PCs- Plots PC1 against PC#
                        % % % for i=2:N
                        % % %     compca1=pcascore(1:400,[1,i])';
                        % % %     compca2=pcascore(401:800,[1,i])';
                        % % %     compca3=pcascore(801:1200,[1,i])';
                        % % %     compca4=pcascore(1201:1600,[1,i])';
                        % % %     figure
                        % % %     hold on
                        % % %     scatter(compca1(1,:),compca1(2,:));
                        % % %     scatter(compca2(1,:),compca2(2,:));
                        % % %     scatter(compca3(1,:),compca3(2,:));
                        % % %     scatter(compca4(1,:),compca4(2,:));
                        % % %     legend('Tilt 1','Tilt 3','Tilt 4','Tilt 6');
                        % % %     hold off
                        % % % end
                        
                        %% Choosing a Hemisphere
                        % Separate Right & Left Hemispheres
                        % hemisphere = input('Choose Right(1) or Left(2) Hemisphere: ');
                        %
                        % if hemisphere == 1
                        %     % When plotting graphs label figures with specific figure name, NOT
                        %     % just Figure(#)
                        %     hemiSide = ' Right Hemisphere '; % Used to label Figure as 'Right Hemisphere'
                        %     for i = 1:16
                        %         for j = 2:5
                        %             if length(allts{j,i}) >= 1 ;
                        %                 spk = [spk,allts(j,i)];
                        %             end
                        %         end
                        %     end
                        %
                        % else % hemisphere == 2
                        %     hemiSide = ' Left Hemisphere '; % Used to label Figure as 'Left Hemisphere'
                        %     for i = 17:32
                        %         for j = 2:5
                        %             if length(allts{j,i}) >= 1 ;
                        %                 spk = [spk,allts(j,i)];
                        %             end
                        %         end
                        %     end
                        % end
                        %
                        % % for spike = 1:length(spk)
                        % %     spk{spike} = cell2mat(spk(spike)) * 1000;
                        % % end
                        %
                        % spiketimes=spk;
                        
                        %% Bin is now 1 ms
                        %bin size is 0.001 seconds (1 ms)
                        dimensions=length(spiketimes);
                        bin = 0.001;
                        edge=0:bin:0.4;
                        total_trials = 100;
                        % classifier bin size
                        binsize=.001;  %seconds
                        bin_size = binsize;
                        % classifier window before time zero
                        pretime=0;   %seconds
                        pre_time = pretime;
                        % classifier window after time zero
                        posttime=.4;   %seconds
                        post_time = posttime;
                        
                        total_bins = (length([-abs(pre_time):bin_size:abs(post_time)]) - 1);
                        
                        neurons = spk;
                        all_events(1) = {event1};
                        all_events(2) = {event3};
                        all_events(3) = {event4};
                        all_events(4) = {event6};
                        %totalrelspikes is the (400 trials)x(Bins*Neurons) matrix which has each event trial for each
                        %neuron with data put into 100 bins (-0.2 : 0.2) seconds.
                        %Binned every 1 ms(see edge above)
                        
                        %This next one is important-has spikes in bins
                        
                        % reltotalspikes = [relspikes1;relspikes3;relspikes4;relspikes6];
                        reltotalspikes = event_spike_times(neurons, all_events, total_trials, total_bins, bin_size, pre_time, post_time);
                       
                        relspikes1 =reltotalspikes(1:100,:);
                        relspikes3 =reltotalspikes(101:200,:);
                        relspikes4 =reltotalspikes(201:300,:);
                        relspikes6 =reltotalspikes(301:400,:);
                        %% Uses graphPSTH to create a single PSTH for a neuron given an event
                        % Uses graphtestPSTH to create a PSTH which represents every trial except one, so that the left out trial can be used to test a decoder.
                        
                        % % neuronnum= input('Enter Neuron Number(1-#neurons):');
                        % % graphPSTH(neuronnum,reltotalspikes,edge);
                        
                        %% Creates a matrix of spikes to perform pca on
                        %sums each column (bin) for every neuron on a per event basis.
                        
                        %count1 = sum(relspikes1);
                        %count2 = sum(relspikes3);
                        %count3 = sum(relspikes4);
                        %count4 = sum(relspikes6);
                        
                        %% Normalize counts
                        
                        %[normcount1] = normcount(count1, edge);
                        %[normcount2] = normcount(count2, edge);
                        %[normcount3] = normcount(count3, edge);
                        %[normcount4] = normcount(count4, edge);
                        
                        % %%  PCA
                        %
                        % pcacount1=[];
                        % pcacount2=[];
                        % pcacount3=[];
                        % pcacount4=[];
                        %
                        % for i=1:dimensions
                        %     pcacount1(1:(length(edge)-1),i)=normcount1((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        %     pcacount2(1:(length(edge)-1),i)=normcount2((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        %     pcacount3(1:(length(edge)-1),i)=normcount3((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        %     pcacount4(1:(length(edge)-1),i)=normcount4((1:length(edge)-1)+((i-1)*(length(edge)-1)));
                        % end
                        %
                        % pcatotalcount=[];
                        % coeff=[]; %PC by column
                        % pcascore=[];
                        % latent=[]; % eigenvalues
                        % tsquared=[]; % Hotelling's T-squared statistic for each observation
                        % explained=[]; % eigenvalue percentage of total variance
                        % mu=[]; % estimated mean of each variable
                        %
                        % %pcatotalcount is a matrix with the counts such that the bins have been
                        % %transposed so that pca can be performed
                        %
                        % pcatotalcount=[pcacount1;pcacount2;pcacount3;pcacount4];
                        % [coeff,pcascore,latent,tsquared,explained,mu]=pca(pcatotalcount);
                        
                        %% Plot PCA from 1 to N
                        %Explained is percentage of total variance in each PC
                        
                        % %  figure;
                        % %  plot(explained) %Figure 1: This plot will show how much variance is explained by each Principal Component (PC)
                        % %  figure
                        % %  for i =1:length(explained)
                        % %      sumexplained(i) = sum(explained(1:i));
                        % %  end
                        % %  plot(sumexplained) %Figure 2: This plot will show what percentage (out of 100) of the variance is explained by by using each PC from 1 to x.
                        
                        %% Graph PCs
                        %Sumexplained shows a cumulative sum of the variance percentage
                        %Choose N for how many PC to graph
                        
                        % % N=input('How Many PCs would you like to graph (1-# of neurons)? ');
                        % % for i=1:N
                        % %     figure
                        % %     hold on
                        % %     pca1=pcascore(1:400,i);
                        % %     pca2=pcascore(401:800,i);
                        % %     pca3=pcascore(801:1200,i);
                        % %     pca4=pcascore(1201:1600,i);
                        % %     plot(pca1);
                        % %     plot(pca2);
                        % %     plot(pca3);
                        % %     plot(pca4);
                        % %     legend('Tilt 1','Tilt 3','Tilt 4','Tilt 6');
                        % %     title([file,' PC',num2str(i)]);
                        % %     hold off
                        % % end
                        
                        %pca(1,2,3,4) will graph the pc score relative to an event(tilt 1,3,4,6)
                        %(1:Right Fast, 3:Right Slow, 4:Left Fast, 6:Left Slow)
                        
                        %Compare PCs- Plots PC1 against PC#
                        % % for i=2:N
                        % %     compca1=pcascore(1:400,[1,i])';
                        % %     compca2=pcascore(401:800,[1,i])';
                        % %     compca3=pcascore(801:1200,[1,i])';
                        % %     compca4=pcascore(1201:1600,[1,i])';
                        % %     figure
                        % %     hold on
                        % %     scatter(compca1(1,:),compca1(2,:));
                        % %     scatter(compca2(1,:),compca2(2,:));
                        % %     scatter(compca3(1,:),compca3(2,:));
                        % %     scatter(compca4(1,:),compca4(2,:));
                        % %     legend('Tilt 1','Tilt 3','Tilt 4','Tilt 6');
                        % %     hold off
                        % % end
                        %%
                        %Added from other code PSTHclassifer_example_code_original.m
                        % create PSTH object
                        % [file,path]=uigetfile('*.plx'); % Keep just commented out for faster
                        
                        %% Commented out these two lines since we will be running new days, not saved .mat files
                        % [file_mat, path_mat] = uigetfile('*.mat');
                        % load([path_mat, file_mat]);
                        plxDir = char(strcat(path, "\"));
                        plxName = file;
                        
                        %% Define strobed_events
                        % if ~isempty(tsevs(17)) % This does not work
                        % if length(tsevs(17)) > 1
                        if length(svStrobed) > 1
                            strobed_events = true;
                        else
                            strobed_events = false; %Switched to true get PRAC to run correctly
                        end
                        %decoderPathold='C:\Users\Adrian & Gloria\Desktop\UC LEADS\Moxon Lab\+NeuroToolbox\+NeuroToolbox';
                        decoderPath=char(strcat(pwd,"\+NeuroToolbox\+NeuroToolbox"));
                        %% Use plx2mat3
                        % load neuron timestamps "neuronsIn" and event timestamps "eventTS"
                        if strobed_events == true
                            [neuronsIn, eventTS] = plx2mat3(plxDir, plxName,decoderPath);
                            Reference=eventTS;
                        else
                            [neuronsIn, eventTS] = plx2mat3(plxDir, plxName);
                            % remove blank event timestamps
                            Reference=eventTS(:,~cellfun(@isempty,eventTS(2,:)))';
                        end
                        
                        % warn user if size of eventTS matrix unusual
                        if size(Reference,1)~=4
                            warning('on')
                            warning(['You have ',num2str(size(eventTS,1)),...
                                ' tilt events when you should have 4'])
                            warning('Open "eventTS" variable and investigate')
                        end
                        
                        % remove neuron waveform data
                        Spikes=neuronsIn(:,[1,2]);
                        totalSpikesLength = length(Spikes);
                        
                        if hemisphere == 1
                            Spikes = Spikes(1:length(spk),:);
                        else
                            Spikes = Spikes((totalSpikesLength-length(spk)+2:end),:);
                        end
                        
                        %% Run PSTH
                        PSTH = NeuroToolbox.PSTHToolbox.PSTH(Spikes,Reference,...
                            'bin_size',binsize,'PSTH_window',[pretime,posttime],...
                            'show_progress',true);
                        
                        % create template from PSTH object
                        template = NeuroToolbox.PSTHToolbox.SU_Classifier(PSTH);
                        
                        % peform classification using template
                        DecoderOutput = template.classify(Spikes,Reference,...
                            'SameDataSet',true);
                        
                        % quantify performance
                        Correct_Trials = cellfun(@strcmp,DecoderOutput.Decision,DecoderOutput.Event);
                        Accuracy=mean(Correct_Trials)
                        
                        %% Confusion Matrix
                        %Row # is relative event #
                        %Column # is Classification of Event
                        ConfusionValues = DecoderOutput.Classification_Parameter;
                        
                        event_key = DecoderOutput.DecoderSpec.Template.TemplateSource.event_key;
                        event_decode = [];
                        
                        min_values = [];
                        
                        for i = 1:length(Correct_Trials)
                            [min_values(i), event_decode(i)] = min(ConfusionValues(i,:));
                        end
                        
                        confusion_numbers = confusionmat(event_key, event_decode');
                        correct1 = confusion_numbers(1,1);
                        correct2 = confusion_numbers(2,2);
                        correct3 = confusion_numbers(3,3);
                        correct4 = confusion_numbers(4,4);
                        confusion_matrix = (confusionmat(event_key, event_decode'))./length(Correct_Trials);
                        
                        %% Get Correct Trials
                        
                        newreltotalspikes = [];
                        for i = 1:length(Correct_Trials)
                            if Correct_Trials(i)
                                newreltotalspikes = [newreltotalspikes ; reltotalspikes(i,:)];
                            end
                        end
                        
                        
                        %% Organize Data
                        % Puts data from reltotalspikes into format for GPFA
                        
                        [rows, columns] = size(newreltotalspikes);
                        % While loop to create structure data with first field 'fieldId'
                        
                        for k = 1:rows
                            dat(k).trialId = k;
                            k = k+1;
                        end
                        
                        % While loop to create structure data with second field 'spikes'
                        timeSlots = (length(edge)-1);
                        
                        
                        % Organize reltotalspikes into an array
                        for x = 1:rows
                            dat(x).spikes = [reshape(newreltotalspikes(x,:), dimensions, timeSlots)];
                        end
                        
                        %%
                        fprintf('\n');
                        fprintf('Basic extraction and plotting of neural trajectories is complete.\n');
                        fprintf('Press any key to start cross-validation...\n');
                        fprintf('[Depending on the dataset, this can take many minutes to hours.]\n');
                        %% ===========================================
                        % 1) Basic extraction of neural trajectories
                        % ===========================================
                        %load('C:\Users\Adrian & Gloria\Desktop\UC LEADS\Moxon Lab\gpfa_v0203\gpfa_v0203\mat_sample\sample_dat');
                        %load('C:\Users\Adrian & Gloria\Desktop\UC LEADS\Moxon Lab\ourData2'); % Our version of Data
                        
                        % Results will be saved in mat_results/runXXX/, where XXX is runIdx.
                        % Use a new runIdx for each dataset.
                        
                        disp('*View Learning Animals Matlab Information.xlsx file to find Run Index Numbers*');
                        %runIdx = input('Enter Run Index Number: ');
                        filesplit = strsplit(file, '.');
                        %filerun grabs the Animal Type ('TNC', 'PRAC', etc,),
                        %animal number, and day
                        filerun = strcat(filesplit(1),filesplit(2),filesplit(4));
                        %runIDx does not need an input using filerun instead.
                        runIdx = filerun;
                        dimsToPlot1 = [1,2,3];
                        dimsToPlot2 = [1,2,3];
                        dimsToPlot3 = [1,2,3];
                        dimsToPlot4 = [1,2,3];
                        % Select method to extraclcct neural trajectories:
                        % 'gpfa' -- Gaussian-process factor analysis
                        % 'fa'   -- Smooth and factor analysis
                        % 'ppca' -- Smooth and probabilistic principal components analysis
                        % 'pca'  -- Smooth and principal components analysis
                        method = 'gpfa';
                        
                        % Select number of latent dimensions
                        xDim = 6;
                        % NOTE: The optimal dimensionality should be found using
                        %       cross-validation (Section 2) below.
                        
                        % If using a two-stage method ('fa', 'ppca', or 'pca'), select
                        % standard deviation (in msec) of Gaussian smoothing kernel.
                        kernSD = 30;
                        % NOTE: The optimal kernel width should be found using
                        %       cross-validation (Section 2) below.
                        
                        % Extract neural trajectories
                        result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,...
                            'kernSDList', kernSD);
                        % NOTE: This function does most of the heavy lifting.
                        
                        % Orthonormalize neural trajectories
                        [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
                        % NOTE: The importance of orthnormalization is described on
                        %       pp.621-622 of Yu et al., J Neurophysiol, 2009.
                        
                        % Original equation
                        for train = 1:length(seqTrain)
                            
                            for normal = 1:length(seqTrain)
                                seqTrain(normal).xorth = (seqTrain(normal).xorth)./ (max(abs(seqTrain(normal).xorth),[],2));
                            end
                        end
                        
                        % Plot neural trajectories in 3D space
                        
                        event1_neural_traj = [];
                        event3_neural_traj = [];
                        event4_neural_traj = [];
                        event6_neural_traj = [];
                        
                        event1_neural_traj = seqTrain(1:correct1); % Event 1 correct trials
                        
                        event3_neural_traj = seqTrain((correct1 + 1):(correct1 + correct2)); % Event 3 correct trials
                        
                        event4_neural_traj = seqTrain((correct1 + correct2 + 1):(correct1 + correct2 + correct3)); % Event 4 correct trials
                        
                        event6_neural_traj = seqTrain((correct1 + correct2 + correct3 + 1):(correct1 + correct2 + correct3 + correct4)); % Event 6 correct trials
                        
                        
                        % Plot & Save 3D graphs
                        % plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3);
                        
                        % Event 1
                        
                        plot3D(event1_neural_traj, 'xorth', 'dimsToPlot', dimsToPlot1);
                        legend({'Neural Trajectory', 'Start of Tilt', 'Decision Made', 'End of Tilt'},'TextColor', 'black');
                        title(legend, 'Event 1');
                        % COPY the below code & PASTE into whatever EVENT# you need to find clusters for
                        hold on
                        
                        [startTiltMeanFactor1Event1, startTiltMeanFactor2Event1, startTiltMeanFactor3Event1, startTiltSDFactor1Event1, startTiltSDFactor2Event1, startTiltSDFactor3Event1, decisionMadeMeanFactor1Event1, decisionMadeMeanFactor2Event1, decisionMadeMeanFactor3Event1, decisionMadeSDFactor1Event1, decisionMadeSDFactor2Event1, decisionMadeSDFactor3Event1, endTiltMeanFactor1Event1, endTiltMeanFactor2Event1, endTiltMeanFactor3Event1, endTiltSDFactor1Event1, endTiltSDFactor2Event1, endTiltSDFactor3Event1] = ellipse_mean_sd(dimsToPlot1, event1_neural_traj, xDim);
                        [decisionClusterX, decisionClusterY, decisionClusterZ] = ellipsoid(decisionMadeMeanFactor1Event1, decisionMadeMeanFactor2Event1, decisionMadeMeanFactor3Event1, decisionMadeSDFactor1Event1, decisionMadeSDFactor2Event1, decisionMadeSDFactor3Event1);
                        surf(decisionClusterX, decisionClusterY, decisionClusterZ, 'DisplayName','Decision Made',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.756862759590149 0.866666674613953 0.776470601558685]);
                        hold off
                        
                        hold on
                        
                        [endClusterX,endClusterY, endClusterZ] = ellipsoid(endTiltMeanFactor1Event1, endTiltMeanFactor2Event1, endTiltMeanFactor3Event1, endTiltSDFactor1Event1, endTiltSDFactor2Event1, endTiltSDFactor3Event1);
                        surf(endClusterX, endClusterY, endClusterZ, 'DisplayName','End Tilt',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.729411780834198 0.831372559070587 0.95686274766922]);
                        hold off
                        % End copying of code
                        msg = '\nUse single or double quotes. \nWrite Description. \nExample: TNC(#)_Day(#)\nFigureName: ';
                        figureName = filerun;
                        newFigureName = char(strcat(figureName, ' Dimensions-',num2str(dimsToPlot1), ', Event 1, ', hemiSide, '-',  method));
                        set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
                        savefig(gcf,char(newFigureName));
                        
                        % Event 3
                        
                        plot3D(event3_neural_traj, 'xorth', 'dimsToPlot', dimsToPlot2);
                        legend({'Neural Trajectory', 'Start of Tilt', 'Decision Made', 'End of Tilt'},'TextColor', 'black');
                        title(legend, 'Event 3');
                        
                        % COPY the below code & PASTE into whatever EVENT# you need to find clusters for
                        hold on
                        [startTiltMeanFactor1Event3, startTiltMeanFactor2Event3, startTiltMeanFactor3Event3, startTiltSDFactor1Event3, startTiltSDFactor2Event3, startTiltSDFactor3Event3, decisionMadeMeanFactor1Event3, decisionMadeMeanFactor2Event3, decisionMadeMeanFactor3Event3, decisionMadeSDFactor1Event3, decisionMadeSDFactor2Event3, decisionMadeSDFactor3Event3, endTiltMeanFactor1Event3, endTiltMeanFactor2Event3, endTiltMeanFactor3Event3, endTiltSDFactor1Event3, endTiltSDFactor2Event3, endTiltSDFactor3Event3] = ellipse_mean_sd(dimsToPlot2, event3_neural_traj, xDim);
                        [decisionClusterX, decisionClusterY, decisionClusterZ] = ellipsoid(decisionMadeMeanFactor1Event3, decisionMadeMeanFactor2Event3, decisionMadeMeanFactor3Event3, decisionMadeSDFactor1Event3, decisionMadeSDFactor2Event3, decisionMadeSDFactor3Event3);
                        surf(decisionClusterX, decisionClusterY, decisionClusterZ, 'DisplayName','Decision Made',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.756862759590149 0.866666674613953 0.776470601558685]);
                        hold off
                        
                        hold on
                        [endClusterX,endClusterY, endClusterZ] = ellipsoid(endTiltMeanFactor1Event3, endTiltMeanFactor2Event3, endTiltMeanFactor3Event3, endTiltSDFactor1Event3, endTiltSDFactor2Event3, endTiltSDFactor3Event3);
                        surf(endClusterX, endClusterY, endClusterZ, 'DisplayName','End Tilt',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.729411780834198 0.831372559070587 0.95686274766922]);
                        hold off
                        % End copying of code
                        
                        newFigureName = char(strcat(figureName, ' Dimensions-',num2str(dimsToPlot2), ', Event 3, ', hemiSide, '-',  method));
                        set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
                        savefig(gcf,char(newFigureName));
                        
                        % Event 4
                        
                        plot3D(event4_neural_traj, 'xorth', 'dimsToPlot', dimsToPlot3);
                        legend({'Neural Trajectory', 'Start of Tilt', 'Decision Made', 'End of Tilt'},'TextColor', 'black');
                        title(legend, 'Event 4');
                        
                        % COPY the below code & PASTE into whatever EVENT# you need to find clusters for
                        hold on
                        [startTiltMeanFactor1Event4, startTiltMeanFactor2Event4, startTiltMeanFactor3Event4, startTiltSDFactor1Event4, startTiltSDFactor2Event4, startTiltSDFactor3Event4, decisionMadeMeanFactor1Event4, decisionMadeMeanFactor2Event4, decisionMadeMeanFactor3Event4, decisionMadeSDFactor1Event4, decisionMadeSDFactor2Event4, decisionMadeSDFactor3Event4, endTiltMeanFactor1Event4, endTiltMeanFactor2Event4, endTiltMeanFactor3Event4, endTiltSDFactor1Event4, endTiltSDFactor2Event4, endTiltSDFactor3Event4] = ellipse_mean_sd(dimsToPlot3, event4_neural_traj, xDim);
                        [decisionClusterX, decisionClusterY, decisionClusterZ] = ellipsoid(decisionMadeMeanFactor1Event4, decisionMadeMeanFactor2Event4, decisionMadeMeanFactor3Event4, decisionMadeSDFactor1Event4, decisionMadeSDFactor2Event4, decisionMadeSDFactor3Event4);
                        surf(decisionClusterX, decisionClusterY, decisionClusterZ, 'DisplayName','Decision Made',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.756862759590149 0.866666674613953 0.776470601558685]);
                        hold off
                        
                        hold on
                        [endClusterX,endClusterY, endClusterZ] = ellipsoid(endTiltMeanFactor1Event4, endTiltMeanFactor2Event4, endTiltMeanFactor3Event4, endTiltSDFactor1Event4, endTiltSDFactor2Event4, endTiltSDFactor3Event4);
                        surf(endClusterX, endClusterY, endClusterZ, 'DisplayName','End Tilt',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.729411780834198 0.831372559070587 0.95686274766922]);
                        hold off
                        % End copying of code
                        
                        newFigureName = char(strcat(figureName, ' Dimensions-',num2str(dimsToPlot3), ', Event 4, ', hemiSide, '-',  method));
                        set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
                        savefig(gcf,char(newFigureName));
                        
                        % Event 6
                        
                        plot3D(event6_neural_traj, 'xorth', 'dimsToPlot', dimsToPlot4);
                        legend({'Neural Trajectory', 'Start of Tilt', 'Decision Made', 'End of Tilt'},'TextColor', 'black');
                        title(legend, 'Event 6');
                        
                        % COPY the below code & PASTE into whatever EVENT# you need to find clusters for
                        hold on
                        [startTiltMeanFactor1Event6, startTiltMeanFactor2Event6, startTiltMeanFactor3Event6, startTiltSDFactor1Event6, startTiltSDFactor2Event6, startTiltSDFactor3Event6, decisionMadeMeanFactor1Event6, decisionMadeMeanFactor2Event6, decisionMadeMeanFactor3Event6, decisionMadeSDFactor1Event6, decisionMadeSDFactor2Event6, decisionMadeSDFactor3Event6, endTiltMeanFactor1Event6, endTiltMeanFactor2Event6, endTiltMeanFactor3Event6, endTiltSDFactor1Event6, endTiltSDFactor2Event6, endTiltSDFactor3Event6] = ellipse_mean_sd(dimsToPlot4, event6_neural_traj, xDim);
                        [decisionClusterX, decisionClusterY, decisionClusterZ] = ellipsoid(decisionMadeMeanFactor1Event6, decisionMadeMeanFactor2Event6, decisionMadeMeanFactor3Event6, decisionMadeSDFactor1Event6, decisionMadeSDFactor2Event6, decisionMadeSDFactor3Event6);
                        surf(decisionClusterX, decisionClusterY, decisionClusterZ, 'DisplayName','Decision Made',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.756862759590149 0.866666674613953 0.776470601558685]);
                        hold off
                        
                        hold on
                        [endClusterX,endClusterY, endClusterZ] = ellipsoid(endTiltMeanFactor1Event6, endTiltMeanFactor2Event6, endTiltMeanFactor3Event6, endTiltSDFactor1Event6, endTiltSDFactor2Event6, endTiltSDFactor3Event6);
                        surf(endClusterX, endClusterY, endClusterZ, 'DisplayName','End Tilt',...
                            'FaceAlpha',0.5,...
                            'LineStyle','none',...
                            'FaceColor',[0.729411780834198 0.831372559070587 0.95686274766922]);
                        hold off
                        % End copying of code
                        
                        newFigureName = char(strcat(figureName, ' Dimensions-',num2str(dimsToPlot4), ', Event 6, ', hemiSide, '-', method));
                        set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
                        savefig(gcf,char(newFigureName));
                        
                        ellipsoid_data = [startTiltMeanFactor1Event1, startTiltMeanFactor2Event1, startTiltMeanFactor3Event1, decisionMadeMeanFactor1Event1, decisionMadeMeanFactor2Event1, decisionMadeMeanFactor3Event1, endTiltMeanFactor1Event1, endTiltMeanFactor2Event1, endTiltMeanFactor3Event1, startTiltSDFactor1Event1, startTiltSDFactor2Event1, startTiltSDFactor3Event1, decisionMadeSDFactor1Event1, decisionMadeSDFactor2Event1, decisionMadeSDFactor3Event1, endTiltSDFactor1Event1, endTiltSDFactor2Event1, endTiltSDFactor3Event1;...
                            startTiltMeanFactor1Event3, startTiltMeanFactor2Event3, startTiltMeanFactor3Event3, decisionMadeMeanFactor1Event3, decisionMadeMeanFactor2Event3, decisionMadeMeanFactor3Event3, endTiltMeanFactor1Event3, endTiltMeanFactor2Event3, endTiltMeanFactor3Event3, startTiltSDFactor1Event3, startTiltSDFactor2Event3, startTiltSDFactor3Event3, decisionMadeSDFactor1Event3, decisionMadeSDFactor2Event3, decisionMadeSDFactor3Event3, endTiltSDFactor1Event3, endTiltSDFactor2Event3, endTiltSDFactor3Event3;...
                            startTiltMeanFactor1Event4, startTiltMeanFactor2Event4, startTiltMeanFactor3Event4, decisionMadeMeanFactor1Event4, decisionMadeMeanFactor2Event4, decisionMadeMeanFactor3Event4, endTiltMeanFactor1Event4, endTiltMeanFactor2Event4, endTiltMeanFactor3Event4, startTiltSDFactor1Event4, startTiltSDFactor2Event4, startTiltSDFactor3Event4, decisionMadeSDFactor1Event4, decisionMadeSDFactor2Event4, decisionMadeSDFactor3Event4, endTiltSDFactor1Event4, endTiltSDFactor2Event4, endTiltSDFactor3Event4;...
                            startTiltMeanFactor1Event6, startTiltMeanFactor2Event6, startTiltMeanFactor3Event6, decisionMadeMeanFactor1Event6, decisionMadeMeanFactor2Event6, decisionMadeMeanFactor3Event6, endTiltMeanFactor1Event6, endTiltMeanFactor2Event6, endTiltMeanFactor3Event6, startTiltSDFactor1Event6, startTiltSDFactor2Event6, startTiltSDFactor3Event6, decisionMadeSDFactor1Event6, decisionMadeSDFactor2Event6, decisionMadeSDFactor3Event6, endTiltSDFactor1Event6, endTiltSDFactor2Event6, endTiltSDFactor3Event6];
                        newfile = strcat(file,".mat");
                        save(newfile,'dat','ellipsoid_data','Accuracy','confusion_numbers');
                        toc(starttime)
                        close all
                    catch ME
                        file_name = animal_name_path(day).name;
                        failed_calculating = {};
                        failed_calculating{end + 1} = file_name;
                        failed_calculating{end, 2} = ME;
                        filename = ['FAILED.', file_name, '.mat'];
                        warning('%s failed to calculate\n', file_name);
                        warning('Error: %s\n', ME.message);
                        matfile = fullfile(original_path, filename);
                        save(matfile, 'failed_calculating');
                    end
                end
            end
        end
    end
end
disp('done')
%cleanup
% close
% clear
% clc

% NOTES:
% - This figure shows the time-evolution of neural population
%   activity on a single-trial basis.  Each trajectory is extracted from
%   the activity of all units on a single trial.
% - This particular example is based on multi-electrode recordings
%   in premotor and motor cortices within a 400 ms period starting 300 ms
%   before movement onset.  The extracted trajectories appear to
%   follow the same general path, but there are clear trial-to-trial
%   differences that can be related to the physical arm movement.
% - Analogous to Figure 8 in Yu et al., J Neurophysiol, 2009.
% WARNING:
% - If the optimal dimensionality (as assessed by cross-validation in
%   Section 2) is greater than 3, then this plot may mask important
%   features of the neural trajectories in the dimensions not plotted.
%   This motivates looking at the next plot, which shows all latent
%   dimensions.

% %% Plot Dimension of Neural Trajectories vs. Time
% % Plot each dimension of neural trajectories versus time
%
% % Event 1
% plotEachDimVsTime(event1_neural_traj, 'xorth', result.binWidth);
% % msg = '\nUse single or double quotes. \nWrite Description. \nExample: TNC(#)_Day(#)\nFigureName: ';
% newFigureName = strcat(figureName, ', Event 1, ', ' Dimension Vs Time ', hemiSide, '-',  method);
% set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
% % set(gcf, 'Name', newFigureName, ...
% %          'NumberTitle', 'off', ...
% %          'Units', 'Pixels',...
% %          'PaperSize', [5 5] ...
% %      );
% %
% % new_figure_name = strcat(newFigureName,'.pdf');
% % saveas(gcf,new_figure_name)
% savefig(newFigureName);
%
% % Event 3
% plotEachDimVsTime(event3_neural_traj, 'xorth', result.binWidth);
% newFigureName = strcat(figureName, ', Event 3, ', ' Dimension Vs Time ',  hemiSide, '-',  method);
% set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
% savefig(newFigureName);
%
% % Event 4
% plotEachDimVsTime(event4_neural_traj, 'xorth', result.binWidth);
% newFigureName = strcat(figureName, ', Event 4, ', ' Dimension Vs Time ',  hemiSide, '-',  method);
% set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
% savefig(newFigureName);
%
% % Event 6
% plotEachDimVsTime(event6_neural_traj, 'xorth', result.binWidth);
% newFigureName = strcat(figureName, ', Event 6, ', ' Dimension Vs Time ',  hemiSide, '-',  method);
% set(gcf, 'Name', newFigureName, 'NumberTitle', 'off');
% savefig(newFigureName);
%
% % NOTES:
% % - These are the same neural trajectories as in the previous figure.
% %   The advantage of this figure is that we can see all latent
% %   dimensions (one per panel), not just three selected dimensions.
% %   As with the previous figure, each trajectory is extracted from the
% %   population activity on a single trial.  The activity of each unit
% %   is some linear combination of each of the panels.  The panels are
% %   ordered, starting with the dimension of greatest covariance
% %   (in the case of 'gpfa' and 'fa') or variance (in the case of
% %   'ppca' and 'pca').
% % - From this figure, we can roughly estimate the optimal
% %   dimensionality by counting the number of top dimensions that have
% %   'meaningful' temporal structure.   In this example, the optimal
% %   dimensionality appears to be about 5.  This can be assessed
% %   quantitatively using cross-validation in Section 2.
% % - Analogous to Figure 7 in Yu et al., J Neurophysiol, 2009.
%
% fprintf('\n');
% fprintf('Basic extraction and plotting of neural trajectories is complete.\n');
% fprintf('Press any key to start cross-validation...\n');
% fprintf('[Depending on the dataset, this can take many minutes to hours.]\n');
% %pause;
%
%
% %% ========================================================
% % 2) Full cross-validation to find:
% %  - optimal state dimensionality for all methods
% %  - optimal smoothing kernel width for two-stage methods
% % ========================================================
% %
% % Select number of cross-validation folds
% numFolds = 4;
% %
% % Perform cross-validation for different state dimensionalities.
% % Results are saved in mat_results/runXXX/, where XXX is runIdx.
% for xDim = [2 5 8]
%   neuralTraj(runIdx, dat, 'method',  'pca', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(runIdx, dat, 'method', 'ppca', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(runIdx, dat, 'method',   'fa', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
% end
% fprintf('\n');
% % NOTES:
% % - These function calls are computationally demanding.  Cross-validation
% %   takes a long time because a separate model has to be fit for each
% %   state dimensionality and each cross-validation fold.
% %
% % Plot prediction error versus state dimensionality.
% % Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
% % kernSD = 30; % select kernSD for two-stage methods
% plotPredErrorVsDim(runIdx, kernSD);
% % NOTES:
% % - Using this figure, we i) compare the performance (i.e,,
% %   predictive ability) of different methods for extracting neural
% %   trajectories, and ii) find the optimal latent dimensionality for
% %   each method.  The optimal dimensionality is that which gives the
% %   lowest prediction error.  For the two-stage methods, the latent
% %   dimensionality and smoothing kernel width must be jointly
% %   optimized, which requires looking at the next figure.
% % - In this particular example, the optimal dimensionality is 5. This
% %   implies that, even though the raw data are evolving in a
% %   53-dimensional space (i.e., there are 53 units), the system
% %   appears to be using only 5 degrees of freedom due to firing rate
% %   correlations across the neural population.
% % - Analogous to Figure 5A in Yu et al., J Neurophysiol, 2009.
% %
% % Plot prediction error versus kernelSD.
% % Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
% xDim = 5; % select state dimensionality
% plotPredErrorVsKernSD(runIdx, xDim);
% % NOTES:
% % - This figure is used to find the optimal smoothing kernel for the
% %   two-stage methods.  The same smoothing kernel is used for all units.
% % - In this particular example, the optimal standard deviation of a
% %   Gaussian smoothing kernel with FA is 30 ms.
% % - Analogous to Figures 5B and 5C in Yu et al., J Neurophysiol, 2009.

%% Functions
%% Find Mean & Standard Deviation Start Tilt
% EDIT code for whatever event is needed

function [startTiltMeanFactor1, startTiltMeanFactor2, startTiltMeanFactor3, startTiltSDFactor1, startTiltSDFactor2, startTiltSDFactor3, decisionMadeMeanFactor1, decisionMadeMeanFactor2, decisionMadeMeanFactor3, decisionMadeSDFactor1, decisionMadeSDFactor2, decisionMadeSDFactor3, endTiltMeanFactor1, endTiltMeanFactor2, endTiltMeanFactor3, endTiltSDFactor1, endTiltSDFactor2, endTiltSDFactor3] = ellipse_mean_sd(dimsToPlot, event3_neural_traj, xDim)

for newTrainEvent3Factor1 = 1:length(event3_neural_traj)
    startTiltEvent3Factor1(newTrainEvent3Factor1) = event3_neural_traj(newTrainEvent3Factor1).xorth(dimsToPlot(1),1);
end

for newTrainEvent3Factor2 = 1:length(event3_neural_traj)
    startTiltEvent3Factor2(newTrainEvent3Factor2) = event3_neural_traj(newTrainEvent3Factor2).xorth(dimsToPlot(2),1);
end

for newTrainEvent3Factor3 = 1:length(event3_neural_traj)
    startTiltEvent3Factor3(newTrainEvent3Factor3) = event3_neural_traj(newTrainEvent3Factor3).xorth(dimsToPlot(3),1);
end

startTiltMeanFactor1 = mean(startTiltEvent3Factor1);
startTiltSDFactor1 = std(startTiltEvent3Factor1);

startTiltMeanFactor2 = mean(startTiltEvent3Factor2);
startTiltSDFactor2 = std(startTiltEvent3Factor2);

startTiltMeanFactor3 = mean(startTiltEvent3Factor3);
startTiltSDFactor3 = std(startTiltEvent3Factor3);

%% Find Mean & Standard Deviation Decision Made
% EDIT code for whatever event is needed
for newTrainEvent3Factor1 = 1:length(event3_neural_traj)
    decisionMadeEvent3Factor1(newTrainEvent3Factor1) = event3_neural_traj(newTrainEvent3Factor1).xorth(dimsToPlot(1),(end/2));
end

for newTrainEvent3Factor2 = 1:length(event3_neural_traj)
    decisionMadeEvent3Factor2(newTrainEvent3Factor2) = event3_neural_traj(newTrainEvent3Factor2).xorth(dimsToPlot(2),(end/2));
end

for newTrainEvent3Factor3 = 1:length(event3_neural_traj)
    decisionMadeEvent3Factor3(newTrainEvent3Factor3) = event3_neural_traj(newTrainEvent3Factor3).xorth(dimsToPlot(3),(end/2));
end

decisionMadeMeanFactor1 = mean(decisionMadeEvent3Factor1);
decisionMadeSDFactor1 = std(decisionMadeEvent3Factor1);

decisionMadeMeanFactor2 = mean(decisionMadeEvent3Factor2);
decisionMadeSDFactor2 = std(decisionMadeEvent3Factor2);

decisionMadeMeanFactor3 = mean(decisionMadeEvent3Factor3);
decisionMadeSDFactor3 = std(decisionMadeEvent3Factor3);


%% Find Mean & Standard Deviation End of Tilt
% EDIT code for whatever event is needed
for newTrainEvent3Factor1 = 1:length(event3_neural_traj)
    endTiltEvent3Factor1(newTrainEvent3Factor1) = event3_neural_traj(newTrainEvent3Factor1).xorth(dimsToPlot(1),end);
end

for newTrainEvent3Factor2 = 1:length(event3_neural_traj)
    endTiltEvent3Factor2(newTrainEvent3Factor2) = event3_neural_traj(newTrainEvent3Factor2).xorth(dimsToPlot(2),end);
end

for newTrainEvent3Factor3 = 1:length(event3_neural_traj)
    endTiltEvent3Factor3(newTrainEvent3Factor3) = event3_neural_traj(newTrainEvent3Factor3).xorth(dimsToPlot(2),end);
end

endTiltMeanFactor1 = mean(endTiltEvent3Factor1);
endTiltSDFactor1 = std(endTiltEvent3Factor1);

endTiltMeanFactor2 = mean(endTiltEvent3Factor2);
endTiltSDFactor2 = std(endTiltEvent3Factor2);

endTiltMeanFactor3 = mean(endTiltEvent3Factor3);
endTiltSDFactor3 = std(endTiltEvent3Factor3);
end