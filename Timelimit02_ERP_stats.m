%% Compute STATS on time-series within/between subject
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: re_created on July 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: 

    OUTPUT: 

    FIXME: **

%}
%=========================================================================%
%% START of the script 

%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear workspace (if needed)
if input('clear all?  (1/0) ... ')
    clearvars; close all;
end

% set paths (if needed)
BT_setpath

% choose subj & go to the right folder
BT_getsubj

clear LevelAnalysis name numlines prompt subj_folders 

%% More specific paths (maybe set this in the start script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

behavioral_folder= [results_Path, '/Behaviour']; % it can be also current_subj_folder
if ~exist(fullfile(behavioral_folder)); mkdir(fullfile(behavioral_folder)); end;

timeseries_folder= [results_Path, '/Timeseries']; % it can be also current_subj_folder
if ~exist(fullfile(timeseries_folder)); mkdir(fullfile(timeseries_folder)); end;
cd(timeseries_folder);

powerspectra_folder= [results_Path, '/Powerspect']; % it can be also current_subj_folder
    if ~exist(fullfile(powerspectra_folder)); mkdir(fullfile(powerspectra_folder)); end;
    cd(powerspectra_folder);
    
correlation_folder= [results_Path, '/Correlations']; % it can be also current_subj_folder
    if ~exist(fullfile(correlation_folder)); mkdir(fullfile(correlation_folder)); end;  
    
regression_folder= [results_Path, '/Regressions']; % it can be also current_subj_folder
    if ~exist(fullfile(regression_folder)); mkdir(fullfile(regression_folder)); end;
    
statistics_folder= [results_Path, '/Statistics']; % it can be also current_subj_folder
    if ~exist(fullfile(statistics_folder)); mkdir(fullfile(statistics_folder)); end;

%% 0) Other tests tried
    
%% 1a) CORRELATION between RP amplitudes (Y) vs conditions (X) before log transformation (raw data)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Set parameters and right path:
% 
% nChans= 60;
% nTimes= 2001; % corresponds to [-3 +1]s
% R= zeros(nChans,nTimes); P= zeros(nChans,nTimes);
% 
% clear data_Path; 
% % set new data_Path
% data_Path= timeseries_folder; % need to fix it in the start up script
% 
% 
% % Load in useful data
% 
%     % Y= Time-locked amplitudes (alias RP avg)
%     cd(timeseries_folder);
%     for subi=1:nSubjs;
% 
%         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
%         pickupER(subi) = load(fname_ER);
% 
%     end
%     
% %      save 'pickupER' pickupER;
% 
%     % X= Conditions (alias Time limits)
%     cd(powerspectra_folder);
%     for subi=1:nSubjs;
%         
%         fname_Cond= sprintf('subj%02d_usefulinfo',subi);
%         pickupCond(subi) = load(fname_Cond);
%         
%     end
%     
% %      save 'pickupCond' pickupCond;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Compute correlations here:
% 
% cd(regression_folder);
% 
%         for subi=1:nSubjs;
% 
%             % Just re-converting Inf values into an integer number
%             pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
%             
%             X= [pickupCond(subi).newcond' ones(length(pickupCond(subi).newcond'),1)];
%             Y= pickupER(subi).avg_trl.trial; %.avg or .trial?
%             if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;
% 
%             profile on
% 
%             R(:)= 0; P(:)=0;
% 
%             for i=1:nChans
% 
%                     for j= 1:nTimes
% 
%                        % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
%                        [tmpR,tmpP] = corrcoef(Y(:,i,j),X(:,1));
%                        R(i,j) = tmpR(1,2); P(i,j) = tmpP(1,2);
% 
%                     end
% 
% 
%                 disp(['We are at channel' num2str(i)]);
% 
%             end
% 
%             
%             filename= [sprintf('subj%02d_Corr_ER1', subi)]; % add one if all trials mixed by condition
%             save(filename,'R','P','-v7.3');
% 
%             profile off;
%             profile viewer;
%             disp(['Subject ' num2str(subi) ' done']);
% 
%         end
% 
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
% % Compute Grandaverages of correlations here:
% 
%             cd(timeseries_folder);
%             
%             load subj22_TimeS_one; % all conditions in one, trials averaged
%             templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
%             
%             % make a matrix with all the subjects
%             cd(regression_folder);
% 
%             R_all= {};
%             for subi=1:nSubjs;
%                 
%                 fname_Corr= sprintf('subj%02d_Corr_ER1',subi);
%                 pickupCorrs(subi) = load(fname_Corr);
%                 
%                 templateER.avg = squeeze(pickupCorrs(subi).R); %squeeze to collapse trials dimension
%                 R_all{subi}= templateER;
%                 
%             end
%             
% 
%             cd(regression_folder);
%             save('R_all_ER1', 'R_all','templateER');
% 
% 
%             GAVG_R= ft_timelockgrandaverage([], R_all{:});
%             
% % Plot correlations:
% 
%             cfg=[];
%             % cfg.preproc.lpfilter='yes';
%             % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
%             cfg.layout = 'eeg_64_NM20884N.lay';
%             cfg.linewidth = 1.5;
%             cfg.showlabels= 'yes';
%             figure
%             ft_multiplotER(cfg,GAVG_R);
            
%% 1b) CORRELATION between RP amplitudes (Y) vs log-transformed conditions (X) 

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
R= zeros(nChans,nTimes); P= zeros(nChans,nTimes); 

% clear data_Path; 
% % set new data_Path
% data_Path= timeseries_folder; % need to fix it in the start up script


% Load in useful data

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    load('pickupER');

    for subi=1:nSubjs;

        fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
        pickupER(subi) = load(fname_ER);

    end    
%     save('pickupER','pickupER','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 


    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    load('pickupCond');
    
%     for subi=1:nSubjs;
%         
%         fname_Cond= sprintf('subj%02d_usefulinfo',subi);
%         pickupCond(subi) = load(fname_Cond);
%         
%     end
%      save 'pickupCond' pickupCond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlations here:
correlation_folder2= [correlation_folder, '/CorrelationsContr']; % it can be also current_subj_folder
    if ~exist(fullfile(correlation_folder2)); mkdir(fullfile(correlation_folder2)); end;  
    
cd(correlation_folder2);

        for subi=1:nSubjs;

            % Just re-converting Inf values into an integer number
            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
            pickupCond(subi).Logcond= log(pickupCond(subi).newcond);
             
            X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)];
            Y= pickupER(subi).avg_trl.trial; %.avg or .trial? /// pickup_LogER
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

            R(:)= 0; P(:)=0;

            for i=1:nChans

                    for j= 1:nTimes

                       % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                       [tmpR,tmpP] = corrcoef(Y(:,i,j),X(:,1));
                       R(i,j) = tmpR(1,2); P(i,j) = tmpP(1,2);

                    end


                disp(['We are at channel' num2str(i)]);

            end

            
            filename= [sprintf('subj%02d_CONTRcorr_ER1Log_semiLogX', subi)]; % add one if all trials mixed by condition
            save(filename,'R','P','-v7.3');

%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Compute Grandaverages of correlations here:

            cd(timeseries_folder);
            
            cd(contrTS_folder);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(correlation_folder2);

            R_all= {};
            for subi=1:nSubjs;
                
                fname_Corr= sprintf('subj%02d_CONTRcorr_ER1Log_semiLogX',subi);
                pickupCorrs(subi) = load(fname_Corr);
                
                templateER.avg = squeeze(pickupCorrs(subi).R); %squeeze to collapse trials dimension
                R_all{subi}= templateER;
                
            end
            

            cd(correlation_folder2);
            save('R_all_ER1Contr_semiLogX', 'R_all','templateER','-v7.3');


            GAVG_R= ft_timelockgrandaverage([], R_all{:});
            
% Plot correlations:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.showlabels= 'yes';
            figure
            ft_multiplotER(cfg,GAVG_R);
            
%% 2a) CORRELATION between RP amplitudes (Y) vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Set parameters and right path:
% 
% nChans= 60;
% nTimes= 2001; % corresponds to [-3 +1]s
% R= zeros(nChans,nTimes); P= zeros(nChans,nTimes);
% 
% clear data_Path; %set right path
% data_Path= timeseries_folder; % need to fix it in the start up script
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Load in useful data:
% 
%     % Y= Time-locked amplitudes (alias RP avg)
%     cd(timeseries_folder);
%     for subi=1:nSubjs;
% 
%         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
%         pickupER(subi) = load(fname_ER);
% 
%     end
%     
%     % X= Response times (alias Waiting Times)
%     cd(behavioral_folder);
%     for subi=1:nSubjs;
% 
%         fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
%         pickupBehav(subi) = load(fname_BehavData);
% 
%     end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Compute correlation here:
% 
% cd(correlation_folder2);
% 
%         for subi=1:nSubjs;
% 
%             X= [pickupBehav(subi).RESPTIMES' ones(length(pickupBehav(subi).RESPTIMES'),1)];
%             Y= pickupER(subi).avg_trl.trial; %.avg or .trial?
%             if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;
% 
%             profile on
% 
%             R(:)= 0; P(:)=0;
% 
%             for i=1:nChans
% 
%                     for j= 1:nTimes
% 
%                        % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
%                        [tmpR,tmpP] = corrcoef(Y(:,i,j),X(:,1));
%                        R(i,j) = tmpR(1,2); P(i,j) = tmpP(1,2);
% 
%                     end
% 
% 
%                 disp(['We are at channel' num2str(i)]);
% 
%             end
% 
% 
%             filename= [sprintf('subj%02d_Corr_ER2', subi)]; % add one if all trials mixed by condition
%             save(filename,'R','P','-v7.3');
% 
%             profile off;
%             profile viewer;
%             disp(['Subject ' num2str(subi) ' done']);
% 
%         end
% 
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
% % Compute Grandaverages of correlations here:
% 
%             cd(data_Path);
%             
%             load subj22_TimeS_one; % all conditions in one, trials averaged
%             templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
%             
%             % make a matrix with all the subjects
%             cd(regression_folder);
% 
%             R_all= {};
%             for subi=1:nSubjs;
% 
%                 fname_Corr= sprintf('subj%02d_Corr_ER2',subi);
%                 pickupCorrs(subi) = load(fname_Corr);
% 
%                 templateER.avg = squeeze(pickupCorrs(subi).R); %squeeze to collapse trials dimension
%                 R_all{subi}= templateER;
% 
%             end
%             
%             
%             cd(regression_folder);
%             save('R_all_ER2', 'R_all','templateER');
% 
% 
%             GAVG_R= ft_timelockgrandaverage([], R_all{:});
%             
% % Plot correlations:
% 
%             cfg=[];
%             % cfg.preproc.lpfilter='yes';
%             % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
%             cfg.layout = 'eeg_64_NM20884N.lay';
%             cfg.linewidth = 1.5;
%             cfg.showlabels= 'yes';
%             figure
%             ft_multiplotER(cfg,GAVG_R);

%% 2b) CORRELATION between RP amplitudes (Y) vs log-transformed response times (X)

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
R= zeros(nChans,nTimes); P= zeros(nChans,nTimes); 

% clear data_Path; 
% % set new data_Path
% data_Path= timeseries_folder; % need to fix it in the start up script


% Load in useful data

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
%     load('pickup_LogER');

%     for subi=1:nSubjs;
% 
%         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
%         pickupER(subi) = load(fname_ER);
% 
%     end    
%     save('pickupER','pickupER','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 


     % X= Response times (alias Waiting Times)
    cd(behavioral_folder);
    load('pickupBehav');
    
%     for subi=1:nSubjs;
% 
%         fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
%         pickupBehav(subi) = load(fname_BehavData);
% 
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlations here:

cd(correlation_folder2);

        for subi=1:nSubjs;

            X= [pickupBehav(subi).LogRESPS' ones(length(pickupBehav(subi).LogRESPS'),1)];
            Y= pickupER(subi).avg_trl.trial; %.avg or .trial?  /// pickup_LogER
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

            R(:)= 0; P(:)=0;

            for i=1:nChans

                    for j= 1:nTimes

                       [tmpR,tmpP] = corrcoef(Y(:,i,j),X(:,1));
                       R(i,j) = tmpR(1,2); P(i,j) = tmpP(1,2);

                    end


                disp(['We are at channel' num2str(i)]);

            end

            
            filename= [sprintf('subj%02d_Corr_ER2Contr_semiLogX', subi)]; % add one if all trials mixed by condition
            save(filename,'R','P','-v7.3');

%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Compute Grandaverages of correlations here:

            cd(timeseries_folder);
            
            cd(contrTS_folder);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(correlation_folder2);

            R_all= {};
            for subi=1:nSubjs;
                
                fname_Corr= sprintf('subj%02d_Corr_ER2Contr_semiLogX',subi);
                pickupCorrs(subi) = load(fname_Corr);
                
                templateER.avg = squeeze(pickupCorrs(subi).R); %squeeze to collapse trials dimension
                R_all{subi}= templateER;
                
            end
            

            cd(correlation_folder2);
            save('R_all_ER2Contr_semiLogX', 'R_all','templateER','-v7.3');


            GAVG_R= ft_timelockgrandaverage([], R_all{:});
            
% Plot correlations:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.showlabels= 'yes';
            figure
            ft_multiplotER(cfg,GAVG_R);
            
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% 1a) REGRESSION between RP amplitudes (Y) vs conditions (X), 'raw'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Set parameters and right path:
% 
% nChans= 60;
% nTimes= 2001; % corresponds to [-3 +1]s
% % B= zeros(2,nChans,nTimes); R= zeros(nChans,nTimes); stats= zeros(nChans,nTimes);
% 
% % clear data_Path; 
% % % set new data_Path
% % data_Path= timeseries_folder; % need to fix it in the start up script
% % 
% 
% % Load in useful data
% 
%     % Y= Time-locked amplitudes (alias RP avg)
%     cd(timeseries_folder);
%     load('pickupER');
%     
% %     for subi=1:nSubjs;
% % 
% %         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
% %         pickupER(subi) = load(fname_ER);
% % 
% %     end
% %     save(filename,'pickupER','pickupER','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 
% % 
% 
%     % X= Conditions (alias Time limits)
%     cd(powerspectra_folder);
%     load('pickupCond');
%     
% %     for subi=1:nSubjs;
% %         
% %         fname_Cond= sprintf('subj%02d_usefulinfo',subi);
% %         pickupCond(subi) = load(fname_Cond);
% %         
% %     end
% %      save 'pickupCond' pickupCond;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Compute regressions here:
% 
% cd(regression_folder);
% 
%         for subi=1:nSubjs;
% 
%             % Just re-converting Inf values into an integer number
%             pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
%              
%             X= [pickupCond(subi).RESPTIMES' ones(length(pickupCond(subi).RESPTIMES'),1)];
%             Y= pickupER(subi).avg_trl.trial; %.avg or .trial?
%             if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;
% 
% %             profile on
% 
%            clear B R stats
% 
%             for i=1:nChans
% 
%                     for j= 1:nTimes
% 
%                        [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X);
% 
%                     end
% 
% 
%                 disp(['We are at channel' num2str(i)]);
% 
%             end
% 
%             
%             filename= [sprintf('subj%02d_Regr_ER1', subi)]; % add one if all trials mixed by condition CHANGE TO Log
%             save(filename,'B','R','stats','-v7.3');
% 
% %             profile off;
% %             profile viewer;
%             disp(['Subject ' num2str(subi) ' done']);
% 
%         end
% 
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
% % Compute Grandaverages of regressions here:
% 
%             cd(timeseries_folder);
%             
%             load subj22_TimeS_one; % all conditions in one, trials averaged
%             templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
%             
%             % make a matrix with all the subjects
%             cd(regression_folder);
% 
%             B_all= {};
%             for subi=1:nSubjs;
%                 
%                 fname_Regr= sprintf('subj%02d_Regr_ER1',subi);
%                 pickupRegr(subi) = load(fname_Regr);
%                 
%                 templateER.avg = squeeze(pickupRegr(subi).B); %squeeze to collapse trials dimension
%                 B_all{subi}= templateER;
%                 
%             end
%             
%             statsB_all= {};
%             for subi=1:nSubjs;
%                 
%                 fname_Regr= sprintf('subj%02d_Regr_ER1',subi);
%                 pickupRegr(subi) = load(fname_Regr);
%                 
%                 templateER.avg = squeeze(pickupRegr(subi).stats); %squeeze to collapse trials dimension
%                 statsB_all{subi}= templateER;
%                 
%             end
%             
% 
%             cd(regression_folder);
%             save('B_all_ER1', 'B_all', 'statsB_all','templateER','-v7.3');
% 
% 
%             GAVG_B= ft_timelockgrandaverage([], B_all{1});
%             
% % Plot regressions:
% 
%             cfg=[];
%             % cfg.preproc.lpfilter='yes';
%             % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
%             cfg.layout = 'eeg_64_NM20884N.lay';
%             cfg.linewidth = 1.5;
%             cfg.showlabels= 'yes';
%             figure
%             ft_multiplotER(cfg,GAVG_B);

%% 1b) REGRESSION between RP amplitudes (Y) vs log-transformed conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
% B= zeros(2,nChans,nTimes); R= zeros(nChans,nTimes); stats= zeros(nChans,nTimes);

% clear data_Path; 
% % set new data_Path
% data_Path= timeseries_folder; % need to fix it in the start up script
% 

% Load in useful data

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    load('pickup_LogER');
    
%     for subi=1:nSubjs;
% 
%         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
%         pickupER(subi) = load(fname_ER);
% 
%     end
%     save(filename,'pickupER','pickupER','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 
% 

    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    load('pickupCond');
    
%     for subi=1:nSubjs;
%         
%         fname_Cond= sprintf('subj%02d_usefulinfo',subi);
%         pickupCond(subi) = load(fname_Cond);
%         
%     end
%      save 'pickupCond' pickupCond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute regressions here:

cd(regression_folder);

        for subi=1:nSubjs;
            
            disp(['We are at subject' num2str(subi)]);

            % Just re-converting Inf values into an integer number
            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
            pickupCond(subi).Logcond= log(pickupCond(subi).newcond);
             
            X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)];
            Y= pickupER(subi).avg_trl.trial; %.avg or .trial? //// pickup_LogER
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

           clear B R stats

            for i=1:nChans

                    for j= 1:nTimes

                       [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X); % Multiple linear regression
%                        LM(:,i,j)= fitlm(X(:,1),Y(:,i,j)) % linear regression


                    end


                disp(['We are at channel' num2str(i)]);

            end

            
            filename= [sprintf('subj%02d_Regr_ER1Log_semiLogX', subi)]; % add one if all trials mixed by condition CHANGE TO Log
            save(filename,'B','R','stats','-v7.3');

%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Compute Grandaverages of regressions here:

            cd(timeseries_folder);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(regression_folder);

            B_all= {};
            for subi=1:nSubjs;
                
                fname_Regr= sprintf('subj%02d_Regr_ER1Log_semiLogX',subi);
                pickupRegr(subi) = load(fname_Regr);
                
                templateER.avg = squeeze(pickupRegr(subi).B(1,:,:)); %squeeze to collapse trials dimension
                B_all{subi}= templateER;
                
            end
            
%             statsB_all= {};
%             for subi=1:nSubjs;
%                 
%                 fname_Regr= sprintf('subj%02d_Regr_ER1Log',subi);
%                 pickupRegr(subi) = load(fname_Regr);
%                 
%                 templateER.avg = squeeze(pickupRegr(subi).stats); %squeeze to collapse trials dimension
%                 statsB_all{subi}= templateER;
%                 
%             end
            


            cd(regression_folder);
            save('B_all_ER1Log_semiLogX', 'B_all','templateER','-v7.3');


            GAVG_B= ft_timelockgrandaverage([], B_all{:});
            save('GAVG_B_ER1_semiLogX', 'GAVG_B','-v7.3');
            
% Plot regressions:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.showlabels= 'yes';
            cfg.comment= 'no'
            cfg.highlight= 'numbers'; %'numbers' %''labels'
            %cfg.highlightchannel=  {'EEG030'};
%             cfg.highlightsymbol= 'x'; %default = 'o')
%             cfg.highlightcolor= highlight marker color (default = [0 0 0] (black))
%             cfg.highlightsize= 48;
            cfg.highlightfontsize= 48;
            cfg.colorbar= 'yes';
            cfg.style= 'straight'; %both' or %'fill'
            cfg.gridscale= 400;
            cfg.xlim= [-0.714 -0.306];
            figure
            ft_topoplotER(cfg,GAVG_B);

%% 2a) REGRESSION between RP amplitudes (Y) vs log-transformed response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Set parameters and right path:
% 
% nChans= 60;
% nTimes= 2001; % corresponds to [-3 +1]s
% 
% % 
% % clear data_Path; %set right path
% % data_Path= timeseries_folder; % need to fix it in the start up script
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Load in useful data:
% 
%     % Y= Time-locked amplitudes (alias RP avg)
%     cd(timeseries_folder);
%     load('pickupER');
%     
% %     for subi=1:nSubjs;
% % 
% %         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
% %         pickupER(subi) = load(fname_ER);
% % 
% %     end
%     
%     % X= Response times (alias Waiting Times)
%     cd(behavioral_folder);
%     load('pickupBehav');
%     
% %     for subi=1:nSubjs;
% % 
% %         fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
% %         pickupBehav(subi) = load(fname_BehavData);
% % 
% %     end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Compute regressions here:
% 
% cd(regression_folder);
% 
%         for subi=1:nSubjs;
% 
%             X= [pickupBehav(subi).RESPTIMES' ones(length(pickupBehav(subi).RESPTIMES'),1)];
%             Y= pickupER(subi).avg_trl.trial; %.avg or .trial?
%             if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;
% 
% %             profile on
% 
%             clear B R stats;
% 
%             for i=1:nChans
% 
%                     for j= 1:nTimes
% 
%                        [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X);
%             
%                     end
% 
% 
%                 disp(['We are at channel' num2str(i)]);
% 
%             end
% 
% 
%             filename= [sprintf('subj%02d_Regr_ER2', subi)]; % add one if all trials mixed by condition
%             save(filename,'B','R','stats','-v7.3');
%             
% %             profile off;
% %             profile viewer;
%             disp(['Subject ' num2str(subi) ' done']);
% 
%         end
% 
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
% % Compute Grandaverages of regressions here:
% 
%             cd(data_Path);
%             
%             load subj22_TimeS_one; % all conditions in one, trials averaged
%             templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
%             
%             % make a matrix with all the subjects
%             cd(regression_folder);
% 
%             R_all= {};
%             for subi=1:nSubjs;
% 
%                 fname_Corr= sprintf('subj%02d_Regr_ER2',subi);
%                 pickupRegr(subi) = load(fname_Corr);
% 
%                 templateER.avg = squeeze(pickupRegr(subi).R); %squeeze to collapse trials dimension
%                 B_all{subi}= templateER;
% 
%             end
%             
%             
%             cd(regression_folder);
%             save('B_all_ER2', 'B_all','templateER');
% 
% 
%             GAVG_R= ft_timelockgrandaverage([], B_all{:});
%             save 'GAVG_R_ER2' GAVG_R;
%             
% % Plot regressions:
% 
%             cfg=[];
%             % cfg.preproc.lpfilter='yes';
%             % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
%             cfg.layout = 'eeg_64_NM20884N.lay';
%             cfg.linewidth = 1.5;
%             cfg.showlabels= 'yes';
%             figure
%             ft_multiplotER(cfg,GAVG_R);  

%% 2b) REGRESSION between RP amplitudes (Y) vs log-transformed response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s

% 
% clear data_Path; %set right path
% data_Path= timeseries_folder; % need to fix it in the start up script
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in useful data:

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    load('pickupER');
    
%     for subi=1:nSubjs;
% 
%         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
%         pickupER(subi) = load(fname_ER);
% 
%     end
    
    % X= Response times (alias Waiting Times)
    cd(behavioral_folder);
    load('pickupBehav');
    
%     for subi=1:nSubjs;
% 
%         fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
%         pickupBehav(subi) = load(fname_BehavData);
% 
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute regressions here:

cd(regression_folder);

        for subi=1:nSubjs;

            X= [pickupBehav(subi).LogRESPS' ones(length(pickupBehav(subi).LogRESPS'),1)];
            Y= pickupER(subi).avg_trl.trial; %.avg or .trial? /// pickup_LogER
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

            clear B R stats;

            for i=1:nChans

                    for j= 1:nTimes

                       [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X);
            
                    end


                disp(['We are at channel' num2str(i)]);

            end


            filename= [sprintf('subj%02d_Regr_ER2Log_semiLogX', subi)]; % add one if all trials mixed by condition
            save(filename,'B','R','stats','-v7.3');
            
%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of regressions here:

            cd(timeseries_folder);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(regression_folder);

            B_all= {};
            for subi=1:nSubjs;

                fname_Regr= sprintf('subj%02d_Regr_ER2Log_semiLogX',subi);
                pickupRegr(subi) = load(fname_Regr);

                templateER.avg = squeeze(pickupRegr(subi).B(1,:,:)); %squeeze to collapse trials dimension
                B_all{subi}= templateER;

            end
            
            
            cd(regression_folder);
            save('B_all_ER2Log_semiLogX', 'B_all','templateER');


            GAVG_B= ft_timelockgrandaverage([], B_all{:});
            save('GAVG_B_ER2_semiLogX', 'GAVG_B','-v7.3');
            
% Plot regressions:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.commet= 'no';
            cfg.showlabels= 'yes';
            cfg.highlight= 'numbers'; %'numbers' %''labels'
            cfg.highlightchannel=  {'EEG030'};
%             cfg.highlightsymbol= 'x'; %default = 'o')
%             cfg.highlightcolor= highlight marker color (default = [0 0 0] (black))
%             cfg.highlightsize= 48;
            cfg.highlightfontsize= 48;
            cfg.colorbar= 'yes';
            cfg.style= 'straight'; %both' or %'fill'
            cfg.gridscale= 400;
            cfg.xlim= [-0.742 -0.26];
            figure
            ft_topoplotER(cfg,GAVG_B);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END
