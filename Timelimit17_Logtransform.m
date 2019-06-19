%% Logtransformation 
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on May 2019; modified in June 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the 

    OUTPUT: R{condi}or P{condi} in 

    FIXME: 

%}
%=========================================================================%
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
    
%% 1a) Behavioral: correlation between response times (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Response times (alias Waiting Times): containing alredy log
    % transformed data
    
    load(pickupBehav);

    % X= Conditions (alias Time limits): do Log transformation in following
    % loop
    for subi=1:nSubjs;

        cd(powerspectra_folder);
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlation here:

cd(regression_folder);
R= zeros(nSubjs,1); P= zeros(nSubjs,1);

    for subi=1:nSubjs;

        % Just re-converting Inf values into an integer number 
        pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32; 
        pickupCond(subi).Logcond= log(pickupCond(subi).newcond);
    
        X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)];
        Y= pickupBehav(subi).LogRESPS';
        if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

    %     profile on

        R(:)= 0; P(:)=0;

                   % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                   [tmpR,tmpP] = corrcoef(Y(:),X(:,1));
                   R = tmpR(1,2); P= tmpP(1,2);

                   filename= [sprintf('subj%02d_Corr_BhvLog', subi)]; % add one if all trials mixed by condition
                   save(filename,'R','P','-v7.3');

    %                profile off;
    %                profile viewer;
                   disp(['Subject ' num2str(subi) ' done']);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reload obtained correlations values in a matrix:

    for subi=1:nSubjs;

        fname_Corr= sprintf('subj%02d_Corr_BhvLog',subi);
        pickupCorrs(subi) = load(fname_Corr);

    end
    
    save 'pickupCorrs_BehavLog' pickupCorrs;

% Plot correlations:

%     a=2; r=2;n=5;
%     s = a*r.^(0:n-1);
%     set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}); 
%     xlabel('Conditions (sec)');
%     ylabel('Waiting Times (sec)');
%     plot(s,pickupCorrs(:).R);
    % 
    % for subi= 1:nSubjs
    %     
    %     plot(pickupBehav(1).RESPTIMES,
    %     
    % end

%% 1a) Behavioral: regression between response times (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Response times (alias Waiting Times): containing alredy log
    % transformed data
    
    load(pickupBehav);

    % X= Conditions (alias Time limits): do Log transformation in following
    % loop
    for subi=1:nSubjs;

        cd(powerspectra_folder);
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute regression here:

cd(regression_folder);
B= zeros(nSubjs,1); R= zeros(nSubjs,1); stats= zeros(nSubjs,1);

    for subi=1:nSubjs;

        % Just re-converting Inf values into an integer number 
        pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32; 
        pickupCond(subi).Logcond= log(pickupCond(subi).newcond);
    
        X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)];
        Y= pickupBehav(subi).LogRESPS';
        if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

    %     profile on

        B(:)= 0; R(:)=0; stats(:)=0;

                   [B,~,R,~,stats]= regress(Y,X);
%                    [tmpR,tmpP] = corrcoef(Y(:),X(:,1));
%                    B= tmpB(1,2); R{:}= tmpR(1,:);

                   filename= [sprintf('subj%02d_Regr_BhvLog', subi)]; % add one if all trials mixed by condition
                   save(filename,'B','R', 'stats','-v7.3');

    %                profile off;
    %                profile viewer;
                   disp(['Subject ' num2str(subi) ' done']);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reload obtained regressions values in a matrix:

    for subi=1:nSubjs;

        fname_Regr= sprintf('subj%02d_Regr_BhvLog',subi);
        pickupRegr(subi) = load(fname_Regr);

    end
    
    save(filename,'pickupRegr_BehavLog', 'pickupRegr','-v7.3');

% Plot correlations:

%     a=2; r=2;n=5;
%     s = a*r.^(0:n-1);
%     set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}); 
%     xlabel('Conditions (sec)');
%     ylabel('Waiting Times (sec)');
%     plot(s,pickupCorrs(:).R);
    % 
    % for subi= 1:nSubjs
    %     
    %     plot(pickupBehav(1).RESPTIMES,
    %     
    % end
%% PREPARE FOR TIME-SERIES
nChans= 60;
nTrls= length(pickupER(subi).avg_trl.trial(:,1)); % you need other lines of code (see below)

% R= zeros(nChans,nTimes); P= zeros(nChans,nTimes);
pickup_posER= [];
pickup_LogER= [];

for subi=1:nSubjs;
    
    for k= 1:length(pickupER(subi).avg_trl.trial(:,1))
        
        for i=1:60
            
            Min_ER{subi}(k,i)= min(pickupER(subi).avg_trl.trial(k,i,:));
            pickup_posER(subi).avg_trl.trial(k,i,:)= pickupER(subi).avg_trl.trial(k,i,:)+ abs(Min_ER{subi}(k,i));
            if isequal(length(pickup_posER(subi).avg_trl.trial(k,i,:)),length(pickupER(subi).avg_trl.trial(k,i,:)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

            pickup_LogER(subi).avg_trl.trial(k,i,:)= log(pickup_posER(subi).avg_trl.trial(k,i,:));
            
        end
        
    end
    
end

save('pickup_LogER', 'pickup_LogER','-v7.3');

%% 2a) Time-series: correlation between RP amplitudes (Y) vs conditions (X)
%%%%%%%%%%%%%%%%pickup_LogER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
R= zeros(nChans,nTimes); P= zeros(nChans,nTimes); 

clear data_Path; 
% set new data_Path
data_Path= timeseries_folder; % need to fix it in the start up script


% Load in useful data

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    for subi=1:nSubjs;

        fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
        pickupER(subi) = load(fname_ER);

    end
    
    save('pickupER','pickupER','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 


    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);
        
    end
    
     save 'pickupCond' pickupCond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlations here:

cd(correlation_folder);

        for subi=1:nSubjs;

            % Just re-converting Inf values into an integer number
            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
            pickupCond(subi).Logcond= log(pickupCond(subi).newcond);
             
            X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)];
            Y= pickupER(subi).avg_trl.trial; %.avg or .trial?
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

            profile on

            R(:)= 0; P(:)=0;

            for i=1:nChans

                    for j= 1:nTimes

                       % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                       [tmpR,tmpP] = corrcoef(Y(:,i,j),X(:,1));
                       R(i,j) = tmpR(1,2); P(i,j) = tmpP(1,2);

                    end


                disp(['We are at channel' num2str(i)]);

            end

            
            filename= [sprintf('subj%02d_Corr_ER1Log', subi)]; % add one if all trials mixed by condition
            save(filename,'R','P','-v7.3');

            profile off;
            profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Compute Grandaverages of correlations here:

            cd(timeseries_folder);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(correlation_folder);

            R_all= {};
            for subi=1:nSubjs;
                
                fname_Corr= sprintf('subj%02d_Corr_ER1Log',subi);
                pickupCorrs(subi) = load(fname_Corr);
                
                templateER.avg = squeeze(pickupCorrs(subi).R); %squeeze to collapse trials dimension
                R_all{subi}= templateER;
                
            end
            

            cd(correlation_folder);
            save('R_all_ER1Log', 'R_all','templateER','-v7.3');


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
            
%% 2a) Time-series: regression between RP amplitudes (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
B= zeros(2,nChans,nTimes); R= zeros(nChans,nTimes); stats= zeros(nChans,nTimes);

clear data_Path; 
% set new data_Path
data_Path= timeseries_folder; % need to fix it in the start up script


% Load in useful data

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    for subi=1:nSubjs;

        fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
        pickupER(subi) = load(fname_ER);

    end
    
    save(filename,'pickupER','pickupER','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 


    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);
        
    end
    
%      save 'pickupCond' pickupCond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute regressions here:

cd(regression_folder);

        for subi=1:nSubjs;

            % Just re-converting Inf values into an integer number
            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
            pickupCond(subi).Logcond= log(pickupCond(subi).newcond);
             
            X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)];
            Y= pickup_LogER(subi).avg_trl.trial; %.avg or .trial?
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

           clear B R stats

            for i=1:nChans

                    for j= 1:nTimes

                       [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X);
%                        [tmpR,tmpP] = corrcoef(Y(:,i,j),X(:,1));
%                        R(i,j) = tmpR(1,2); P(i,j) = tmpP(1,2);

                    end


                disp(['We are at channel' num2str(i)]);

            end

            
            filename= [sprintf('subj%02d_Regr_ER1Log', subi)]; % add one if all trials mixed by condition CHANGE TO Log
            save(filename,'B','R','stats','-v7.3');

%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Compute Grandaverages of correlations here:

            cd(timeseries_folder);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(regression_folder);

            B_all= {};
            for subi=1:nSubjs;
                
                fname_Regr= sprintf('subj%02d_Regr_ER1',subi);
                pickupRegr(subi) = load(fname_Regr);
                
                templateER.avg = squeeze(pickupRegr(subi).B); %squeeze to collapse trials dimension
                B_all{subi}= templateER;
                
            end
            
            statsB_all= {};
            for subi=1:nSubjs;
                
                fname_Regr= sprintf('subj%02d_Regr_ER1',subi);
                pickupRegr(subi) = load(fname_Regr);
                
                templateER.avg = squeeze(pickupRegr(subi).stats); %squeeze to collapse trials dimension
                statsB_all{subi}= templateER;
                
            end
            

            cd(regression_folder);
            save('B_all_ER1', 'B_all', 'statsB_all','templateER','-v7.3');


            GAVG_B= ft_timelockgrandaverage([], B_all{1});
            
% Plot correlations:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.showlabels= 'yes';
            figure
            ft_multiplotER(cfg,GAVG_R);
            
%% 2b) Time-series: regression between RP amplitudes (Y) vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
% R= zeros(nChans,nTimes); P= zeros(nChans,nTimes);

clear data_Path; %set right path
data_Path= timeseries_folder; % need to fix it in the start up script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in useful data:

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    for subi=1:nSubjs;

        fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
        pickupER(subi) = load(fname_ER);

    end
    
    % X= Response times (alias Waiting Times)
    cd(behavioral_folder);
    for subi=1:nSubjs;

        fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
        pickupBehav(subi) = load(fname_BehavData);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlation here:

cd(regression_folder);

        for subi=1:nSubjs;

            X= [pickupBehav(subi).LogRESPS' ones(length(pickupBehav(subi).LogRESPS'),1)];
            Y= pickup_LogER(subi).avg_trl.trial; %.avg or .trial?
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

            clear B R stats;

            for i=1:nChans

                    for j= 1:nTimes

                       [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X);
            
                    end


                disp(['We are at channel' num2str(i)]);

            end


            filename= [sprintf('subj%02d_Regr_ER2Log', subi)]; % add one if all trials mixed by condition
            save(filename,'B','R','stats','-v7.3');
            
%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of correlations here:

            cd(data_Path);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(regression_folder);

            R_all= {};
            for subi=1:nSubjs;

                fname_Corr= sprintf('subj%02d_Regr_ER2Log',subi);
                pickupRegr(subi) = load(fname_Corr);

                templateER.avg = squeeze(pickupRegr(subi).R); %squeeze to collapse trials dimension
                B_all{subi}= templateER;

            end
            
            
            cd(regression_folder);
            save('B_all_ER2Log', 'B_all','templateER');


            GAVG_R= ft_timelockgrandaverage([], B_all{:});
            save 'GAVG_R_ER2Log' GAVG_R;
            
% Plot correlations:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.showlabels= 'yes';
            figure
            ft_multiplotER(cfg,GAVG_R);  
            
            