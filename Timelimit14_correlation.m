%% CORRELATIONS
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on April/May 2019
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

    % Y= Response times (alias Waiting Times)
    for subi=1:nSubjs;

        cd(behavioral_folder);
        fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
        pickupBehav(subi) = load(fname_BehavData);

    end
    
    
    % X= Conditions (alias Time limits)
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
    
        X= [pickupCond(subi).newcond' ones(length(pickupCond(subi).newcond'),1)];
        Y= pickupBehav(subi).RESPTIMES';
        if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

    %     profile on

        R(:)= 0; P(:)=0;

                   % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                   [tmpR,tmpP] = corrcoef(Y(:),X(:,1));
                   R = tmpR(1,2); P= tmpP(1,2);

                   filename= [sprintf('subj%02d_Corr_Bhv', subi)]; % add one if all trials mixed by condition
                   save(filename,'R','P','-v7.3');

    %                profile off;
    %                profile viewer;
                   disp(['Subject ' num2str(subi) ' done']);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reload obtained correlations values in a matrix:

    for subi=1:nSubjs;

        fname_Corr= sprintf('subj%02d_Corr_Bhv',subi);
        pickupCorrs(subi) = load(fname_Corr);

    end
    
    save 'pickupCorrs_Behav' pickupCorrs;

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

    

%% 2a) Time-series: correlation between RP amplitudes (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%      save 'pickupER' pickupER;

    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);
        
    end
    
%      save 'pickupCond' pickupCond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlations here:

cd(regression_folder);

        for subi=1:nSubjs;

            % Just re-converting Inf values into an integer number
            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
            
            X= [pickupCond(subi).newcond' ones(length(pickupCond(subi).newcond'),1)];
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

            
            filename= [sprintf('subj%02d_Corr_ER1', subi)]; % add one if all trials mixed by condition
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
            cd(regression_folder);

            R_all= {};
            for subi=1:nSubjs;
                
                fname_Corr= sprintf('subj%02d_Corr_ER1',subi);
                pickupCorrs(subi) = load(fname_Corr);
                
                templateER.avg = squeeze(pickupCorrs(subi).R); %squeeze to collapse trials dimension
                R_all{subi}= templateER;
                
            end
            

            cd(regression_folder);
            save('R_all_ER1', 'R_all','templateER');


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
      
%% 2b) Time-series: correlation between RP amplitudes (Y) vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
R= zeros(nChans,nTimes); P= zeros(nChans,nTimes);

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

            X= [pickupBehav(subi).RESPTIMES' ones(length(pickupBehav(subi).RESPTIMES'),1)];
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


            filename= [sprintf('subj%02d_Corr_ER2', subi)]; % add one if all trials mixed by condition
            save(filename,'R','P','-v7.3');

            profile off;
            profile viewer;
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

                fname_Corr= sprintf('subj%02d_Corr_ER2',subi);
                pickupCorrs(subi) = load(fname_Corr);

                templateER.avg = squeeze(pickupCorrs(subi).R); %squeeze to collapse trials dimension
                R_all{subi}= templateER;

            end
            
            
            cd(regression_folder);
            save('R_all_ER2', 'R_all','templateER');


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


%% 3a) Time-frequency: correlation between power-spectra (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nFreqs= 40;
nTimes= 61;
R= zeros(nChans,nFreqs,nTimes); P= zeros(nChans,nFreqs,nTimes);

data_Path= powerspectra_folder; % need to fix it in the start up script
cd(data_Path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Time-frequencies (alias TFR)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_TFR= sprintf('subj%02d_TFR_bytrial',subi); 
        pickupTFR(subi) = load(fname_TFR);
      
    end
    
%     save pickupTFRtrl pickupTFR;
    
    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlation here:

        
cd(regression_folder);
        
        for subi=1:nSubjs;

            profile on

             R(:)= 0; P(:)=0;

            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;  
            X= [pickupCond(subi).newcond' ones(length(pickupCond(subi).newcond'),1)];
            Y= pickupTFR(subi).TFR_trl.powspctrm;

                if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

                for i=1:nChans
                    for k= 1:nFreqs
                        for j= 1:nTimes

                            % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                            [tmpR,tmpP] = corrcoef(Y(:,i,k,j),X(:,1));
                            R(i,k,j)= tmpR(1,2);   P(i,k,j) = tmpP(1,2);

                        end
                    end

                    disp(['We are at channel' num2str(i)]);

                % save stuff
              
                filename= [sprintf('subj%02d_Corr_TFR1', subi)]; % add one if all trials mixed by condition
                save(filename,'R','P','-v7.3');

            end%     
            
            profile off;
            profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of correlations here:

            cd(data_Path);

            load subj22_TFR_oneCorrect; % all conditions in one, trials averaged
            templateTFR= TFR_one; % fake Fieldtrip structure for inserting correlation values

            % make a matrix with all the subjects
            cd(regression_folder);

            R_all= {};
            for subi=1:nSubjs;

                fname_Corr= sprintf('subj%02d_Corr_TFR1',subi);
                pickupCorrs(subi) = load(fname_Corr);

                templateTFR.powspctrm = squeeze(pickupCorrs(subi).R);
                R_all{subi}= templateTFR;


            end
            
            
            cd(regression_folder);
            save('R_all_TFR1', 'R_all','templateTFR');

            GAVG_R= ft_freqgrandaverage([], R_all{:});
            
% Plot correlations:

            cfg=[];
            % cfg.baseline= [-2 -1.5];
            % cfg.baselinetype= 'relative';
            cfg.zlim= [-0.05 -0.01];
%             cfg.xlim= [-2.3 -1.3];
            cfg.xlim= [-1.3 -0.7];
%             cfg.ylim= [11 31];
            cfg.ylim= [7 15];
            cfg.layout= 'eeg_64_NM20884N.lay';
            cfg.linewidth = 2;
            cfg.showlabels= 'yes';
            figure
%             ft_multiplotTFR(cfg,GAVG_R);
            ft_topoplotTFR(cfg,GAVG_R);
            
%% 3b) Time-frequency: correlation between power-spectra (Y) vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nFreqs= 40;
nTimes= 61;
R= zeros(nChans,nFreqs,nTimes); P= zeros(nChans,nFreqs,nTimes);

data_Path= powerspectra_folder; % need to fix it in the start up script
cd(data_Path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Time-frequencies (alias TFR)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_TFR= sprintf('subj%02d_TFR_bytrial',subi); 
        pickupTFR(subi) = load(fname_TFR);
        
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

            profile on

             R(:)= 0; P(:)=0;

          X= [pickupBehav(subi).RESPTIMES' ones(length(pickupBehav(subi).RESPTIMES'),1)];
            Y= pickupTFR(subi).TFR_trl.powspctrm;
            
                if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

                for i=1:nChans
                    for k= 1:nFreqs
                        for j= 1:nTimes

                            % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                            [tmpR,tmpP] = corrcoef(Y(:,i,k,j),X(:,1));
                            R(i,k,j)= tmpR(1,2);   P(i,k,j) = tmpP(1,2);
                            
                        end
                    end

                    disp(['We are at channel' num2str(i)]);

                % save stuff
              
                filename= [sprintf('subj%02d_Corr_TFR2', subi)]; % add one if all trials mixed by condition
                save(filename,'R','P','-v7.3');

                end
                
            profile off;
            profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of correlations here:

            cd(data_Path);

            load subj22_TFR_oneCorrect; % all conditions in one, trials averaged
            templateTFR= TFR_one; % fake Fieldtrip structure for inserting correlation values

            % make a matrix with all the subjects
            cd(regression_folder);

            R_all= {};
            for subi=1:nSubjs;

                fname_Corr= sprintf('subj%02d_Corr_TFR2',subi);
                pickupCorrs(subi) = load(fname_Corr);

                templateTFR.powspctrm = squeeze(pickupCorrs(subi).R);
                R_all{subi}= templateTFR;


            end
            
            
            cd(regression_folder);
            save('R_all_TFR2', 'R_all','templateTFR');
            

            GAVG_R= ft_freqgrandaverage([], R_all{:});
            
% Plot correlations:

            cfg=[];
            % cfg.baseline= [-2 -1.5];
            % cfg.baselinetype= 'relative';
            cfg.zlim= [-0.06 -0.02];
            cfg.xlim= [-2.75 -2.2];
            cfg.ylim= [12 31];
            % cfg.ylim= [15 30];
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 2;
            cfg.showlabels= 'yes';
            figure
%             ft_multiplotTFR(cfg,GAVG_R);
            ft_topoplotTFR(cfg,GAVG_R);

            
%% 3c) Time-frequency: correlation between power-spectra (Y) sorted by conditions vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nFreqs= 40;
nTimes= 61;
R= zeros(nChans,nFreqs,nTimes); P= zeros(nChans,nFreqs,nTimes);

data_Path= powerspectra_folder; % need to fix it in the start up script
cd(data_Path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Time-frequencies (alias TFR)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_TFR= sprintf('subj%02d_TFR_cond_trl',subi);
        pickupTFR(subi) = load(fname_TFR);
        
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

% sorted by condition
uconds=  [2     4     8    16   Inf];

        
        for subi=1:nSubjs;

            profile on

            R={}; P={};

            for condi = 1:length(uconds)

                X{condi}= [pickupBehav(subi).good_resps_cond{condi}' ones(length(pickupBehav(subi).good_resps_cond{condi}'),1)];
                Y{condi}= [pickupTFR(subi).TFR_cond{condi}.powspctrm];

                if isequal(length(Y{condi}(:,1)),length(X{condi}(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

                for i=1:nChans
                    for k= 1:nFreqs
                        for j= 1:nTimes

                            [tmpR,tmpP] = corrcoef(Y{condi}(:,i,k,j),X{condi}(:,1));
                            R{condi}(i,k,j) = tmpR(1,2);   P{condi}(i,k,j) = tmpP(1,2);

                        end
                    end

                    disp(['We are at channel' num2str(i)]);

                end

                % save stuff
                filename= [sprintf('subj%02d_Corr_TFR3', subi)]; % add one if all trials mixed by condition
                save(filename,'R','P','-v7.3');

                disp(['Condition ' num2str(condi) ' done']);

            end%

            profile off;
            profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanity Check

for condi = 1:length(uconds)
    for i=1:nChans
        for k= 1:nFreqs
            for j= 1:nTimes
                if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at nTime n° ' num2str(j) ' of nFreq n°' num2str(j) ' nChan n° ' num2str(i) ' condi n° ' num2str(condi) ]); end;
            end
            if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at nFreq n° ' num2str(k) ' nChan n° ' num2str(i) ' condi n° ' num2str(condi)]); end;
        end
        if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at nChan n° ' num2str(i) ' condi n° ' num2str(condi)]); end;
    end
    if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at condi n° ' num2str(condi)]); end;
end
        
% isequal(template{1}.powspctrm,template{2}.powspctrm,template{3}.powspctrm,template{4}.powspctrm,template{5}.powspctrm)
% isequal(pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of correlations here:

            cd(data_Path);

            load subj22_TFR_cond; % all conditions in one, trials averaged
            templateTFR_cond= TFR_cond; % fake Fieldtrip structure for inserting correlation values

            % make a matrix with all the subjects
            cd(regression_folder);

           
            R_all= {};
            for subi=1:nSubjs;
              
                fname_Corr= sprintf('subj%02d_Corr_TFR3',subi);
                pickupCorrs(subi) = load(fname_Corr);
                
                for condi = 1:length(uconds)
                    
                    templateTFR_cond{condi}.powspctrm= pickupCorrs(subi).R{condi}; 
                    R_all{subi,condi}= templateTFR_cond{condi};
                    
                end
                
            end

            save 'R_all_TFR3' 'R_all';

%Grandaverage

            GAVG_R= {};
            for condi = 1:length(uconds)

                GAVG_R{condi}= ft_freqgrandaverage([], R_all{:,condi});

            end
            
            save GAVG_R_TFR3 GAVG_R;

% Normalize grandavg values

            for condi = 1:length(uconds)

                GAVG_R{condi}.powspctrm = normalize(GAVG_R{condi}.powspctrm);

            end
            
            save GAVG_R_TFR3_norm GAVG_R;


% Plot correlations:

            for condi = 1:length(uconds)

                cfg=[];

                % cfg.baseline= [-2 -1.5];
                % cfg.baselinetype= 'relative';
                % cfg.zlim= [-3e-25 3e-25];
                cfg.xlim= [-2.75 0];
                cfg.ylim= [0 30];
                cfg.zlim= [-0.2 0.2];
                % cfg.ylim= [15 30];
                cfg.layout = 'eeg_64_NM20884N.lay';
                cfg.linewidth = 2;
                cfg.showlabels= 'yes';
                figure

                ft_multiplotTFR(cfg,GAVG_R{condi});

            end




