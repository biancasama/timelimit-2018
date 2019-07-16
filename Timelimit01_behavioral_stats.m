%% STATS on behavioural data
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on July 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: 
    OUTPUT: 

    FIXME: 

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

statistics_folder= [results_Path, '/Statistics']; % it can be also current_subj_folder
    if ~exist(fullfile(statistics_folder)); mkdir(fullfile(statistics_folder)); end;
    
%% Load behavioral data 

cd(behavioral_folder);
load 'DescriptiveStats'; load 'pickupBehav';
cd(statistics_folder);

%% Kolgoromov test for normality 
% [h,p,k,c] = kstest(x,'Tail','larger')

% copy code from Windows computer

%% Kruskalwallis test =non_parametric ANOVA) + multcompare 

% on medians
X= behavStats.mdWT;
GROUP= {'2s','4s','8s','16s','Inf'};
[P,ANOVATAB,STATS] = kruskalwallis(X,GROUP);

[c,m,h,nms] = multcompare(STATS);
save WTstats1 P c;

% on standard deviations
clear X STATS;

X= behavStats.stdWT;
GROUP= {'2s','4s','8s','16s','Inf'};
[P,ANOVATAB,STATS] = kruskalwallis(X,GROUP);

[c,m,h,nms] = multcompare(STATS);
save WTstats2 P c;

%% Correlation between mean and standard deviation
clear X;

X= behavStats.mdWT; Y= behavStats.stdWT; 
[r,p] = corrcoef(X,Y);
save WTcorr1 r p;

%% CORRELATION between mean/std response times (Y) vs conditions (X)

% cfr. Logtransform code and correlation code for trial-by-trial basis.
% here we do on average
clear X;

A = 1:5; 
% or non-log: 
a=2; r=2;n=5;
A = a*r.^(0:n-1); 

X = repmat(A,22,1); 
% mean 
Y= behavStats.mdWT;
% std
Y= behavStats.stdWT;

[r,p] = corrcoef(X,Y);
save WTcorr2 r p;
save WTcorr3 r p;

%% Correlation between SINGLE-TRIAL response times (Y) vs conditions (X)
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
    
save pickupCorrs_Behav pickupCorrs;

%% If semi-LOG transformed 

% Compute correlation here:

cd(correlation_folder);
R= zeros(nSubjs,1); P= zeros(nSubjs,1);

    for subi=1:nSubjs;

        % Just re-converting Inf values into an integer number 
        pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32; 
        pickupCond(subi).Logcond= log(pickupCond(subi).newcond); % LOG HERE
    
        X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)]; % LOG HERE
        Y= pickupBehav(subi).RESPTIMES'; % LogRESPS
        if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

    %     profile on

        R(:)= 0; P(:)=0;

                   % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                   [tmpR,tmpP] = corrcoef(Y(:),X(:,1));
                   R = tmpR(1,2); P= tmpP(1,2);

                   filename= [sprintf('subj%02d_Corr_Bhv_semiLog', subi)]; % add one if all trials mixed by condition
                   save(filename,'R','P','-v7.3');

    %                profile off;
    %                profile viewer;
                   disp(['Subject ' num2str(subi) ' done']);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reload obtained correlations values in a matrix:

    for subi=1:nSubjs;

        fname_Corr= sprintf('subj%02d_Corr_Bhv_semiLog',subi);
        pickupCorrs(subi) = load(fname_Corr);

    end
    
save 'pickupCorrs_Behav_semiLog' pickupCorrs;
    
%% END
