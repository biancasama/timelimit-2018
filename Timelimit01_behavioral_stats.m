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

%% END


