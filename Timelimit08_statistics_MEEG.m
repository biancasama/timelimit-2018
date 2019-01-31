%% Statistics: Non-parametric tests 

%% Wilcoxon signed rank test
%=========================================================================%
% This routine computes Wilcoxon paired tests between pairs of conditions
% from the Timelimit experiment (2018) and does FDR correction for multiple
% comparisons.

% AUTHOR: Bianca Trov?
% DATE: 27th July 2018; modified: 21st September 2018, 15 October 2018.

%{
Syntax
p = signrank(x)
p = signrank(x,y)
p = signrank(x,y,Name,Value)
[p,h] = signrank(___)
[p,h,stats] = signrank(___)
[___] = signrank(x,m)
[___] = signrank(x,m,Name,Value)

'method','approximate' if N>15
%}

%=========================================================================%
%% Set path, LOAD DATA

clearvars; close all;

% Define some paths 
if strcmp(computer, 'MACI64')
    script_Path= '/Volumes/USB DISK/TIMELIMIT_backup/SCRIPTS_ANALYSES/MEEG'; % here you find all the scripts for preprocessing/analysing MEEG data
    parent_folder = '/Volumes/USB DISK/TIMELIMIT_backup/MEEG_fif_files'; % = parent_folder: all the raw data are stored here.
    ft_Path = '//Users/bt_neurospin/matlab/FIELDTRIP/fieldtrip-20170405'; % Fieldtrip tools 
    tool_Path= '/Users/bt_neurospin/matlab/matlab_internal'; % other useful functions for matlab (written by Aaron)
end

% add some paths
addpath(genpath(script_Path)); % general pre-proc path
addpath(genpath(tool_Path));
addpath([ft_Path '/engine']);

resultsfolder = fullfile(parent_folder,'/Results');
cd(resultsfolder);

%% MEAN RP AMPLITUDE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel 20
clearvars;load('mean_premov_amp_ch20');

% Alpha= 0.05,  tail 'left'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(mean_premov_amp.ch20(:,1),mean_premov_amp.ch20(:,5),'tail','left') % vs Inf
[p2,h2,stats] =signrank(mean_premov_amp.ch20(:,1),mean_premov_amp.ch20(:,4),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(mean_premov_amp.ch20(:,1),mean_premov_amp.ch20(:,3),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(mean_premov_amp.ch20(:,1),mean_premov_amp.ch20(:,2),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(mean_premov_amp.ch20(:,2),mean_premov_amp.ch20(:,5),'tail','left') % vs Inf
[p6,h6,stats] =signrank(mean_premov_amp.ch20(:,2),mean_premov_amp.ch20(:,4),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(mean_premov_amp.ch20(:,2),mean_premov_amp.ch20(:,3),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(mean_premov_amp.ch20(:,3),mean_premov_amp.ch20(:,5),'tail','left') % vs Inf
[p9,h9,stats] =signrank(mean_premov_amp.ch20(:,3),mean_premov_amp.ch20(:,4),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(mean_premov_amp.ch20(:,4),mean_premov_amp.ch20(:,5),'tail','left') % vs Inf

% FDR
pvals_ch20_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch20_05,0.05,'pdep','yes');

%=================================================================================================
% channel 28
clearvars;load('mean_premov_amp_ch28');

% Alpha= 0.05,  tail 'left'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(mean_premov_amp.ch28(:,1),mean_premov_amp.ch28(:,5),'tail','left') % vs Inf
[p2,h2,stats] =signrank(mean_premov_amp.ch28(:,1),mean_premov_amp.ch28(:,4),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(mean_premov_amp.ch28(:,1),mean_premov_amp.ch28(:,3),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(mean_premov_amp.ch28(:,1),mean_premov_amp.ch28(:,2),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(mean_premov_amp.ch28(:,2),mean_premov_amp.ch28(:,5),'tail','left') % vs Inf
[p6,h6,stats] =signrank(mean_premov_amp.ch28(:,2),mean_premov_amp.ch28(:,4),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(mean_premov_amp.ch28(:,2),mean_premov_amp.ch28(:,3),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(mean_premov_amp.ch28(:,3),mean_premov_amp.ch28(:,5),'tail','left') % vs Inf
[p9,h9,stats] =signrank(mean_premov_amp.ch28(:,3),mean_premov_amp.ch28(:,4),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(mean_premov_amp.ch28(:,4),mean_premov_amp.ch28(:,5),'tail','left') % vs Inf

% FDR
pvals_ch28_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch28_05,0.05,'pdep','yes');

%=================================================================================================
% channel 30
clearvars;load('mean_premov_amp_ch30');

% Alpha= 0.05,  tail 'left'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(mean_premov_amp.ch30(:,1),mean_premov_amp.ch30(:,5),'tail','left') % vs Inf
[p2,h2,stats] =signrank(mean_premov_amp.ch30(:,1),mean_premov_amp.ch30(:,4),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(mean_premov_amp.ch30(:,1),mean_premov_amp.ch30(:,3),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(mean_premov_amp.ch30(:,1),mean_premov_amp.ch30(:,2),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(mean_premov_amp.ch30(:,2),mean_premov_amp.ch30(:,5),'tail','left') % vs Inf
[p6,h6,stats] =signrank(mean_premov_amp.ch30(:,2),mean_premov_amp.ch30(:,4),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(mean_premov_amp.ch30(:,2),mean_premov_amp.ch30(:,3),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(mean_premov_amp.ch30(:,3),mean_premov_amp.ch30(:,5),'tail','left') % vs Inf
[p9,h9,stats] =signrank(mean_premov_amp.ch30(:,3),mean_premov_amp.ch30(:,4),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(mean_premov_amp.ch30(:,4),mean_premov_amp.ch30(:,5),'tail','left') % vs Inf

% FDR
pvals_ch30_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch30_05,0.05,'pdep','yes');

%% SLOPE RP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel 20
clearvars;load('Slope_all');

% Alpha= 0.05,  tail 'left'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{5}(:),'tail','left') % vs Inf
[p2,h2,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{4}(:),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{3}(:),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{2}(:),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Slope_all.ch20{2}(:),Slope_all.ch20{5}(:),'tail','left') % vs Inf
[p6,h6,stats] =signrank(Slope_all.ch20{2}(:),Slope_all.ch20{4}(:),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(Slope_all.ch20{2}(:),Slope_all.ch20{3}(:),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Slope_all.ch20{3}(:),Slope_all.ch20{5}(:),'tail','left') % vs Inf
[p9,h9,stats] =signrank(Slope_all.ch20{3}(:),Slope_all.ch20{4}(:),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Slope_all.ch20{4}(:),Slope_all.ch20{5}(:),'tail','left') % vs Inf

% FDR
pvals_ch20_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch20_05,0.05,'pdep','yes');

%=================================================================================================
clearvars;load('Slope_all');

% Alpha= 0.05,  tail 'right'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{5}(:),'tail','right') % vs Inf
[p2,h2,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{4}(:),'tail','right') % vs 16sec
[p3,h3,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{3}(:),'tail','right') % vs 8sec
[p4,h4,stats] =signrank(Slope_all.ch20{1}(:),Slope_all.ch20{2}(:),'tail','right') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Slope_all.ch20{2}(:),Slope_all.ch20{5}(:),'tail','right') % vs Inf
[p6,h6,stats] =signrank(Slope_all.ch20{2}(:),Slope_all.ch20{4}(:),'tail','right') % vs 16sec
[p7,h7,stats] =signrank(Slope_all.ch20{2}(:),Slope_all.ch20{3}(:),'tail','right') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Slope_all.ch20{3}(:),Slope_all.ch20{5}(:),'tail','right') % vs Inf
[p9,h9,stats] =signrank(Slope_all.ch20{3}(:),Slope_all.ch20{4}(:),'tail','right') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Slope_all.ch20{4}(:),Slope_all.ch20{5}(:),'tail','right') % vs Inf

% FDR
pvals_ch20_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch20_05,0.05,'pdep','yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel 28
clearvars;load('Slope_all');

% Alpha= 0.05,  tail 'left'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{5}(:),'tail','left') % vs Inf
[p2,h2,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{4}(:),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{3}(:),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{2}(:),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Slope_all.ch28{2}(:),Slope_all.ch28{5}(:),'tail','left') % vs Inf
[p6,h6,stats] =signrank(Slope_all.ch28{2}(:),Slope_all.ch28{4}(:),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(Slope_all.ch28{2}(:),Slope_all.ch28{3}(:),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Slope_all.ch28{3}(:),Slope_all.ch28{5}(:),'tail','left') % vs Inf
[p9,h9,stats] =signrank(Slope_all.ch28{3}(:),Slope_all.ch28{4}(:),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Slope_all.ch28{4}(:),Slope_all.ch28{5}(:),'tail','left') % vs Inf

% FDR
pvals_ch28_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch28_05,0.05,'pdep','yes');

%=================================================================================================
% channel 28
clearvars;load('Slope_all');

% Alpha= 0.05,  tail 'right'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{5}(:),'tail','right') % vs Inf
[p2,h2,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{4}(:),'tail','right') % vs 16sec
[p3,h3,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{3}(:),'tail','right') % vs 8sec
[p4,h4,stats] =signrank(Slope_all.ch28{1}(:),Slope_all.ch28{2}(:),'tail','right') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Slope_all.ch28{2}(:),Slope_all.ch28{5}(:),'tail','right') % vs Inf
[p6,h6,stats] =signrank(Slope_all.ch28{2}(:),Slope_all.ch28{4}(:),'tail','right') % vs 16sec
[p7,h7,stats] =signrank(Slope_all.ch28{2}(:),Slope_all.ch28{3}(:),'tail','right') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Slope_all.ch28{3}(:),Slope_all.ch28{5}(:),'tail','right') % vs Inf
[p9,h9,stats] =signrank(Slope_all.ch28{3}(:),Slope_all.ch28{4}(:),'tail','right') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Slope_all.ch28{4}(:),Slope_all.ch28{5}(:),'tail','right') % vs Inf

% FDR
pvals_ch28_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch28_05,0.05,'pdep','yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel 30
clearvars;load('Slope_all');

% Alpha= 0.05,  tail 'left'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{5}(:),'tail','left') % vs Inf
[p2,h2,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{4}(:),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{3}(:),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{2}(:),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Slope_all.ch30{2}(:),Slope_all.ch30{5}(:),'tail','left') % vs Inf
[p6,h6,stats] =signrank(Slope_all.ch30{2}(:),Slope_all.ch30{4}(:),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(Slope_all.ch30{2}(:),Slope_all.ch30{3}(:),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Slope_all.ch30{3}(:),Slope_all.ch30{5}(:),'tail','left') % vs Inf
[p9,h9,stats] =signrank(Slope_all.ch30{3}(:),Slope_all.ch30{4}(:),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Slope_all.ch30{4}(:),Slope_all.ch30{5}(:),'tail','left') % vs Inf

% FDR
pvals_ch30_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch30_05,0.05,'pdep','yes');

%=================================================================================================
% channel 30
clearvars;load('Slope_all');

% Alpha= 0.05,  tail 'right'

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{5}(:),'tail','right') % vs Inf
[p2,h2,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{4}(:),'tail','right') % vs 16sec
[p3,h3,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{3}(:),'tail','right') % vs 8sec
[p4,h4,stats] =signrank(Slope_all.ch30{1}(:),Slope_all.ch30{2}(:),'tail','right') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Slope_all.ch30{2}(:),Slope_all.ch30{5}(:),'tail','right') % vs Inf
[p6,h6,stats] =signrank(Slope_all.ch30{2}(:),Slope_all.ch30{4}(:),'tail','right') % vs 16sec
[p7,h7,stats] =signrank(Slope_all.ch30{2}(:),Slope_all.ch30{3}(:),'tail','right') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Slope_all.ch30{3}(:),Slope_all.ch30{5}(:),'tail','right') % vs Inf
[p9,h9,stats] =signrank(Slope_all.ch30{3}(:),Slope_all.ch30{4}(:),'tail','right') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Slope_all.ch30{4}(:),Slope_all.ch30{5}(:),'tail','right') % vs Inf

% FDR
pvals_ch30_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch30_05,0.05,'pdep','yes');

%% INTERCEPT RP (added 21/09)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel 20
clearvars;load('Intercept_all');

% Alpha= 0.05

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Intercept_all.ch20{1}(:),Intercept_all.ch20{5}(:),'tail','left') % vs Inf
[p2,h2,stats] =signrank(Intercept_all.ch20{1}(:),Intercept_all.ch20{4}(:),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(Intercept_all.ch20{1}(:),Intercept_all.ch20{3}(:),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(Intercept_all.ch20{1}(:),Intercept_all.ch20{2}(:),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Intercept_all.ch20{2}(:),Intercept_all.ch20{5}(:),'tail','left') % vs Inf
[p6,h6,stats] =signrank(Intercept_all.ch20{2}(:),Intercept_all.ch20{4}(:),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(Intercept_all.ch20{2}(:),Intercept_all.ch20{3}(:),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Intercept_all.ch20{3}(:),Intercept_all.ch20{5}(:),'tail','left') % vs Inf
[p9,h9,stats] =signrank(Intercept_all.ch20{3}(:),Intercept_all.ch20{4}(:),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Intercept_all.ch20{4}(:),Intercept_all.ch20{5}(:),'tail','left') % vs Inf

% FDR
pvals_ch20_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch20_05,0.05,'pdep','yes');

%===================================================================================================
% channel 28 - added 25.09.18
clearvars;load('Intercept_all');

% Alpha= 0.05

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Intercept_all.ch28{1}(:),Intercept_all.ch28{5}(:),'tail','left') % vs Inf
[p2,h2,stats] =signrank(Intercept_all.ch28{1}(:),Intercept_all.ch28{4}(:),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(Intercept_all.ch28{1}(:),Intercept_all.ch28{3}(:),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(Intercept_all.ch28{1}(:),Intercept_all.ch28{2}(:),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Intercept_all.ch28{2}(:),Intercept_all.ch28{5}(:),'tail','left') % vs Inf
[p6,h6,stats] =signrank(Intercept_all.ch28{2}(:),Intercept_all.ch28{4}(:),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(Intercept_all.ch28{2}(:),Intercept_all.ch28{3}(:),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Intercept_all.ch28{3}(:),Intercept_all.ch28{5}(:),'tail','left') % vs Inf
[p9,h9,stats] =signrank(Intercept_all.ch28{3}(:),Intercept_all.ch28{4}(:),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Intercept_all.ch28{4}(:),Intercept_all.ch28{5}(:),'tail','left') % vs Inf

% FDR
pvals_ch28_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch28_05,0.05,'pdep','yes');

%=====================================================================================================
% channel 30

clearvars;load('Intercept_all');
% Alpha= 0.05, default

% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(Intercept_all.ch30{1}(:),Intercept_all.ch30{5}(:),'tail','left') % vs Inf
[p2,h2,stats] =signrank(Intercept_all.ch30{1}(:),Intercept_all.ch30{4}(:),'tail','left') % vs 16sec
[p3,h3,stats] =signrank(Intercept_all.ch30{1}(:),Intercept_all.ch30{3}(:),'tail','left') % vs 8sec
[p4,h4,stats] =signrank(Intercept_all.ch30{1}(:),Intercept_all.ch30{2}(:),'tail','left') % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(Intercept_all.ch30{2}(:),Intercept_all.ch30{5}(:),'tail','left') % vs Inf
[p6,h6,stats] =signrank(Intercept_all.ch30{2}(:),Intercept_all.ch30{4}(:),'tail','left') % vs 16sec
[p7,h7,stats] =signrank(Intercept_all.ch30{2}(:),Intercept_all.ch30{3}(:),'tail','left') % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(Intercept_all.ch30{3}(:),Intercept_all.ch30{5}(:),'tail','left') % vs Inf
[p9,h9,stats] =signrank(Intercept_all.ch30{3}(:),Intercept_all.ch30{4}(:),'tail','left') % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(Intercept_all.ch30{4}(:),Intercept_all.ch30{5}(:),'tail','left') % vs Inf

% FDR
pvals_ch30_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_ch30_05,0.05,'pdep','yes');

%% Kruskal-Wallis test

%=====================================================================================================
% channel 20
clearvars;load('mean_premov_amp_ch20.mat')
[p,tbl,stats] = kruskalwallis(mean_premov_amp_ch20);
[c,m,h,nms] = multcompare(stats);

clearvars;load('median_premov_amp_ch20.mat')
[p,tbl,stats] = kruskalwallis(median_premov_amp_ch20);
[c,m,h,nms] = multcompare(stats);

%=====================================================================================================
% channel 28
clearvars;load('mean_premov_amp_ch28.mat');
[p,tbl,stats] = kruskalwallis(mean_premov_amp_ch28);
[c,m,h,nms] = multcompare(stats);

clearvars;load('median_premov_amp_ch28.mat');
[p,tbl,stats]= kruskalwallis(median_premov_amp_ch28);
[c,m,h,nms] = multcompare(stats);

%=====================================================================================================
% channel 30
clearvars;load('mean_premov_amp_ch30.mat');
[p,tbl,stats] = kruskalwallis(mean_premov_amp_ch30);
[c,m,h,nms] = multcompare(stats);

clearvars;load('median_premov_amp_ch30.mat');
[p,tbl,stats]= kruskalwallis(median_premov_amp_ch30);
[c,m,h,nms] = multcompare(stats);


% FDR
pvals_chans_mean= [p11,p13,p15];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_chans_mean,0.05,'pdep','yes');

pvals_chans_median= [p12,p14,p16];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_chans_median,0.05,'pdep','yes');

% Dunn test 
% x= reshape(mean_premov_amp_ch20,1,[]);
% g= [ones(1,11) repmat(2,1,11) repmat(3,1,11) repmat(4,1,11) repmat(5,1,11)];
% dunn(x,g) % there is no control group


%% Friedman test + Multiple comparison test
%=====================================================================================================
% channel 20
clearvars;load('mean_premov_amp_ch20.mat');
[p,tbl,stats] = friedman(mean_premov_amp_ch20);
[c,m,h,nms] = multcompare(stats);

clearvars;load('median_premov_amp_ch20.mat')
[p,tbl,stats] = friedman(median_premov_amp_ch20);
[c,m,h,nms] = multcompare(stats);

%=====================================================================================================
% channel 28
clearvars;load('mean_premov_amp_ch28.mat');
[p,tbl,stats] = friedman(mean_premov_amp_ch28);
[c,m,h,nms] = multcompare(stats);

clearvars;load('median_premov_amp_ch28.mat');
[p,tbl,stats] = friedman(median_premov_amp_ch28);
[c,m,h,nms] = multcompare(stats);

%=====================================================================================================
% channel 30
clearvars;load('mean_premov_amp_ch30.mat');
[p,tbl,stats] = friedman(mean_premov_amp_ch30);
[c,m,h,nms] = multcompare(stats);

clearvars;load('median_premov_amp_ch30.mat');
[p,tbl,stats] = friedman(median_premov_amp_ch30);
c,m,h,nms] = multcompare(stats);

%% END
%%%%%%%%%%%%%%%%%%%
