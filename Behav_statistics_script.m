%% Statistics on behavioural data
% AUTHOR: Bianca
% DATE: 20th of September. 

%{
    This script takes the results from behavioural analysis and computes
    non-parametric statistics (signrank) between conditions
%}

%% START: initialize path
%% Load data

parent_folder='/Volumes/USB_DISK/SCIENCE/TIMELIMIT_backup/Behavioural_MATLAB_output'; %'/Users/bt_neurospin/Desktop/TimeLim_Behav'
% subj_folders = dir(fullfile(parent_folder, 'Subject*'))
% nSubjects= length(subj_folders);
resultsfolder = fullfile(parent_folder,'/Results');
cd(resultsfolder);
% load('conditions_all');
load('condRT_allsubj.mat');

%% SIGNRANK: Wilcoxon test
% [p,h,stats] =signrank(X,Y,'tail','left'): direction for X<Y
% [p,h,stats] =signrank(X,Y,'tail','right'): direction for Y>X

% Loop across conditions and then take unique values to not have doubles,
% then store for all the possible combinations and do correction for
% multiple comparisons.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha= 0.05

% Condition n.1 vs all the others
[p1,h1,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,5),'tail','left')
[p2,h2,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,4),'tail','left')
[p3,h3,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,3),'tail','left')
[p4,h4,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,2),'tail','left')

% Condition n.2 vs all the others
[p5,h5,stats] =signrank(condRT_allsubj.mRT(:,2),condRT_allsubj.mRT(:,5),'tail','left')
[p6,h6,stats] =signrank(condRT_allsubj.mRT(:,2),condRT_allsubj.mRT(:,4),'tail','left')
[p7,h7,stats] =signrank(condRT_allsubj.mRT(:,2),condRT_allsubj.mRT(:,3),'tail','left')

% Condition n.3 vs all the others
[p8,h8,stats8] =signrank(condRT_allsubj.mRT(:,3),condRT_allsubj.mRT(:,5),'tail','left')
[p9,h9,stats] =signrank(condRT_allsubj.mRT(:,3),condRT_allsubj.mRT(:,4),'tail','left')

% Condition n.4 vs 5
[p10,h10,stats] =signrank(condRT_allsubj.mRT(:,4),condRT_allsubj.mRT(:,5),'tail','left')


pvals_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha= 0.001

% Condition n.1 vs all the others
[p11,h11,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,5),'tail','left','alpha',0.001)
[p12,h12,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,4),'tail','left','alpha',0.001)
[p13,h13,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,3),'tail','left','alpha',0.001)
[p14,h14,stats] =signrank(condRT_allsubj.mRT(:,1),condRT_allsubj.mRT(:,2),'tail','left','alpha',0.001)

% Condition n.2 vs all the others
[p15,h15,stats] =signrank(condRT_allsubj.mRT(:,2),condRT_allsubj.mRT(:,5),'tail','left','alpha',0.001)
[p16,h16,stats] =signrank(condRT_allsubj.mRT(:,2),condRT_allsubj.mRT(:,4),'tail','left','alpha',0.001)
[p17,h17,stats] =signrank(condRT_allsubj.mRT(:,2),condRT_allsubj.mRT(:,3),'tail','left','alpha',0.001)

% Condition n.3 vs all the others
[p18,h18,stats] =signrank(condRT_allsubj.mRT(:,3),condRT_allsubj.mRT(:,5),'tail','left','alpha',0.001)
[p19,h19,stats] =signrank(condRT_allsubj.mRT(:,3),condRT_allsubj.mRT(:,4),'tail','left','alpha',0.001)

% Condition n.4 vs 5
[p20,h20,stats] =signrank(condRT_allsubj.mRT(:,4),condRT_allsubj.mRT(:,5),'tail','left','alpha',0.001)


pvals_001= [p11,p12,p13,p14,p15,p16,p17,p18,p19,p20];


%% FDR correction
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);

%{
        Required Input:
        pvals - A vector or matrix (two dimensions or more) containing the
        p-value of each individual test in a family of tests.

        Optional Inputs:
        q       - The desired false discovery rate. {default: 0.05}
        method  - ['pdep' or 'dep']
        report  - ['yes' or 'no']
%}

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_05,0.05,'pdep','yes');
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_05,0.01,'pdep','yes');
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_05,0.005,'pdep','yes');

% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_001,0.05,'pdep','yes');
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_001,0.01,'pdep','yes');
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_001,0.005,'pdep','yes');

%% End (for now)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha= 0.05

% Condition n.1 vs all the others
[p1,h1,stats] =signrank(condRT_allsubj.std(:,1),condRT_allsubj.std(:,5),'tail','left')
[p2,h2,stats] =signrank(condRT_allsubj.std(:,1),condRT_allsubj.std(:,4),'tail','left')
[p3,h3,stats] =signrank(condRT_allsubj.std(:,1),condRT_allsubj.std(:,3),'tail','left')
[p4,h4,stats] =signrank(condRT_allsubj.std(:,1),condRT_allsubj.std(:,2),'tail','left')

% Condition n.2 vs all the others
[p5,h5,stats] =signrank(condRT_allsubj.std(:,2),condRT_allsubj.std(:,5),'tail','left')
[p6,h6,stats] =signrank(condRT_allsubj.std(:,2),condRT_allsubj.std(:,4),'tail','left')
[p7,h7,stats] =signrank(condRT_allsubj.std(:,2),condRT_allsubj.std(:,3),'tail','left')

% Condition n.3 vs all the others
[p8,h8,stats8] =signrank(condRT_allsubj.std(:,3),condRT_allsubj.std(:,5),'tail','left')
[p9,h9,stats] =signrank(condRT_allsubj.std(:,3),condRT_allsubj.std(:,4),'tail','left')

% Condition n.4 vs 5
[p10,h10,stats] =signrank(condRT_allsubj.std(:,4),condRT_allsubj.std(:,5),'tail','left')



