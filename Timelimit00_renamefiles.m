%% RENAME (specific)
% function rename_many_files(filespec,findPatt,replPatt,newPrefix,[actuallyDoIt])
% PUT ZERO first and then ONE

% subj01_Regr_ER2Logsubj*_Corr_ER2Log.mat > subj01_Corr_ER2Log.mat

rename_many_files('*R_all_ER2Log.mat','R_all_ER2Log','R_all_ER2mix',[],1);

% subj*_TimeS_cond > subj*_TimeS_condTrl

rename_many_files('subj*_TimeS_cond.mat','cond','condTrl',[],1);

%

rename_many_files('subj*_CONTRcorr*','CONTRcorr_ER1Log','Corr_ER1Contr',[],1);

