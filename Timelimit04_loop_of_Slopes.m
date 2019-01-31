%% Compute SLOPEs for Timelimit2018
%=========================================================================%
% This routine computes slopes for individual participants from the Timelimit
% experiment 2018 and all its conditions and plots them on mean RP amplitudes.

% Selected channels for this analysis: 20= FC1; 28= C3; 30= CZ.

% Output: LateRP{condi}(j,:), RP amplitudes from the window of interest;
%         SlopeRP{condi}(:,j), containing slopes and intercepts;
%         Y{condi}(l,k), corresponding to Y= bX + a *;
%         *formula for regression from Howell's statistics book (2013). 

% AUTHOR: Bianca Trov?
% DATE: 23 July 2018; modified: 24 July, 10 August, 25 September 2018.


% HOW:
%{ 

1)Choose the subject you want to analyse
2)Load subject's mean RP amplitude file ('avg_EEG.mat'). 
3)Loop across the 5 conditions and across channels to put the RP mean averages
  in a 3x500 matrix (3 channels x 500 samples= 1 sec of window of
  interest).
4)Loop across the 5 conditions and across channels to compute the slope.
5)Loop across the 5 conditions and across channels to fit the line.
6)Save results in folder called 'Slopes'.
7)Plot the slopes vs the RP amplitudes.

%}

%% Set paths (27.09)
clearvars;
% Define some paths 
if strcmp(computer, 'MACI64')% on my laptop
    script_Path= '/Volumes/USB_DISK/TIMELIMIT_backup/SCRIPTS_ANALYSES/MEEG'; % here you find all the scripts for preprocessing/analysing MEEG data
    data_Path = '/Volumes/USB_DISK/TIMELIMIT_backup/MEEG_fif_files'; % = parent_folder: all the raw data are stored here.
    ft_Path = '//Users/bt_neurospin/matlab/FIELDTRIP/fieldtrip-20170405'; % Fieldtrip tools 
    tool_Path= '/Users/bt_neurospin/matlab'; % other useful functions for matlab (written by Aaron)
end

% add some paths
addpath(genpath(script_Path)); % general pre-proc path
addpath(genpath(tool_Path));
addpath(ft_Path); ft_defaults; % start up fieldtrip [NEW]
addpath([ft_Path '/engine']); % start up FT engines [NEW] NOTE: what is this actually doing? 

%% Choose the subject

prompt={'Which participant do you want to look at?'};
name='Subject number';
numlines=1;
answer=inputdlg(prompt,name,numlines);
subjnum= str2double(answer);

%% Move to the right folder/path

parent_folder= data_Path; 
subj_folders = dir(fullfile(parent_folder, 'subj*'));
current_subj_folder = fullfile(parent_folder, subj_folders(subjnum).name);
% cd(current_subj_folder);
% cd(fullfile(current_subj_folder,'Timeseries/Amplitudes/Baseline_Haggard'));
cd(fullfile(current_subj_folder,'Timeseries/Amplitudes/No_Baseline'));

%% Load data from current subject/folder
% clear ans answer parent_folder name numlines prompt subj_folders 
load('/timeseries_EEG.mat') % modified 11th November

%% Parameters loop

conditions_all= [1 2 3 4 5];
n_conditions= length(conditions_all);
channels= [20 28 30]; %20= FC1; 28= C3; 30= CZ (EEG cap 60 electrodes)

%% Loop starts: are you sure about the window? 

for condi= 1: n_conditions
    
    for j= 1:numel(channels)
        
        % Put together the RP amplitudes from the *window of interest*
        % [-1sec:0sec]or 1001:1500 corresponding to 1 sec or 500 samples - sampling
        % rate 500Hz) and from the 3 channels of interest in a matrix
        % 3x500.
        LateRP{condi} = zeros(3,500);
%         LateRP{condi}(j,:)=(avg_EEG{condi}.avg(channels(j),1001:1500));
        LateRP{condi}(j,:)=(avg{condi}.avg(channels(j),1001:1500));
%         LateRPrs{condi} = zeros(3,500); % NO NEEDED
        % Reshape matrix for computing regression: "rows corresponding to observations
        % and columns to predictor variables".
%       
        LateRPrs{condi}(:,j)= LateRP{condi}(j,:)';
%         
%         SlopeRP{condi} = zeros(3,500); % NO NEEDED
        % COMPUTE THE SLOPE HERE
        % Check help regress: B = regress(Y,X). X is an n-by-p design matrix, with rows
        % corresponding to observations and columns to predictor variables.  Y is
        % an n-by-1 vector of response observations. X should include a column of ones so that the model contains a constant
        % term.
        
        SlopeRP{condi}(:,j) = regress(LateRPrs{condi}(:,j),[[1:500]' ones(500,1)]);
        
    end
end

%% added 11th November

% compute percentage non-negative RP slopes

% for condi= 1: n_conditions
%     
%     store{condi}=(SlopeRP{condi}(1,:)>=0)
% %     percentage= (ratio*100)/15;
% %     disp(ratio)
%     tot{condi}=sum(store{condi})
%     if percentage>
% end
% 
% end

%% Following Howell's book of Statistical Methods for Psychology, p. 262
% (chapter 9, Correlation and Regression):
% Y= bX + a (equation to fit a line): b= slope, 1st value on output of SlopeRP, a=intercept, 2nd value on output of SlopeRP
for condi= 1: n_conditions
    for l= 1:3 % because 3 channels so 3 outputs
        for k = 1:500 % 500 samples for 1 sec of time window [-1s 0s]
            Y{condi}(l,k) = SlopeRP{condi}(1,l)*k + SlopeRP{condi}(2,l); % Y= bX + a
        end
    end
    
end

%% Save stuff in the right folder (correcred: 10th December Baseline_Haggard)

% Create the folder if it doesn't exist already.
[status, msg] = mkdir(current_subj_folder,'Timeseries/Slopes/No_Baseline');
timeseries_folder= fullfile(current_subj_folder,'Timeseries');
cd(fullfile(timeseries_folder,'/Slopes/No_Baseline'));
save SlopeRP SlopeRP; % removed "save LateRP LateRP; save LateRPrs LateRPrs";
save Y Y;

%% Now plot it again with the line fit

% Total (18.07.18)
figfolder = fullfile(timeseries_folder,'/Slopes/No_Baseline');
cd(figfolder);

 h=figure;
for condi= 1: n_conditions
    
    for j= 1:numel(channels)
       
        %subplot?
        plot(avg{condi}.time(:),avg{condi}.avg(channels(j),:)) %this time use samples on the time axes
        hold on
        plot([-1:(1/500):-0.0020],Y{condi}(j,:),'LineWidth',2)
        hold off
        title(['Readiness Potential for condition ' int2str(condi) ', Channel ' int2str(channels(j)),', Subj ' num2str(subjnum) ])
        xlabel('Time (s)')
        ylabel('Amplitude (\muV)')
        %     str= {'Slope:' int2str(SlopeRP(1,j))} % NOT WORKING
        %     text(2,2,str) % NOT WORKING
        hold on
 %       filename= sprintf('SlopesRP_cond%s.png',condi);
        filename= ['SlopesRP_cond_' int2str(condi) '_chan' int2str(channels(j)) '.png'];

        saveas(h,filename) % save more than one figure in one unique file 
        hold on
        hold off
%         exportfig(h, filename) % not really working
    end
end

%% END 

cd(parent_folder)
fprintf('Subject number %d = DONE.\n',subjnum);
%%%%

