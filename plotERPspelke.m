cd(pwd);
for subi=1:nSubjs;
    
    fname_TimeS= sprintf('subj%02d_TimeS_cond',subi);
    pickupTimeS(subi) = load(fname_TimeS);
   
end
% %
TimeSmatrix=[];
for subi= 1:nSubjs; %nGoodSubjects
    
    for k= 1:5
        TimeSmatrix{subi,k}= pickupTimeS(subi).avg_cond{k}; %pickupSub(i).avg_EEG{k}
    end
end
% 
save TimeSmatrix TimeSmatrix
% 
Grand_ER=[]; 
for k= 1:5
    
    cfg=[];
%     cfg.keepindividual= 'yes';
%     cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
    Grand_ER{k}= ft_timelockgrandaverage(cfg,TimeSmatrix{:,k});
%     Grand_Freq_ROI{k}= mean(Grand_Freq{k}.powspctrm,1);
    
end
save Grand_ER_Ind Grand_ER
%% EELKE's way (http://eelkespaak.nl/blog/customizing-common-m-eeg-plots-part-1-the-event-related-potential-field-erp-f/)
figure('color','white');
for condi=1:5
    
    tl = Grand_ER{condi};
    % Find the index of the channels of interest in the data
    chaninds = match_str(tl.label, {'EEG030'});
    
    % Compute average over channels
    erp = tl.avg(chaninds,:); % add squeeze if indiv subjs | tl.individual
    
    % Plot
%     figure();
hold on;
    plot(tl.time, erp);
    
    micro = 1e-6;
    erp_microV = erp ./ micro;
    
%     figure();
%     colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]};
    Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]};
    plot(tl.time, erp_microV,'LineWidth',2.5); % 'Color',Colors{condi}
%     xlim([-3.0 1.0]);
    xlabel('Time (s)');
    ylabel('Amplitude electric potential (\mu V)');
    
    %grid
    
    % Change axis
%     ax = gca(); % this Gets the Current Axis so we can set properties
%     ax.XAxisLocation = 'origin';
%     ax.YAxisLocation = 'origin';
%     %ax.TickDir = 'out';
%     
%     % Remove the box around the plot, while we're at it:
%     box off;
%     
%     % And move the x-axis label to underneath the axis:
%     ax.XLabel.Position(2) = -30;
    clear tl erp
end
%%
addpath(genpath('/Users/bt_neurospin/matlab/kakearney-boundedline-pkg-50f7e4b'));
conds= [2 4 8 16 Inf];
mycolormap= [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.4940 0.1840 0.5560]; %linespec
channels= [20 28 30]; %20= FC1; 28= C3; 30= CZ (EEG cap 60 electrodes)
chan_names= {'FC1','C3','Cz'};
ROI = [20 21 29 30 31 39 40];
chaninds= [30 20 28 {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'}];
if chaninds == 30
    chan_names= {'Cz'}
elseif chaninds == 20
    chan_names= {'FC1'}
elseif chaninds == 28
    chan_names== {'C3'}
else chaninds== {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'}; %Undefined operator '==' for input arguments of type 'cell'.
    chan_names= {'ROI'}
end

    f1=figure('color','white');
    for condi=1:length(conds);
        %     figure('color','white');
        tl = Grand_ER{condi};
%         colors= Colors{condi};
        % Find the index of the channels of interest in the data
        %chaninds = match_str(tl.label, {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'});
        chaninds = match_str(tl.label, {'EEG030'});
        % Compute average over channels
        erp = squeeze(mean(tl.individual(:,chaninds,:),1)); % add squeeze if indiv subjs | tl.individual
        micro = 1e-6;
        erp_microV = erp ./ micro;
        % First, get the individual trial data, averaged over our channels of
        % interest:
        trialdata = squeeze(mean(tl.individual(:,chaninds,:), 2)) ./ micro;
        % Use the standard deviation over trials as error bounds:
        bounds = sem(trialdata, 1); %sem or std?
        % boundedline will replace the call to plot():
        b=boundedline(tl.time, erp_microV, bounds, 'cmap', mycolormap(condi,:),'alpha'); % alpha makes bounds transparent
        b.LineWidth= 2;
        %stattimes= [-1.0320 -0.316];
%         stattimes1= linspace(-0.742,-0.26,60);
% %         stattimes2= linspace(-1.398,-1.098,40);
%         y1= linspace(1,1,60);
% %         y2= linspace(1,1,40);
%         txt1 = {'p < .05'};
%         t1=text(-0.75,2,txt1,'FontSize',25);
% %         txt2 = {'p < .05'};
% %         t2=text(-1.4,2,txt2,'FontSize',25);
%         l1=line(stattimes1,y1); %,'Color','black','Marker','*'
%         l1.Color= 'black'; 
%         l1.LineWidth= 3;
%         l1.Marker= '*';
%         l1.MarkerSize= 10;
%         darkGrey1  = [0.4 0.4 0.4]; darkGrey2  = [0.2 0.2 0.2]; darkGrey3  = [0.1 0.1 0.1]; 
% %         l2=line(stattimes2,y2); %,'Color','black','Marker','*'
% %         l2.Color= darkGrey1;
% %         l2.LineWidth= 2;
% %         l2.Marker= '*';
% %         l2.MarkerSize= 10;
        
        
        % Add all our previous improvements:
        xlabel('Time (s)','FontSize',34);
        ylabel('ERP amplitude (\muV)','FontSize',34);
        title(cellstr(chan_names),'FontSize',34);
        set(gca,'FontSize',34);
        %     legend(subset,'FontSize',12)
        %     ax = gca();
        %     ax.XAxisLocation = 'origin';
        %     ax.YAxisLocation = 'origin';
        %     ax.TickDir = 'out';
        %     box off;
        %     ax.XLabel.Position(2) = -60;
        hold on;
        hold off;
    end
    
    f2=figure('color','white');
    for condi=1:length(conds);
        %     figure('color','white');
        tl = Grand_ER{condi};
%         colors= Colors{condi};
        % Find the index of the channels of interest in the data
        chaninds = match_str(tl.label, {'EEG030'});
        % Compute average over channels
        erp = squeeze(mean(tl.individual(:,chaninds,:),1)); % add squeeze if indiv subjs | tl.individual
        micro = 1e-6;
        erp_microV = erp ./ micro;
        % First, get the individual trial data, averaged over our channels of
        % interest:
        trialdata = squeeze(mean(tl.individual(:,chaninds,:), 2)) ./ micro;
        % Use the standard deviation over trials as error bounds:
        bounds = sem(trialdata, 1); %sem or std?
        % boundedline will replace the call to plot():
        b=boundedline(tl.time, erp_microV, bounds, 'cmap', mycolormap(condi,:),'alpha'); % alpha makes bounds transparent
        b.LineWidth= 2;
        stattimes1= linspace(-0.742,-0.26,60);%linspace(-0.722,-0.448,40);
        %stattimes2= linspace(-1.868,-1.63,20);
        y1= linspace(1,1,60);
%         y2= linspace(1,1,20);
        txt1 = {'p < .05'};
        t1=text(-0.75,2,txt1,'FontSize',25);
%         txt2 = {'p < .05'};
%         t2=text(-1.87,4,txt2,'FontSize',25);
        l1=line(stattimes1,y1); %,'Color','black','Marker','*'
        l1.Color= 'black'; 
        l1.LineWidth= 3;
        l1.Marker= '*';
        l1.MarkerSize= 10;
        darkGrey1  = [0.4 0.4 0.4]; darkGrey2  = [0.2 0.2 0.2]; darkGrey3  = [0.1 0.1 0.1]; 
%         l2=line(stattimes2,y2); %,'Color','black','Marker','*'
%         l2.Color= darkGrey2;
%         l2.LineWidth= 2;
%         l2.Marker= '*';
%         l2.MarkerSize= 10;
        
        
        % Add all our previous improvements:
        xlabel('Time (s)','FontSize',34);
        ylabel('ERP amplitude (\muV)','FontSize',34);
        title(cellstr(chan_names),'FontSize',34);
        set(gca,'FontSize',34);
%         lgd= legend(,'FontSize',12)
        %     ax = gca();
        %     ax.XAxisLocation = 'origin';
        %     ax.YAxisLocation = 'origin';
        %     ax.TickDir = 'out';
        %     box off;
        %     ax.XLabel.Position(2) = -60;
        hold on;
        hold off;
    end
    
%% per individual subjects (code not working)
nSubjs=22;
condi=5;
    for subi=1:length(nSubjs);
        figure('color','white');
        tl = Grand_ER{condi};
%         colors= Colors{condi};
        % Find the index of the channels of interest in the data
        %chaninds = match_str(tl.label, {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'});
        chaninds = match_str(tl.label, {'EEG030'});
        % Compute average over channels
        erp = squeeze(tl.individual(subi,chaninds,:)); % add squeeze if indiv subjs | tl.individual
        micro = 1e-6;
        erp_microV = erp ./ micro;
        % First, get the individual trial data, averaged over our channels of
        % interest:
        trialdata = (tl.individual(subi,:) ./ micro);
        % Use the standard deviation over trials as error bounds:
        bounds = std(trialdata, 1); %sem or std?
        % boundedline will replace the call to plot():
        b=boundedline(tl.time, erp_microV, bounds, 'cmap', mycolormap(condi,:),'alpha'); % alpha makes bounds transparent
        b.LineWidth= 2;
        
        % Add all our previous improvements:
        xlabel('Time (s)','FontSize',34);
        ylabel('ERP amplitude (\muV)','FontSize',34);
        title(['subj ' num2str(subi) cellstr(chan_names)],'FontSize',34);
        set(gca,'FontSize',34);
        %     legend(subset,'FontSize',12)
        %     ax = gca();
        %     ax.XAxisLocation = 'origin';
        %     ax.YAxisLocation = 'origin';
        %     ax.TickDir = 'out';
        %     box off;
        %     ax.XLabel.Position(2) = -60;
%         hold on;
        %filename= ['RPprofile_subj_' num2str(subi) '.png'];
%         hold off;
    end
    