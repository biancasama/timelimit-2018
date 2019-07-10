%% Raincloud plots
% --------------------------
% Example: h = raincloud('X', myData, 'box_on', 1, 'color', [0.5 0.5 0.5]) 
% Parameters 
Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]};
Cond_names= [2 4 8 16 Inf];
x_limits= {[0 2],[0 4],[0 8],[0 16],[0 Inf]};
ylims= [
%'box_dodge' ,1,'box_dodge_amount',.15,'dot_dodge_amount',.35,'box_col_match',1

%% Medians

figure('color','white');
for condi= 1:5
    %figure('color','white');
    ax(condi)=subplot(2,5,condi);
    subplot(ax(condi));
%     axis([ax(condi).XLim 0 3]);
    
    myData= behavStats.mdWT(:,condi);
    mycolors= Colors{condi};
    h=raincloud_plot('X',myData,'color',mycolors,'alpha',.7,'line_width',1,'box_on', 1); 
    
    title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',20);
    set(ax(condi),'fontsize', 20);
    xlabel('Waiting times (sec)', 'fontsize', 18);
%     xticks([0:.5:2]);
    box off;
    hold on;
end
    

% all trials
    for condi= 1:length(Posiz1)
        %     figure;
        ax(condi)=subplot(2,5,condi);
        subplot(ax(condi));
        title(['Condition ' num2str(Cond_names(condi)) ' sec']);
        axis([ax(condi).XLim 0 0.6]);
        %     set(gcf,'color','w');
        figure();
        for subi= 1:22
            myData= pickupBehav(subi).good_resps_cond{condi};
            mycolors= Colors{condi};
            h=raincloud_plot('X',myData,'color',mycolors,'FaceAlpha',0.6); %

            hold on;
        end
        hold on;
    end