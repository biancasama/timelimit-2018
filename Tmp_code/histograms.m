%% Histograms

% averages

Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]};
Cond_names= [2 4 8 16 Inf];
Posiz1= 1:5;
Posiz2= 6:10;

figure('color','white');
    for condi= 1:length(Posiz1)   
    figure('color','white');
        ax(condi)=subplot(2,5,condi);
        ax(condi).FontSize=34;
        subplot(ax(condi));
        xlabel('Waiting times (sec)','FontSize',34);
        ylabel('Probability','FontSize',34);
        title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',34);
        axis([ax(condi).XLim 0 0.6]);
    %     set(gcf,'color','w');

            h1(condi)=histogram(behavStats.mdWT(:,condi),'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.6); %'Normalization','count','BinMethod','auto'
            h1(condi).NumBins= 10;
        box off;
        hold on;
    end   

    Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]};
Cond_names= [2 4 8 16 Inf];
Cond_values= [2 4 8 16 36];
Posiz1= 1:5;
Posiz2= 6:10;

figure('color','white');
for condi= 1:length(Posiz1)   
    %     figure;
    ax(condi)=subplot(2,5,condi);
    subplot(ax(condi));
    
    for subi= 1:22
        
        mydata= pickupBehav(subi).good_resps_cond{condi};
        
        h1(condi)=histogram(mydata,'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.6); %'Normalization','count','BinMethod','auto'
%                   h1(condi).NumBins= 20;
%                 h1(condi).BinWidth= 0.5;
        %         NBins= morebins(h(condi));
        %         histogram(pickupBehav(subi).normRESPCond{condi});
        xlabel('Waiting times (sec)','FontSize',18);
        ylabel('Probability','FontSize',18);
        title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',18);
        set(ax(condi),'fontsize', 18);
        ax(condi).XLim= [ax(condi).XLim(1) Cond_values(condi)];
        if condi== 1
            ax(condi).YLim== [ax(condi).YLim(1) 0.7];
        end
        if condi== 1
            xticks([0:1:2]);
        elseif condi== 2
            xticks([0:2:4]);
        elseif condi== 3
            xticks([0:4:8]);
        elseif condi== 4
            xticks([0:4:16]);
        else condi== 5
            xticks([0:8:36]);
        end
    
        hold on;
        
    end
    
%         dataforfit= pickupBehav(22).good_resps_cond{condi};
%         distribution = 'Lognormal';
%         pd = fitdist(dataforfit',distribution)
%         x_values = linspace(0, median(pickupBehav(subi).good_resps_cond{condi}), 100);
%         y = pdf(pd,x_values);
%         % subplot(2, 2, 3);
%         hold on;
%         plot(x_values,y,'LineWidth',2)
      
    hold on;
end   

%%%%%
figure('color','white');
for condi= 1:5
    %figure('color','white');
    ax(condi)=subplot(2,5,condi);
    subplot(ax(condi));
%     axis([ax(condi).XLim 0 3]);
    
    myData= behavStats.mdWT(:,condi);
    mycolors= Colors{condi};
    
    h3(condi)=histogram(myData,'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.7);
    
    xlabel('Waiting times (sec)','FontSize',18);
    ylabel('Probability','FontSize',18);
    title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',18);
    set(ax(condi),'fontsize', 18);
    ax(condi).XLim= [ax(condi).XLim(1) Cond_values(condi)];
%     if condi== 1
%         ax(condi).YLim== [ax(condi).YLim(1) 0.7];
%     end
%     if condi== 1
%         xticks([0:1:max(behavStats.mdWT(:,1))]);
%     elseif condi== 2
%         xticks([0:2:max(behavStats.mdWT(:,2))]);
%     elseif condi== 3
%         xticks([0:4:max(behavStats.mdWT(:,3))]);
%     elseif condi== 4
%         xticks([0:4:max(behavStats.mdWT(:,4))]);
%     else condi== 5
%         xticks([0:8:max(behavStats.mdWT(:,5))]);
%     end

    % fitting not working
%     distribution = 'Lognormal';
%     pd = fitdist(myData,distribution)
%     x_values = linspace(0, max(myData),100);
%     y = pdf(pd,x_values);
%     % subplot(2, 2, 3);
%     hold on;
%     plot(x_values,y,'LineWidth',2);
%    
    box off;
    hold on;
    
end
