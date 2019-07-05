%% Run Signrank test with cubic matrix and plot

RPSubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
nRP= length(RPSubjs);
R_CUBE= zeros(60,2001,nRP);

for subi=1:RPSubjs;

    R_CUBE(:,:,subi) = myft_ftavg2dataCube(R_all{subi});

end

TAIL= 'both';
for i=1:nChans

            for j= 1:nTimes

[P(i,j)] = signrank(squeeze(R_CUBE(i,j,:)));

            end

end

SR=500;
t0= 3;
tAx = [0:2000]./SR - 3;
figure; plot(tAx, P(28,:), 'o','Linewidth',2);
title('Correlations between RP average amplitudes from -1s to -0.2s and participants Waiting Times (N=22)');
Xlabel('Time(s)'); Ylabel('P values');%legend('Correlations p-values','P-value= 0.05','Location','Best');
hold on; plot([tAx(1) tAx(end)],[0.05 0.05],'r','Linewidth',2);
txt = 'P-value= 0.05';
text(-2, 0.07,txt,'FontSize',14);
plot([0 0],[0 1],'k--','Linewidth',1.5);

%%%%%%%%%%%%%%%%%

