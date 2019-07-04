%% Eelke Spaak's way
addpath(genpath('/Users/bt_neurospin/matlab/DrosteEffect-BrewerMap-5b84f95')); % or wherever you put brewermap
ROI= [20,21,29,30,31,39,40];
    for condi= 1: 5
        cfg = [];
        cfg.frequency = [5 30];
        cfg.channel = 30;
        cfg.latency = [-2 0];
        freq = ft_selectdata(cfg, Grand_TFR{condi})
        meanpow = squeeze(freq.powspctrm);
        % The finer time and frequency axes:
        tim_interp = linspace(-2.75, 0, 512);
        freq_interp = linspace(5, 30, 512);
        % We need to make a full time/frequency grid of both the original and
        % interpolated coordinates. Matlab's meshgrid() does this for us:
        [tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
        [tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);
        % And interpolate:
        pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,...
            tim_grid_interp, freq_grid_interp, 'spline');
        figure();
        imagesc(tim_interp, freq_interp, pow_interp);
        xlim([-2 0]);
        axis xy;
        xlabel('Time (s)','FontSize',34);
        ylabel('Frequency (Hz)','FontSize',34);
        title(['Time-frequency in Cz, condition ' num2str(condi)],'FontSize',34);
        set(gca,'FontSize',34);
        %clim = max(abs(meanpow(:)));
        caxis([-4.5*10^-12 4.5*10^-12]);
        hold on;
        h = colorbar();
        ylabel(h, 'Power vs baseline (dB)');
        colormap(brewermap(256, '*RdYlBu'));
        % Let's also add a line at t = 0s for clarity:
        hold on;
        plot(zeros(size(freq_interp)), freq_interp, 'k:');
    end
    
time_of_interest = 0.3;
freq_of_interest = 50;
timind = nearest(tim_interp, time_of_interest);
freqind = nearest(freq_interp, freq_of_interest);
pow_at_toi = pow_interp(:,timind);
pow_at_foi = pow_interp(freqind,:);
Collapse




