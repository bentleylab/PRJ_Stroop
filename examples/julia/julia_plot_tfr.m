function plot_tfr(subject)

% function plot_fft(subject)
%
%       Plots the mean FFT averaged across trials for a particular electrode 
%       for 6 conditions (edastan, edaodd, idastan, idaodd, switch2eda, switch2ida). 
%       Data already baseline-corrected, and only contains post-stim window
%       (tones: 0-1s, switch: 0.5-3s)
%
%       subject - subject ID, i.e.,'CP15'. taken in as string.
%
%       Usage: plot_fft('CP15')
%       Output: average TFR per electrode for one condition;
%       20 electrodes per figure


taskdir = ['/home/knight/julia/ECoG_DirCog/' subject '/analysis/task/'];
TFAdir  = [taskdir 'TFA/']; 
 
%% Load Data

load([TFAdir 'tfr_byCondition.mat'])


%% Plot

% ~~ select latency of interest
cfg = [];
cfg.latency  = [-0.5 1];
[tfr_eda]       = ft_selectdata(cfg, tfr_eda);
[tfr_ida]       = ft_selectdata(cfg, tfr_ida);
[tfr_stan]      = ft_selectdata(cfg, tfr_stan);
[tfr_odd]       = ft_selectdata(cfg, tfr_odd);
[tfr_edastan]   = ft_selectdata(cfg, tfr_edastan);
[tfr_edaodd]    = ft_selectdata(cfg, tfr_edaodd);
[tfr_idastan]   = ft_selectdata(cfg, tfr_idastan);
[tfr_idaodd]    = ft_selectdata(cfg, tfr_idaodd);


if length(tfr_edastan.label) <= 20 %fft_tones_lf_edastan
    
    tfrplots=figure;
    for chan = 1:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
         
    
elseif length(tfr_edastan.label) > 20 && length(tfr_edastan.label) <= 40
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
       
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);    

    tfrplots=figure;
    for chan = 21:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
 elseif length(tfr_edastan.label) > 40 && length(tfr_edastan.label) <= 60
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
       
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

     tfrplots=figure;
    for chan = 41:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

elseif length(tfr_edastan.label) > 60 && length(tfr_edastan.label) <= 80
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
       
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
  
     tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
 
 tfrplots=figure;
    for chan = 61:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
  
elseif length(tfr_edastan.label) > 80 && length(tfr_edastan.label) <= 100
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
       
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
  
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
   
     tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

 tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
     tfrplots=figure;
    for chan = 81:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 81:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 81:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 81:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

elseif length(tfr_edastan.label) > 100 && length(tfr_edastan.label) <= 120
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs'; 
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
       
    tfrplots=figure;
    for chan = 1:20; 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs1-20.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
   
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 21:40;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-20);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs21-40.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

     tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 41:60;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-40);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs41-60.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
  
 tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 61:80;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-60);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs61-80.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

     tfrplots=figure;
    for chan = 81:100;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 81:100;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 81:100;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 81:100;
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-80);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs81-100.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);

     tfrplots=figure;
    for chan = 101:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-100);
        ft_singleplotTFR(cfg, tfr_edastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edastan_elecs101-120.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 101:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-100);
        ft_singleplotTFR(cfg, tfr_edaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_edaodd_elecs101-120.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 101:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-100);
        ft_singleplotTFR(cfg, tfr_idastan); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idastan_elecs101-120.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
    
    tfrplots=figure;
    for chan = 101:length(tfr_edastan.label); 
        cfg = [];
        cfg.channel = chan;
        cfg.zlim = 'maxabs';
        subplot(4,5,chan-100);
        ft_singleplotTFR(cfg, tfr_idaodd); 
        title(tfr_edastan.label{chan})
    end
    filename = strcat(subject, '_', 'TFR_idaodd_elecs101-120.png');
    print(tfrplots, [TFAdir filename], '-dpng');
    close(tfrplots);
      
    
end