 % compute the fractal and original spectra
    cfg               = [];
    cfg.foilim        = [1 60];
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_peaks';
    fractal = ft_freqanalysis(cfg, roi_trl);
    
    original = ft_freqanalysis(cfg, roi_trl);

    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2./x1';  % equivalent to 10^(log10(x2)-log10(x1))
    oscillatory_alt = ft_math(cfg, fractal, original);

    % display the spectra on a log-log scale
    figure();
    subplot(1,2,1); hold on;
    plot(log(original.freq), log(original.powspctrm),'k');
    plot(log(fractal.freq), log(fractal.powspctrm));
    plot(log(fractal.freq), log(oscillatory.powspctrm));
    xlabel('log-freq'); ylabel('log-power'); grid on;
    legend({'original','fractal','oscillatory = spectrum-fractal'},'location','southwest');
    if F~=0 && O==0
      title('pure fractal signal');
    elseif F==0 && O~=0
      title('pure oscillatory signal');
    elseif F~=0 && O~=0
      title('mixed signal');
    end
    subplot(1,2,2); hold on;
    plot(log(original.freq), log(original.powspctrm),'k');
    plot(log(fractal.freq), log(fractal.powspctrm));
    plot(log(oscillatory_alt.freq), log(oscillatory_alt.powspctrm));
    xlabel('log-freq'); ylabel('log-power'); grid on;
    legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
    title('oscillatory = spectrum / fractal');