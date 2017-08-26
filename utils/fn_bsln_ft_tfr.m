function bslnd_tfr = fn_bsln_ft_tfr(tfr, bsln_tfr, bsln_type)
%% Baseline correct one TFR based on another TFR, both from ft_freqanalysis

if ~strcmp(tfr.dimord,'rpt_chan_freq_time')
    error('Check dimord to be sure trial dimension is first!')
end

bslnd_tfr = tfr;
for t = 1:size(tfr.powspctrm,1)
    for ch = 1:size(tfr.powspctrm,2)
        for f = 1:size(tfr.powspctrm,3)
            trl  = tfr.powspctrm(t,ch,f,:);
            bsln = bsln_tfr.powspctrm(t,ch,f,:);
            switch bsln_type
                case 'zscore'
                    bslnd_tfr.powspctrm(t,ch,f,:) = (trl-nanmean(bsln))/nanstd(bsln);
                case 'demean'
                    bslnd_tfr.powspctrm(t,ch,f,:) = trl-nanmean(bsln);
                case 'my_relchange'
                    bslnd_tfr.powspctrm(t,ch,f,:) = (trl-nanmean(bsln))/nanmean(bsln);
                otherwise
                    error('Why are you using this if not for zscore? Adeen isnt implemented yet!')
            end
        end
    end
end

end