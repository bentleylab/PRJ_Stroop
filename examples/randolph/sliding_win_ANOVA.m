% equivalent to what Siegel et al. (2015) used in the Science paper
% w2 is omega2 (debiased effect size; multiply with 100 to get PEV)

% prepare the file in FT format
trange = [hfa.time(1) + 0.25 hfa.time(end) - 0.25];
twin = round(0.05 * hfa.fsample); % define smoothing window

w2.time = hfa.time(nearest(hfa.time,trange(1)):nearest(hfa.time,trange(2)));
w2.label = hfa.label;
w2.trial = zeros([3 length(hfa.label) length(w2.time)]);
w2.pval = w2.trial;
w2.dimord = 'rpt_chan_time';
w2.trialinfo = {'Type', 'Prediction', 'ACC'};

% start the loop
for CH = 1:length(hfa.label)
    
    t = 1;
    for tidx = nearest(hfa.time,trange(1)):nearest(hfa.time,trange(2));
        
        % do the test
        [p table] = anovan(squeeze(nanmean(hfa.trial(:,CH,tidx-twin:tidx+twin),3)), ...
            {hfa.trialinfo(:,1) hfa.trialinfo(:,13) hfa.trialinfo(:,9)}, ...
            'model', 'linear', 'sstype',2, ...
            'varnames', strvcat('Type', 'Prediction', 'ACC'), 'display', 'off');
        
        % calc w2
        w2.trial(1,CH,t) = (table{2,2} - (table{2,3} * table{5,5})) / (table{6,2} + table{5,5});
        w2.trial(2,CH,t) = (table{3,2} - (table{3,3} * table{5,5})) / (table{6,2} + table{5,5});
        w2.trial(3,CH,t) = (table{4,2} - (table{4,3} * table{5,5})) / (table{6,2} + table{5,5});
        
        % get significance
        w2.pval(:,CH,t) = p;
        
        % trial keeping
        clear p table
        t = t+1;
        
    end
end
