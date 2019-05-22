%% Merge data
clearvars;
close('all');
clc

dropboxfolder = '/Users/nick/Dropbox (Attention Group)/Attention Group Team Folder/Myers';
dropboxfolder = 'E:\nmyers\Dropbox (Attention Group)\Attention Group Team Folder\myers';
studypath     = [dropboxfolder '/TemplateSwitching/ECoG'];
toolboxpath   = [dropboxfolder '/toolbox'];
cd(studypath)

fs = filesep;
datapath  = sprintf('%s%sdata',studypath,fs);

fieldtrippath = '/Work/Matlab/toolbox/fieldtrip';
fieldtrippath = '/Users/nick/OneDrive - Nexus365/Matlab/toolbox/fieldtrip/';
fieldtrippath = 'C:\Matlab\fieldtrip-20190419';
addpath(fieldtrippath);
ft_defaults;

addpath([studypath '/scripts'])
addpath([toolboxpath])
addpath([toolboxpath '/decoding_functions'])
addpath([toolboxpath '/CircStat'])
addpath([toolboxpath '/jv10'])
addpath([toolboxpath '/stat'])
addpath(genpath([dropboxfolder '/toolbox/sigproc/chi/HOLO']))
addpath(genpath([dropboxfolder '/toolbox/sigproc/emd']))
addpath(genpath([toolboxpath '/plt']))

sublist = {'CP24a' 'CP24b' 'IR67' 'IR69' 'IR75' 'IR82' 'IR84' 'IR85' 'IR87'};
nsubs   = length(sublist);
%%
for isub = 1:nsubs %subject loop
    substrg  = sublist{isub};
    subpath  = sprintf('%s%s%s',datapath,fs,substrg);
    filename = sprintf('%s%sTF%s%s_MultitaperTF_hifreq.mat',datapath,fs,fs,substrg);
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename);
    ntrials = size(datatfr.powspctrm,1);
    
    fprintf('\n');
    % get rejected trials    
    rejfile = sprintf('%s%s%s_semiautomaticAR.mat',subpath,fs,substrg);
    if ~exist(rejfile,'file')
        warning(sprintf('%s does not exist!',rejfile));
        rejsemiautomatic = [];
    else
        load(rejfile,'artefacts');
        rejsemiautomatic = artefacts.rejsemiautomatic;
    end
    
    % get behavior file
    behavfile = sprintf('%s%s%s_behav_final.mat',subpath,fs,substrg);
    if ~exist(behavfile,'file')
        error(sprintf('%s not found!',behavfile));
    else
        load(behavfile,'behav');
    end
    
    switch substrg
        case 'IR84'
            skiptrl = 1:144;
            bfield = fieldnames(behav);
            bfield(ismember(bfield,'ntrial')) = [];
            for ifield = 1:length(bfield)
                behav.(bfield{ifield})(skiptrl,:,:,:) = [];
            end
            behav.ntrial = size(behav.trial,1);
    end     
    
    artefact                   = zeros(ntrials,1);
    artefact(rejsemiautomatic) = 1;
    usable = ~artefact;
    usable = usable & ~isinf(behav.time) & ~isnan(behav.time) & behav.responded;
    fprintf('\n\n%s: %03d usable trials (%1.2f percent).\n',substrg,nnz(usable),100*nnz(usable)/length(usable));    
    
    cue     = behav.cue(usable);
    stims   = behav.stim(usable);
    stimset = behav.stimset(usable);
    relrule = behav.relrule(usable,:);
    rule    = behav.rule(usable,:);
    rulecol = behav.rulecolor(usable,:);
    block   = behav.block(usable);
    tnum    = behav.tnum(usable);
    targ    = behav.target(usable);
    dstr    = behav.distr(usable);
    targtype              = double(targ);
    targtype(dstr==1)     = 2;
    targtype(~dstr&~targ) = 3;
    cuecol  = rulecol;
    cuecol(cue==2,:) = cuecol(cue==2,[2 1]);
    
    acc   = behav.corr(usable);
    rt    = behav.time(usable);
    rsp   = behav.resp(usable);
    ispractice = behav.ispractice(usable);
    
    time  = datatfr.time;
    freq  = datatfr.freq;
    nfreq = length(freq);
    dat   = permute(datatfr.powspctrm(usable,:,:,:),[1 2 4 3]);
    
    chanlist = [1:length(datatfr.label)]';
    nchans   = length(chanlist);
    
    % Basic setup for final preprocessing
    do_log      = 0;
    base_single = 0;
    base_zscore = 1;
    do_zscore   = 0;
    avg_freq    = 1;
    if do_log
       dat = 10*log10(dat); 
    end
    if base_single
       bwin        = [-1 -0.5];
       ib          = time>=bwin(1) & time<=bwin(2);
       dat = bsxfun(@minus,dat,mean(dat(:,:,ib,:),3)); 
    end
    if base_zscore
       bwin = [-0.75 -0.25];
       ib   = time>=bwin(1) & time<=bwin(2);
       %stddat = std(dat(:,:,ib,:),[],3);
       stddat = repmat(mean(std(dat(:,:,ib,:),[],1),3),[size(dat,1) 1 1 1]);
       dat  = bsxfun(@minus,dat,mean(dat(:,:,ib,:),3));
       dat  = bsxfun(@rdivide,dat,stddat);
       clear stddat
    end
    if avg_freq
        freqwin = [70 150];
        freqwin = [50 150];
        ifr     = freq>=freqwin(1) & freq<=freqwin(2);
        dat     = mean(dat(:,:,:,ifr),4);
        freq    = 1;
        nfreq   = 1;
    end
    %% do stim identity ANOVA on each channel
    do_stim_anova = 1;
    if do_stim_anova
        % get data
        ichan   = 1:size(dat,2);
        fprintf('\n Running Channelwise ANOVA:');  
        itrl    = ~ispractice; %use all trials for this
        nitrl   = nnz(itrl);
        tmpdat  = dat(itrl,ichan,:,:);
        tmpstm  = stims(itrl,1);
        tmprep  = tnum(itrl,1)>1;
        
        grp_col = {targtype(itrl,1) stims(itrl) tnum(itrl)>1 relrule(itrl,1) (rt(itrl,1)>median(rt(itrl)))};
        ncond   = length(grp_col);
        tic
        w2 = struct;
        w2.cond   = {'Targ','Stim','Rept','Rule', 'RTmed'};
        w2.time   = time;
        w2.label  = datatfr.label;
        w2.dimord = 'rpt_chan_time';
        groups    = w2.cond;
        w2.trial  = massANOVA(tmpdat,grp_col);
        w2.design = grp_col;
        
        %get randomized values
        nperm = 500;
        w2.nperm = nperm;
        b = '';
        for iperm = 1:nperm
           m = sprintf(' permutation %d/%d', iperm,nperm);
           fprintf([b m]); b = repmat('\b',[1 length(m)]);
           for icond = 1:ncond
            grp_col{icond} = grp_col{icond}(randperm(nitrl));
           end
           w2.boot(:,:,:,iperm)  = massANOVA(tmpdat,grp_col);
        end
        anovatime = toc;
        fprintf(' - done (it took %1.1f secs).\n',anovatime);
        
        %get zscore
        w2.pval     = sum(bsxfun(@ge,w2.boot,w2.trial),4)./nperm;
        w2.zscore   = norminv(1-w2.pval,0,1);
        w2.bootmean = mean(w2.boot,4);
        w2.bootstd  = std(w2.boot,[],4);
        w2 = rmfield(w2,'boot');
        w2.zscore(isinf(w2.zscore)) = norminv(1-1/nperm/2,0,1);
    end
    
    do_save = 1;
    if do_save
        resultsdir = sprintf('%s/results/anova',studypath); if ~exist(resultsdir,'dir'), mkdir(resultsdir); end
        fresult    = sprintf('%s/%s_ANOVA_50-150Hz.mat',resultsdir,substrg);        
        %fresult = sprintf('%s_%s.mat',fresult,datestr(now,'yyyymmdd'));
        save(fresult,'w2');
    end    
end
disp('Done');
%%
close all
clc
SetupPlotting();
figure
colormap(cbrewer('seq','Reds',64))
%colormap(cbrewer('div','BuRd',64))
imagesc(w2.time,[],squeeze(w2.zscore(4,:,:)),[0 3])
colorbar

nchan = size(w2.zscore,2);
ylims = [0.5 nchan+0.5];
hold on
plot([0 0],ylims,'k-')
plot([0 0]+1.1,ylims,'k-')
%%
close all
clc
figure

d = datatfr.powspctrm;
t = datatfr.time;
f = datatfr.freq;
for ifr = 1:length(f), flab{ifr} = sprintf('%d',f(ifr)); end
ib = t>=-0.75&t<=-0.25;
d = bsxfun(@minus,d,mean(d(:,:,:,ib),4));
d = bsxfun(@rdivide,d,mean(std(d(:,:,:,ib),[],1),4));
for ichan =1:61
%for ichan = 80
    subplot(8,8,ichan)
    %imagesc(time,[],squeeze(mean(mean(d(:,ichan,:,:),1),2)),[-1 1]/4)
    imagesc(time,[],squeeze(w2.trial(1,ichan,:,:)))
    %colorbar
    hold on
    axis xy
    set(gca,'ytick',[1:5:35 35],'yticklabel',flab([1:5:35 35]))
    set(gca,'xtick',[],'ytick',[]);
    plot([0 0],[0.5 f(end)+0.5],'k-')
    plot([0 0]+1.1,[0.5 f(end)+0.5],'k-');
    title(sprintf('%s %d',datatfr.label{ichan},ichan))
end
%%
do_save = 0;
if do_save
    resultsdir = sprintf('%s/results',studypath); if ~exist(resultsdir,'dir'), mkdir(resultsdir); end
    fresult    = sprintf('%s/DV_Decoding_PreTargBaseline_ERP_Pilot_09112017.mat',resultsdir);
    save(fresult,'time','Cproj');
end
%% basic plotting
close all
clc
figure
set(gcf,'color','white');
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
%set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPosition',[0 0 8 2]);
for iplot = 1%:6
    %subplot(2,3,iplot)
    set(gca,'fontname','Helvetica')
    if iplot == 1
        plotsubs = [1:7];
    else
        plotsubs = iplot-1;
    end
%plotsubs  = 1:3;
plotrules = [1 2 ];
plotstat  = 1;
plotmeas  = 1;
C = -Cproj(:,:,:,:,plotmeas);
if plotmeas == 1    
    C(:,:,:,1,:) = C(:,:,:,1) * 100;
else
    if plotstat == 3
    	C(:,:,:,1,:) = C(:,:,:,1) * 100;
    else
        C(:,:,:,1,:) = C(:,:,:,1) * 25000;
    end
end
C(:,:,:,3,:) = C(:,:,:,1,:)./C(:,:,:,2,:);

do_filt = 1;
if do_filt    
    ftype   = 'gaussian';
    fsize   = 0.016*fsample;
    C = filtfast(C,2,[],ftype,fsize);
    smoothsuffix = sprintf('_%s%dms',ftype,round(fsize/fsample*1000));
else
    smoothsuffix = '';
end

plotcolor = [0.2 0.8 0.2;
             0.8 0.2 0.2];
         
%plotcolor = [0.1 0.1 0.1];
plotcolor = [0.1 0.1 0.8;
             0.7 0.7 0.7;
             0   0   0   ];

%plotcolor = [0.8 0.2 0.2;
%             0.6 0.3 0.3];

hold on
ylims = [-.50 1.0]*75;
xlims = [-0.2 1.0];
plot([0 0],ylims,'k-','color',[1 1 1]*0.65);
plot([0.8 0.8],ylims,'k-','color',[1 1 1]*0.65);
plot(xlims,[0 0],'k-','color',[1 1 1]*0.65);
%plot([1 1]*.8,ylims,'k--','color',[1 1 1]*0.65);
cfg = [];
if iplot > 1
    cfg.linewidth = 2;
end
plotdata = C(plotsubs,:,plotrules,plotstat);
hplot = plotpatch(plotdata,time,plotcolor,cfg);
%if iplot == 4, legend(hplot,{'Active Rule' 'Inactive Rule'},'location','northeast'); legend boxoff; end
%set(gca,'ytick',sort([0 ylims]))
set(gca,'xtick',[-1:.5:.5 .8 1:.5:2])
set(gca,'tickdir','out','box','off')
set(gca,'xticklabel',{'' '' 'Cue' '.5' 'Target' '1'})
xlabel('Time (s)')
if iplot == 1, ylabel('Decoding (a.u.)'); end
xlim(xlims)
ylim(ylims)
if iplot == 1
    title(sprintf('Group (N=%d)',length(plotsubs)));
else
    title(sprintf('P%d',iplot-1));
end


plot_sig = 1;
if size(plotdata,3)>1,
    plot_diff = 0;
else
    plot_diff = false;    
end
if plot_sig
    if plot_diff
    plotdata(:,:,3) = diff(plotdata,[],3); end
    [h p c s] = ttest(plotdata);
    nsig = size(plotdata,3);
    if plot_diff
    plotdata(:,:,3) = []; end
    for isig = 1:nsig
       barcol = plotcolor(isig,:);
       pthresh = 0.05;
       yheight = min(ylims) + range(ylims)*0.025*(isig+1);
       b      = bwconncomp(p(:,:,isig)<pthresh);
       for ib = 1:b.NumObjects
           currid = b.PixelIdxList{ib};
           if nnz(currid)>5
            plot(time(currid),ones(nnz(currid),1)*yheight,'k-','linewidth',4,'color',barcol);
           end
       end
    end
end

end

do_print = 0;
if do_print
   figpath =  [studypath '/figures']; if ~exist(figpath,'dir'), mkdir(figpath); end
   fname = sprintf('%s/StimDecoding_PreTargBaseline_30Jan2018',figpath);
   print(gcf,'-dpng',fname);
   print(gcf,'-depsc2',fname);
end
%%
plotreg  = [1 2];
plotsubs = [1:30];
plotsess = [1 2];
%plotbeta = squeeze(mean(betamat(plotsubs,:,plotreg,plotsess),4));
plotbeta = squeeze(mean(atanh(groupcorr(plotsubs,:,plotreg,plotsess)),4));
ib = time>=-0.200 & time <=-0.05;
%plotbeta = bsxfun(@minus,plotbeta,mean(plotbeta(:,ib,:),2));

do_filt = 0;
if do_filt    
    ftype   = 'gaussian';
    fsize   = 0.012*fsample;
    plotbeta = filtfast(plotbeta,2,[],ftype,fsize);
    smoothsuffix = sprintf('_%s%dms',ftype,round(fsize/fsample*1000));
else
    smoothsuffix = '';
end

close all
clc
figure
set(gcf,'color','white')

ylims = [-1 +1]*.1;
xlims = time([1 end]);
%xlims = [-.1 0.5];

hold on
plotcol = [0.1 0.1 0.8;
           0.1 0.8 0.1];

plotcol = cat(1,plotcol,[0 0 0]);
%plotcol = colormap(jet(length(plotelem)));
       
plot_sig = 1;
if size(plotbeta,3)>1,
    plot_diff = 0;
else
    plot_diff = false;    
end
if plot_sig
    if plot_diff
    plotbeta(:,:,3) = diff(plotbeta,[],3); end
    [h p c s] = ttest(plotbeta);
    nsig = size(plotbeta,3);
    if plot_diff
    plotbeta(:,:,3) = []; end
    for isig = 1:nsig
       barcol = plotcol(isig,:);
       pthresh = 0.05;
       yheight = min(ylims) + range(ylims)*0.025*(isig+1);
       b      = bwconncomp(p(:,:,isig)<pthresh);
       for ib = 1:b.NumObjects
           currid = b.PixelIdxList{ib};
           plot(time(currid),ones(nnz(currid),1)*yheight,'k-','linewidth',4,'color',barcol);
       end
    end
end


hplot = plotpatch(plotbeta,time,plotcol);
legend(hplot,{'Abs Probe Offs.' 'Sign'},'location','northwest'), legend boxoff
hold on
plot([0 0],ylims,'k-')
plot(xlims,[0 0],'k-')
xlim(xlims)
ylim(ylims)
%set(gca,'tickdir','out','box','off')

do_print = 0;
if do_print
   figpath =  [studypath '/figures']; if ~exist(figpath,'dir'), mkdir(figpath); end
   fname = sprintf('%s/PV_Probe_2AFC_IEMCorrelationAbsSign%s_N%02d_7Aug2017',figpath,smoothsuffix,nsubs);
   print(gcf,'-dpng',fname);
   print(gcf,'-depsc2',fname);
end
%%
clc
itcorr  = time>=0 & time <=1;
corrtim = time(itcorr);
corrdat = squeeze(plotbeta(:,itcorr,2));
nperm   = 1000;
[pcorr] = ClusterCorrection2(corrdat,nperm,0.05);
pclustu = unique(pcorr);
npclust = nnz(pclustu < 0.5);
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');
for ipclust = 1:npclust
   currind = pcorr == pclustu(ipclust);
   currt   = corrtim(currind);
   fprintf('Cluster %02d (%02d timepoints, %d ms): p=%g, %d-%dms.\n',ipclust,nnz(currind),round(range(currt)*1000),pclustu(ipclust),round(min(currt)*1000),round(max(currt)*1000));
end
%%
clc
cfg = [];
cfg.channel = 1:60;
if ~exist('dataavg')
    dataavg = ft_timelockanalysis(cfg, data);
end

if length(dataavg.time)>length(time)
    dataavg.avg = dataavg.avg(:,1:ntimes);
    dataavg.var = dataavg.var(:,1:ntimes);
    dataavg.dof = dataavg.dof(:,1:ntimes);
    dataavg.time = time;
end

xlims = [-0.2 1.0];
ylims = [-5 +5];

plotitems = [1];
plotsubs = [1:30];
plotsess = 1:2;
plotelem = [2];
plottask = 1:2;        

plotbeta = squeeze(mean(mean(betamat(plotsubs,:,:,plotelem,plottask,plotitems,2,plotsess),8),9));

%
close all
figure
set(gcf,'color','white')
hold on

interactive_mode  = 0;
plotallconditions = 1;
plotallconditionsplusdiff = 1;
if interactive_mode
    [h p c s] = ttest(plotbeta);
    cfg = [];
    plot_t = 1;
    if plot_t
        plotdat = squeeze(s.tstat);
        cfg.clim = [-5 5];
    else
        plotdat = squeeze(mean(plotbeta,1));
        cfg.zlim = 'maxabs';
    end
    
    cfg.layout = 'EEG1010.lay';
    
    nplots = size(plotdat,3);
    plotstruct = cell(1,nplots);
    for iplot = 1:nplots    
    plotstruct{iplot}     = dataavg;
    plotstruct{iplot}.avg = squeeze(plotdat(:,:,iplot));
    end
    
    ft_multiplotER(cfg,plotstruct{:});
elseif plotallconditions
    [h p c s] = ttest(plotbeta);
    cfg = [];
    plot_t = 0;
    if plot_t
        plotbeta = squeeze(s.tstat);
        cfg.zlim = [-5 5];
    else
        plotbeta = squeeze(mean(plotbeta,1));
        cfg.zlim = [-2 2];
    end
    
    plotstruct = cell(1,2);
    plotstruct{1}     = dataavg;
    plotstruct{1}.avg = squeeze(plotbeta(:,:,1));
    plotstruct{2}     = dataavg;
    plotstruct{2}.avg = squeeze(plotbeta(:,:,2));
    
    plotnames = {'ForcedChoice','FreeRecall'};
    for iplot = 1:2
        hold on
        subplot(1,2,iplot)
        cfg = [];    
        %cfg.layout = 'EEG1010.lay'; 
        cfg.layout = 'elec1010.lay'; 
        cfg.interactive = 'yes';
        cfg.interactive = 'no';
        cfg.xlim     = [.36 .415];        
        %cfg.xlim     = [.10 .20];
        cfg.xlim     = [.3 0.650];
        cfg.comment  = 'no';
        cfg.colorbar = 'no';
        ft_topoplotER(cfg,plotstruct{iplot});
        title(sprintf('%s\n%d-%dms',plotnames{iplot},round(cfg.xlim*1000)));
    end
else
    plotstruct = cell(1,1);
    plotstruct{1}     = dataavg;
    plotstat = squeeze(diff(pltobeta,[],3));
    [h p c s] = ttest(plotstat);
    show_t_map = 0;
    if show_t_map
        plotstruct{1}.avg = squeeze(s.tstat);
        tflag = 'Tstat';
    else
        plotstruct{1}.avg = squeeze(mean(plotstat,1));
        tflag = 'Mean';
    end
    
    cfg = [];    
    cfg.layout = 'EEG1010.lay'; 
    cfg.interactive = 'yes';
    cfg.interactive = 'no';
    cfg.xlim = [2.25 2.95];
     if show_t_map
    cfg.zlim = [-5 5];
     else
    cfg.zlim = 'maxabs';
     end
    cfg.colorbar = 'yes';
    cfg.comment = 'no';
    ft_topoplotER(cfg,plotstruct{1});
    title(sprintf('Free Recall - Forced Choice\n%d-%dms',round(cfg.xlim*1000)));
end

do_print = 0;
if do_print
   figpath =  [studypath '/figures']; if ~exist(figpath,'dir'), mkdir(figpath); end
   fname = sprintf('%s/NoiseAnalysis_TaskComparison_Topography_N%02d_%d-%dms',figpath,nsubs,round(cfg.xlim*1000));
   print(gcf,'-dpng',fname);
   print(gcf,'-depsc2','-painters',fname);
end
%% single subject EMG plotting - summary
clc
close all
figure
ylims = [-10 10];
xlims = [-.1 1];
rtrange = [.3 .4; .4 .5; .5 .6];
rtrange = [.2 .25; .25 .3; .3 .35; .35 .4; .4 .45; .45 .5; .5 .6; .6 .7];

for irt = 1:size(rtrange,1)
    for pltrsp = 1:2
        subplot(size(rtrange,1),2,pltrsp + 2*(irt-1))

        plottrial = rsp==pltrsp & rt>=rtrange(irt,1) & rt<=rtrange(irt,2);

        d = squeeze(mean(dataerp_elem(plottrial,62:63,:),1));
        %d = zscore(d,[],2);
        hp = plot(time,d);
        %legend(hp,data.label(64:65)), legend boxoff
        hold on
        r = quantile(rt(plottrial),[.25 .5 .75]);
        for ir = 1:length(r)
            if ir == 2
                plot([1 1]*r(ir),ylims,'k-');
            else
                plot([1 1]*r(ir),ylims,'k--');
            end
        end
        xlim(xlims)
        ylim(ylims)
        resplabel = {'L' 'R'};
        if length(plottrial)==1
            title(sprintf('Resp: %s',resplabel{rsp(plottrial)}))
        else
            title(sprintf('Resp: %s, %dms',resplabel{pltrsp},round(r(2)*1000)))
        end
    end
end
%% single trial EMG plotting - single subject
clc
close all
figure
ylims = [-30 30];
xlims = [-.1 1];

plottrial = 1215;

d = squeeze(mean(dataerp_elem(plottrial,62:63,:),1));
hp = plot(time,d);
legend(hp,data.label(64:65)), legend boxoff
hold on
r = quantile(rt(plottrial),[.25 .5 .75]);
for ir = 2
    plot([1 1]*r(ir),ylims,'k-','linewidth',3);                
end
xlim(xlims)
ylim(ylims)
resplabel = {'L' 'R'};
title(sprintf('Resp: %s, %dms',resplabel{rsp(plottrial)},round(r(2)*1000)));