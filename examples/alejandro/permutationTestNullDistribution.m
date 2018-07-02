function [stats] = permutationTestNullDistribution(data,cfg)

%{
Function to compute the Permutation test vs null distribution for ECoG/SEEG
data matrices. This test allows to assess the degree of statistical
significance of samples against a null surrogate distribution.  
Surrogate distribution is build from shifted trials in time or phase(fft
phase)

Multiplecomparison problems are considered in a cluster level statistics of
each channel.

% need extra files from fieldtrip and spm to run
% addpath('/data/processing_codes/private/')

INPUT:
    data: Ch x samples x trials
    cfg.mode: Generate surrogates by shifting
           'phase' in the fft transform of every trial
           'time'  in every trial
    cfg.permutations: the number of permutations that one wants to perfom.
            By default, 5000.
    cfg.alpha: the significance threshold. By default, 0.05
    cfg.clusteralpha: the significance threshold to build clusters. 
            By default, 0.05
    cfg.neighbours : electrode neighbour infrmation or AKA:Adjacency matrix
            (ch x ch).
            Use eye(ch) if spatial clustering is not desired.
            [], all channels are considered individually - Not working by now. 
 
OUTPUT:
    stats.h: test decision for the null hypothesis. Ch x samples
    stats.pvalue: cluster coorrected pvalue. Ch x samples
    stats.stat:  cluster cummulative sum statistic. Ch x samples
    stats.clusInd: cluster index. Ch x samples 

Based on Maris & Oostenveld. Nonparametric statistical testing of 
EEG- and MEG-data. Journal of neuroscience methods. 2007


TESTCODE:
run permutationTestNullDistribution_test_script

Code with some vestigeus pieces from JesÃºs Monge-Ã?lvarez

A Blenkmann 2017

%}
%% Configuration parameters:
% Intializatons:
if ~isfield(cfg,'mode');                cfg.mode='time';                end
if ~isfield(cfg,'permutations');        cfg.permutations=5000;          end
if ~isfield(cfg,'alpha');               cfg.alpha = 0.05;               end
if ~isfield(cfg,'clusteralpha');        cfg.clusteralpha = 0.05;        end
if ~isfield(cfg,'neighbours');          cfg.neighbours = [];            end


% The number of permutation must allow to apply the significance threshold,
% that is to say, the number of permutations must be greater than a
% recommended minimum:
control = ~((1/cfg.permutations) > cfg.alpha);
assert(control,'The number of permutations is not valid for the selected significance threshold.');

[nCh,nS,~]=size(data); % nCh x nSamples x nTrials
h=zeros([nCh,nS]);
clusPvalue = ones([nCh,nS]);
clusTstat = zeros([nCh,nS]);
clusInd = zeros([nCh,nS]);

%% alocate some memory 

% We will save surrogates results here
statsSurrogateSup=zeros(nCh,nS);
hSurrogateSup=zeros(nCh,nS);
statsSurrogateInf=zeros(nCh,nS);
hSurrogateInf=zeros(nCh,nS);

% we will save cluster statistics values here
if ~isempty(cfg.neighbours)
    % all channels together
    sumClusterDistributionSup = zeros(1,cfg.permutations);
    sumClusterDistributionInf = zeros(1,cfg.permutations);
    numClusterDistributionSup = zeros(1,cfg.permutations);
    numClusterDistributionInf = zeros(1,cfg.permutations);
else
    % channel by channel clustering
    sumClusterDistributionSup = zeros(nCh,cfg.permutations);
    sumClusterDistributionInf = zeros(nCh,cfg.permutations);
    numClusterDistributionSup = zeros(nCh,cfg.permutations);
    numClusterDistributionInf = zeros(nCh,cfg.permutations);
end
%% build surrogate data
for i = 1:cfg.permutations
    parfor ch=1:nCh
        chData = squeeze(data(ch,:,:));  % nS x nT
        chData = chData-mean(chData(:)); % remove mean
    
        % surrogate data:
        chSurrData = buildSurrogate(chData,cfg.mode); %  nS x nT
        
        % build statistics for that surrogate (2 tails)
        [surrChHSup,~,~,surrChTstatSup] = ttest(chSurrData',0,'Alpha',cfg.clusteralpha/2,'Tail','right'); %  1 x nS
        [surrChHInf,~,~,surrChTstatInf] = ttest(chSurrData',0,'Alpha',cfg.clusteralpha/2,'Tail','left'); %  1 x nS
        
        % Get t statistics and H. Data is overwrite on each permutation 
        statsSurrogateSup(ch,:) = surrChTstatSup.tstat;
        statsSurrogateInf(ch,:) = surrChTstatInf.tstat;
        hSurrogateSup(ch,:) = surrChHSup;
        hSurrogateInf(ch,:) = surrChHInf;
    end
    
    %% build surrogate cluster distribution
    
    if ~isempty(cfg.neighbours)
    % clustering all channeles together
        % fildtrip-spm function
        [clusterSup, num] = findcluster(hSurrogateSup, cfg.neighbours, 1);
        tempSum=zeros(1,num);
        tempNum=zeros(1,num);
        parfor j=1:num
            % calculate cluster sum and number
            tempSum(j)=sum(statsSurrogateSup(clusterSup==j));
            tempNum(j)=numel(statsSurrogateSup(clusterSup==j));    
        end
        % keep only the maximum cluster sum and size
        sumClusterDistributionSup(i)=max([tempSum 0]);
        numClusterDistributionSup(i)=max([tempNum 0]);
        
        [clusterSup, num] = findcluster(hSurrogateInf, cfg.neighbours, 1);
        tempSum=zeros(1,num);
        tempNum=zeros(1,num);
        parfor j=1:num
            % calculate cluster sum and number
            tempSum(j)=sum(statsSurrogateInf(clusterSup==j));
            tempNum(j)=numel(statsSurrogateInf(clusterSup==j));    
        end
        % keep only the maximum cluster sum and size
        sumClusterDistributionInf(i)=min([tempSum 0]);
        numClusterDistributionInf(i)=max([tempNum 0]);
        
    else
    % channel by channel clustering
        
        parfor ch=1:nCh
            
            % build the cluster level statistics. Sum of tvalues in the bigest cluster
            
            % positive clusters - only max cluster is considered
            [cSum, cNumElem,~] = clusterSum(hSurrogateSup(ch,:),statsSurrogateSup(ch,:));
            
            sumClusterDistributionSup(ch,i) = max(cSum);
            numClusterDistributionSup(ch,i) = max(cNumElem);
            
            % negegative clusters - only min cluster is considered
            [cSum, cNumElem,~] = clusterSum(hSurrogateInf(ch,:),statsSurrogateInf(ch,:));
            sumClusterDistributionInf(ch,i) = min(cSum);
            numClusterDistributionInf(ch,i) = max(cNumElem);
        end
    end

    if sum(i == floor([1:cfg.permutations/10:cfg.permutations+1]-1))
        fprintf('%.0f out %.0f permutations done \n',i,cfg.permutations);
    end
end
clear statsSurrogateSup statsSurrogateInf surrChHSup surrChHInf
%% compare with observed measure
if ~isempty(cfg.neighbours)
    % clustering all channeles together
    statsSup=zeros(nCh,nS);
    statsInf=zeros(nCh,nS);
    hSup=zeros(nCh,nS);
    hInf=zeros(nCh,nS);
    parfor ch=1:nCh
        
        chData=squeeze(data(ch,:,:));  % nS x nT
        chData=chData-mean(chData(:)); % remove mean
       
        [hS,~,~,statsS] =ttest(chData',0,'Alpha',cfg.clusteralpha/2,'Tail','right'); %  1 x nS
        statsSup(ch,:) = statsS.tstat;
        hSup(ch,:) = hS;
        
        % build statistics for channel data - Left side
        [hI,~,~,statsI] =ttest(chData',0,'Alpha',cfg.clusteralpha/2,'Tail','left'); %  1 x nS
        statsInf(ch,:) = statsI.tstat;
        hInf(ch,:) = hI;
    end
    
    % clustering
    [clusterSup, numSup] = findcluster(hSup, cfg.neighbours, 1);
    cSumSup=zeros(1,numSup);
    cNumSup=zeros(1,numSup);
    parfor j=1:numSup
        % calculate cluster sum and number
        cSumSup(j)=sum(statsSup(clusterSup==j));
        cNumSup(j)=numel(statsSup(clusterSup==j));
    end
    
    [clusterInf, numInf] = findcluster(hInf, cfg.neighbours, 1);
    cSumInf=zeros(1,numInf);
    cNumInf=zeros(1,numInf);
    parfor j=1:numInf
        % calculate cluster sum and number
        cSumInf(j)=sum(statsInf(clusterInf==j));
        cNumInf(j)=numel(statsInf(clusterInf==j));
    end
    
    % thresholding
    thrClusSup=prctile(sumClusterDistributionSup,(1-cfg.alpha/2)*100);
    thrClusInf=prctile(sumClusterDistributionInf,(cfg.alpha/2)*100);

    % Reject clusters smaller than thrClusSup and thrClusNeg !!
    rejClusSup=cSumSup<thrClusSup;
    rejClusInf=cSumInf>thrClusInf;
    
    k=1;
    for i=1:length(cSumSup)
        if ~rejClusSup(i)
            clusPvalue(clusterSup==i) = (sum(cSumSup(i) < sumClusterDistributionSup)+1) / cfg.permutations; %proportion
            clusTstat(clusterSup==i)  = cSumSup(i);
            clusInd (clusterSup==i)   = k;
            h(clusterSup==i)          = 1;
            k=k+1;
        end
    end
    
    for i=1:length(cSumInf)
        if ~rejClusInf(i)
            clusPvalue(clusterInf==i) = (sum(cSumInf(i) > sumClusterDistributionInf)+1) / cfg.permutations; %proportion
            clusTstat (clusterInf==i) = cSumInf(i);
            clusInd (clusterInf==i)   = k;
            h(clusterInf==i)          = 1;
            k=k+1;
        end
    end
     
     
else
    % channel by channel clustering
    for ch=1:nCh
        
        chData=squeeze(data(ch,:,:));  % nS x nT
        chData=chData-mean(chData(:)); % remove mean
        
        %% statistics on the empirical data
        % get percentiles from cluster distribution
        thrClusSup=prctile(sumClusterDistributionSup(ch,:),(1-cfg.alpha/2)*100);
        thrClusInf=prctile(sumClusterDistributionInf(ch,:),(cfg.alpha/2)*100);
        
        %     hist(sumClusterDistributionSup(ch,:),20); hold on; hist(sumClusterDistributionInf(ch,:),20)
        %     line([thrClusInf,thrClusInf],[0 cfg.permutations]);
        %     line([thrClusSup,thrClusSup],[0 cfg.permutations]);
        
        % build statistics for channel data - Positive clusters
        [hSup,~,~,statsS] =ttest(chData',0,'Alpha',cfg.alpha/2,'Tail','right'); %  1 x nS
        tstat=statsS.tstat;
        % positive clusters
        [cSumSup, cNumElemSup,hIndSup]=clusterSum(hSup,tstat);
        
        % build statistics for channel data - Left side
        [hInf,~,~,statsI] =ttest(chData',0,'Alpha',cfg.alpha/2,'Tail','left'); %  1 x nS
        tstat=statsI.tstat;
        % positive clusters
        [cSumInf, cNumElemInf,hIndInf]=clusterSum(hInf,tstat);
        
        
        % Reject clusters smaller than thrClusSup and thrClusNeg !!
        rejClusSup=cSumSup<thrClusSup;
        rejClusInf=cSumInf>thrClusInf;
        
        k=1;
        for i=1:length(cSumSup)
            if ~rejClusSup(i)
                clusPvalue(ch,hIndSup{i}) = (sum(cSumSup(i) < sumClusterDistributionSup(ch,:))+1) / cfg.permutations; %proportion
                clusTstat(ch,hIndSup{i})  = cSumSup(i);
                clusInd (ch,hIndSup{i})   = k;
                h(ch,hIndSup{i})          = 1;
                k=k+1;
            end
        end
        
        for i=1:length(cSumInf)
            if ~rejClusInf(i)
                clusPvalue(ch,hIndInf{i}) = (sum(cSumInf(i) > sumClusterDistributionInf(ch,:))+1) / cfg.permutations; %proportion
                clusTstat (ch,hIndInf{i}) = cSumInf(i);
                clusInd (ch,hIndInf{i})   = k;
                h(ch,hIndInf{i})          = 1;
                k=k+1;
            end
        end
        
    end
end
stats.h=h;
stats.hSup=hSup;
stats.hInf=hInf;
stats.clusterSup=clusterSup;
stats.clusterInf=clusterInf;
stats.statsSup=statsSup;
stats.statsInf=statsInf;
stats.thrClusSup=thrClusSup;
stats.thrClusInf=thrClusInf;
stats.clusPvalue=clusPvalue;
stats.clusTstat=clusTstat;
stats.clusInd=clusInd;
