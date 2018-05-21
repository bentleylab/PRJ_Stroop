function check_trial_info(SBJ, data_id)
% compare trial types found by A02 to _log.txt; print results

% Script Parameters
SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
log_dir  = strcat(SBJ_dir,'00_raw/');
proc_dir = strcat(SBJ_dir,'04_proc/');
log_filename  = strcat(log_dir,SBJ,'_strooptask_log.txt');
proc_filename = strcat(proc_dir,data_id,'.mat');

% Load data
load(proc_filename,'trial_info');
log = importdata(log_filename,'\t',5);
trial_type = log.textdata(6:end,5);
trial_type(trial_info.ignore_trials) = [];
block_type = log.textdata(6:end,6);
block_type(trial_info.ignore_trials) = [];
% dash_concat = repmat(let,size(trial_type));
cond_lab = {};
for ix = 1:length(trial_type)
    cond_lab{ix} = strcat(trial_type(ix), '-', block_type(ix));
end

% Convert labels to numbers
cond_num_txt = convert_trial_lab2num(cond_lab);

% Logic comparison
lab_compare = cond_num_txt==trial_info.condition_n;

% Print results
if sum(lab_compare)==length(lab_compare)
    fprintf('All %i trials match! Good to go! :)\n',length(lab_compare));
else
    mismatch_cnt = length(lab_compare)-sum(lab_compare);
    fprintf('WARNING!!! Mismatch on %i trials! Check event scripts!\n',mismatch_cnt);
    block_nums = unique(trial_info.block_n);
    for b_ix = 1:length(block_nums)
        block_type = unique(trial_info.blocktype(trial_info.block_n==b_ix));
        block_trials = trial_info.trialtype(trial_info.block_n==b_ix);
        num_inc = strfind(block_trials,'inc');
        num_neu = strfind(block_trials,'neu');
        num_con = strfind(block_trials,'con');
        fprintf('block %i (%s) %i total = %4.2d con, %4.2d neu, %4.2d inc\n',b_ix,block_type{1},length(block_trials),...
            sum([num_inc{:}])/length(block_trials), sum([num_neu{:}])/length(block_trials), sum([num_con{:}])/length(block_trials));
    end
end

% Printing both to compare
% both_str = {};
% for r = 1:size(both,1)
%     for c = 1:size(both,2)
%         both_str{r,c} = trial_info.condition_types(both(r,c));
%     end
% end

end