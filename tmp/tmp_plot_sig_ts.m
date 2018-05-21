%% Plot significance times within subject across electrodes
max_len=0;
for i = 1:length(sig_times)
    max_len=max([max_len sig_times{i}]);
end
% for i = 1:length(ir39.sig_times);max_len=max([max_len ir39.sig_times{i}]);end

sig_ts = zeros([1 max_len]);
for i = 1:length(sig_times)
    sig_ts(sig_times{i})=sig_ts(sig_times{i})+1;
end
% sig_ts39 = zeros([1 max_len]);
% for i = 1:length(ir39.sig_times);sig_ts39(ir39.sig_times{i})=sig_ts39(ir39.sig_times{i})+1;end
figure;hold on;
plot(sig_ts,'r');
% plot(sig_ts39,'b');