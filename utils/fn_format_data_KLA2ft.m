function [data_ft] = fn_format_data_KLA2ft(data,hdr)
%% Save data in KLA data+header format to Fieldtrip format

data_ft.trial   = {data};
data_ft.time    = {[1:size(data,2)]./hdr.sample_rate};
data_ft.fsample = hdr.sample_rate;
data_ft.label   = hdr.channel_labels;

end