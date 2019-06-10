% Parameters for HFA actvation vs. baseline
st.model_lab = 'actv';
st.groups    = {'actv'};    % for compatibility
st.trial_cond  = {'all'};
st.ep_lab   = 'D';
st.evnt_lab = 'S';
st.stat_lim = [-0.1 0];
st.lim_adj  = {'min(RT)', 'RT'};
st.cust_win = 1;            % custom windows per trials
st.min_rt   = 0;

st.actv_win = 0;  % minimum window to be different than baseline
st.alpha    = 0.05;
