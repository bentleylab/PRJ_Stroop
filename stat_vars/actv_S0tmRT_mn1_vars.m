% Parameters for HFA actvation vs. baseline
st.model_lab = 'actv';
st.groups    = {'actv'};    % for compatibility
st.evnt_lab = 'S';
st.stat_lim = [0 0];
st.lim_adj  = {'', 'min(RT)'};
st.cust_win = 0;            % custom windows per trials
st.actv_win = 0.1;  % minimum window to be different than baseline
st.alpha    = 0.05;
