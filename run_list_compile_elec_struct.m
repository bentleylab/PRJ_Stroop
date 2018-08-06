SBJs = {};%,
for s = 1:numel(SBJs)
    fprintf('============================================================================\n');
    fprintf(['\t' SBJs{s} ' pat\n']);
    fprintf('============================================================================\n');
    fn_compile_elec_struct(SBJs{s},'main_ft','pat','');
    fprintf('============================================================================\n');
    fprintf(['\t' SBJs{s} ' mni_v\n']);
    fprintf('============================================================================\n');
    fn_compile_elec_struct(SBJs{s},'main_ft','mni','v');
    fprintf('============================================================================\n');
    fprintf(['\t' SBJs{s} ' mni_s\n']);
    fprintf('============================================================================\n');
    fn_compile_elec_struct(SBJs{s},'main_ft','mni','s');
end


% IN PROGRESS:
TBD_SBJs = {'IR26','IR32','IR37'};
% 26- needs all processing
% 32- figure out how to deal with interhemispheric L/R (same chanpos, so deleted as duplicates in ft_appendsens
% 37- needs elec update from Ricahrd

% DONE:
done_SBJs = {'IR21','IR27','IR31','IR35','IR39','IR41','IR48','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};
