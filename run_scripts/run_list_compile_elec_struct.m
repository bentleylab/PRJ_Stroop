SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};
for s = 1:numel(SBJs)
    fprintf('============================================================================\n');
    fprintf(['\t' SBJs{s} ' pat\n']);
    fprintf('============================================================================\n');
    fn_compile_elec_struct(SBJs{s},'main_ft','pat','',0);
%     fn_compile_elec_struct(SBJs{s},'main_ft','pat','',1);
%     fprintf('============================================================================\n');
%     fprintf(['\t' SBJs{s} ' mni_v\n']);
%     fprintf('============================================================================\n');
%     fn_compile_elec_struct(SBJs{s},'main_ft','mni','v');
%     fprintf('============================================================================\n');
%     fprintf(['\t' SBJs{s} ' mni_s\n']);
%     fprintf('============================================================================\n');
%     fn_compile_elec_struct(SBJs{s},'main_ft','mni','s');
end


% IN PROGRESS:
TBD_SBJs = {'IR26'};
% 26- needs all processing
% 37- needs elec update from Ricahrd
% never: IR37 (lesion)

% DONE:
done_SBJs = {'IR21','IR27','IR31','IR32','IR35','IR39','IR41','IR48','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};
