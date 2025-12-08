function run_all_modalities()
% RUN_ALL_MODALITIES
% Convenience driver to run the pipeline for multiple FC modalities.
%
% Adapt the lists below to your setup. For each modality, you can:
%   - give it a label (used for subfolder under out/)
%   - point to an edges CSV
%   - optionally override which levels to run

% % ----- USER EDIT SECTION ----------------------------------------------
% fc_labels = {
%     'BOLD',       ...
%     'MEG_alpha',  ...
%     'MEG_beta'    ...
% };
% 
% edges_files = {
%     'edges_BOLD.csv',      ... % in cfg.data_dir
%     'edges_MEG_alpha.csv', ...
%     'edges_MEG_beta.csv'   ...
% };
% 
% % Optional: customize levels per modality (or leave [] to use defaults)
% levels_main_per_mod = {
%     {'global','rsn_pairs','nodewise'}, ...
%     {'global','rsn_pairs'},            ...
%     {'global'}                         ...
% };
% 
% levels_binned_per_mod = {
%     {'global','rsn_pairs','nodewise'}, ...
%     {'global','rsn_pairs'},            ...
%     {'global'}                         ...
% };
% ----- USER EDIT SECTION ----------------------------------------------

fc_labels = {
    'BOLDin'       
};

edges_files = {
    'edges_fc_BOLDin.csv'
};

% Optional: customize levels per modality (or leave [] to use defaults)
levels_main_per_mod = {
    {'global','rsn_pairs','nodewise'}
};

levels_binned_per_mod = {
    {'global'}                         
};

% ---------------------------------------------------------------------

% Sanity check
n_mod = numel(fc_labels);
assert(numel(edges_files) == n_mod, 'edges_files must match fc_labels in length.');
if ~isempty(levels_main_per_mod)
    assert(numel(levels_main_per_mod) == n_mod, 'levels_main_per_mod length mismatch.');
end
if ~isempty(levels_binned_per_mod)
    assert(numel(levels_binned_per_mod) == n_mod, 'levels_binned_per_mod length mismatch.');
end

% Ensure MATLAB sees the code
addpath(genpath('matlab'));

for k = 1:n_mod
    fprintf('\n=============================================\n');
    fprintf(' Modality %d / %d: %s\n', k, n_mod, fc_labels{k});
    fprintf('=============================================\n');

    % Start from defaults
    cfg = config.default_config();

    % Modality-specific overrides
    cfg.fc_label   = fc_labels{k};
    cfg.edges_file = edges_files{k};

    if ~isempty(levels_main_per_mod)
        cfg.levels_main = levels_main_per_mod{k};
    end
    if ~isempty(levels_binned_per_mod)
        cfg.levels_binned = levels_binned_per_mod{k};
    end

    % Run pipeline for this modality
    main(cfg);
end

fprintf('\nAll modalities finished.\n');
end
