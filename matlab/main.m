function main(cfg)
% MAIN  Entry point for Myelin–FC Coupling analysis (MATLAB)

% addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
% cfg = config.default_config();

if nargin < 1
    cfg = config.default_config();
else
    % window for user to set FC modality
    cfg_default = config.default_config();
    fn = fieldnames(cfg_default);
    for i = 1:numel(fn)
        if ~isfield(cfg, fn{i}) || isempty(cfg.(fn{i}))
            cfg.(fn{i}) = cfg_default.(fn{i});
        end
    end
end

% FC modality-specific output directory
fc_subdir = util.sanitize_label(cfg.fc_label);
cfg.out_dir = fullfile(cfg.out_dir, fc_subdir);
util.ensure_dir(cfg.out_dir);

fprintf('\n=== Running modality: %s ===\n', cfg.fc_label);
fprintf('Output dir: %s\n', cfg.out_dir);


% 1) --- Load ---
[edges, nodes] = io.load_inputs(cfg);


% 2) --- Build design tables ---
D_global    = prep.build_design_table(edges, nodes, 'global', cfg);
D_rsn_pairs = prep.build_design_table(edges, nodes, 'rsn_pairs', cfg);
D_nodewise  = prep.build_design_table(edges, nodes, 'nodewise', cfg);

% Bundle design tables by name
D_all = struct();
D_all.global    = D_global;
D_all.rsn_pairs = D_rsn_pairs;
D_all.nodewise  = D_nodewise;

% Build lists for main and binned analyses based on cfg
D_list_main   = cellfun(@(nm) D_all.(nm), cfg.levels_main,   'UniformOutput', false);
D_list_binned = cellfun(@(nm) D_all.(nm), cfg.levels_binned, 'UniformOutput', false);


% 3) --- Models: with & without interactions (only for levels_main) ---
res_full    = models.fit_models(D_list_main,   cfg, true);
res_reduced = models.fit_models(D_list_main,   cfg, false);


% 4) --- Myelin binned analyses (only for levels_binned) ---
res_bins = models.fit_models_binned(D_list_binned, cfg);
bin_corr = metrics.compute_bin_correlations(D_list_binned, cfg);


% 5) --- Save ---
util.ensure_dir(cfg.out_dir);
save(fullfile(cfg.out_dir,'results_all.mat'), 'res_full','res_reduced','res_bins','bin_corr','-v7');

% full and reduced model results
util.save_struct_by_level(res_full,    cfg.levels_main,   cfg.out_dir, 'full');
util.save_struct_by_level(res_reduced, cfg.levels_main,   cfg.out_dir, 'reduced');

% binned model results
util.save_struct_by_level(res_bins,    cfg.levels_binned, cfg.out_dir, 'bins');

% binned correlations
util.save_struct_by_level(bin_corr,    cfg.levels_binned, cfg.out_dir, 'bin_corr');

fprintf('\nDone. Results in %s\n', cfg.out_dir);
end
