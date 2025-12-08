function cfg = default_config()
% DEFAULT_CONFIG  Centralized parameters

cfg = struct();

% I/O
cfg.root      = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg.data_dir  = fullfile(cfg.root,'data');

% Base output directory (root for all FC modalities)
cfg.out_dir   = fullfile(cfg.root,'out');

% Label for the FC modality (used as subdirectory name under out_dir)
cfg.fc_label  = 'BOLDin';                   % e.g., 'BOLDin', 'MEG_alpha', 'MEG_beta'

cfg.edges_file= 'edges.csv';                % or .mat with variable 'edges'
cfg.nodes_file= 'nodes.csv';                % optional

% Columns
cfg.col_i      = 'i';
cfg.col_j      = 'j';
cfg.col_FC     = 'FC';
cfg.col_cal    = 'caliber';
cfg.col_myelin = 'myelin';
cfg.col_len    = 'length';
cfg.col_rsn_i  = 'rsn_i';
cfg.col_rsn_j  = 'rsn_j';

% RSN / grouping
cfg.rsn_name_list = {'Visual','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
cfg.include_within_rsn_pairs = true;

% Modeling
cfg.standardize_predictors = true;
cfg.standardize_y          = true;
cfg.robust = false; % fitlm RobustOpts

% Which levels to run in the main regression
% Options: 'global', 'rsn_pairs', 'nodewise'
cfg.levels_main   = {'global','rsn_pairs','nodewise'};

% Which levels to run in the myelin-binned analyses
cfg.levels_binned = {'global'};

% Binning
cfg.myelin_num_bins = 5;
cfg.myelin_bin_mode = 'quantile'; % or 'equalwidth'

rng(1,'twister');
end
