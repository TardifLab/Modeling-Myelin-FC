function D = build_design_table(edges, nodes, level, cfg) %#ok<INUSD>
% BUILD_DESIGN_TABLE  Table of predictors for a given level
switch lower(level)
    case 'global',    D = local_global(edges, cfg);
    case 'rsn_pairs', D = local_rsn_pairs(edges, cfg);
    case 'nodewise',  D = local_nodewise(edges, cfg);
    otherwise, error('Unknown level: %s', level);
end
end

function D = local_global(edges, cfg)
D       = edges(:, {cfg.col_FC,cfg.col_cal,cfg.col_myelin,cfg.col_len,cfg.col_rsn_i,cfg.col_rsn_j});
D.level = repmat("global", height(D), 1);
D.key   = ones(height(D),1);
end

function D = local_rsn_pairs(edges, cfg)
ri      = string(edges.(cfg.col_rsn_i)); rj = string(edges.(cfg.col_rsn_j));
pair    = arrayfun(@(a,b) strjoin(sort([a,b]),'–'), ri, rj);
if ~cfg.include_within_rsn_pairs, keep = ri~=rj; else, keep = true(size(pair)); end
D       = edges(keep, {cfg.col_FC,cfg.col_cal,cfg.col_myelin,cfg.col_len});
D.level = repmat("rsn_pairs", height(D), 1);
D.key   = string(pair(keep));
end

function D = local_nodewise(edges, cfg)
Ti      = edges(:, {cfg.col_i, cfg.col_FC,cfg.col_cal,cfg.col_myelin,cfg.col_len});
Tj      = edges(:, {cfg.col_j, cfg.col_FC,cfg.col_cal,cfg.col_myelin,cfg.col_len});
Ti.Properties.VariableNames{1} = 'node';
Tj.Properties.VariableNames{1} = 'node';
D       = [Ti; Tj];
D.level = repmat("nodewise", height(D), 1);
D.key   = D.node;
end
