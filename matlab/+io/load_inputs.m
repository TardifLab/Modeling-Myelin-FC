function [edges, nodes] = load_inputs(cfg)
% LOAD_INPUTS  Load edges and optional nodes; ensure RSN columns on edges.

% Load edges
epath = fullfile(cfg.data_dir, cfg.edges_file);
[~,~,ext] = fileparts(epath);
switch lower(ext)
    case '.csv', edges = readtable(epath);
    case '.mat', S = load(epath); assert(isfield(S,'edges')); edges = S.edges;
    otherwise, error('Unsupported edges extension: %s',ext);
end

% Load nodes (optional)
nodes = table();
if ~isempty(cfg.nodes_file) && isfile(fullfile(cfg.data_dir,cfg.nodes_file))
    npath = fullfile(cfg.data_dir,cfg.nodes_file);
    [~,~,ext2] = fileparts(npath);
    switch lower(ext2)
        case '.csv', nodes = readtable(npath);
        case '.mat', S2 = load(npath); assert(isfield(S2,'nodes')); nodes = S2.nodes;
        otherwise, error('Unsupported nodes extension: %s',ext2);
    end
end

% Ensure RSN columns
has_i = ismember(cfg.col_rsn_i, edges.Properties.VariableNames);
has_j = ismember(cfg.col_rsn_j, edges.Properties.VariableNames);
if ~(has_i && has_j)
    assert(~isempty(nodes), 'RSN labels missing on edges; provide nodes.csv with node_id & rsn');
    assert(all(ismember({'node_id','rsn'}, nodes.Properties.VariableNames)), 'nodes: need node_id, rsn');
    T1 = nodes(:,{'node_id','rsn'}); T1.Properties.VariableNames = {cfg.col_i, cfg.col_rsn_i};
    T2 = nodes(:,{'node_id','rsn'}); T2.Properties.VariableNames = {cfg.col_j, cfg.col_rsn_j};
    edges = outerjoin(edges, T1, 'Keys', cfg.col_i, 'MergeKeys', true);
    edges = outerjoin(edges, T2, 'Keys', cfg.col_j, 'MergeKeys', true);
end

need = {cfg.col_i,cfg.col_j,cfg.col_FC,cfg.col_cal,cfg.col_myelin,cfg.col_len,cfg.col_rsn_i,cfg.col_rsn_j};
assert(all(ismember(need, edges.Properties.VariableNames)), 'Edges table missing columns');

end
