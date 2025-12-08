function save_struct_by_level(S, levels, out_dir, suffix)
% SAVE_STRUCT_BY_LEVEL  Save tables from a struct by level to CSV files.
%
% S      : struct with fields for each level (e.g., .global, .rsn_pairs)
% levels : cellstr of level names to consider
% out_dir: output directory
% suffix : char/string appended to filename, e.g., 'full', 'reduced'

for i = 1:numel(levels)
    lev = levels{i};
    if isfield(S, lev)
        fname = fullfile(out_dir, sprintf('%s_%s.csv', lev, suffix));
        writetable(S.(lev), fname);
    end
end
end
