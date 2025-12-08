function res = compute_bin_correlations(D_list, cfg)
% COMPUTE_BIN_CORRELATIONS  Pearson correlations within myelin bins

res = struct();
for li = 1:numel(D_list)
    D = D_list{li};

    levelName = util.get_level_name(D, 'fit_models');

    my = D.(cfg.col_myelin);
    nb = cfg.myelin_num_bins;
    edges = quantile(my, linspace(0,1,nb+1));
    edges(1) = -inf; edges(end) = inf;
    [~,bin] = histc(my, edges);
    rows = [];

    for b = 1:nb
        sub = D(bin==b, :);

        if isempty(sub); continue; end

        R = table;
        R.bin = b;

        vars = {cfg.col_FC,cfg.col_cal,cfg.col_myelin,cfg.col_len};

        for v = 2:numel(vars)                           % correlate FC with each predictor
            x = sub.(vars{v}); y = sub.(cfg.col_FC);
            ok = isfinite(x) & isfinite(y);

            rName = [vars{v} '_r'];                     % char, safe for table var name
            nName = [vars{v} '_n'];

            if any(ok)
                R.(rName) = corr(x(ok), y(ok), 'Type','Pearson');
                R.(nName) = sum(ok);
            else
                R.(rName) = NaN;
                R.(nName) = 0;
            end
        end
        rows = [rows; R];   %#ok<AGROW>
    end
    res.(levelName) = rows;
end
end
