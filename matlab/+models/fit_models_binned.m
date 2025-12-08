function res = fit_models_binned(D_list, cfg)
% FIT_MODELS_BINNED  Run models within myelin bins (quantiles by default)

res = struct();
for li = 1:numel(D_list)
    D = D_list{li};

    levelName = util.get_level_name(D, 'fit_models');

    my = D.(cfg.col_myelin);
    nb = cfg.myelin_num_bins;
    edges = quantile(my, linspace(0,1,nb+1));
    % protect uniqueness
    edges(1) = -inf; edges(end) = inf;
    [~,bin] = histc(my, edges);

    rows = [];

    for b = 1:nb
        sub = D(bin==b, :);
        if isempty(sub); continue; end

        % Run model on this subset (single level only)
        S = models.fit_models({sub}, cfg, false);   % without interactions
        fns = fieldnames(S);
        assert(numel(fns) == 1, 'Expected a single level from fit_models in binned analysis.');

        T = S.(fns{1});
        T.bin = repmat(b, height(T), 1);
        rows = [rows; T]; %#ok<AGROW>
    end
    res.(levelName) = rows;
end
end
