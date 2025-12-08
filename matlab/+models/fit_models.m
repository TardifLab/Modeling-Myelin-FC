function res = fit_models(D_list, cfg, with_interactions)
% FIT_MODELS  Run regressions at global, rsn_pairs, nodewise

nLevels     = numel(D_list);
res         = struct();
for li = 1:nLevels
    D = D_list{li};

    levelName = util.get_level_name(D, 'fit_models');

  % Fit model
    T = aggregate_and_fit(D, cfg, with_interactions);
    res.(levelName) = T;
end
end

function T = aggregate_and_fit(D, cfg, with_interactions)
keys = unique(D.key);
rows = [];
for k = 1:numel(keys)
    sub = D(D.key==keys(k), :);

    X = sub(:, {cfg.col_cal,cfg.col_myelin,cfg.col_len});
    y = sub.(cfg.col_FC);


    % Optional: standardize dependent variable (FC) per group
    if cfg.standardize_y
        mu_y = mean(y,'omitnan');
        sd_y = std(y,'omitnan');
        % protect against zero-variance edge case
        if sd_y > 0
            y = (y - mu_y) ./ sd_y;
        else
            % if all FC identical in this group, leave as-is
            % (model would be degenerate anyway)
        end
    end

    % Standardize predictors per group
    if cfg.standardize_predictors
        varNames = X.Properties.VariableNames;
        for ii = 1:numel(varNames)
            nm = varNames{ii};  % this is a char row vector
            X.(nm) = (X.(nm) - mean(X.(nm),'omitnan')) ./ std(X.(nm),'omitnan');
        end
    end

    % Build interaction columns explicitly (for dominance analysis)
    if with_interactions
        X.int_my_cal = X.(cfg.col_myelin) .* X.(cfg.col_cal);
        X.int_my_len = X.(cfg.col_myelin) .* X.(cfg.col_len);
    end

    % Model specification
    if with_interactions
        spec = sprintf('%s ~ 1 + %s + %s + %s + int_my_cal + int_my_len', ...
            cfg.col_FC, cfg.col_cal, cfg.col_myelin, cfg.col_len);
    else
        spec = sprintf('%s ~ 1 + %s + %s + %s', cfg.col_FC, cfg.col_cal, cfg.col_myelin, cfg.col_len);
    end
    tbl = [table(y,'VariableNames',{cfg.col_FC}) X];
    M = fitlm(tbl, spec);

    % Dominance analysis (custom method; supports interactions via mdlspec)
    % Map predictors to x1..xN for dominance_interaction
    Xmat = [X.(cfg.col_cal), X.(cfg.col_myelin), X.(cfg.col_len)];
    if with_interactions
        mdlspec = 'y~1+x1+x2+x3+x2*x1+x2*x3';
    else
        mdlspec = 'y~1+x1+x2+x3';
    end
    [tab_lo, tab_hi] = models.dominance_interaction(Xmat, y, mdlspec); %#ok<ASGLU>
    if with_interactions
        Tdom = tab_hi;
    else
        Tdom = tab_lo;
    end
% Extract summary
    r = table;
    r.key   = keys(k);
    r.N     = height(sub);
    r.R2    = M.Rsquared.Ordinary;
    r.R2adj = M.Rsquared.Adjusted;
    r.RMSE  = M.RMSE;
    % Coefficients
    C = M.Coefficients;
    for ii = 1:height(C)
        nm = matlab.lang.makeValidName(C.Properties.RowNames{ii});
        r.(nm + "_B") = C.Estimate(ii);
        r.(nm + "_p") = C.pValue(ii);
    end
    % Add dominance metrics (prefixed dom_ and domfrac_)
    for jj = 1:height(Tdom)
        nm = matlab.lang.makeValidName(Tdom.Variable_name(jj));
        r.("dom_" + nm)     = Tdom.Dominance_value(jj);
        r.("domfrac_" + nm) = Tdom.Dominance_percentage(jj);
    end
    rows = [rows; r]; %#ok<AGROW>
end
T = rows;
end
