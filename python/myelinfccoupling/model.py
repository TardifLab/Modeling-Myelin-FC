import pandas as pd, numpy as np
import statsmodels.formula.api as smf
def _dominance_interaction(df: pd.DataFrame, with_interactions: bool):
    """
    Dominance analysis mirroring user's MATLAB method.
    Uses variables x1=caliber, x2=myelin, x3=length; interactions: x2*x1 and x2*x3 when enabled.
    Returns dicts (dom_lo, domfrac_lo, dom_hi, domfrac_hi).
    """
    import itertools
    from math import factorial
    # Build X matrix and drop NA rows
    cols = ['caliber','myelin','length','FC']
    D = df[cols].dropna().copy()
    X = D[['caliber','myelin','length']].rename(columns={'caliber':'x1','myelin':'x2','length':'x3'})
    y = D['FC']

    # Lower-order only dominance (tab_lo)
    # R2_full_lo: y ~ 1 + x1 + x2 + x3
    m = smf.ols('FC ~ 1 + x1 + x2 + x3', data=pd.concat([y, X], axis=1)).fit()
    R2_full_lo = m.rsquared
    # R2_indi_lo: each Xi alone
    R2_indi_lo = {}
    for k in ['x1','x2','x3']:
        m = smf.ols(f'FC ~ 1 + {k}', data=pd.concat([y, X[[k]]], axis=1)).fit()
        R2_indi_lo[k] = m.rsquared
    # R2_rest_lo: model excluding Xi
    R2_rest_lo = {}
    for k in ['x1','x2','x3']:
        others = [t for t in ['x1','x2','x3'] if t != k]
        m = smf.ols(f"FC ~ 1 + {' + '.join(others)}", data=pd.concat([y, X], axis=1)).fit()
        R2_rest_lo[k] = m.rsquared
    R2_marg_lo = {k: R2_full_lo - R2_rest_lo[k] for k in ['x1','x2','x3']}
    dom_lo = {k: (R2_marg_lo[k] + R2_indi_lo[k])/2 for k in ['x1','x2','x3']}
    s = sum(dom_lo.values()) or float('nan')
    domfrac_lo = {k: (dom_lo[k]/s if s==s else float('nan')) for k in dom_lo}

    # High-order with interactions (tab_hi)
    if with_interactions:
        # full model: y ~ 1 + x1 + x2 + x3 + x2*x1 + x2*x3
        data = pd.concat([y, X], axis=1)
        m = smf.ols('FC ~ 1 + x1 + x2 + x3 + x2*x1 + x2*x3', data=data).fit()
        R2_full_hi = m.rsquared

        # Identify interaction terms list and their variables
        inter_terms = [('x2:x1', ['x2','x1']), ('x2:x3', ['x2','x3'])]

        # Individual contributions
        R2_indi_hi = {}
        # lower-order terms: same as lo
        for k in ['x1','x2','x3']:
            R2_indi_hi[k] = R2_indi_lo[k]
        # interactions: difference between models with and without interaction among its terms
        for name, vars_ in inter_terms:
            df_sub = data[['FC'] + vars_]
            m_lo = smf.ols(f"FC ~ 1 + {' + '.join(vars_)}", data=df_sub).fit()
            m_hi = smf.ols(f"FC ~ 1 + {'*'.join(vars_)}", data=df_sub).fit()
            R2_indi_hi[name] = m_hi.rsquared - m_lo.rsquared

        # R2_rest_hi: remove Xi and all its interactions; for interaction term, remove only that term
        R2_rest_hi = {}
        # for lower-order predictors
        for k in ['x1','x2','x3']:
            # remaining terms: drop k and any inter containing k
            keep_inters = [it for it,_ in inter_terms if k not in it.split(':')]
            rhs_parts = [t for t in ['x1','x2','x3'] if t != k]
            rhs = ' + '.join(rhs_parts)
            if keep_inters:
                rhs += ' + ' + ' + '.join([it.replace(':','*') for it in keep_inters])
            m = smf.ols(f'FC ~ 1 + {rhs}', data=data).fit()
            R2_rest_hi[k] = m.rsquared
        # for interaction terms: remove just that interaction
        for name, vars_ in inter_terms:
            keep_inters = [it for it,_ in inter_terms if it != name]
            rhs = 'x1 + x2 + x3'
            if keep_inters:
                rhs += ' + ' + ' + '.join([it.replace(':','*') for it in keep_inters])
            m = smf.ols(f'FC ~ 1 + {rhs}', data=data).fit()
            R2_rest_hi[name] = m.rsquared

        R2_marg_hi = {k: R2_full_hi - R2_rest_hi[k] for k in list(R2_rest_hi.keys())}
        dom_hi = {k: (R2_marg_hi[k] + R2_indi_hi[k])/2 for k in R2_marg_hi}
        s = sum(dom_hi.values()) or float('nan')
        domfrac_hi = {k: (dom_hi[k]/s if s==s else float('nan')) for k in dom_hi}
    else:
        dom_hi = dom_lo
        domfrac_hi = domfrac_lo

    return dom_lo, domfrac_lo, dom_hi, domfrac_hi

def _dominance_lmg(df: pd.DataFrame, predictors, target='FC'):
    """LMG/Shapley dominance analysis for small p.
    Returns dict {name: (raw, frac)}
    """
    import itertools
    # Drop rows with NaNs in any used column
    cols = list(predictors) + [target]
    D = df.dropna(subset=cols)
    p = len(predictors)
    if p == 0:
        return {}
    # Precompute R2 for all subsets
    subsets = []
    R2 = {}
    for k in range(p+1):
        for comb in itertools.combinations(range(p), k):
            subsets.append(comb)
            if k == 0:
                R2[comb] = 0.0
            else:
                colsX = [predictors[i] for i in comb]
                f = target + ' ~ 1 + ' + ' + '.join(colsX)
                m = smf.ols(f, data=D).fit()
                R2[comb] = m.rsquared
    # Shapley values
    import math
    shap = np.zeros(p)
    for j in range(p):
        for comb in subsets:
            if j not in comb:
                withj = tuple(sorted(list(comb) + [j]))
                s = len(comb)
                w = math.factorial(s) * math.factorial(p - s - 1) / math.factorial(p)
                shap[j] += w * (R2[withj] - R2[comb])
    R2_full = R2[tuple(range(p))]
    dom_raw = shap
    dom_frac = dom_raw / dom_raw.sum() if dom_raw.sum() > 0 else np.full_like(dom_raw, np.nan)
    return {predictors[i]: (float(dom_raw[i]), float(dom_frac[i])) for i in range(p)}


def _standardize(df, cols):
    return (df[cols] - df[cols].mean())/df[cols].std(ddof=0)

def _fit(df, with_interactions=True, standardize=True):
    cols = dict(FC='FC', cal='caliber', my='myelin', ln='length')
    if standardize:
        df[['caliber','myelin','length']] = _standardize(df, ['caliber','myelin','length'])
    if with_interactions:
        formula = f"{cols['FC']} ~ 1 + {cols['my']}*{cols['cal']} + {cols['my']}*{cols['ln']}"
    else:
        formula = f"{cols['FC']} ~ 1 + {cols['cal']} + {cols['my']} + {cols['ln']}"
    m = smf.ols(formula, data=df).fit()
    dom_lo, domfrac_lo, dom_hi, domfrac_hi = _dominance_interaction(df, with_interactions=with_interactions)
    dom = dom_hi if with_interactions else dom_lo
    domfrac = domfrac_hi if with_interactions else domfrac_lo
    out = {'N': len(df), 'R2': m.rsquared, 'R2adj': m.rsquared_adj, 'AIC': m.aic, 'BIC': m.bic}
    for name, coef in m.params.items():
        out[f"{name}_B"] = coef
        out[f"{name}_p"] = m.pvalues.get(name, np.nan)
    name_map = {'x1':'caliber','x2':'myelin','x3':'length','x2:x1':'int_my_cal','x2:x3':'int_my_len'}
    for k,v in dom.items():
        out[f"dom_{name_map.get(k,k)}"] = v
    for k,v in domfrac.items():
        out[f"domfrac_{name_map.get(k,k)}"] = v
    return pd.Series(out)

def run(edges: pd.DataFrame, level: str, with_interactions=True):
    rows = []
    if level=='global':
        rows.append(_fit(edges, with_interactions))
    elif level=='rsn_pairs':
        pair = pd.DataFrame({'p': np.where(edges['rsn_i']<edges['rsn_j'],
                                           edges['rsn_i']+'–'+edges['rsn_j'],
                                           edges['rsn_j']+'–'+edges['rsn_i'])})
        edges2 = edges.join(pair)
        for k, sub in edges2.groupby('p'):
            s = _fit(sub, with_interactions); s['key'] = k; rows.append(s)
    elif level=='nodewise':
        for k, sub in pd.concat([edges.rename(columns={'i':'node'}), edges.rename(columns={'j':'node'})]).groupby('node'):
            s = _fit(sub, with_interactions); s['key'] = k; rows.append(s)
    else:
        raise ValueError("level must be 'global','rsn_pairs','nodewise'")
    return pd.DataFrame(rows)
