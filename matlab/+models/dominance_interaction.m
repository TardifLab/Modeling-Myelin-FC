function [tab_lo,tab_hi] = dominance_interaction(X,y,mdlspec)
%% function for conducitng domiance analysis based on multiple regression
% --- INPUT:
% X       : matrix with columns holding the explanatory variables
% y       : vector holding the response variable 
% mdlspec : model terms in Wilkinson Notation (should have interactions)
% --- OUTPUT
% tab_lo  : table sumarizing the results with lower order terms only
% tab_hi  : table sumarizing the results incl interactions
%
% --- Extended by Mark C Nelson, 2024, MNI
% + mdlspec input & tab_hi output to handle interaction terms
%
% ======= METHOD
%
% (1) R2_full computed from the full model (incl interactions)
% (2) R2_indi (each predictors individual contribution) is computed:
%   - lower order terms: Y ~ 1 + Xi 
%   - interaction terms: (Y ~ 1 + Xi + Xj + Xi:Xj) - (Y ~ 1 + Xi + Xj)
% (3) R2_rest (the absence of each predictor) is computed:
%   - lower order term Xi by leaving out Xi & all of its interaction terms 
%   - interaction terms by leaving them out
% (4) R2_marg = R2_full - R2_rest
% (5) Dominance_value = (R2_marg + R2_indi) / 2
%
%
%--------------------------------------------------------------------------

% remove missing data
  N = size(X,2);
  data = rmmissing([X y],1);
  X = data(:,1:N);
  y = data(:,N+1);

%% ==============    Lower order terms only
% It is assumed that all columns in X are included in the model as lower
% order terms

% R2
  mdl = fitlm(X,y);
  R2_full_lo = mdl.Rsquared.Ordinary;

% calculte R2 for individual linear models (y = X_i, i = 1:N)
  R2_indi_lo = zeros(N,1);
  for i = 1:N
      mdl = fitlm(X(:,i),y);
      R2_indi_lo(i) = mdl.Rsquared.Ordinary;
  end

% calculate R2 for marginal linear models (y = X_k+...,+X_{N-1}|Xi, i~=k)
  R2_rest_lo = zeros(N,1);
  for i = 1:N
      X_rest = X;
      X_rest(:,i) = [];
      mdl = fitlm(X_rest,y);
      R2_rest_lo(i) = mdl.Rsquared.Ordinary;
  end
  R2_marg_lo = R2_full_lo - R2_rest_lo;

% calculate the general dominance (Shapley) values
  Dominance_value = (R2_marg_lo + R2_indi_lo)/2;

% calculate the standardized domiance (Shapley) values in percentage
  Dominance_percentage = Dominance_value./sum(Dominance_value)*100;

% rank the numbers
  [~,p] = sort(Dominance_percentage,'descend');
  Rank = 1:length(Dominance_percentage);
  Rank = Rank';
  Rank(p) = Rank;

% Wrap up results in a table
  Variable_name = cell(N,1);
  for i = 1:N
      Variable_name{i} = strcat('X',num2str(i));
  end
  tab_lo = table(Variable_name,Dominance_value, Dominance_percentage,Rank);


%% ==============    Including interactions

% R2
  mdl = fitlm(X,y,mdlspec);
  R2_full_hi = mdl.Rsquared.Ordinary;

% --- Get info about interactions
  mdlterms=mdl.Coefficients.Properties.RowNames(2:end);
  N_mt=length(mdlterms);
  i_int=cellstrfind(mdlterms,':');                                          % indices for interaction terms
  str_intterms=mdlterms(i_int);
  mdlspec_intterms=strrep(str_intterms,':','*');                            % used to test model instances with interaction

% --- Get anything that preceeds the model terms
  t_mdlspec=mdlspec;
  t_mdlspec_split=strsplit(t_mdlspec,'+');
  t_mdlspec_start=t_mdlspec_split{1};

% --- Full version of model spec (avoids Wilkinson notation)
  mdlspec_full=[t_mdlspec_start '+' strjoin(mdlterms','+')];

% ---- calculte R2 for individual linear models (y = X_i, i = 1:N)
  R2_indi_hi=zeros(N_mt,1);
  R2_indi_hi(1:N)=R2_indi_lo(1:N);                                          % Lower order terms copied over

% Individual models for interaction terms
% Method 1: Difference of R2; lower order terms subtracted from model with
% both lower order & interactions (R2_indi_i = R2_hi_i - R2_lo_i)
% R2_hi from model: y~1+x1*...,*xM; M is # terms for given interaction
% R2_lo from model: y~1+x1+...,+xM; i.e., lower order terms only
for i = i_int'
  % Get indices for the x terms in this interaction
    t_modelspec=mdlspec_intterms{i-N};
    t_terms=strsplit(t_modelspec,'*');                                      % all terms for this interaction
    t_inds_terms=str2num(cell2mat(strrep(t_terms,'x',' ')));                % extract indices of x terms
    t_N_terms=length(t_inds_terms);
  % Create new model spec for these terms with interaction
    use_intterms_hi='x1'; 
    for xx=2:t_N_terms; use_intterms_hi=[use_intterms_hi '*x' num2str(xx)]; end
    use_mdlspec_hi=['y~1+' use_intterms_hi]; 
  % Create new model spec for these terms WITHOUT interaction
    use_mdlspec_lo=strrep(use_mdlspec_hi,'*','+');
  % Run models for these terms with & without interaction
    mdl = fitlm(X(:,t_inds_terms),y,use_mdlspec_lo);
    t_R2_lo=mdl.Rsquared.Ordinary;
    mdl = fitlm(X(:,t_inds_terms),y,use_mdlspec_hi);
    t_R2_hi=mdl.Rsquared.Ordinary;
    R2_indi_hi(i)=t_R2_hi-t_R2_lo;                                         
end


% --- R2 for marginal linear models (y = X_k+...,+X_{N-1}|Xi, i~=k)
% Method 1: Xi & interaction terms removed
R2_rest_hi = zeros(N_mt,1);
for i = 1:N
    t_str_rmX=['x' num2str(i)];
  % Identify all lower order terms & interactions for this predictor
    t_mdlspec=mdlspec_full;
    t_mdlspec_split=strsplit(t_mdlspec,'+');
    t_i_rmintterms=cellstrfind(t_mdlspec_split,t_str_rmX);                  % indices for model terms to be excluded
  % Run model sans term & all its interactions
    t_mdlspec_split(t_i_rmintterms)=[];
    t_mdlspec=strjoin(t_mdlspec_split,'+');
    mdl = fitlm(X,y,t_mdlspec);
    R2_rest_hi(i) = mdl.Rsquared.Ordinary;
end


% % % % --- R2 for marginal linear models (y = X_k+...,+X_{N-1}|Xi, i~=k)
% % % % Method 2: interaction terms added back in
% % % % Excessively penalizes variables with interaction terms
% % % % Can yield dominance < 0 if lower order term has strong interactions
% % % R2_rest_hi = zeros(N_mt,1);
% % % for i = 1:N
% % %     t_str_rmX=['x' num2str(i)];
% % %   % Get info about any interactions involving this model term
% % %     t_mdlterms=mdlterms';
% % %     t_i_rmintterms=N+cellstrfind(str_intterms,t_str_rmX);                   % indices for interaction terms to be excluded
% % %     t_R2_intterms=R2_indi_hi(t_i_rmintterms);                                 % get their R2_indi_hi
% % %   % Remove terms
% % %     t_mdlterms(t_i_rmintterms)=[];                                          % interactions
% % %     t_i_rmloterm=cellstrfind(t_mdlterms,t_str_rmX);
% % %     t_mdlterms(t_i_rmloterm)=[];                                            % lower order term
% % %   % Run model sans term & all of its interactions
% % %     t_mdlspec=strjoin([t_mdlspec_start t_mdlterms],'+');
% % %     mdl=fitlm(X,y,t_mdlspec);
% % %     R2_rest_hi(i)=mdl.Rsquared.Ordinary + sum(t_R2_intterms(:));            % add R2 from interactions 
% % % end


% --- Repeat for interaction terms
% models of form: y~1+x1+...,+xM; M=# of terms in interaction 
% note: interaction is excluded
for i = i_int'
      t_mdlterms=mdlterms';
      t_mdlterms(i)=[];                                                     % Remove this interaction term
      t_mdlspec=strjoin([t_mdlspec_start t_mdlterms],'+');
      mdl=fitlm(X,y,t_mdlspec);
      R2_rest_hi(i)=mdl.Rsquared.Ordinary;
end

% compute R2 marginal
R2_marg_hi = R2_full_hi - R2_rest_hi;
%calculate the general dominance (Shapley) values
Dominance_value = (R2_marg_hi + R2_indi_hi)/2;
%calculate the standardized domiance (Shapley) values in percentage
Dominance_percentage = Dominance_value./sum(Dominance_value)*100;
%rank the numbers
[~,p] = sort(Dominance_percentage,'descend');
Rank = 1:length(Dominance_percentage);
Rank = Rank';
Rank(p) = Rank;
%Wrap up results in a table
Variable_name = cell(N_mt,1);
for i = 1:N_mt
    Variable_name{i} = mdlterms{i};
end
tab_hi = table(Variable_name,Dominance_value, Dominance_percentage,Rank); 
end