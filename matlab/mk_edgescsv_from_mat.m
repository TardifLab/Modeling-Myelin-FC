function mk_edgescsv_from_mat()
% MK_EDGESCSV_FROM_MAT 
% Converts NxN connectivity matrices into required .csv file by FC modality
%
% addpath(genpath('matlab'))
% mk_edgescsv_from_mat()


% Setup
  ROOTDIR       = fileparts(fileparts(mfilename('fullpath')));
  DATADIR       = [ROOTDIR '/data'];                                        
  PARCDIR       = [ROOTDIR '/parcellations'];
  addpath(genpath(ROOTDIR))

% ----- USER EDIT SECTION ----------------------------------------------
  fc_fn         = 'fc_BOLDin';                                              % *** TOGGLE FC NETWORK HERE ***

  myelin_fn     = 'mtsat';                                                  % edge myelin
  caliber_fn    = 'commit';                                                 % edge caliber;
  length_fn     = 'edgelength';                                             % edge length
  parc          = 'schaefer-400';                                           % {'schaefer-200' 'schaefer-400'}
% ----------------------------------------------------------------------

% Load data
  loaddir           = [DATADIR '/' parc '_group/'];
  load([loaddir caliber_fn '.mat'],'Dtg');                  caliber_mat=Dtg;
  load([loaddir myelin_fn '.mat'],'Dtg');                   myelin_mat=Dtg;
  load([loaddir length_fn '.mat'],'Dtg');                   length_mat=Dtg;         % Edge Length
  load([loaddir fc_fn '.mat'],'Dtg');                       FC_mat=Dtg;                 % FC
  clear Dtg

% Load SAaxis (not actually used here)
  % SA_map  = readNPY([loaddir 'sa_axis.npy']);

% Set missing values to nan for removal (do for all in case they are different)
  caliber_mat(caliber_mat==0)=nan;
  myelin_mat(myelin_mat==0)=nan;
  length_mat(length_mat==0)=nan;

% load parc info
  pinfo     = conn_getParcInfo({parc},PARCDIR,'cor');                       % structure with parcellation info
  Nstx      = pinfo.nsub;
  Nctx      = pinfo.ncor;
  Nnode     = Nctx + Nstx;
  Nntwk     = length(pinfo.clabels);

% Yeo RSNs
  rsn_idx    = nan(Nnode,1); 
  rsn_labels = cell(Nnode,1);                                               % Nx1 string/cellstr of RSN per node (optional)
  for ii=1:Nntwk 
      rsn_idx(pinfo.cis{ii},1)  = ii; 
      rsn_labels(pinfo.cis{ii}) = pinfo.clabels(ii);
  end 

% --- Build CSV ---
% lower-triangle edge list (exclude self-edges)
  % N             = size(FC_mat,1);
  N             = Nnode;
  [i,j]         = find(tril(true(N), -1));

  edges         = table();
  edges.i       = i;
  edges.j       = j;
  edges.FC      = FC_mat(sub2ind([N N], i, j));                             % NxN functional connectivity
  edges.caliber = caliber_mat(sub2ind([N N], i, j));                        % NxN edge caliber
  edges.myelin  = myelin_mat(sub2ind([N N], i, j));                         % NxN edge myelin
  edges.length  = length_mat(sub2ind([N N], i, j));                         % NxN edge length

% Optional RSN columns (recommended if you’ll run RSN-pairs)
  if exist('rsn_labels','var') && ~isempty(rsn_labels)
        edges.rsn_i = string(rsn_labels(i));
        edges.rsn_j = string(rsn_labels(j));
  end

% remove rows with any missing values
  mask = all(~isnan(edges{:, {'FC','caliber','myelin','length'}}),2);
  edges = edges(mask,:);

% Save
  writetable(edges, fullfile(DATADIR,['edges_' fc_fn '.csv']));
end
