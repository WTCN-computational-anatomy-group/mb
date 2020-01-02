function varargout = spm_mb_param(varargin)
%__________________________________________________________________________
%
% Functions for settings and parameters related.
%
% FORMAT [sett,given] = spm_mb_param('ConditionalSettings',model,sett,ndat)
% FORMAT sett         = spm_mb_param('DefaultSettings')
% FORMAT sz           = spm_mb_param('ZoomSettings',d, Mmu, v_settings, mu_settings, n)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_param
    error('Not enough argument. Type ''help spm_mb_param'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id                   
    case 'ConditionalSettings'
        [varargout{1:nargout}] = ConditionalSettings(varargin{:});        
    case 'DefaultSettings'
        [varargout{1:nargout}] = DefaultSettings(varargin{:});   
    case 'ZoomSettings'
        [varargout{1:nargout}] = ZoomSettings(varargin{:});           
    otherwise
        help spm_mb_param
        error('Unknown function %s. Type ''help spm_mb_param'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% SetFit()
function [sett,given] = ConditionalSettings(model,sett,ndat)
% Change settings based on input model and number of subjects.

% TODO: What if we want to use the provided model to initialise, but still 
%       want to optimise it (like some sort of pretraining)?
%       Instead, we could use SetDefault here, so that if the user 
%       explicitely asked to update something, it is indeed updated.

given          = struct;
given.template = (isfield(model,'shape') && isfield(model.shape,'template'));
given.subspace = (isfield(model,'shape') && isfield(model.shape,'subspace'));
given.appear   = isfield(model,'appear');

if ndat == 1 || given.template && given.appear
    % Fit learned shape and appearance model
    sett.model.groupwise  = false;
    sett.do.updt_int      = false;
    sett.do.updt_template = false;
elseif given.template
    % Fit learned shape model, learn appearance model
    sett.model.groupwise  = false;
    sett.do.updt_int      = true;
    sett.do.updt_template = false;
elseif given.appear
    % Fit learned appearance model, learn shape model
    sett.model.groupwise  = true;
    sett.do.updt_int      = false;
    sett.do.updt_template = true;
else
    % Learn both shape and appearance model
    sett.model.groupwise  = true;
    sett.do.updt_int      = true;
    sett.do.updt_template = true;
end
end
%==========================================================================

%==========================================================================
% Settings()
function s = DefaultSettings(s)

if nargin < 1, s = struct; end

%------------------
% .appear (appearance model)
%------------------

s = SetDefault(s, 'appear.tol_gmm', 1e-4); % Tolerance when fitting GMM

%------------------
% .bf (bias field)
%------------------

s = SetDefault(s, 'bf.fwhm', 60);   % Cutoff on bias field frequency
s = SetDefault(s, 'bf.reg',  1e5);  % Bending energy regularisation

%------------------
% .pca (shape modelling)
%------------------

s = SetDefault(s, 'pca.npc',          10);    % Number of principal components
s = SetDefault(s, 'pca.res_prior',    14);    % Gamma prior on residuals precision - expected value
s = SetDefault(s, 'pca.res_df',       1e-3);  % Gamma prior on residuals precision - degrees of freedom
s = SetDefault(s, 'pca.latent_prior', 1);     % Wishart prior on latent precision  - expected value
s = SetDefault(s, 'pca.latent_df',    1e-3);  % Wishart rior on latent precision   - degrees of freedom
% For df: % 0=ML, (>0)=prior, Inf=fixed/dirac

%------------------
% .clean_z (resp clean-up related, spm_mb_output)
%------------------

s = SetDefault(s, 'clean_z.mrf',       false);  % Markov Random Field cleanup?
s = SetDefault(s, 'clean_z.nit_mrf',   10);     % Number of MRF iterations
s = SetDefault(s, 'clean_z.gwc_tix',   []);     % Grey-White cleanup: tissue id
s = SetDefault(s, 'clean_z.gwc_level', 1);      % ?

%------------------
% .do (enable/disable functionality)
%------------------

s = SetDefault(s, 'do.gmm',                 true);  % Gaussian mixture
s = SetDefault(s, 'do.pca',                 false); % Shape prior (PCA)
s = SetDefault(s, 'do.infer',               1);     % 0,1,2 ??
s = SetDefault(s, 'do.updt_aff',            true);  % Affine registration
s = SetDefault(s, 'do.updt_bf',             true);  % Bias field correction
s = SetDefault(s, 'do.updt_int',            true);  % ?
s = SetDefault(s, 'do.updt_template',       true);  % Template
s = SetDefault(s, 'do.updt_subspace',       true);  % Principal subspace (PCA)
s = SetDefault(s, 'do.updt_latent_prior',   true);  % Latent precision matrix (PCA)
s = SetDefault(s, 'do.updt_res_prior',      true);  % Residual precision (PCA)
s = SetDefault(s, 'do.updt_vel',            true);  % Velocities (diffeomorphic registration)
s = SetDefault(s, 'do.zoom',                true);  % Multiscale fit
% TODO: . can't we just set nz=1 instead of do.zoom = false?
%       . maybe 'updt' should be 'fit'
%       . Do we prefer single level (do.updt_bf) or multilevel (do.updt.bf)?

% Depending on what is disabled, disable other stuff too
if ~s.do.gmm, s.do.updt_bf           = false; end
if ~s.do.pca, s.do.updt_subspace     = false;
              s.do.updt_latent_prior = false;
              s.do.updt_res_prior    = false; end

%------------------
% .gen (general)
%------------------

s = SetDefault(s, 'gen.accel',    0.8); % 0 <= accel <= 1: 0 -> slow & stable; 1 -> fast & unstable
s = SetDefault(s, 'gen.run2d',    0);   % 0, 1, 2, 3 ?
s = SetDefault(s, 'gen.samp',     3);   % Sampling distance (mm)
s = SetDefault(s, 'gen.max_mem',  1);   % Size (GB) above which the template/subspace is stored on disk
% TODO: Maybe accel should move to .optim?

%------------------
% .labels (manual labels)
%------------------

s = SetDefault(s, 'labels.use',    false);  % Use manual labels
s = SetDefault(s, 'labels.w',      0.99);   % Label confidence (0<>1)
% TODO: . Something more telling than 'w'
%       . labels.w = 0 <=> labels.use = false ?

%------------------
% .model
%------------------

s = SetDefault(s, 'model.groupwise',   false); % Groupwise vs Independent
s = SetDefault(s, 'model.init_mu_dm',  8);     % Initial template dimension
s = SetDefault(s, 'model.K',           5);     % Number of template classes
s = SetDefault(s, 'model.vx',          1);     % Final template voxel size
s = SetDefault(s, 'model.ix_init_pop', 1);     % Index of population to use for initialising
s = SetDefault(s, 'model.mg_ix',       1);     % Number of Gaussians per tissue
% TODO: . Do we really need a 'groupwise' option?
%         Can't it be inferred from 'do'?
%         Or should groupwise == false, induce do.<lots> = false?
%       . I'd prefer 'nclass' or 'ntissue' to 'K'

%------------------
% .nit (iterations)
%------------------

s = SetDefault(s, 'nit.appear',   2);   % ?
s = SetDefault(s, 'nit.bf',       1);   % Bias field: gauss-newton iterations
s = SetDefault(s, 'nit.gmm',      32);  % GMM: EM (gauss+resp) iterations
s = SetDefault(s, 'nit.init',     6);   % Initial: EM (template+rigid) iterations
s = SetDefault(s, 'nit.init_mu',  2);   % Template: gauss-newton iterations
s = SetDefault(s, 'nit.gmm_miss', 32);  % GMM: sub EM (missing data) iterations
s = SetDefault(s, 'nit.zm',       3);   % Number of iterations at each zoom level
                                        % The final zoom level uses sett.nit.zm iterations, 
                                        % while, earlier zoom levels use sett.nit.zm + zoom_level.
% TODO: . I like 'iter' better than 'nit'
                                       
%------------------
% .optim (optimisation)
%------------------

s = SetDefault(s, 'optim.nls_bf',  1); % Bias field: number of line searches
s = SetDefault(s, 'optim.scal_q',  1); % Affine:     scaling of GN updates
s = SetDefault(s, 'optim.scal_v',  1); % Velocities: scaling of GN updates
% TODO: . Should nls_bf move to nit?

%------------------
% .show (visualisation)
%------------------

s = SetDefault(s, 'show.axis_3d',                3);   % Axis to select 2d slice from
s = SetDefault(s, 'show.channel',                1);   % Channel for multicontrast data
s = SetDefault(s, 'show.figs',                   {});  % {'model','normalised','segmentations','intensity','parameters'}
s = SetDefault(s, 'show.figname_bf',             '(spm_mb) Bias fields');
s = SetDefault(s, 'show.figname_imtemplatepace', '(spm_mb) Template space data');
s = SetDefault(s, 'show.figname_int',            '(spm_mb) Intensity model');
s = SetDefault(s, 'show.figname_model',          '(spm_mb) Template model');
s = SetDefault(s, 'show.figname_parameters',     '(spm_mb) Parameters');
s = SetDefault(s, 'show.figname_subjects',       '(spm_mb) Segmentations');
s = SetDefault(s, 'show.mx_subjects',            2);
s = SetDefault(s, 'show.print2screen',           true);
% TODO: . print2screen -> verbose?

%------------------
% .shoot (diffeomorphic shooting)
%------------------

s = SetDefault(s, 'shoot.args',       3);      % Timesteps for shooting integration
s = SetDefault(s, 'shoot.s_settings', [3 2]);  % Full multigrid settings
% TODO: . timestep or nb of Euler steps?
%       . why are FMG parameters stored here?

%------------------
% .var (variable - parameters that vary during algorithm progress)
%------------------

s = SetDefault(s, 'var.mu_settings',  [1e-3 0.2 0]); % Template: regularisation
if s.do.updt_aff, s = SetDefault(s, 'var.v_settings', [0 0 0.2 0.05 0.2]*2^4);  % Velocities: regularisation
else,             s = SetDefault(s, 'var.v_settings', [1e-4 0 0.2 0.05 0.2]*2^4); end
% TODO: . I think this structure should be empty (and v_, mu_ settings 
%         moved somewhere else) and only filled during initialisation.

%------------------
% .write (writing to disk)
%------------------

s = SetDefault(s, 'write.bf',        false);         % Bias field
s = SetDefault(s, 'write.clean_def', false);         % Clean warp files
s = SetDefault(s, 'write.clean_vel', false);         % Clean velocity files
s = SetDefault(s, 'write.df',        false(1,2));    % Warps [forward inverse]
s = SetDefault(s, 'write.model',     true);          % Model parts (template/intensity/pca)
s = SetDefault(s, 'write.mu',        [true false]);  % Template [log softmax]
s = SetDefault(s, 'write.dir_res',   '.');           % Results directory
s = SetDefault(s, 'write.im',        false(1,4));    % Input: [native corrected warped warped&corrected] 
s = SetDefault(s, 'write.tc',        true(1,3));     % Tissue classes: [native warped warped&modulated]
s = SetDefault(s, 'write.workspace', false);         % Save workspace variables (.mat)

% Make directories (if does not exist)
if ~(exist(s.write.dir_res,'dir') == 7) 
    ok = mkdir(s.write.dir_res);
    if ~ok, error('Could not create output directory. Check permissions.'); end
end
s.write.dir_res = AbsoluteDir(s.write.dir_res);

end
%==========================================================================

%==========================================================================
% ZoomSettings()
function sz = ZoomSettings(d, Mmu, v_settings, mu_settings, n)
[dz{1:n}] = deal(d);
sz        = struct('Mmu',Mmu,'d',dz,...
                   'v_settings', v_settings,...
                   'mu_settings',mu_settings);

% I'm still not entirely sure how best to deal with regularisation
% when dealing with different voxel sizes.
scale = 1/abs(det(Mmu(1:3,1:3)));
for i=1:n
    sz(i).d           = ceil(d/(2^(i-1)));
    z                 = d./sz(i).d;
    sz(i).Mmu         = Mmu*[diag(z), (1-z(:))*0.5; 0 0 0 1];
    vx                = sqrt(sum(sz(i).Mmu(1:3,1:3).^2));
    sz(i).v_settings  = [vx v_settings *(scale*abs(det(sz(i).Mmu(1:3,1:3))))];
    sz(i).mu_settings = [vx mu_settings*(scale*abs(det(sz(i).Mmu(1:3,1:3))))];
end
end
%==========================================================================


%==========================================================================
function o = SetDefault(o, field, value, exists)
% Set default value in option structure. The value is only set if none 
% existed previously.
%
% FORMAT opt = SetDefault(opt, field, value)
% opt   - Structure
% field - Hierarchy of fields (cell, or '.' separated)
% value - Default value
%
% EXAMPLE
% >> opt.my.field1 = 1;
% >> opt = SetDefault(opt, 'my.field1', 2);
% >> opt.my.field1
% ans = 1
% >> opt = SetDefault(opt, 'my.field2', 3);
% >> opt.my.field2
% ans = 3
if nargin < 4, exists = false; end
if ~iscell(field), field = strsplit(field, '.'); end

if isempty(field)
    if ~exists, o = value; end
else
    exists = isfield(o, field{1});
    if ~exists, o.(field{1}) = []; end
    o.(field{1}) = SetDefault(o.(field{1}), field(2:end), value, exists);
end
end
%==========================================================================

% === AbsoluteDir =========================================================
function path = AbsoluteDir(path)
% Return the absolute path to a directory
path = what(path);
path = path.path;
end