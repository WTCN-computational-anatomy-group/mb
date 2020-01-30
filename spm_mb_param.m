function varargout = spm_mb_param(varargin)
%__________________________________________________________________________
%
% Functions for settings and parameters related.
%
% FORMAT [template_given,appear_given,sett] = spm_mb_param('SetFit',model,sett)
% FORMAT sett                               = spm_mb_param('Settings')
% FORMAT sz                                 = spm_mb_param('ZoomSettings',d, Mmu, v_settings, mu_settings, n)
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
    case 'SetFit'
        [varargout{1:nargout}] = SetFit(varargin{:});
    case 'Settings'
        [varargout{1:nargout}] = Settings(varargin{:});
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
function [template_given,appear_given,sett] = SetFit(model,sett)

template_given = (isfield(model,'shape') && isfield(model.shape,'template'));
appear_given   = isfield(model,'appear');

if template_given && appear_given
    % Fit learned shape and appearance model
    sett.model.groupwise  = false;
    sett.do.updt_int      = false;
    sett.do.updt_template = false;
elseif template_given
    % Fit learned shape model, learn appearance model
    sett.model.groupwise  = false;
    sett.do.updt_int      = true;
    sett.do.updt_template = false;
elseif appear_given
    % Fit learned appearance model, learn shape model
    sett.model.groupwise  = true;
    sett.do.updt_int      = false;
    sett.do.updt_template = true;
end
end
%==========================================================================

%==========================================================================
% Settings()
function sett = Settings(sett)

%------------------
% .appear (appearance model related)
%------------------

if ~isfield(sett,'appear')
    sett.appear = struct;
end
if ~isfield(sett.appear,'tol')
    sett.appear.tol_gmm = 1e-4;
end

%------------------
% .bf (bias field related)
%------------------

if ~isfield(sett,'bf')
    sett.bf = struct;
end
if ~isfield(sett.bf,'fwhm')
    sett.bf.fwhm = 60;
end
if ~isfield(sett.bf,'reg')
    sett.bf.reg = 1e5;
end

%------------------
% .clean_z (resp clean-up related, spm_mb_output)
%------------------

if ~isfield(sett,'clean_z')
    sett.clean_z = struct;
end
if ~isfield(sett.clean_z,'mrf')
    sett.clean_z.mrf = 0; % 1
end
if ~isfield(sett.clean_z,'nit_mrf')
    sett.clean_z.nit_mrf = 10;
end
if ~isfield(sett.clean_z,'gwc_tix')
    sett.clean_z.gwc_tix = [];
end
if ~isfield(sett.clean_z,'gwc_level')
    sett.clean_z.gwc_level = 1; % 1
end

%------------------
% .do (enable/disable functionality)
%------------------

if ~isfield(sett,'do')
    sett.do = struct;
end
if ~isfield(sett.do,'gmm')
    sett.do.gmm = true;
end
if ~isfield(sett.do,'infer')
    sett.do.infer = 1; % 0, 1, 2
end
if ~isfield(sett.do,'mu_bg')
    % Should voxels outside the FOV of the template, when warped to subject
    % space, be replaced? (only when using an already learned template)
    sett.do.mu_bg = false;
end
if ~isfield(sett.do,'updt_aff')
    sett.do.updt_aff = true;
end
if ~isfield(sett.do,'updt_bf')
    sett.do.updt_bf = true;
end
if ~isfield(sett.do,'updt_int')
    sett.do.updt_int = true;
end
if ~sett.do.gmm
    sett.do.updt_int = false;
end
if ~isfield(sett.do,'updt_prop')
    sett.do.updt_prop = true;
end
if ~isfield(sett.do,'updt_template')
    sett.do.updt_template = true;
end
if ~isfield(sett.do,'updt_vel')
    sett.do.updt_vel = true;
end
if ~isfield(sett.do,'zoom')
    sett.do.zoom = true;
end

% Depending on what is disabled, disable other stuff too
if ~sett.do.gmm
    sett.do.updt_bf = false;
end

%------------------
% .gen (general)
%------------------

if ~isfield(sett,'gen')
    sett.gen = struct;
end
if ~isfield(sett.gen,'accel')
    % 0 <= accel <= 1: 0 -> slow & stable; 1 -> fast & unstable
    sett.gen.accel = 0.8;
end
if ~isfield(sett.gen,'run2d')
    sett.gen.run2d = 0; % 0, 1, 2, 3
end
if ~isfield(sett.gen,'samp')
    sett.gen.samp = 3;
end
if ~isfield(sett.gen,'samp_min')
    sett.gen.samp_min = 1;
end

%------------------
% .labels (label related)
%------------------

if ~isfield(sett,'labels')
    sett.labels = struct;
end
if ~isfield(sett.labels,'use')
    sett.labels.use = false;
end
if ~isfield(sett.labels,'w')
    sett.labels.w  = 0.99;
end
if ~isfield(sett.labels,'use_initgmm')
    sett.labels.use_initgmm  = true;
end

%------------------
% .model (model related)
%------------------

if ~isfield(sett,'model')
    sett.model = struct;
end
if ~isfield(sett.model,'crop_mu')
    sett.model.crop_mu = false;
end
if ~isfield(sett.model,'groupwise')
    sett.model.groupwise = false;
end
if ~isfield(sett.model,'init_mu_dm')
    % Minimum dimensions of template
    sett.model.init_mu_dm = 16;
end
if ~isfield(sett.model,'ix_init_pop')
    % Index of population to use for initialising, then more than one
    % population
    sett.model.ix_init_pop = 1;
end
if ~isfield(sett.model,'K')
    sett.model.K = 6;
end
if ~isfield(sett.model,'mg_ix')
    % For using multiple Gaussians per tissue (as in spm_preproc8)
    sett.model.mg_ix = 1;
end
sett.model.mg_ix_intro = sett.model.mg_ix;
sett.model.mg_ix       = 1;
if ~isfield(sett.model,'mu_bg')
    % For dealing with template FOV being smaller than subject image's. If
    % empty pullc/pushc are used when warping template, which is the
    % scenario for learning a template. If a template is given, then this
    % parameter is set using spm_mb_shape('MuValOutsideFOV',mu,sett). Its
    % values are then used to fill in the missing FOV.
    sett.model.mu_bg = [];
end
if ~isfield(sett.model,'vx')
    sett.model.vx = 1;
end
if ~isfield(sett.model,'tol')
    sett.model.tol = 1e-4;
end
if ~isfield(sett.model,'appear_ix')
    sett.model.appear_ix = 1;
end

%------------------
% .nit (iteration related)
%------------------

if ~isfield(sett,'nit')
    sett.nit = struct;
end
if ~isfield(sett.nit,'appear')
    sett.nit.appear = 2;
end
if ~isfield(sett.nit,'bf')
    sett.nit.bf = 1;
end
if ~isfield(sett.nit,'gmm')
    sett.nit.gmm = 32;
end
if ~isfield(sett.nit,'init')
    % The number of iterations, for init rigid alignment.
    sett.nit.init = 16;
end
if ~isfield(sett.nit,'init_mu')
    % The number of template update iterations
    sett.nit.init_mu = 3;
end
if ~isfield(sett.nit,'miss')
    sett.nit.gmm_miss = 32;
end
if ~isfield(sett.nit,'zm')
    % The number of iterations, for updating all model parameters, at each zoom
    % level. The final zoom level uses sett.nit.zm iterations, while
    % earlier zoom levels use sett.nit.zm + zoom_level.
    sett.nit.zm = 3;
end

%------------------
% .optim (optimisation related)
%------------------

if ~isfield(sett,'optim')
    sett.optim = struct;
end
if ~isfield(sett.optim,'nls_bf')
    sett.optim.nls_bf = 1;
end
if ~isfield(sett.optim,'scal_bf')
     % Scaling of bf GN updates
    sett.optim.scal_bf = 1.0;
end
if ~isfield(sett.optim,'scal_q')
     % Scaling of q GN updates
    sett.optim.scal_q = 1.0;
end
if ~isfield(sett.optim,'scal_v')
     % Scaling of v GN updates
    sett.optim.scal_v = 1.0;
end

%------------------
% .registr (registration related)
%------------------

if ~isfield(sett,'registr')
    sett.registr = struct;
end

%------------------
% .show (visualisation related)
%------------------

if ~isfield(sett,'show')
    sett.show = struct;
end
if ~isfield(sett.show,'axis_3d')
    sett.show.axis_3d = 3; % 1, 2, 3
end
if ~isfield(sett.show,'channel')
    sett.show.channel = 1; % 1, ..., C
end
if ~isfield(sett.show,'figs')
    sett.show.figs = {}; % {'model','normalised','segmentations','intensity','parameters','InitGMM'}
end
if ~isfield(sett.show,'figname_bf')
    sett.show.figname_bf = '(spm_mb) Bias fields';
end
if ~isfield(sett.show,'figname_int')
    sett.show.figname_int = '(spm_mb) Intensity model';
end
if ~isfield(sett.show,'figname_gmm')
    sett.show.figname_gmm = '(spm_mb) InitGMM';
end
if ~isfield(sett.show,'figname_model')
    sett.show.figname_model = '(spm_mb) Template model';
end
if ~isfield(sett.show,'figname_parameters')
    sett.show.figname_parameters = '(spm_mb) Parameters';
end
if ~isfield(sett.show,'figname_subjects')
    sett.show.figname_subjects = '(spm_mb) Segmentations';
end
if ~isfield(sett.show,'figname_imtemplatepace')
    sett.show.figname_imtemplatepace = '(spm_mb) Template space data';
end
if ~isfield(sett.show,'print2screen')
    sett.show.print2screen = true;
end
if ~isfield(sett.show,'mx_subjects')
    sett.show.mx_subjects = 2;
end

%------------------
% .shoot (shooting related)
%------------------

if ~isfield(sett,'shoot')
    sett.shoot = struct;
end
if ~isfield(sett.shoot,'args')
    % Timesteps for shooting integration
    sett.shoot.args = 8;
end
if ~isfield(sett.shoot,'s_settings')
    % Full multigrid settings
    sett.shoot.s_settings = [3 2];
end

%------------------
% .var (variable - parameters that vary during algorithm progress)
%------------------

if ~isfield(sett,'var')
    sett.var = struct;
end
if ~isfield(sett.var,'mu_settings')
    sett.var.mu_settings = [1e-3 0.2 0];
end
if ~isfield(sett.var,'v_settings')
    if sett.do.updt_aff
        sett.var.v_settings = [0 0 0.2 0.05 0.2]*2^2;
    else
        sett.var.v_settings = [1e-4 0 0.2 0.05 0.2]*2^2;
    end
end
%------------------
% .write (writing to disk)
%------------------

if ~isfield(sett,'write')
    sett.write = struct;
end
if ~isfield(sett.write,'bf')
    sett.write.bf = false; % field
end
if ~isfield(sett.write,'clean_def')
    % Remove nifti file containing deformation after algorithm finishes
    sett.write.clean_def = true;
end
if ~isfield(sett.write,'clean_vel')
    % Remove nifti file containing velocities after algorithm finishes
    sett.write.clean_vel = true;
end
if ~isfield(sett.write,'df')
    sett.write.df = false(1,2); % forward, inverse
end
if ~isfield(sett.write,'labels')
    % Write labels
    sett.write.labels = [false false]; % native, template
end
if ~isfield(sett.write,'model')
    sett.write.model = true;
end
if ~isfield(sett.write,'mu')
    sett.write.mu = [true false]; % log, softmax
end
if ~isfield(sett.write,'dir_res')
    sett.write.dir_res = '';
end
if ~isfield(sett.write,'im')
    sett.write.im = false(1,4); % image, corrected, warped, warped corrected
end
if ~isfield(sett.write,'vel')
    sett.write.vel = false;
end
if ~isfield(sett.write,'affine')
    sett.write.affine = false;
end
if ~isfield(sett.write,'workspace')
    sett.write.workspace = false;
end
if ~isfield(sett.write,'tc')
    sett.write.tc = false(1,3); % native, warped, warped-mod
end

% Make directories (if does not exist)
if ~isempty(sett.write.dir_res) && ~(exist(sett.write.dir_res,'dir') == 7)
    mkdir(sett.write.dir_res);
end
s                  = what(sett.write.dir_res); % Get absolute path
sett.write.dir_res = s.path;

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
