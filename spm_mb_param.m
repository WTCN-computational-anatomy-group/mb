function varargout = spm_mb_param(varargin)
%__________________________________________________________________________
%
% Functions for settings and parameters related.
%
% FORMAT [template_given,appear_given,sett] = spm_mb_param('SetFit',model,sett)
% FORMAT sett                               = spm_mb_param('Settings')
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
sdef = struct('appear',  struct('tol_gmm',1e-4),...
              'bf',      struct('fwhm', 60, 'reg', 1e5),...
              'clean_z', struct('mrf', 0, 'nit_mrf', 10, 'gwc_tix', [], 'gwc_level', 1),...
              'do',      struct('gmm', true, 'infer',true, 'mu_bg', 0, 'updt_aff', true,...
                                'updt_bf', true, 'updt_prop', true, 'updt_template', true,...
                                'updt_vel', true, 'zoom', true),...
              'gen',     struct('accel', 0.8, 'run2d', false, 'samp', 3, 'samp_min', 1, 'num_workers', 0),...
              'labels',  struct('use', false, 'w', 0.99, 'use_initgmm', true),...
              'model',   struct('crop_mu', false, 'groupwise', false, 'init_mu_dm', 16, 'ix_init_pop', true,...
                                'K', 6, 'mg_ix', true, 'mg_ix_intro', true, 'mu_bg', [], 'vx', 1, 'tol_aff', 1e-5,...
                                'tol_diffeo', 0, 'appear_ix', true, 'appear_chn', []),...
              'nit',     struct('appear', 2, 'bf', 1, 'gmm', 32, 'init', 32, 'init_mu', 3, 'gmm_miss', 32, 'zm', 3),...
              'optim',   struct('nls_bf', 1, 'scal_bf', 1, 'scal_q', 1, 'scal_v', 1),...
              'registr', struct(),...
              'show',    struct('axis_3d', 3, 'channel', 1, 'figs', {{}}, 'figname_bf', '(spm_mb) Bias fields',...
                                'figname_int', '(spm_mb) Intensity model', 'figname_gmm', '(spm_mb) InitGMM',...
                                'figname_model', '(spm_mb) Template model', 'figname_parameters', '(spm_mb) Parameters',...
                                'figname_subjects', '(spm_mb) Segmentations', 'print2screen', true, 'mx_subjects', 8, 'dir_vis', 'tmp-vis'),...
              'var',     struct('mu_settings', [1.0000e-02 0.5000 0], 'v_settings',  [4e-4 0 0.8000 0.2000 0.8000]),...
              'write',   struct('bf', false, 'clean_def', true, 'clean_vel', true, 'df', true, 'model', true,...
                                'mu', [true false], 'dir_res', '', 'im', [false false false false], 'vel', false,...
                                'affine', false, 'workspace', false, 'tc', [false false false]));
 
if nargin==0, sett = struct; end
sett = copy_fields(sdef,sett);

% If no affine registration, include a penalty on absolute deformations
if numel(sett.var.v_settings)>1
    if sett.do.updt_aff
        sett.var.v_settings(1) = 0;
    else
        sett.var.v_settings(1) = 4e-4;
    end
end


% Make directories (if does not exist)
if ~isempty(sett.write.dir_res) && ~(exist(sett.write.dir_res,'dir') == 7)
    mkdir(sett.write.dir_res);
end
s                  = what(sett.write.dir_res); % Get absolute path
sett.write.dir_res = s.path;
sett.show.dir_vis  = fullfile(sett.write.dir_res,sett.show.dir_vis);
if any(strcmp(sett.show.figs,'segmentations')) || any(strcmp(sett.show.figs,'parameters'))
    if exist(sett.show.dir_vis,'dir') == 7, rmdir(sett.show.dir_vis,'s'); end
    mkdir(sett.show.dir_vis)
    s                 = what(sett.show.dir_vis); % Get absolute path
    sett.show.dir_vis = s.path;
end
end
%==========================================================================

%==========================================================================
function sett = copy_fields(sdef,sett)
if nargin==1 || isempty(sett)
    sett = sdef;
    return
end
fields = fieldnames(sdef);
for i = 1:numel(fields)
    field = fields{i};
    if ~isfield(sett,field)
        sett.(field) = sdef.(field);
    elseif isa(sdef.(field),'struct') && isa(sett.(field),'struct')
        sett.(field) = copy_fields(sett.(field),sdef.(field));
    else
        sett.(field) = sdef.(field);
    end
end
end
%==========================================================================

