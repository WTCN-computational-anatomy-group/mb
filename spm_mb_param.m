function varargout = spm_mb_param(varargin)
%__________________________________________________________________________
%
% Functions for settings and parameters related.
%
% FORMAT nw                                 = spm_mb_param('GetNumWork',sett)
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
    case 'GetNumWork'
        [varargout{1:nargout}] = GetNumWork(varargin{:});    
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
% GetNumWork()
function nw = GetNumWork(sett)
% Estimate number of parfor workers from available system RAM
% (if sett.gen.num_workers = -1)

% Parse function settings
MemMax       = sett.gen.memmx;   % max memory usage (in MB)
NumWork      = sett.gen.num_workers;
dm           = sett.var.d;       % current template dimensions
K            = sett.model.K;     % template classes
print2screen = sett.show.print2screen;

if NumWork >= 0
    % Only estimates number of workers if sett.gen.num_workers = -1
    nw = NumWork;
    return
end

if MemMax == 0 % default
    try
        % Get memory info automatically (in MB)
        if ispc
            % Windows
            [~,meminfo] = memory;
            MemMax      = meminfo.PhysicalMemory.Available;
        else         
            % UNIX
            [~,meminfo] = system('free --mega');
            meminfo     = string(meminfo);
            meminfo     = strsplit(meminfo,' ');
            MemMax      = str2double(meminfo{14}); % field that holds available RAM (MATLAB 2018a)
        end
    catch
        MemMax = 8096;
    end
end

% Get memory requirement (with current template size)
K1             = K + 1;
Multiplier     = 8;                              % A multiplier that has been estimated empirically based on algorithm memory requirement
NumFloats      = Multiplier*prod(dm(1:3))*K1;    % times two..we also keep images, etc in memory (rough)
FloatSizeBytes = 4;                              % One float is four bytes (single precision)
MemReq         = (NumFloats*FloatSizeBytes)/1e6; % to MB

nw = floor(MemMax/MemReq) - 1; % Number of parfor workers to use (minus one..for main thread)

if print2screen, printf('Current amount of available RAM allows parfor to use %i workers\n',nw); end
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
    sett.write.model      = false;
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
              'do',      struct('gmm', true, 'infer',1, 'mu_bg', 0, 'updt_aff', true,...
                                'updt_bf', true, 'updt_prop', true, 'updt_template', true,...
                                'updt_vel', true, 'zoom', true, 'updt_int', true),...
              'gen',     struct('accel', 0.8, 'run2d', 0, 'samp', 3, 'samp_min', 1, 'num_workers', 0, 's_settings', [3 2],'args', 8, 'memmx', 0),...
              'labels',  struct('use', false, 'w', 0.99, 'use_initgmm', true),...
              'model',   struct('crop_mu', false, 'groupwise', false, 'init_mu_dm', 16, 'ix_init_pop', 1,...
                                'K', 6, 'mg_ix', 1, 'mg_ix_intro', true, 'mu_bg', [], 'vx', 1, 'tol_aff', 1e-5,...
                                'tol_diffeo', 0, 'appear_ix', 1, 'appear_chn', []),...
              'nit',     struct('appear', 2, 'bf', 1, 'gmm', 32, 'init', 32, 'init_mu', 3, 'gmm_miss', 32, 'zm', 3),...
              'optim',   struct('nls_bf', 1, 'scal_bf', 1, 'scal_q', 1, 'scal_v', 1),...              
              'show',    struct('axis_3d', 3, 'channel', 1, 'figs', {{}}, 'figname_bf', '(spm_mb) Bias fields',...
                                'figname_int', '(spm_mb) Intensity model', 'figname_init', '(spm_mb) Init',...
                                'figname_model', '(spm_mb) Template model', 'figname_parameters', '(spm_mb) Parameters',...
                                'figname_subjects', '(spm_mb) Segmentations', 'print2screen', true, 'mx_subjects', 8, 'dir_vis', 'tmp-vis'),...
              'var',     struct('mu_settings', [1.0000e-02 0.5000 0], 'v_settings',  [4e-4 0 0.8000 0.2000 0.8000]*0.5),...
              'write',   struct('bf', false, 'clean_def', true, 'clean_vel', true, 'df', true, 'model', true,...
                                'mu', [true false], 'dir_res', '', 'im', [false false false false], 'vel', false,...
                                'affine', false, 'workspace', false, 'tc', [false false false]));
 
if nargin==0, sett = struct; end
sett = copy_fields(sdef,sett);

% Depending on what is disabled, disable other stuff too
if ~sett.do.gmm
    sett.do.updt_int = false;
    sett.do.updt_bf  = false;
end

% For including multiple Gaussians per tissue
sett.model.mg_ix_intro = sett.model.mg_ix;
sett.model.mg_ix       = 1;

% If no affine registration, include a penalty on absolute deformations
if numel(sett.var.v_settings)>1
    if sett.do.updt_aff
        sett.var.v_settings(1) = 0;
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

