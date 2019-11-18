function varargout = spm_multireg_par(varargin)
%__________________________________________________________________________
%
% Parameter/settings functions for spm_multireg.
%
% FORMAT B    = spm_multireg_par('AffineBases',code)
% FORMAT sett = spm_multireg_par('Settings')
% FORMAT sz   = spm_multireg_par('ZoomSettings',d, Mmu, v_settings, mu_settings, n)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_par
    error('Not enough argument. Type ''help spm_multireg_par'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'AffineBases'
        [varargout{1:nargout}] = AffineBases(varargin{:});                      
    case 'Settings'
        [varargout{1:nargout}] = Settings(varargin{:});   
    case 'ZoomSettings'
        [varargout{1:nargout}] = ZoomSettings(varargin{:});           
    otherwise
        help spm_multireg_par
        error('Unknown function %s. Type ''help spm_multireg_par'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% AffineBases()
function B = AffineBases(code)
% This should probably be re-done to come up with a more systematic way of defining the
% groups.
g     = regexpi(code,'(?<code>\w*)\((?<dim>\d*)\)','names');
g.dim = str2num(g.dim);
if numel(g.dim)~=1 || (g.dim ~=0 && g.dim~=2 && g.dim~=3)
    error('Can not use size');
end
if g.dim==0
    B        = zeros(4,4,0);
elseif g.dim==2
    switch g.code
    case 'T' 
        B        = zeros(4,4,2);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
    case 'SO'
        B        = zeros(4,4,1);
        B(1,2,1) =  1;
        B(2,1,1) = -1;
    case 'SE' 
        B        = zeros(4,4,3);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(1,2,3) =  1;
        B(2,1,3) = -1;
    otherwise
        error('Unknown group.');
    end
elseif g.dim==3
    switch g.code
    case 'T' 
        B        = zeros(4,4,3);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(3,4,3) =  1;
    case 'SO' 
        B        = zeros(4,4,3);
        B(1,2,1) =  1;
        B(2,1,1) = -1;
        B(1,3,2) =  1;
        B(3,1,2) = -1;
        B(2,3,3) =  1;
        B(3,2,3) = -1;
    case 'SE' 
        B        = zeros(4,4,6);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(3,4,3) =  1;
        B(1,2,4) =  1;
        B(2,1,4) = -1;
        B(1,3,5) =  1;
        B(3,1,5) = -1;
        B(2,3,6) =  1;
        B(3,2,6) = -1;
    otherwise
        error('Unknown group.');
    end
end
end
%==========================================================================

%==========================================================================
% Settings()
function sett = Settings(sett)

%------------------
% .bf (bias field related)
%------------------

if ~isfield(sett,'bf')
    sett.bf = struct;
end
if ~isfield(sett.bf,'reg')
    sett.bf.reg = 1e4;
end
if ~isfield(sett.bf,'fwhm')
    sett.bf.fwhm = 60;
end

%------------------
% .do (enable/disable functionality)
%------------------

if ~isfield(sett,'do')
    sett.do = struct;
end
if ~isfield(sett.do,'bf_norm')
    sett.do.bf_norm = true;
end
if ~isfield(sett.do,'gmm')
    sett.do.gmm = true;
end
if ~isfield(sett.do,'updt_aff')
    sett.do.updt_aff = true;
end
if ~isfield(sett.do,'updt_bf')
    sett.do.updt_bf = true;
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
if ~sett.do.updt_bf
    sett.do.bf_norm = false;
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
if ~isfield(sett.gen,'samp_gmm')
    sett.gen.samp_gmm = 3;
end
if ~isfield(sett.gen,'threads')
    sett.gen.threads = Inf;
end

%------------------
% .model (model related)
%------------------

if ~isfield(sett,'model')
    sett.model = struct;
end
if ~isfield(sett.model,'groupwise')
    sett.model.groupwise = false;
end
if ~isfield(sett.model,'K')
    sett.model.K = 3;
end
if ~isfield(sett.model,'vx')    
    sett.model.vx = 1;    
end

%------------------
% .nit (iteration related)
%------------------

if ~isfield(sett,'nit')
    sett.nit = struct;
end
if ~isfield(sett.nit,'bf')   
    sett.nit.bf = 1;
end
if ~isfield(sett.nit,'gmm')   
    sett.nit.gmm = 10;
end
if ~isfield(sett.nit,'init')
    % The number of iterations, at largest zoom level.
    sett.nit.init = 3;
end
if ~isfield(sett.nit,'init_mu')
    % The number of template update iterations
    sett.nit.init_mu = 2;
end
if ~isfield(sett.nit,'miss')   
    sett.nit.gmm_miss = 10;
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
     % Scaling of q GN updates
    sett.optim.nls_bf = 6;
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
if ~isfield(sett.registr,'B')
    if sett.gen.run2d
        sett.registr.B = AffineBases('SE(2)');
    else
        sett.registr.B = AffineBases('SE(3)');
    end
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
if ~isfield(sett.show,'figname_bf')
    sett.show.figname_bf = '(spm_multireg) Bias fields';
end
if ~isfield(sett.show,'figname_int')
    sett.show.figname_int = '(spm_multireg) Intensity model';
end
if ~isfield(sett.show,'figname_model')
    sett.show.figname_model = '(spm_multireg) Template model';
end
if ~isfield(sett.show,'figname_parameters')
    sett.show.figname_parameters = '(spm_multireg) Parameters';
end
if ~isfield(sett.show,'figname_subjects')
    sett.show.figname_subjects = '(spm_multireg) Segmentations';
end
if ~isfield(sett.show,'level')
    sett.show.level = 1;
end
if ~isfield(sett.show,'mx_subjects')
    sett.show.mx_subjects = 8;
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
        sett.var.v_settings = [0 0 0.2 0.05 0.2]*4;
    else
        sett.var.v_settings = [1e-4 0 0.2 0.05 0.2]*4;
    end
end
%------------------
% .write (writing to disk)
%------------------

if ~isfield(sett,'write')
    sett.write = struct;
end
if ~isfield(sett.write,'dir_res')
    sett.write.dir_res = '.';
end

% Make directories (if does not exist)
if ~(exist(sett.write.dir_res,'dir') == 7)  
    mkdir(sett.write.dir_res);  
end

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
scale     = 1/abs(det(Mmu(1:3,1:3)));
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