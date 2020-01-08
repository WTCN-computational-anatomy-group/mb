function varargout = spm_mb_animate(varargin)
%__________________________________________________________________________
%
% Functions for creating nice visualisations.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_animate
    error('Not enough argument. Type ''help spm_mb_shape'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'AnimatePC'
        [varargout{1:nargout}] = AnimatePC(varargin{:});
    case 'Frames2Movie'
        [varargout{1:nargout}] = Frames2Movie(varargin{:});
    case 'LoadModel'
        [varargout{1:nargout}] = LoadModel(varargin{:});
    case 'SampleBrain'
        [varargout{1:nargout}] = SampleBrain(varargin{:});
    case 'SampleWarp'
        [varargout{1:nargout}] = SampleWarp(varargin{:});
    otherwise
        help spm_mb_animate
        error('Unknown function %s. Type ''help spm_mb_animate'' for help.', id)
end
end
%==========================================================================

%==========================================================================
function anim_model = LoadModel(shape_model,ftemplate,fsubspace)
% FORMAT model = LoadModel(model,[ftemplate],[fsubspace])
% INPUT:  model generated using spm_mb_fit (model.shape).
% OUTPUT: structure that can be used with the animate functions.

if nargin < 3 || isempty(fsubspace)
    fsubspace = shape_model.subspace;
end
if nargin < 2 || isempty(ftemplate)
    ftemplate = shape_model.template;
end

template = nifti(ftemplate);
if template.dat.dim(4) == 1, template.dat.dim(4) = []; end
mat = template.mat;

subspace = nifti(fsubspace);
if subspace.dat.dim(4) == 1, subspace.dat.dim(4) = []; end

vs  = sqrt(sum(mat(1:3,1:3).^2));
reg = shape_model.v_settings;

anim_model     = struct;
anim_model.A   = shape_model.latent_prior;
anim_model.U   = subspace.dat;
anim_model.mu  = template.dat;
anim_model.mat = mat;
anim_model.prm = [vs reg];
anim_model.is_anim = true;

end
%==========================================================================

%==========================================================================
function frames = AnimatePC(model, sl, ax, pc, nsigma, nframes)
% FORMAT frames = AnimatePC(model, sl, ax, pc, nsigma, nframes)
if nargin < 6 || ~isfinite(nframes), nframes = 65;   end
if nargin < 5 || ~isfinite(nsigma),  nsigma  = 3;    end
if nargin < 4 || ~isfinite(pc),      pc      = 1;    end
if nargin < 6 || ~isfinite(ax),      ax      = 3;    end
if nargin < 2,                       sl       = NaN;  end
if ~isfield(model, 'is_anim'), model = LoadModel(model); end

sigma = (diag(inv(model.A))).^(0.5);
sigma = sigma(pc);

z   = linspace(-nsigma*sigma,+nsigma*sigma,nframes);
u   = single(model.U(:,:,:,:,pc));
mu0 = single(model.mu());

if ~isfinite(sl), sl = ceil(size(mu0,ax)/2); end

spm_diffeo('boundary', 0);
f = spm_mb_show('SetFigure', '(spm_mb) Animate PC');
frames(numel(z)) = struct('cdata',[],'colormap',[]);
    
for i=1:numel(z)
    y  = spm_shoot3d(1E10*z(i)*u, model.prm, Inf);
    switch ax
        case 1
            y = y(sl,:,:,:);
        case 2
            y = y(:,sl,:,:);
        case 3
            y = y(:,:,sl,:);
        otherwise
            error('Axis must be 1, 2 or 3.')
    end
    mu = spm_diffeo('pull', mu0, y);
    mu = spm_padarray(mu, [0 0 0 1], 0, 'post');
    mu = spm_mb_shape('Softmax', mu, 4);
    p  = spm_mb_show('ShowCat', mu, ax);
    frames(i) = getframe(p.Parent);
end
end
%==========================================================================

%==========================================================================
function [mu,z] = SampleBrain(model)
% FORMAT [wmu,z] = SampleBrain(model)

spm_diffeo('boundary', 0);

% -------------------------------------------------------------------------
% Sample warp
[y,~,z] = SampleWarp(model.subspace, model.prm, model.A);

% -------------------------------------------------------------------------
% Deform template
mu = single(model.mu());
mu = spm_diffeo('pull', mu, y);
mu = spm_padarray(mu, [0 0 0 1], 0, 'post');
mu = spm_mb_shape('Softmax', mu, 4);

end
%==========================================================================

%==========================================================================
function [y,iy,z] = SampleWarp(subspace, prm, AorZ)
% FORMAT [iy,y] = SampleWarp(U, prm, [A])
% FORMAT [iy,y] = SampleWarp(U, prm, z)
%   U   - [Nx Ny Nz 3 M] Subspace (array or file_array)
%   prm - prm(1:3) = voxel size / prm(4:8) = regularisation settings
%   A   - [M M] Latent precision matrix [default: eye]
%   z   - [M 1] Latent code.
%   y   - [Nx Ny Nz 3] Forward deformation (used to forward deform volumes)
%   iy  - [Nx Ny Nz 3] Inverse deformation (used to forward deform meshes)
%   z   - [M 1] Sampled latent code

% -------------------------------------------------------------------------
% Load model
dim = [size(subspace) 1 1 1 1];
if dim(4) == 1
    dim(4) = [];
    subspace = reshape(subspace, dim);  % Remove (empty) 4th dimension
end
M   = dim(5);                           % Number of principal components       
dim = dim(1:3);                         % Dimensions of the lattice

% -------------------------------------------------------------------------
% Sample latent code
if nargin < 3
    AorZ = eye(M);
end
if size(AorZ,1) == size(AorZ,2)
    A = AorZ;
    [U,S] = svd(inv(A));
    S = S^(0.5);
    z = randn(M,1);
    z = U*S*z;
else
    z = AorZ;
end

% -------------------------------------------------------------------------
% Sample initial velocity
v = zeros([dim 3], 'single');
for m=1:M
    if z(m)
        v = v + z(m) * single(subspace(:,:,:,:,m));
    end
end

% -------------------------------------------------------------------------
% Shoot deformations
spm_diffeo('boundary', 0);
if nargout == 1
    y = spm_shoot3d(v, prm, Inf);
else
    [y,~,~,iy] = spm_shoot3d(v, prm, Inf);
end
end
%==========================================================================

%==========================================================================
function Frames2Movie(frames, fname, duration, rebound, middle)
% FORMAT Frames2Movie(frames, fname, duration, loop, middle)
% frames   - Frame array built using `Animate*`
% fname    - Ouptut filename. Accepeted extensions: mp4, avi, gif
% duration - Total duration in seconds                          [10]
% rebound  - Rebound animation (plays once in each direction)   [false]
% middle   - Start from the middle frame                        [false]
if nargin < 5 || isempty(middle)
    middle = false;
    if nargin < 4 || isempty(rebound)
        rebound = false;
        if nargin < 3 || isempty(duration)
            duration = 10;
        end
    end
end

n = numel(frames);

if      rebound &&  middle, iter = [ceil(n/2):n n:-1:1 1:ceil(n/2)];
elseif  rebound && ~middle, iter = [1:n n:-1:1];
elseif ~rebound &&  middle, iter = [ceil(n/2):n n:-1:ceil(n/2)];
elseif ~rebound && ~middle, iter = 1:n;
end
framerate = numel(iter)/duration;
if endsWith(fname, 'mp4'), vidopt = {'MPEG-4'};
else,                      vidopt = {};
end
if endsWith(fname, 'avi') || endsWith(fname, 'mp4')
    vid = VideoWriter(fname, vidopt{:});
    vid.FrameRate = framerate;
    open(vid);
end
first_frame = true;
for i=iter
    if endsWith(fname, 'gif')
        [imind, cm] = rgb2ind(frame2im(frames(i)), 256);
        if first_frame
            imwrite(imind, cm, fname, 'gif', 'Loopcount', inf, 'DelayTime', 1/framerate); 
            first_frame = false;
        else
            imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append', 'DelayTime', 1/framerate); 
        end
    elseif endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
        writeVideo(vid, frames(i));
    end
end
if endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
    close(vid);
end
end
%==========================================================================
