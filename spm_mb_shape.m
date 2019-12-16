function varargout = spm_mb_shape(varargin)
%__________________________________________________________________________
%
% Functions for shape model related.
%
% FORMAT psi0          = spm_mb_shape('Affine',d,Mat)
% FORMAT psi           = spm_mb_shape('Compose',psi1,psi0)
% FORMAT id            = spm_mb_shape('Identity',d)
% FORMAT dat           = spm_mb_shape('Init',dat,sett)
% FORMAT [dat,mu,sett] = spm_mb_shape('InitMu',dat,K,samp_mu,sett)
% FORMAT l             = spm_mb_shape('LSE',mu,dr)
% FORMAT a1            = spm_mb_shape('Pull1',a0,psi,r)
% FORMAT [f1,w1]       = spm_mb_shape('Push1',f,psi,d,r)
% FORMAT mu1           = spm_mb_shape('ShrinkTemplate',mu,oMmu,sett)
% FORMAT P             = spm_mb_shape('Softmax',mu,dr)
% FORMAT [Mmu,d]       = spm_mb_shape('SpecifyMean',dat,vx)
% FORMAT E             = spm_mb_shape('TemplateEnergy',mu,sett)
% FORMAT dat           = spm_mb_shape('UpdateAffines',dat,mu,sett)
% FORMAT [mu,dat]      = spm_mb_shape('UpdateMean',dat, mu, sett)
% FORMAT dat           = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett)
% FORMAT [mu,dat]      = spm_mb_shape('UpdateSimpleMean',dat, mu, sett)
% FORMAT dat           = spm_mb_shape('UpdateVelocities',dat,mu,sett)
% FORMAT dat           = spm_mb_shape('UpdateWarps',dat,sett)
% FORMAT dat           = spm_mb_shape('VelocityEnergy',dat,sett)
% FORMAT [dat,mu]      = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_shape
    error('Not enough argument. Type ''help spm_mb_shape'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id    
    case 'Affine'
        [varargout{1:nargout}] = Affine(varargin{:});             
    case 'Compose'
        [varargout{1:nargout}] = Compose(varargin{:});
    case 'Identity'
        [varargout{1:nargout}] = Identity(varargin{:});
    case 'Init'
        [varargout{1:nargout}] = Init(varargin{:});                
    case 'InitMu'
        [varargout{1:nargout}] = InitMu(varargin{:});         
    case 'LSE'
        [varargout{1:nargout}] = LSE(varargin{:});
    case 'Pull1'
        [varargout{1:nargout}] = Pull1(varargin{:});
    case 'Push1'
        [varargout{1:nargout}] = Push1(varargin{:});    
    case 'ShrinkTemplate'
        [varargout{1:nargout}] = ShrinkTemplate(varargin{:});
    case 'Softmax'
        [varargout{1:nargout}] = Softmax(varargin{:});       
    case 'SpecifyMean'
        [varargout{1:nargout}] = SpecifyMean(varargin{:});       
    case 'TemplateEnergy'
        [varargout{1:nargout}] = TemplateEnergy(varargin{:});          
    case 'UpdateAffines'
        [varargout{1:nargout}] = UpdateAffines(varargin{:});                   
    case 'UpdateMean'
        [varargout{1:nargout}] = UpdateMean(varargin{:});
    case 'UpdateSimpleAffines'
        [varargout{1:nargout}] = UpdateSimpleAffines(varargin{:});
    case 'UpdateSimpleMean'
        [varargout{1:nargout}] = UpdateSimpleMean(varargin{:});
    case 'UpdateVelocities'
        [varargout{1:nargout}] = UpdateVelocities(varargin{:});
    case 'UpdateWarps'
        [varargout{1:nargout}] = UpdateWarps(varargin{:});  
    case 'ZoomVolumes'
        [varargout{1:nargout}] = ZoomVolumes(varargin{:});        
    case 'VelocityEnergy'
        [varargout{1:nargout}] = VelocityEnergy(varargin{:});              
    otherwise
        help spm_mb_shape
        error('Unknown function %s. Type ''help spm_mb_shape'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% Affine()
function psi0 = Affine(d,Mat)
id    = Identity(d);
psi0  = reshape(reshape(id,[prod(d) 3])*Mat(1:3,1:3)' + Mat(1:3,4)',[d 3]);
end
%==========================================================================

%==========================================================================
% Compose()
function psi = Compose(psi1,psi0)
psi = spm_diffeo('comp',psi1,psi0);
end
%==========================================================================

%==========================================================================
% Identity()
function id = Identity(d)
id = zeros([d(:)',3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
end
%==========================================================================

%==========================================================================
% Init()
function dat = Init(dat,sett)

% Parse function settings
B       = sett.registr.B;
d       = sett.var.d;
dir_res = sett.write.dir_res;
Mmu     = sett.var.Mmu;

v    = zeros([d,3],'single');
psi1 = spm_mb_shape('Identity',d);
for n=1:numel(dat)
    dat(n).q = zeros(size(B,3),1);
    if isnumeric(dat(n).f)
        dat(n).v   = v;
        dat(n).psi = psi1;
    else
        if isa(dat(n).f,'nifti')
            [~,nam,~] = fileparts(dat(n).f(1).dat.fname);
            vname    = fullfile(dir_res,['v_' nam '.nii']);
            pname    = fullfile(dir_res,['psi_' nam '.nii']);
            
            fa       = file_array(vname,[d(1:3) 1 3],'float32',0);
            nii      = nifti;
            nii.dat  = fa;
            nii.mat  = Mmu;
            nii.mat0 = Mmu;
            nii.descrip = 'Velocity';
            create(nii);
            nii.dat(:,:,:,:) = v;
            dat(n).v    = nii;

            nii.dat.fname = pname;
            nii.descrip = 'Deformation (WIP)';
            create(nii);
            nii.dat(:,:,:,:) = psi1;
            dat(n).psi  = nii;
        end
    end
end
end
%==========================================================================

%==========================================================================
% InitMu()
function [dat,mu,sett] = InitMu(dat,K,samp_mu,sett)
% Make 'quick' initial estimates of GMM posteriors and template on very coarse
% scale

% Parse function settings
do_gmm = sett.do.gmm;

% Uniform template
mu = zeros([sett.var.d K],'single');

if ~do_gmm, return; end

% Change some settings
do_updt_bf      = sett.do.updt_bf;
samp            = sett.gen.samp;
nit_appear      = sett.nit.appear;
nit_gmm         = sett.nit.gmm;
sett.nit.gmm    = 200;
sett.gen.samp   = min(samp_mu,5);
sett.nit.appear = 1;
sett.do.updt_bf = false;

% Get population indices
p_ix       = spm_mb_appearance('GetPopulationIdx',dat);
Npop       = numel(p_ix);
% first_subj = true;
for p=2:Npop
    % To make the algorithm more robust when using multiple populations,
    % set posterior and prior means (m) of GMMs of all but the first population to
    % uniform  
    avg_po = 0;
    avg_pr = 0;
    for n=p_ix{p}
        avg_po = avg_po + dat(n).mog.po.m;
        avg_pr = avg_pr + dat(n).mog.pr.m;
    end
    avg_po = avg_po./numel(p_ix{p});
    avg_po = mean(avg_po,2);
    avg_pr = avg_pr./numel(p_ix{p});
    avg_pr = mean(avg_pr,2);
    
    for n=p_ix{p}
        dat(n).mog.pr.m = repmat(avg_pr,[1 K + 1]);        
%         if first_subj            
%             first_subj = false;
%             continue
%         end        
        dat(n).mog.po.m = repmat(avg_po,[1 K + 1]);        
    end
end

% Update template based on only first subject..
% [mu,dat(1)]       = spm_mb_shape('UpdateSimpleMean',dat(1), mu, sett);
% ..then propagate to all other subjects in populations..
[mu,dat(p_ix{1})] = spm_mb_shape('UpdateSimpleMean',dat(p_ix{1}), mu, sett);
if Npop > 1
    % ..if more than one population, use template learned on first
    % population to initialise other populations' subjects
    [mu,dat]      = spm_mb_shape('UpdateSimpleMean',dat,    mu, sett);    
end

% % Update template based on only first population
% [mu,dat(p_ix{1})] = spm_mb_shape('UpdateSimpleMean',dat(p_ix{1}), mu, sett);
% if Npop > 1
%     % Use template learned from first population to initialise all other
%     % populations' subjects
%     [mu,dat]      = spm_mb_shape('UpdateSimpleMean',dat,    mu, sett);    
% end

% Restore settings
sett.do.updt_bf = do_updt_bf;
sett.gen.samp   = samp;
sett.nit.appear = nit_appear;
sett.nit.gmm    = nit_gmm;
end
%==========================================================================

%==========================================================================
% LSE()
function l = LSE(mu,dr)
mx = max(mu,[],dr);
l  = log(exp(-mx) + sum(exp(mu - mx),dr)) + mx;
end
%==========================================================================

%==========================================================================
% Pull1()
function a1 = Pull1(a0,psi,r)
% Resample an image or set of images
% FORMAT a1 = Pull1(a0,psi,r)
%
% a0  - Input image(s)
% psi - Deformation
% r   - subsampling density in each dimension (default: [1 1 1])
%
% a1  - Output image(s)
%
% There are also a couple of Push1 and Pull1 functions, which might be of
% interest.  The idea is to reduce aliasing effects in the pushed images,
% which might also be useful for dealing with thick-sliced images.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$
if nargin<3, r=[1 1 1]; end

if isempty(a0)
    a1 = a0;
elseif isempty(psi)
    a1 = a0;
else
    if r==1
        a1 = spm_diffeo('pullc',a0,psi);
        return
    end
    d  = [size(a0) 1 1];
    if d(3)>1, zrange = Range(r(3)); else, zrange = 0; end
    if d(2)>1, yrange = Range(r(2)); else, yrange = 0; end
    if d(1)>1, xrange = Range(r(1)); else, xrange = 0; end
    dp = size(psi); 
    id = Identity(dp(1:3));     
    a1 = zeros([size(psi,1),size(psi,2),size(psi,3),size(a0,4)],'single');
    for l=1:d(4)
        tmp = single(0);
        al  = single(a0(:,:,:,l));
        for dz=zrange
            for dy=yrange
                for dx=xrange
                    ids  = id  + cat(4,dx,dy,dz);
                    psi1 = spm_diffeo('pull',psi-id,    ids)+ids;
                    as   = spm_diffeo('pull',al,psi1);
                   %ids  = id  - cat(4,dx,dy,dz);
                    tmp  = tmp + spm_diffeo('push',as,  ids);
                end
            end
        end
        a1(:,:,:,l) = tmp/(numel(zrange)*numel(yrange)*numel(xrange));
    end
end
end
%==========================================================================

%==========================================================================
% Push1()
function [f1,w1] = Push1(f,psi,d,r)
% Push an image (or set of images) accorging to a spatial transform
% FORMAT [f1,w1] = Push1(f,psi,d,r)
%
% f   - Image (3D or 4D)
% psi - Spatial transform
% d   - dimensions of output (default: size of f)
% r   - subsampling density in each dimension (default: [1 1 1])
%
% f1  - "Pushed" image
%
% There are also a couple of Push1 and Pull1 functions, which might be of
% interest.  The idea is to reduce aliasing effects in the pushed images,
% which might also be useful for dealing with thick-sliced images.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$
if nargin<4, r = [1 1 1]; end
if nargin<3, d = [size(f,1) size(f,2) size(f,3)]; end

%msk    = isfinite(f);
%f(msk) = 0;
if ~isempty(psi)
    if r==1
        if nargout==1
            f1      = spm_diffeo('pushc',single(f),psi,d);
        else
            [f1,w1] = spm_diffeo('pushc',single(f),psi,d);
        end
        return
    end

    if d(3)>1, zrange = Range(r(3)); else, zrange = 0; end
    if d(2)>1, yrange = Range(r(2)); else, yrange = 0; end
    if d(1)>1, xrange = Range(r(1)); else, xrange = 0; end
    
    dp    = size(psi);
    id    = Identity(dp(1:3));
    f1    = single(0);
    w1    = single(0);
    for dz=zrange
        for dy=yrange
            for dx=xrange
                ids       = id + cat(4,dx,dy,dz);
                psi1      = spm_diffeo('pull',psi-id,    ids)+ids;
                fs        = spm_diffeo('pull',single(f), ids);
               %fs=single(f);
                if nargout==1
                    fs        = spm_diffeo('push',fs,        psi1,d);
                    f1        = f1  + fs;
                else
                    [fs,ws]   = spm_diffeo('push',fs,        psi1,d);
                    f1        = f1  + fs;
                    w1        = w1  + ws;
                end
            end
        end
    end
    scale = 1/(numel(zrange)*numel(yrange)*numel(xrange));
    f1    = f1*scale;
    w1    = w1*scale;
else
    msk      = isfinite(f);
    f1       = f;
    f1(~msk) = 0;
    w1       = single(all(msk,4));
end
end
%==========================================================================

%==========================================================================
% ShrinkTemplate()
function mu1 = ShrinkTemplate(mu,oMmu,sett)

% Parse function settings
d   = sett.var.d;
Mmu = sett.var.Mmu;

d0      = [size(mu,1) size(mu,2) size(mu,3)];
Mzoom   = Mmu\oMmu;
if norm(Mzoom-eye(4))<1e-4 && all(d0==d)
    mu1 = mu;
else
    y       = reshape(reshape(Identity(d0),[prod(d0),3])*Mzoom(1:3,1:3)'+Mzoom(1:3,4)',[d0 3]);
    [mu1,c] = Push1(mu,y,d);
    mu1     = mu1./(c+eps);
end
end
%==========================================================================

%==========================================================================
% Softmax()
function P = Softmax(mu,dr)
mx  = max(mu,[],dr);
E   = exp(mu-mx);
den = sum(E,dr)+exp(-mx);
P   = E./den;
end
%==========================================================================

%==========================================================================
% SpecifyMean()
function [Mmu,d] = SpecifyMean(dat,vx)
dims = zeros(numel(dat),3);
Mat0 = zeros(4,4,numel(dat));
for n=1:numel(dat)
    dims(n,:)   = spm_mb_io('GetSize',dat(n).f)';    
    Mat0(:,:,n) = dat(n).Mat;
end

[Mmu,d] = ComputeAvgMat(Mat0,dims);

% Adjust voxel size
if numel(vx) == 1
    vx = vx.*ones([1 3]);
end
vxmu = sqrt(sum(Mmu(1:3,1:3).^2));
samp = vxmu./vx;   
D    = diag([samp 1]);
Mmu  = Mmu/D;
d    = floor(D(1:3,1:3)*d')';

if unique(dims(:,3),'rows') == 1
    % 2D
    bb      = [-Inf Inf ; -Inf Inf ; 1 1]';
    bb      = round(bb);
    bb      = sort(bb);
    bb(1,:) = max(bb(1,:),[1 1 1]);
    bb(2,:) = min(bb(2,:),d(1:3));
    
    d(1:3) = diff(bb) + 1;
    Mmu    = Mmu*spm_matrix((bb(1,:)-1));
end
end
%==========================================================================

%==========================================================================
% TemplateEnergy()
function E = TemplateEnergy(mu,sett)

% Parse function settings
do_updt_template = sett.do.updt_template;
mu_settings      = sett.var.mu_settings;

if do_updt_template
    g = spm_field('vel2mom', mu, mu_settings);
    E = 0.5*mu(:)'*g(:);
else
    E = 0;
end
end
%==========================================================================

%==========================================================================
% UpdateAffines()
function dat = UpdateAffines(dat,mu,sett)

% Parse function settings
B           = sett.registr.B;
do_updt_aff = sett.do.updt_aff;
groupwise   = sett.model.groupwise;

if ~do_updt_aff, return; end

% Update the affine parameters
if ~isempty(B)
    for n=1:numel(dat)
        dat(n) = UpdateAffinesSub(dat(n),mu,sett);
    end

    if groupwise
        % Zero-mean the affine parameters
        mq = sum(cat(2,dat(:).q),2)/numel(dat);
        for n=1:numel(dat)
            dat(n).q = dat(n).q - mq;
        end
    end
end
end
%==========================================================================

%==========================================================================
% UpdateMean()
function [mu,dat] = UpdateMean(dat, mu, sett)

% Parse function settings
accel            = sett.gen.accel;
do_updt_template = sett.do.updt_template;
mu_settings      = sett.var.mu_settings;
s_settings       = sett.shoot.s_settings;

if ~do_updt_template, return; end

g  = spm_field('vel2mom', mu, mu_settings);
M  = size(mu,4);
H  = zeros([sett.var.d M*(M+1)/2],'single');
H0 = AppearanceHessian(mu,accel);
for n=1:numel(dat)
    [gn,Hn,dat(n)] = UpdateMeanSub(dat(n),mu,H0,sett);
    g              = g + gn;
    H              = H + Hn;
end
mu = mu - spm_field(H, g, [mu_settings s_settings]);  
end
%==========================================================================

%==========================================================================
% UpdateSimpleAffines()
function dat = UpdateSimpleAffines(dat,mu,sett)

% Parse function settings
accel       = sett.gen.accel;
B           = sett.registr.B;
do_updt_aff = sett.do.updt_aff;
groupwise   = sett.model.groupwise;

if ~do_updt_aff, return; end

% Update the affine parameters
G  = spm_diffeo('grad',mu);
H0 = VelocityHessian(mu,G,accel);

if ~isempty(B)
    for n=1:numel(dat)
        dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
    end

    if groupwise
        % Zero-mean the affine parameters
        mq = sum(cat(2,dat(:).q),2)/numel(dat);
        for n=1:numel(dat)
            dat(n).q = dat(n).q - mq;
        end
    end
end
end
%==========================================================================

%==========================================================================
% UpdateSimpleMean()
function [mu,dat] = UpdateSimpleMean(dat, mu, sett)

% Parse function settings
accel       = sett.gen.accel;
mu_settings = sett.var.mu_settings;
s_settings  = sett.shoot.s_settings;

w  = zeros(sett.var.d,'single');
gf = zeros(size(mu),'single');
for n=1:numel(dat)
    [gn,wn,dat(n)] = UpdateSimpleMeanSub(dat(n),mu,sett);
    gf             = gf + gn;
    w              = w  + wn;
end
for it=1:ceil(4+2*log2(numel(dat)))
    H  = w.*AppearanceHessian(mu,accel);
    g  = w.*Softmax(mu,4) - gf;
    g  = g  + spm_field('vel2mom', mu, mu_settings);
    mu = mu - spm_field(H, g, [mu_settings s_settings]);
end
end
%==========================================================================

%==========================================================================
% UpdateVelocities()
function dat = UpdateVelocities(dat,mu,sett)

% Parse function settings
accel       = sett.gen.accel;
do_updt_vel = sett.do.updt_vel;

if ~do_updt_vel, return; end

G  = spm_diffeo('grad',mu);
H0 = VelocityHessian(mu,G,accel);
if size(G,3) == 1
    % Data is 2D -> add some regularisation
    H0(:,:,:,3) = H0(:,:,:,3) + mean(reshape(H0(:,:,:,[1 2]),[],1));
end
for n=1:numel(dat)
    dat(n) = UpdateVelocitiesSub(dat(n),mu,G,H0,sett);
end
end
%==========================================================================

%==========================================================================
% UpdateWarps()
function dat = UpdateWarps(dat,sett)

% Parse function settings
groupwise  = sett.model.groupwise;
v_settings = sett.var.v_settings;

if groupwise
    % Total initial velocity should be zero (Khan & Beg)
    avg_v = single(0);
    for n=1:numel(dat)
        avg_v = avg_v + spm_mb_io('GetData',dat(n).v); % For mean correcting initial velocities
    end
    avg_v = avg_v/numel(dat);
    d     = [size(avg_v,1) size(avg_v,2) size(avg_v,3)];
else
    avg_v = [];
    d     = spm_mb_io('GetSize',dat(1).v);
end
kernel = Shoot(d,v_settings);
for n=1:numel(dat)
    dat(n) = UpdateWarpsSub(dat(n),avg_v,sett,kernel);
end
end
%==========================================================================

%==========================================================================
% ZoomVolumes()
function [dat,mu] = ZoomVolumes(dat,mu,sett,oMmu)

% Parse function settings
d       = sett.var.d;
Mmu     = sett.var.Mmu;

d0    = [size(mu,1) size(mu,2) size(mu,3)];
z     = single(reshape(d./d0,[1 1 1 3]));
Mzoom = oMmu\Mmu;
y     = reshape(reshape(Identity(d),[prod(d),3])*Mzoom(1:3,1:3)' + Mzoom(1:3,4)',[d 3]);
mu    = spm_diffeo('pullc',mu,y);
for n=1:numel(dat)
    v          = spm_mb_io('GetData',dat(n).v);
    v          = spm_diffeo('pullc',v,y).*z;
    dat(n).v   = ResizeFile(dat(n).v  ,d,Mmu);
    dat(n).psi = ResizeFile(dat(n).psi,d,Mmu);
    dat(n).v   = spm_mb_io('SetData',dat(n).v,v);
end
end
%==========================================================================

%==========================================================================
% VelocityEnergy()
function dat = VelocityEnergy(dat,sett)

% Parse function settings
v_settings = sett.var.v_settings;

for n=1:numel(dat)
    v           = spm_mb_io('GetData',dat(n).v);
    u0          = spm_diffeo('vel2mom', v, v_settings); % Initial momentum
    dat(n).E(2) = 0.5*sum(u0(:).*v(:));                 % Prior term
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%========================================================================== 
% AffineHessian()
function [H,g] = AffineHessian(mu,G,a,w,accel)
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = Horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(1:d(1),1:d(2),1,1);
for i=1:d(3)
    x{3} = x{3}*0+i;
    gv   = reshape(sum(a(:,:,i,:).*G(:,:,i,:,:),4),[d(1:2) 1 3]);
    Hv   = w(:,:,i).*VelocityHessian(mu(:,:,i,:),G(:,:,i,:,:),accel);
    for i1=1:12
        k1g   = rem(i1-1,3)+1;
        k1x   = floor((i1-1)/3)+1;
        g(i1) = g(i1) + sum(sum(sum(x{k1x}.*gv(:,:,:,k1g))));
        for i2=1:12
            k2g      = rem(i2-1,3)+1;
            k2x      = floor((i2-1)/3)+1;
            H(i1,i2) = H(i1,i2) + sum(sum(sum(x{k1x}.*Hv(:,:,:,I(k1g,k2g)).*x{k2x})));
        end
    end
end
end
%========================================================================== 

%==========================================================================
% AppearanceHessian()
function H = AppearanceHessian(mu,accel)
M  = size(mu,4);
d  = [size(mu,1) size(mu,2) size(mu,3)];
if accel>0, s  = spm_mb_shape('Softmax',mu,4); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
I  = Horder(M);
H  = zeros([d (M*(M+1))/2],'single');
for m1=1:M
    for m2=m1:M
        if accel==0
            tmp = Ab(m1,m2)*ones(d,'single');
        else
            if m2~=m1
                tmp = accel*(-s(:,:,:,m1).*s(:,:,:,m2))           + (1-accel)*Ab(m1,m2);
            else
                tmp = accel*(max(s(:,:,:,m1).*(1-s(:,:,:,m1)),0)) + (1-accel)*Ab(m1,m2);
            end
        end
        H(:,:,:,I(m1,m2)) = tmp;
    end
end
end
%==========================================================================

%==========================================================================
% ComputeAvgMat()
function [M_avg,d] = ComputeAvgMat(Mat0,dims)
% Compute an average voxel-to-world mapping and suitable dimensions
% FORMAT [M_avg,d] = spm_compute_avg_mat(Mat0,dims)
% Mat0  - array of matrices (4x4xN)
% dims  - image dimensions (Nx3)
% M_avg - voxel-to-world mapping
% d     - dimensions for average image
%
%__________________________________________________________________________
% Copyright (C) 2012-2019 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


% Rigid-body matrices computed from exp(p(1)*B(:,:,1)+p(2)+B(:,:,2)...)
%--------------------------------------------------------------------------
B = spm_mb_param('AffineBases','SE(3)');

% Find combination of 90 degree rotations and flips that brings all
% the matrices closest to axial
%--------------------------------------------------------------------------
Matrices = Mat0;
pmatrix  = [1,2,3; 2,1,3; 3,1,2; 3,2,1; 1,3,2; 2,3,1];
for i=1:size(Matrices,3)
    vx    = sqrt(sum(Matrices(1:3,1:3,i).^2));
    tmp   = Matrices(:,:,i)/diag([vx 1]);
    R     = tmp(1:3,1:3);
    minss = Inf;
    minR  = eye(3);
    for i1=1:6
        R1 = zeros(3);
        R1(pmatrix(i1,1),1)=1;
        R1(pmatrix(i1,2),2)=1;
        R1(pmatrix(i1,3),3)=1;
        for i2=0:7
            F  = diag([bitand(i2,1)*2-1, bitand(i2,2)-1, bitand(i2,4)/2-1]);
            R2 = F*R1;
            ss = sum(sum((R/R2-eye(3)).^2));
            if ss<minss
                minss = ss;
                minR  = R2;
            end
        end
    end
    rdim = abs(minR*dims(i,:)');
    R2   = inv(minR);
    minR = [R2 R2*((sum(R2,1)'-1)/2.*(rdim+1)); 0 0 0 1];
    Matrices(:,:,i) = Matrices(:,:,i)*minR;
end

% Average of these matrices
%--------------------------------------------------------------------------
M_avg = spm_meanm(Matrices);

% If average involves shears, then find the closest matrix that does not
% require them
%--------------------------------------------------------------------------
p = spm_imatrix(M_avg);
if sum(p(10:12).^2)>1e-8

    % Zooms computed from exp(p(7)*B2(:,:,1)+p(8)*B2(:,:,2)+p(9)*B2(:,:,3))
    %----------------------------------------------------------------------
    B2        = zeros(4,4,3);
    B2(1,1,1) = 1;
    B2(2,2,2) = 1;
    B2(3,3,3) = 1;

    p      = zeros(9,1); % Parameters
    for it=1:10000
        [R,dR] = spm_dexpm(p(1:6),B);  % Rotations + Translations
        [Z,dZ] = spm_dexpm(p(7:9),B2); % Zooms

        M  = R*Z; % Voxel-to-world estimate
        dM = zeros(4,4,6);
        for i=1:6, dM(:,:,i)   = dR(:,:,i)*Z; end
        for i=1:3, dM(:,:,i+6) = R*dZ(:,:,i); end
        dM = reshape(dM,[16,9]);

        d   = M(:)-M_avg(:); % Difference
        gr  = dM'*d;         % Gradient
        Hes = dM'*dM;        % Hessian
        p   = p - Hes\gr;    % Gauss-Newton update
        if sum(gr.^2)<1e-8, break; end
    end
    M_avg = M;
end

% Ensure that the FoV covers all images, with a few voxels to spare
%--------------------------------------------------------------------------
mn    =  Inf*ones(3,1);
mx    = -Inf*ones(3,1);
for i=1:size(Mat0,3)
    dm      = [dims(i,:) 1 1];
    corners = [
        1 dm(1)    1  dm(1)   1  dm(1)    1  dm(1)
        1    1  dm(2) dm(2)   1     1  dm(2) dm(2)
        1    1     1     1 dm(3) dm(3) dm(3) dm(3)
        1    1     1     1    1     1     1     1];
    M  = M_avg\Mat0(:,:,i);
    vx = M(1:3,:)*corners;
    mx = max(mx,max(vx,[],2));
    mn = min(mn,min(vx,[],2));
end
mx    = ceil(mx);
mn    = floor(mn);
o     = 3;
d     = (mx-mn+(2*o+1))';
M_avg = M_avg * [eye(3) mn-(o+1); 0 0 0 1];
end
%==========================================================================

%==========================================================================
% Horder()
function I = Horder(d)
I = diag(1:d);
l = d;
for i1=1:d
    for i2=(i1+1):d
        l = l + 1;
        I(i1,i2) = l;
        I(i2,i1) = l;
    end
end
end
%==========================================================================

%==========================================================================
% Mask()
function f = Mask(f,msk)
f(~isfinite(f)) = 0;
f = f.*msk;
end
%==========================================================================

%==========================================================================
% Range()
function r = Range(n)
r = (-floor((n-1)/2):ceil((n-1)/2))/n;
end
%==========================================================================

%==========================================================================
% ResizeFile()
function fin = ResizeFile(fin,d,Mat)
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    for m=1:numel(fin)
        fin(m).dat.dim(1:3) = d(1:3);
        fin(m).mat  = Mat;
        fin(m).mat0 = Mat;
        create(fin(m));
    end
end
end
%==========================================================================

%==========================================================================
% Shoot()
function varargout = Shoot(v0,kernel,args)
% Geodesic shooting
% FORMAT psi = Shoot(v0,kernel,args)
%
% v0       - Initial velocity field n1*n2*n3*3 (single prec. float)
% kernel   - structure created previously
% args     - Integration parameters
%            - [1] Num time steps
%
% psi      - Inverse deformation field n1*n2*n3*3 (single prec. float)
%
% FORMAT kernel = Shoot(d,v_settings)
% d          - dimensions of velocity fields
% v_settings - 8 settings
%              - [1][2][3] Voxel sizes
%              - [4][5][6][7][8] Regularisation settings.
%              Regularisation uses the sum of
%              - [4] - absolute displacements
%              - [5] - laplacian
%              - [6] - bending energy
%              - [7] - linear elasticity mu
%              - [8] - linear elasticity lambda
%
% kernel     - structure encoding Greens function
% 
% This code generates inverse deformations from
% initial velocity fields by gedesic shooting.  See the work of Miller,
% Younes and others.
%
% LDDMM (Beg et al) uses the following evolution equation:
%     d\phi/dt = v_t(\phi_t)
% where a variational procedure is used to find the stationary solution
% for the time varying velocity field.
% In principle though, once the initial velocity is known, then the
% velocity at subsequent time points can be computed.  This requires
% initial momentum (u_0), computed (using differential operator L) by:
%     u_0 = L v_0
% Then (Ad_{\phi_t})^* m_0 is computed:
%     u_t = |d \phi_t| (d\phi_t)^T u_0(\phi_t)
% The velocity field at this time point is then obtained by using
% multigrid to solve:
%     v_t = L^{-1} u_t
%
% These equations can be found in:
% Younes (2007). "Jacobi fields in groups of diffeomorphisms and
% applications". Quarterly of Applied Mathematics, vol LXV,
% number 1, pages 113-134 (2007).
%
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2013-2017)

% John Ashburner
% $Id$

if nargin==2
    if numel(v0)>5
        d = [size(v0) 1];
        d = v0(1:3);
    else
        d = v0;
    end
    v_settings = kernel;
    spm_diffeo('boundary',0);
    F   = spm_shoot_greens('kernel',d,v_settings);
    varargout{1} = struct('d',d, 'v_settings',v_settings, 'F', F);
    return;
end

if isempty(v0)
    varargout{1} = [];
    varargout{2} = [];
    return;
end

args0 = 8;
if nargin<3
    args = args0;
else
    if numel(args)<numel(args0)
        args = [args args0((numel(args)+1):end)];
    end
end

T     = args(1);   % # Time steps
d     = size(v0);
d     = d(1:3);
id    = Identity(d);

if sum(v0(:).^2)==0
    varargout{1} = id;
    varargout{2} = v0;
end

if ~isfinite(T)
    % Number of time steps from an educated guess about how far to move
    T = double(floor(sqrt(max(max(max(v0(:,:,:,1).^2+v0(:,:,:,2).^2+v0(:,:,:,3).^2)))))+1);
end

spm_diffeo('boundary',0);
v   = v0;
u   = spm_diffeo('vel2mom',v,kernel.v_settings); % Initial momentum (u_0 = L v_0)
psi = id - v/T;

for t=2:abs(T)
    % The update of u_t is not exactly as described in the paper, but describing this might be a bit
    % tricky. The approach here was the most stable one I could find - although it does lose some
    % energy as < v_t, u_t> decreases over time steps.
    Jdp         = spm_diffeo('jacobian',id-v/T);
    u1          = zeros(size(u),'single');
    u1(:,:,:,1) = Jdp(:,:,:,1,1).*u(:,:,:,1) + Jdp(:,:,:,2,1).*u(:,:,:,2) + Jdp(:,:,:,3,1).*u(:,:,:,3);
    u1(:,:,:,2) = Jdp(:,:,:,1,2).*u(:,:,:,1) + Jdp(:,:,:,2,2).*u(:,:,:,2) + Jdp(:,:,:,3,2).*u(:,:,:,3);
    u1(:,:,:,3) = Jdp(:,:,:,1,3).*u(:,:,:,1) + Jdp(:,:,:,2,3).*u(:,:,:,2) + Jdp(:,:,:,3,3).*u(:,:,:,3);
    Jdp = [];
    u           = spm_diffeo('pushc',u1,id+v/T);

    % v_t \gets L^g u_t
    v            = spm_shoot_greens(u,kernel.F,kernel.v_settings); % Convolve with Greens function of L

    if size(v,3)==1, v(:,:,:,3) = 0; end

    % $\psi \gets \psi \circ (id - \tfrac{1}{T} v)$
    % I found that simply using $\psi \gets \psi - \tfrac{1}{T} (D \psi) v$ was not so stable.
    psi          = spm_diffeo('comp',psi,id-v/T);

    if size(v,3)==1, psi(:,:,:,3) = 1; end
end
varargout{1} = psi;
varargout{2} = v;
end
%==========================================================================

%==========================================================================
% SimpleAffineHessian()
function [H,g] = SimpleAffineHessian(mu,G,H0,a,w)
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = Horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(1:d(1),1:d(2),1,1);
for i=1:d(3)
    x{3} = x{3}*0+i;
    gv   = reshape(sum(a(:,:,i,:).*G(:,:,i,:,:),4),[d(1:2) 1 3]);
    Hv   = w(:,:,i).*H0(:,:,i,:);
    for i1=1:12
        k1g   = rem(i1-1,3)+1;
        k1x   = floor((i1-1)/3)+1;
        g(i1) = g(i1) + sum(sum(sum(x{k1x}.*gv(:,:,:,k1g))));
        for i2=1:12
            k2g      = rem(i2-1,3)+1;
            k2x      = floor((i2-1)/3)+1;
            H(i1,i2) = H(i1,i2) + sum(sum(sum(x{k1x}.*Hv(:,:,:,I(k1g,k2g)).*x{k2x})));
        end
    end
end
end
%==========================================================================

%==========================================================================
% UpdateAffinesSub()
function datn = UpdateAffinesSub(datn,mu,sett)
% This could be made more efficient.

% Parse function settings
accel = sett.gen.accel;
B     = sett.registr.B;
d     = sett.var.d;
Mmu   = sett.var.Mmu;
scal  = sett.optim.scal_q;

df   = spm_mb_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,B);
dM   = zeros(12,size(B,3));
for m=1:size(B,3)
    tmp     = Mmu\dM3(:,:,m)*Mn;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi1 = spm_mb_io('GetData',datn.psi);
psi0 = Affine(df,Mmu\Mr*Mn);
J    = spm_diffeo('jacobian',psi1);
J    = reshape(Pull1(reshape(J,[d 3*3]),psi0),[df 3 3]);
psi  = Compose(psi1,psi0);
psi0 = []; psi1 = [];  

mu1  = Pull1(mu,psi);
[f,datn] = spm_mb_io('GetClasses',datn,mu1,sett);
M    = size(mu,4);
G    = zeros([df M 3],'single');
for m=1:M
    [~,Gm{1},Gm{2},Gm{3}] = spm_diffeo('bsplins',mu(:,:,:,m),psi,[1 1 1  0 0 0]);
    for i1=1:3
        tmp = single(0);
        for j1=1:3
            tmp = tmp + J(:,:,:,j1,i1).*Gm{j1};
        end
        tmp(~isfinite(tmp)) = 0;
        G(:,:,:,m,i1) = tmp;
    end
    Gm = [];
end
J = []; mu = [];

msk       = all(isfinite(f),4);
a         = Mask(f - Softmax(mu1,4),msk);
[H,g]     = AffineHessian(mu1,G,a,single(msk),accel);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-6 + 0.1);
q         = q + scal*(H\g);
datn.q    = q;
end
%==========================================================================

%==========================================================================
% UpdateMeanSub()
function [g,H,datn] = UpdateMeanSub(datn,mu,H0,sett)

% Parse function settings
accel = sett.gen.accel;
B     = sett.registr.B;
d     = sett.var.d;
Mmu   = sett.var.Mmu;

df  = spm_mb_io('GetSize',datn.f);
q   = double(datn.q);
Mn  = datn.Mat;
psi = Compose(spm_mb_io('GetData',datn.psi),Affine(df, Mmu\spm_dexpm(q,B)*Mn));
mu  = Pull1(mu,psi);
[f,datn] = spm_mb_io('GetClasses',datn,mu,sett);
% if isempty(H0)
%     g     = Push1(Softmax(mu,4) - f,psi,d);
%     H     = Push1(AppearanceHessian(mu,accel),psi,d);
% else
    % Faster approximation - but might be unstable
    % If there are problems, then revert to the slow
    % way.
    [g,w] = Push1(Softmax(mu,4) - f,psi,d);
    H     = w.*H0;
% end
end
%==========================================================================

%==========================================================================
% UpdateSimpleAffinesSub()
function datn = UpdateSimpleAffinesSub(datn,mu,G,H0,sett)

% Parse function settings
B    = sett.registr.B;
d    = sett.var.d;
Mmu  = sett.var.Mmu;
scal = sett.optim.scal_q;

df   = spm_mb_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,B);
dM   = zeros(12,size(B,3));
for m=1:size(B,3)
    tmp     = Mmu\dM3(:,:,m)*Mmu;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi      = Affine(df,Mmu\Mr*Mn);
mu1      = Pull1(mu,psi);
[f,datn] = spm_mb_io('GetClasses',datn,mu1,sett);

[a,w]     = Push1(f - Softmax(mu1,4),psi,d);
mu1 = []; psi = []; f = [];

[H,g]     = SimpleAffineHessian(mu,G,H0,a,w);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-6 + 0.1);
q         = q + scal*(H\g);
datn.q    = q;
end
%==========================================================================

%==========================================================================
% UpdateSimpleMeanSub()
function [g,w,datn] = UpdateSimpleMeanSub(datn,mu,sett)

% Parse function settings
B   = sett.registr.B;
d   = sett.var.d;
Mmu = sett.var.Mmu;

df    = spm_mb_io('GetSize',datn.f);
q     = double(datn.q);
Mn    = datn.Mat;
psi   = Compose(spm_mb_io('GetData',datn.psi),Affine(df,Mmu\spm_dexpm(q,B)*Mn));
mu    = Pull1(mu,psi);
[f,datn] = spm_mb_io('GetClasses',datn,mu,sett);
[g,w] = Push1(f,psi,d);
end
%==========================================================================

%==========================================================================
% UpdateVelocitiesSub()
function datn = UpdateVelocitiesSub(datn,mu,G,H0,sett)

% Parse function settings
B          = sett.registr.B;
d          = sett.var.d;
Mmu        = sett.var.Mmu;
s_settings = sett.shoot.s_settings;
scal       = sett.optim.scal_v;
v_settings = sett.var.v_settings;

v         = spm_mb_io('GetData',datn.v);
q         = datn.q;
Mn        = datn.Mat;
Mr        = spm_dexpm(q,B);
Mat       = Mmu\Mr*Mn;
df        = spm_mb_io('GetSize',datn.f);
psi       = Compose(spm_mb_io('GetData',datn.psi),Affine(df,Mat));
mu        = Pull1(mu,psi);
[f,datn]  = spm_mb_io('GetClasses',datn,mu,sett);
[a,w]     = Push1(f - Softmax(mu,4),psi,d);
psi = []; f = []; mu = [];

g         = reshape(sum(a.*G,4),[d 3]);
H         = w.*H0;
a = []; w = [];

u0        = spm_diffeo('vel2mom', v, v_settings);                          % Initial momentum
datn.E(2) = 0.5*sum(u0(:).*v(:));                                          % Prior term
v         = v - scal*spm_diffeo('fmg',H, g + u0, [v_settings s_settings]); % Gauss-Newton update

if d(3)==1, v(:,:,:,3) = 0; end % If 2D
if v_settings(1)==0             % Mean displacement should be 0
    avg = mean(mean(mean(v,1),2),3);
    v   = v - avg;
end
datn.v = spm_mb_io('SetData',datn.v,v);
end
%==========================================================================

%==========================================================================
% UpdateWarpsSub()
function datn = UpdateWarpsSub(datn,avg_v,sett,kernel)
v        = spm_mb_io('GetData',datn.v);
if ~isempty(avg_v)
    v    = v - avg_v;
end
datn.v   = spm_mb_io('SetData',datn.v,v);
psi1     = Shoot(v, kernel, sett.shoot.args); % Geodesic shooting
datn.psi = spm_mb_io('SetData',datn.psi,psi1);
end
%==========================================================================

%==========================================================================
% VelocityHessian()
function H = VelocityHessian(mu,G,accel)
d  = [size(mu,1),size(mu,2),size(mu,3)];
M  = size(mu,4);
if accel>0, s  = spm_mb_shape('Softmax',mu,4); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
H  = zeros([d 6],'single');
for m1=1:M
    for m2=1:M
        if accel==0
            tmp = Ab(m1,m2);
        else
            if m2~=m1
                tmp = (-s(:,:,:,m1).*s(:,:,:,m2))*accel           + (1-accel)*Ab(m1,m2);
            else
                tmp = (max(s(:,:,:,m1).*(1-s(:,:,:,m1)),0))*accel + (1-accel)*Ab(m1,m2);
            end
        end
        H(:,:,:,1) = H(:,:,:,1) + tmp.*G(:,:,:,m1,1).*G(:,:,:,m2,1);
        H(:,:,:,2) = H(:,:,:,2) + tmp.*G(:,:,:,m1,2).*G(:,:,:,m2,2);
        H(:,:,:,3) = H(:,:,:,3) + tmp.*G(:,:,:,m1,3).*G(:,:,:,m2,3);
        H(:,:,:,4) = H(:,:,:,4) + tmp.*G(:,:,:,m1,1).*G(:,:,:,m2,2);
        H(:,:,:,5) = H(:,:,:,5) + tmp.*G(:,:,:,m1,1).*G(:,:,:,m2,3);
        H(:,:,:,6) = H(:,:,:,6) + tmp.*G(:,:,:,m1,2).*G(:,:,:,m2,3);
    end
end
end
%==========================================================================