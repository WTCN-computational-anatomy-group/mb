function varargout = spm_multireg_util(varargin)
%__________________________________________________________________________
%
% Utility functions for spm_multireg.
%
% FORMAT psi0      = spm_multireg_util('Affine',d,Mat)
% FORMAT psi       = spm_multireg_util('Compose',psi1,psi0)
% FORMAT id        = spm_multireg_util('Identity',d)
% FORMAT l         = spm_multireg_util('lse',mu,dr)
% FORMAT f         = spm_multireg_util('MaskF',f)
% FORMAT a1        = spm_multireg_util('Pull1',a0,psi,r)
% FORMAT [f1,w1]   = spm_multireg_util('Push1',f,psi,d,r)
% FORMAT             spm_multireg_util('SetBoundCond')
% FORMAT             spm_multireg_util('SetPath')
% FORMAT varargout = spm_multireg_util('Shoot',v0,kernel,args)
% FORMAT mu1       = spm_multireg_util('ShrinkTemplate',mu,oMmu,sett)
% FORMAT P         = spm_multireg_util('softmax',mu,dr)
% FORMAT P         = spm_multireg_util('softmaxmu',mu,dr)
% FORMAT [Mmu,d]   = spm_multireg_util('SpecifyMean',dat,vx)
% FORMAT varargout = spm_multireg_util('SubSample',samp,Mat,d0,varargin)
% FORMAT [dat,mu]  = spm_multireg_util('ZoomVolumes',dat,mu,sett,oMmu)
% FORMAT             spm_multireg_util('WriteNii',f,img,Mmu,descrip);
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_util
    error('Not enough argument. Type ''help spm_multireg_util'' for help.');
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
    case 'lse'
        [varargout{1:nargout}] = lse(varargin{:});
    case 'MaskF'
        [varargout{1:nargout}] = MaskF(varargin{:});    
    case 'Pull1'
        [varargout{1:nargout}] = Pull1(varargin{:});
    case 'Push1'
        [varargout{1:nargout}] = Push1(varargin{:});
    case 'SetBoundCond'
        [varargout{1:nargout}] = SetBoundCond(varargin{:});
    case 'SetPath'
        [varargout{1:nargout}] = SetPath(varargin{:});           
    case 'Shoot'
        [varargout{1:nargout}] = Shoot(varargin{:});
    case 'ShrinkTemplate'
        [varargout{1:nargout}] = ShrinkTemplate(varargin{:});
    case 'softmax'
        [varargout{1:nargout}] = softmax(varargin{:}); 
    case 'softmaxmu'
        [varargout{1:nargout}] = softmaxmu(varargin{:});         
    case 'SpecifyMean'
        [varargout{1:nargout}] = SpecifyMean(varargin{:});        
    case 'SubSample'
        [varargout{1:nargout}] = SubSample(varargin{:});               
    case 'ZoomVolumes'
        [varargout{1:nargout}] = ZoomVolumes(varargin{:});        
    case 'WriteNii'
        [varargout{1:nargout}] = WriteNii(varargin{:});            
    otherwise
        help spm_multireg_util
        error('Unknown function %s. Type ''help spm_multireg_util'' for help.', id)
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
SetBoundCond;
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
% lse()
function l = lse(mu,dr)
mx = max(mu,[],dr);
l  = log(exp(-mx) + sum(exp(mu - mx),dr)) + mx;
end
%==========================================================================

%==========================================================================
% MaskF()
function fn = MaskF(fn)
C = size(fn,2);
for c=1:C
    fn(:,c) = Mask(fn(:,c));
end
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

SetBoundCond;

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
    id = Identity(size(psi));
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

SetBoundCond;

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

    id    = Identity(size(psi));
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
% SetBoundCond()
function SetBoundCond
spm_diffeo('bound',0);
spm_field('bound',1);
end
%==========================================================================

%==========================================================================
% SetPath()
function SetPath
pth=fileparts(which('spm'));
addpath(pth);
addpath(fullfile(pth,'toolbox','Longitudinal'));
addpath(fullfile(pth,'toolbox','Shoot'));
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
    clear Jdp
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
% softmax()
function P = softmax(mu,dr)
mx  = max(mu,[],dr);
E   = exp(mu-mx);
den = sum(E,dr)+exp(-mx);
P   = E./den;
end
%==========================================================================

%==========================================================================
% softmaxmu()
function P = softmaxmu(mu,dr)
d  = size(mu);
d  = [d 1 1];
K  = d(4);
d  = d(1:3);

mx = max(mu,[],dr);
e  = exp(mu - mx);
on = exp(-mx);
se = sum(e,dr);
se = on + se;

P            = zeros([d K + 1],'like',mu);
P(:,:,:,1:K) = e./se;  % The first K classes softmax is: exp(a_i)/(1+sum exp(a_j))
P(:,:,:,end) = on./se; % And the last class is: 1/(1+sum exp(a_j))
end
%==========================================================================

%==========================================================================
% SpecifyMean()
function [Mmu,d] = SpecifyMean(dat,vx)
dims = zeros(numel(dat),3);
Mat0 = zeros(4,4,numel(dat));
for n=1:numel(dat)
    dims(n,:)   = spm_multireg_io('GetSize',dat(n).f)';    
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
% SubSample()
function varargout = SubSample(samp,Mat,d0,varargin)
vx        = sqrt(sum(Mat(1:3,1:3).^2));
samp      = max([1 1 1],round(samp*[1 1 1]./vx));
N         = numel(varargin);
varargout = cell(1,2*N + 1);
cnt       = 1;
for n=1:N
    f              = varargin{n};    
    varargout{cnt} = f; 
    cnt            = cnt + 1;
    
    if d0 > 0
        C = size(f,2);    
        f = reshape(f,[d0(1:3) C]);    
    end
    
    f = f(1:samp:end,1:samp:end,1:samp:end,:);
    d = size(f); 
    d = [d 1];
    
    if d0 > 0
        varargout{cnt} = reshape(f,[prod(d(1:3)) C]); 
    else
        varargout{cnt} = f; 
    end
    cnt = cnt + 1;
end

if d0 > 0
    % For weighting data parts of lowerbound with factor based on amount of
    % downsampling  
    W              = prod(d0(1:3))/prod(d(1:3));
    varargout{end} = W;
end
end
%==========================================================================

%==========================================================================
% ZoomVolumes()
function [dat,mu] = ZoomVolumes(dat,mu,sett,oMmu)

% Parse function settings
d       = sett.var.d;
Mmu     = sett.var.Mmu;
threads = sett.gen.threads;

SetBoundCond;

d0    = [size(mu,1) size(mu,2) size(mu,3)];
z     = single(reshape(d./d0,[1 1 1 3]));
Mzoom = oMmu\Mmu;
y     = reshape(reshape(Identity(d),[prod(d),3])*Mzoom(1:3,1:3)' + Mzoom(1:3,4)',[d 3]);
mu    = spm_diffeo('pullc',mu,y);

if threads>1 && numel(dat)>1
    % Should attempt to change number of threads accrding to how much memory
    % each one requires
    % Memory = 4*(prod(d)*(3+3)+prod(d0)*(3+3));
    parfor(n=1:numel(dat),threads) % PARFOR
        v          = spm_multireg_io('GetData',dat(n).v);
        v          = spm_diffeo('pullc',v.*z,y);
        dat(n).v   = ResizeFile(dat(n).v  ,d,Mmu);
        dat(n).psi = ResizeFile(dat(n).psi,d,Mmu);
        dat(n).v   = spm_multireg_io('SetData',dat(n).v,v);
    end
else
    for n=1:numel(dat)
        v          = spm_multireg_io('GetData',dat(n).v);
        v          = spm_diffeo('pullc',v,y).*z;
        dat(n).v   = ResizeFile(dat(n).v  ,d,Mmu);
        dat(n).psi = ResizeFile(dat(n).psi,d,Mmu);
        dat(n).v   = spm_multireg_io('SetData',dat(n).v,v);
    end
end
end
%==========================================================================

%==========================================================================
% WriteNii()
function WriteNii(f,img,M,descrip)
fa       = file_array(f,size(img),'float32',0);
Nii      = nifti;
Nii.dat  = fa;
Nii.mat  = M;
Nii.mat0 = M;
Nii.descrip = descrip;
create(Nii);
Nii.dat(:,:,:,:) = img;
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
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
B = spm_multireg_par('AffineBases','SE(3)');

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
% Mask()
function f = Mask(f)
f(~isfinite(f) | f == 0) = NaN;
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