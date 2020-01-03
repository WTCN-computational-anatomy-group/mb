function varargout = spm_mb_shape(varargin)
%__________________________________________________________________________
%
% Functions for shape model related.
%
%--------------------------------------------------------------------------
% INITIALISE MODEL
%
% FORMAT dat           = spm_mb_shape('InitDat',dat,sett)
% FORMAT model         = spm_mb_shape('InitModel',sett)
% FORMAT model         = spm_mb_shape('InitSubspace',model,sett)
% FORMAT [dat,mu]      = spm_mb_shape('InitTemplate',dat,K,sett)
%
%--------------------------------------------------------------------------
% UPDATE MODEL
%
% FORMAT dat           = spm_mb_shape('UpdateAffines',dat,mu,sett)
% FORMAT [model,dat]   = spm_mb_shape('UpdateMean',dat, model, sett)
% FORMAT dat           = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett)
% FORMAT [mu,dat]      = spm_mb_shape('UpdateSimpleMean',dat, mu, sett)
% FORMAT dat           = spm_mb_shape('UpdateVelocities',dat,model,sett)
% FORMAT dat           = spm_mb_shape('UpdateWarps',dat,sett)
% FORMAT [dat,model]   = spm_mb_shape('UpdateLatent',dat,model,sett)
% FORMAT model         = spm_mb_shape('UpdateSubspace',dat,model,sett)
% FORMAT model         = spm_mb_shape('UpdateLatentPrecision',model,sett)
% FORMAT model         = spm_mb_shape('UpdateResidualPrecision',model,sett)
% FORMAT [dat,model]   = spm_mb_shape('OrthoSubspace',dat,model,sett)
%
%--------------------------------------------------------------------------
% ENERGY
%
% FORMAT model         = spm_mb_shape('SuffStatVelocities',dat,model,sett)
% FORMAT model         = spm_mb_shape('SuffStatTemplate',dat,model,sett)
% FORMAT E             = spm_mb_shape('ShapeEnergy',model,sett)
% FORMAT [dat,model]   = spm_mb_shape('OrthoSubspace',dat,model,sett)
%
%--------------------------------------------------------------------------
% UTILS
%
% FORMAT psi0          = spm_mb_shape('Affine',d,Mat)
% FORMAT B             = spm_mb_shape('AffineBases',code)
% FORMAT psi           = spm_mb_shape('Compose',psi1,psi0)
% FORMAT id            = spm_mb_shape('Identity',d)
% FORMAT l             = spm_mb_shape('LSE',mu,dr)
% FORMAT a1            = spm_mb_shape('Pull1',a0,psi,r)
% FORMAT [f1,w1]       = spm_mb_shape('Push1',f,psi,d,r)
% FORMAT sd            = spm_mb_shape('SampDens',Mmu,Mn)
% FORMAT P             = spm_mb_shape('Softmax',mu,dr)
% FORMAT mun           = spm_mb_shape('TemplateK1',mun,ax)
% FORMAT mu1           = spm_mb_shape('Shrink',mu,oMmu,sett)
% FORMAT [Mmu,d]       = spm_mb_shape('SpecifyMean',dat,vx)
% FORMAT [dat,model]   = spm_mb_shape('ZoomVolumes',dat,model,sett,oMmu)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

%__________________________________________________________________________
% 
% STUFF TO DISCUSS
% ----------------
% . Function name: Mean or Template (need to be consistant)
% . In th code: math (mu, A, psi) or english (template, latent_prec, warp)
% . In the model structure: same question
% . Clean the public API: do all the utils need to be exposed?
% . It is a bit annoying that (private) subfunctions are very far from
%   (public) main functions in the file.
% . Function arguments: model structure or individual variables?
% . IO: could be a made a bit cleaner using file_arrays.
%__________________________________________________________________________

if nargin == 0
    help spm_mb_shape
    error('Not enough argument. Type ''help spm_mb_shape'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id    
    case 'Affine'
        [varargout{1:nargout}] = Affine(varargin{:});     
    case 'AffineBases'
        [varargout{1:nargout}] = AffineBases(varargin{:});           
    case 'Compose'
        [varargout{1:nargout}] = Compose(varargin{:});
    case 'Identity'
        [varargout{1:nargout}] = Identity(varargin{:});
    case 'InitDat'
        [varargout{1:nargout}] = InitDat(varargin{:});          
    case 'InitModel'
        [varargout{1:nargout}] = InitModel(varargin{:});         
    case 'InitSubspace'
        [varargout{1:nargout}] = InitSubspace(varargin{:});              
    case 'InitTemplate'
        [varargout{1:nargout}] = InitTemplate(varargin{:});          
    case 'LSE'
        [varargout{1:nargout}] = LSE(varargin{:});
    case 'OrthoSubspace'
        [varargout{1:nargout}] = OrthoSubspace(varargin{:});  
    case 'Pull1'
        [varargout{1:nargout}] = Pull1(varargin{:});
    case 'Push1'
        [varargout{1:nargout}] = Push1(varargin{:});   
    case 'ShapeEnergy'
        [varargout{1:nargout}] = ShapeEnergy(varargin{:});   
    case 'SampDens'
        [varargout{1:nargout}] = SampDens(varargin{:});
    case 'ShrinkTemplate'
        [varargout{1:nargout}] = Shrink(varargin{:});
    case 'Softmax'
        [varargout{1:nargout}] = Softmax(varargin{:});       
    case 'SpecifyMean'
        [varargout{1:nargout}] = SpecifyMean(varargin{:});   
    case 'SuffStatVelocities'
        [varargout{1:nargout}] = SuffStatVelocities(varargin{:});    
    case 'SuffStatTemplate'
        [varargout{1:nargout}] = SuffStatTemplate(varargin{:});              
    case 'TemplateEnergy'
        [varargout{1:nargout}] = TemplateEnergy(varargin{:});    
    case 'TemplateK1'
        [varargout{1:nargout}] = TemplateK1(varargin{:});          
    case 'UpdateAffines'
        [varargout{1:nargout}] = UpdateAffines(varargin{:});  
    case 'UpdateLatent'
        [varargout{1:nargout}] = UpdateLatent(varargin{:});
    case 'UpdateLatentPrecision'
        [varargout{1:nargout}] = UpdateLatentPrecision(varargin{:});                 
    case 'UpdateMean'
        [varargout{1:nargout}] = UpdateMean(varargin{:});
    case 'UpdateResidualPrecision'
        [varargout{1:nargout}] = UpdateResidualPrecision(varargin{:});
    case 'UpdateSimpleAffines'
        [varargout{1:nargout}] = UpdateSimpleAffines(varargin{:});
    case 'UpdateSimpleMean'
        [varargout{1:nargout}] = UpdateSimpleMean(varargin{:});
    case 'UpdateSubspace'
        [varargout{1:nargout}] = UpdateSubspace(varargin{:});
    case 'UpdateVelocities'
        [varargout{1:nargout}] = UpdateVelocities(varargin{:});
    case 'UpdateWarps'
        [varargout{1:nargout}] = UpdateWarps(varargin{:});      
    case 'ZoomVolumes'
        [varargout{1:nargout}] = ZoomVolumes(varargin{:});              
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
if d(3) == 1, psi0(:,:,:,3) = 1; end
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
function dat = InitDat(dat,sett)

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
            [pth,nam,~] = fileparts(dat(n).f(1).dat.fname);
            if isempty(dir_res), dir_res = pth; end
            vname       = fullfile(dir_res,['v_' nam '.nii']);
            pname       = fullfile(dir_res,['psi_' nam '.nii']);
            
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
    dat(n).ss.trLSv = 0;
end
end
%==========================================================================


%==========================================================================
% InitModel()
function [model,dat,sett] = InitModel(dat,sett)
% Initialise model structure:
%   . if sett.pca.do: A, nA, lam, nlam, Z, ZZ, Sz, dat.z
%   . ss.trLVV, ss.trLSV


model = struct;

if sett.do.pca
    npc = sett.pca.npc;
    
    % ---------------------------------------------------------------------
    % Latent precision
    if isscalar(sett.pca.latent_prior)
        if ~isfinite(sett.pca.latent_prior)
            sett.pca.latent_prior = 1;
        end
        sett.pca.latent_prior = sett.pca.latent_prior * eye(npc);
    elseif size(sett.pca.latent_prior,1) ~= npc
        npc0 = size(sett.pca.latent_prior,1);
        warning(['Latent prior and number of principal component not ' ...
                 'consistant: %d vs %d'], npc0, npc);
        A0 = sett.pca.latent_prior;
        sett.pca.latent_prior = eye(sett.pca.npc);
        npc0 = min(npc,npc0);
        sett.pca.npc(1:npc0,1:npc0) = A0(1:npc0,1:npc0);
        sett.pca.latent_prior = A0;
    end
    model.A = sett.pca.latent_prior;
    if isfinite(sett.pca.latent_df) && sett.pca.latent_df > 0
        sett.pca.latent_df = max(sett.pca.latent_df, ...
                                 sett.pca.npc - 1 + 1E-3);
        model.nA = sett.pca.latent_df;
    else
        model.nA = Inf;
    end
    
    % ---------------------------------------------------------------------
    % Residual precision
    if ~isfinite(sett.pca.res_prior)
        sett.pca.res_prior = 10;
    end
    model.lam = sett.pca.res_prior;
    if isfinite(sett.pca.res_df) && sett.pca.res_df > 0
        model.nlam = sett.pca.res_df;
    else
        model.nlam = Inf;
    end
    
    % ---------------------------------------------------------------------
    % Latent variable (random but orthogonal initialisation)
    model.Z  = zeros(npc,numel(dat));
    model.Sz = eye(npc)/numel(dat);
    model.ZZ = 0;
    for n=1:numel(dat)
        z            = randn(npc,1);
        dat(n).z     = z;
        model.Z(:,n) = z;
        model.ZZ     = model.ZZ + (z*z');
    end
    [U,S] = svd(model.ZZ);
    Q = S^(-0.5)*U';
    model.ZZ = Q * (model.ZZ / Q);
    model.ZZ = (model.ZZ + model.ZZ')/2;
    model.ZZ = model.ZZ + eye(npc);
    model.Z  = Q * model.Z;
    for n=1:numel(dat)
        dat(n).z     = Q * dat(n).z;
    end
end
    
% -------------------------------------------------------------------------
% Sufficient statistics
model.ss.trLVV = 0;
model.ss.trLSV = 0;
end


%==========================================================================
% InitSubspace()
function model = InitSubspace(model,sett)
% Initialise subspace and relatd sufficient statistics

if ~sett.do.pca, return; end

updt_subspace = sett.do.updt_subspace;
dir_res       = sett.write.dir_res;
max_mem       = sett.gen.max_mem;
lat           = [sett.var.d 1];
lat           = lat(1:3);
npc           = sett.pca.npc;
Mmu           = sett.var.Mmu;
v_setting     = sett.var.v_settings;

if ~isfield(model, 'U')
    if prod(lat)*3*npc*4 > max_mem
        % Create file
        fname       = fullfile(dir_res,'subspace.nii');
        fa          = file_array(fname,[lat 1 3 npc],'float32');
        nii         = nifti;
        nii.dat     = fa;
        nii.mat     = Mmu;
        nii.mat0    = Mmu;
        nii.descrip = 'Subspace';
        create(nii);
        model.U        = nii.dat;
        model.U(end)   = 0;
        model.U.dim(4) = [];
    else
        % Keep in memory
        model.U = zeros([lat 3 npc], 'single');
    end
    model.Su  = eye(npc);
    model.ULU = model.Su * 3*prod(lat);
else
    D    = [size(model.U) 1 1 1];
    D    = prod(D(1:4));             % Number of voxels * 3
    L    = size(model.U,5);          % nb PCs

    % posterior precision
    if updt_subspace
        model.Su = eye(npc);
    else
        model.Su = zeros(npc);
    end

    % 2nd order moment
    ULU = zeros(npc);
    for l=1:L
        U1  = single(model.U(:,:,:,:,l));
        LU1 = spm_diffeo('vel2mom', U1, v_setting);
        ULU(l,l) = sum(LU1(:) .* U1(:), 'double');
        for k=(l+1):L
            U2       = single(model.U(:,:,:,:,k));
            ULU(l,k) = sum(LU1(:) .* U2(:), 'double');
            ULU(k,l) = ULU(l,k);
        end
    end
    clear U1 LU1 U2
    model.ULU = ULU + D * model.Su;  
end
end
%==========================================================================


%==========================================================================
% InitTemplate()
function [dat,mu] = InitTemplate(dat,K,sett)
% Make 'quick' initial estimates of GMM posteriors and template on very coarse
% scale

% Parse function settings
do_gmm = sett.do.gmm;

% Uniform template
mu = zeros([sett.var.d K],'single');

if ~do_gmm, return; end

% % Change some settings
% do_updt_bf      = sett.do.updt_bf;
ix_init         = sett.model.ix_init_pop;
mg_ix           = sett.model.mg_ix;
nit_init_mu     = sett.nit.init_mu;
sett.do.updt_bf = false;
sett.gen.samp   = 5;
sett.nit.appear = 1;
sett.nit.gmm    = 100;

% Parameters
K1  = K + 1;
Kmg = numel(mg_ix);

% Get population indices
p_ix       = spm_mb_appearance('GetPopulationIdx',dat);
Npop       = numel(p_ix);
first_subj = true; % For when not using InitGMM, make posterior means of all subjects uninformative, except first subject in population sett.model.ix_init_pop
pop_rng    = 1:Npop;

% Was InitGMM run on initialising population?
[~,C]       = spm_mb_io('GetSize',dat(p_ix{ix_init}(1)).f);
has_ct      = any(dat(p_ix{ix_init}(1)).is_ct == true);
use_initgmm = C == 1 && ~has_ct;

for p=[ix_init pop_rng(pop_rng ~= ix_init)] % loop over populations (starting index defined by sett.model.ix_init_pop)
    % To make the algorithm more robust when using multiple populations,
    % set posterior and prior means (m) of GMMs of all but the first population to
    % uniform  
    
    if use_initgmm && p == ix_init, continue; end
    
    C = size(dat(p_ix{p}(1)).mog.po.m,1);
   
    avg_m_po  = 0;
    avg_m_pr  = 0;
    avg_vr_po = 0;
    avg_vr_pr = 0;
    for n=p_ix{p}
        % mean
        avg_m_po = avg_m_po + dat(n).mog.po.m;
        avg_m_pr = avg_m_pr + dat(n).mog.pr.m;
        
        % variance
        vr_po = reshape(dat(n).mog.po.n,[1 1 Kmg]).*dat(n).mog.po.W;
        vr_pr = reshape(dat(n).mog.pr.n,[1 1 Kmg]).*dat(n).mog.pr.W;
        for k=1:Kmg
            vr_po(:,:,k) = inv(vr_po(:,:,k));
            vr_pr(:,:,k) = inv(vr_pr(:,:,k));
        end
        avg_vr_po = avg_vr_po + vr_po;
        avg_vr_pr = avg_vr_pr + vr_pr;
    end
    
    % Average means
    avg_m_po = avg_m_po./numel(p_ix{p});
    avg_m_po = mean(avg_m_po,2);
    avg_m_pr = avg_m_pr./numel(p_ix{p});
    avg_m_pr = mean(avg_m_pr,2);
        
    % Average variances
    avg_vr_po = avg_vr_po./numel(p_ix{p});
    avg_vr_po = mean(avg_vr_po,3);
    avg_vr_pr = avg_vr_pr./numel(p_ix{p});
    avg_vr_pr = mean(avg_vr_pr,3);
        
    % Add a bit of random noise to prior    
    mpo = zeros(C,Kmg);
    mpr = zeros(C,Kmg);
    for k=1:K1
        kk = sum(mg_ix == k);
        w  = 1./(1 + exp(-(kk - 1)*0.25)) - 0.5;
                
        rng(1);
        mn                = avg_m_po;
        vr                = avg_vr_po;                
        mpo(:,mg_ix == k) = sqrtm(vr)*sort(randn(C,kk),2)*w + repmat(mn,[1 kk]);
        
        rng(1);
        mn                = avg_m_pr;
        vr                = avg_vr_pr;                
        mpr(:,mg_ix == k) = sqrtm(vr)*sort(randn(C,kk),2)*w + repmat(mn,[1 kk]);
    end
    
    % Assign
    for n=p_ix{p}
        dat(n).mog.pr.m = mpr; % prior
        if ~use_initgmm && first_subj            
            first_subj = false;
            continue
        end        
        dat(n).mog.po.m = mpo; % posterior
    end
    
    if 0
        spm_gmm_lib('plot','gaussprior',{mpo,dat(n).mog.po.b,dat(n).mog.po.W,dat(n).mog.po.n},[],'InitMu');
        spm_gmm_lib('plot','gaussprior',{mpr,dat(n).mog.pr.b,dat(n).mog.pr.W,dat(n).mog.pr.n},[],'InitMu');
    end
end

if ~use_initgmm
    % Update template based on only first subject of population sett.model.ix_init_pop
    for it=1:nit_init_mu
        [mu,dat(p_ix{ix_init}(1))] = spm_mb_shape('UpdateSimpleMean',dat(p_ix{ix_init}(1)), mu, sett);
    end
end
% Update template using all subjects from population sett.model.ix_init_pop
for it=1:nit_init_mu
    [mu,dat(p_ix{ix_init})] = spm_mb_shape('UpdateSimpleMean',dat(p_ix{ix_init}), mu, sett);
end
if Npop > 1
    % If more than one population, use template learned on sett.model.ix_init_pop
    % population to initialise other populations' GMM parameters
    [mu,dat] = spm_mb_shape('UpdateSimpleMean',dat,    mu, sett);    
end
end
%==========================================================================

%==========================================================================
% LSE()
function l = LSE(mu,dr)
% Log-Sum-Exp
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
% SampDens()
function sd = SampDens(Mmu,Mn)
vx_mu = sqrt(sum(Mmu(1:3,1:3).^2,1));
vx_f  = sqrt(sum( Mn(1:3,1:3).^2,1));
sd    = max(round(2.0*vx_f./vx_mu),1);
end
%==========================================================================

%==========================================================================
% Shrink()
function mu1 = Shrink(mu,oMmu,sett)

% Parse function settings
d       = sett.var.d;
Mmu     = sett.var.Mmu;
max_mem = sett.gen.max_mem*1E9; % (GB)
dir_res = sett.write.dir_res;

d0      = [size(mu,1) size(mu,2) size(mu,3)];
Mzoom   = Mmu\oMmu;
if norm(Mzoom-eye(4))<1e-4 && all(d0==d)
    mu1 = mu;
else
    y = reshape(reshape(Identity(d0),[prod(d0),3])*Mzoom(1:3,1:3)'+Mzoom(1:3,4)',[d0 3]);
    if isa(mu, 'file_array') ...
            && prod(d)*size(mu,4)*size(mu,5)*4 > max_mem
        mu1 = nifti(mu.fname);
        [~,fname,ext] = mu1.dat.fname;
        mu1.dat.fname = fullfile(dir_res, [fname '_shrink' ext]);
        mu1.dat.dim(1:3) = d;
        mu1.mat = Mmu;
        create(mu1);
        mu1 = mu1.dat;
        [mu11,c] = Push1(mu(:,:,:,1,1),y,d);
        mu11     = mu11./(c+eps);
        mu1(:,:,:,1,1) = mu11;
        for k=1:size(mu,4)
            for l=1:size(mu,5)
                mu11 = Push1(mu(:,:,:,1,1),y,d);
                mu1(:,:,:,1,1) = mu11./(c+eps);
            end
        end
    else
        [mu1,c] = Push1(mu,y,d);
        mu1     = mu1./(c+eps);
    end
end
end
%==========================================================================

%==========================================================================
% Softmax()
function P = Softmax(mu,dr)
% Safe softmax
mx  = max(mu,[],dr);
E   = exp(mu-mx);
den = sum(E,dr)+exp(-mx);
P   = E./den;
end
%==========================================================================

%==========================================================================
% SpecifyMean()
function [Mmu,d] = SpecifyMean(dat,vx)
% Compute orientation and dimensions of the mean/template space from the
% individual observed volumes.
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
% TemplateK1()
function mun = TemplateK1(mun,ax)
% Add (K+1)-th class to template. This is done in a safe way using the
% log-sum-exp trick.
mx  = max(max(mun,[],ax),0);
lse = mx + log(sum(exp(mun - mx),ax) + exp(-mx));
mun = cat(ax,mun - lse, -lse);
end
%==========================================================================

%==========================================================================
% UpdateAffines()
function dat = UpdateAffines(dat,model,sett)

% Parse function settings
B                = sett.registr.B;
do_updt_aff      = sett.do.updt_aff;
do_updt_template = sett.do.updt_template;

if ~do_updt_aff, return; end

% Update the affine parameters
if ~isempty(B)
    for n=1:numel(dat)
        dat(n) = UpdateAffinesSub(dat(n),model.mu,sett);
    end

    if do_updt_template
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
function [model,dat] = UpdateMean(dat, model, sett)

% Parse function settings
accel            = sett.gen.accel;
do_updt_template = sett.do.updt_template;
mu_settings      = sett.var.mu_settings;
s_settings       = sett.shoot.s_settings;

if ~do_updt_template, return; end

g  = spm_field('vel2mom', model.mu, mu_settings);
M  = size(model.mu,4);
H  = zeros([sett.var.d M*(M+1)/2],'single');
H0 = AppearanceHessian(model.mu,accel);
for n=1:numel(dat)
    [gn,Hn,dat(n)] = UpdateMeanSub(dat(n),model.mu,H0,sett);
    g              = g + gn;
    H              = H + Hn;
end
clear H0 gn Hn
model.mu = model.mu - spm_field(H, g, [mu_settings s_settings]);  
end
%==========================================================================

%==========================================================================
% UpdateSimpleAffines()
function dat = UpdateSimpleAffines(dat,mu,sett)

% Parse function settings
accel            = sett.gen.accel;
B                = sett.registr.B;
do_updt_aff      = sett.do.updt_aff;
do_updt_template = sett.do.updt_template;

if ~do_updt_aff, return; end

% Update the affine parameters
G  = spm_diffeo('grad',mu);
H0 = VelocityHessian(mu,G,accel);

if ~isempty(B)
    for n=1:numel(dat)
        dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
    end

    if do_updt_template
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
function dat = UpdateVelocities(dat,model,sett)

% Parse function settings
accel       = sett.gen.accel;
do_updt_vel = sett.do.updt_vel;

if ~do_updt_vel, return; end

G  = spm_diffeo('grad',model.mu);
H0 = VelocityHessian(model.mu,G,accel);
if size(G,3) == 1
    % Data is 2D -> add some regularisation
    H0(:,:,:,3) = H0(:,:,:,3) + mean(reshape(H0(:,:,:,[1 2]),[],1));
end
for n=1:numel(dat)
    dat(n) = UpdateVelocitiesSub(dat(n),model,G,H0,sett);
end
end
%==========================================================================

%==========================================================================
% UpdateLatent()
function [dat,model] = UpdateLatent(dat,model,sett)
% Posterior update of a multivariate Gaussian distribution.

if ~sett.do.pca, return; end

% Dependencies
N   = numel(dat);         % Number of subjects
L   = size(model.U,5);    % Number of principal components
A   = model.A;            % Precision matrix (posterior expected value)
lam = model.lam;          % Residual precision (posterior expected value)
ULU = model.ULU;          % Subspace 2nd order moment: E[U'LU]

% Posterior variance (common to all subjects)
model.Sz = inv(A + lam * ULU);

% Posterior mean
model.Z  = zeros(L,N);  % Mean for each subject
model.ZZ = zeros(L);    % 2nd order moment: E[Z*Z']
for n=1:N
    dat(n)       = UpdateLatentSub(dat(n), model, sett);
    model.Z(:,n) = dat(n).z;
    model.ZZ     = model.ZZ + dat(n).z*dat(n).z';
end
model.ZZ = model.ZZ + N*model.Sz;

end
%==========================================================================

%==========================================================================
% UpdateLatentSub()
function datn = UpdateLatentSub(datn,model,sett)

v_settings = sett.var.v_settings;

L   = size(model.U,5);    % Number of principal components
A   = model.A;            % Precision matrix (posterior expected value)
lam = model.lam;          % Residual precision (posterior expected value)

m = spm_mb_io('GetData',datn.v);
m = spm_diffeo('vel2mom',m,v_settings); % Compute momentum: L * v
m = m(:);

datn.z = zeros(L,1);
for j=1:L
    U1        = single(model.U(:,:,:,:,j)); % Extract j-th mode
    datn.z(j) = sum(U1(:).*m,'double');     % Project on subspace
end
clear U1

datn.z  = model.lam * model.Sz * datn.z;

end
%==========================================================================

%==========================================================================
% UpdateLatentPrecision()
function model = UpdateLatentPrecision(model,sett)
% Posterior update of a Wishart distribution.

% n0 == 0:   Maximum-likelihood
% n0 >  0:   Posterior
% n0 == Inf: Fixed value

if ~sett.do.pca, return; end
if ~sett.do.updt_latent_prior, return; end
if ~isfinite(sett.pca.latent_df), return; end % Fixed value

ndat = size(model.Z,2);         % Numbe of subjects (ss0)
ZZ   = model.ZZ;                % 2nd order moment  (ss2)
A0   = sett.pca.latent_prior;   % Prior expected value
n0   = sett.pca.latent_df;      % Prior deg. freedom

n = n0 + ndat;
if n0 > 0,  A = (n*A0)/(A0*ZZ +  n0*eye(size(ZZ)));
else,       A = n*inv(ZZ + 1E-3*eye(size(ZZ))); % slight regularisation
end

model.A = A;
if n0 > 0, model.nA = n;
else,      model.nA = Inf;
end

end
%==========================================================================

%==========================================================================
% UpdateResidualPrecision()
function model = UpdateResidualPrecision(model,sett)
% Posterior update of a Gamma distribution.

% n0 == 0:   Maximum-likelihood
% n0 >  0:   Posterior
% n0 == Inf: Fixed value

if ~sett.do.pca, return; end
if ~sett.do.updt_res_prior, return; end
if ~isfinite(sett.pca.res_df), return; end % Fixed value

ndat  = size(model.Z,2);         % Numbe of subjects (ss0)
lam0  = sett.pca.res_prior;      % Prior expected value
n0    = sett.pca.res_df;         % Prior deg. freedom
trLVV = model.ss.trLVV;          % suff stat: tr(L*(V-WZ)*(V-WZ)')
trLSV = model.ss.trLSV;          % suff stat: tr(L*Cov[V])
Su    = model.Su;                % subspace posterior covariance (one vox)
Sz    = model.Sz;                % latent posterior covariance (one subj)
D     = [size(model.U) 1 1 1];
D     = prod(D(1:4));            % Number of voxels * 3

n   = n0 + ndat;
ss2 = trLVV/D + ndat*trace(Su*Sz) + trLSV/D;
if n0 > 0,  lam = (n*lam0)/(lam0*ss2 +  n0);
else,       lam = n/(ss2 + eps); % slight regularisation
end

model.lam = lam;
if n0 > 0, model.nlam = n;
else,      model.nlam = Inf;
end

end
%==========================================================================

%==========================================================================
% UpdateSubspace()
function model = UpdateSubspace(dat,model,sett)
% Posterior update of a Gaussian Matrix distribution.

% n0 == 0:   Maximum-likelihood
% n0 >  0:   Posterior
% n0 == Inf: Fixed value

if ~sett.do.pca, return; end
if ~sett.do.updt_subspace, return; end

v_setting = sett.var.v_settings;
Z    = model.Z;                  % expected latent   (per subject)
ZZ   = model.ZZ;                 % 2nd order moment  (ss2)
D    = [size(model.U) 1 1 1];
D    = prod(D(1:4));             % Number of voxels * 3
lam  = model.lam;                % Residual precision
L    = size(model.U,5);          % nb PCs

% posterior precision
model.Su  = inv(eye(size(ZZ)) + lam * ZZ);

% posterior mean
ZS = lam * Z' * model.Su;
for l=1:L
    U1 = 0;
    for n=1:numel(dat)
        v  = spm_mb_io('GetData',dat(n).v);
        U1 = U1 + v * ZS(n,l);
    end
    model.U(:,:,:,:,l) = U1;
end
clear v U1

% 2nd order moment
ULU = zeros(size(ZZ));
for l=1:L
    U1  = single(model.U(:,:,:,:,l));
    LU1 = spm_diffeo('vel2mom', U1, v_setting);
    ULU(l,l) = sum(LU1(:) .* U1(:), 'double');
    for k=(l+1):L
        U2       = single(model.U(:,:,:,:,k));
        ULU(l,k) = sum(LU1(:) .* U2(:), 'double');
        ULU(k,l) = ULU(l,k);
    end
end
clear U1 LU1 U2
model.ULU = ULU + D * model.Su;

end
%==========================================================================

%==========================================================================
% UpdateWarps()
function dat = UpdateWarps(dat,sett)

% Parse function settings
do_updt_template = sett.do.updt_template;
v_settings       = sett.var.v_settings;

if do_updt_template
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
% OrthoSubspace()
function [dat,model] = OrthoSubspace(dat,model,sett)
% Posterior update of a Gaussian Matrix distribution.

% n0 == 0:   Maximum-likelihood
% n0 >  0:   Posterior
% n0 == Inf: Fixed value

if ~sett.do.pca, return; end
if ~sett.do.updt_subspace, return; end    % Need subspace not fixed
if sett.pca.latent_df == Inf, return; end % Need latent precision not fixed

lat = [size(model.U) 1 1 1];
lat = lat(1:4);
L   = size(model.U,5);         % nb PCs
N   = numel(dat);              % Number of subjects
ULU = model.ULU;               % subspace 2nd order moment 
ZZ  = model.ZZ;                % latent 2nd order moment
A0  = sett.pca.latent_prior;   % Prior expected value
n0  = sett.pca.latent_df;      % Prior deg. freedom


% -------------------------------------------------------------------------
% Orthogonalise (= joint diagonalisation of ULU and ZZ)
[T, iT] = JointDiag(ZZ,ULU);

ULU = iT' * ULU * iT;  % < now a diagonal matrix
ZZ  = T   * ZZ  * T';  % < now a diagonal matrix

% -------------------------------------------------------------------------
% Rescale (= optimise additional diagoal loading)
[Q, iQ] = ScaleSubspace(ULU, ZZ, A0, n0, N, prod(lat));
Q  = Q  * T;
iQ = iT * iQ;

% -------------------------------------------------------------------------
% Rotate subspace (slice-wise for memory)
slice_lat    = lat;
slice_lat(3) = 1;
for z=1:lat(3)
    Uz = single(model.U(:,:,z,:,:));
    Uz = reshape(Uz, [], L) * iQ;
    Uz = reshape(Uz, [slice_lat L]);
    model.U(:,:,z,:,:) = Uz;
end
clear Uz

% -------------------------------------------------------------------------
% Rotate sufficient statistics
model.ULU = iQ' * model.ULU * iQ;
model.Su  = iQ' * model.Su * iQ;
model.ZZ  = Q   * model.ZZ  * Q';
model.Sz  = Q   * model.Sz  * Q';
model.Z   = Q   * model.Z;

% -------------------------------------------------------------------------
% Rotate subjects
for n=1:N
    dat(n).z = Q * dat(n).z;
end

% -------------------------------------------------------------------------
% Update latent prior
model = UpdateLatentPrecision(model,sett);

end
%==========================================================================

%==========================================================================
% ZoomVolumes()
function [dat,model] = ZoomVolumes(dat,model,sett,oMmu)

% Parse function settings
d       = sett.var.d;
Mmu     = sett.var.Mmu;
do_pca  = sett.do.pca;
max_mem = sett.gen.max_mem*1E9; % (GB)
dir_res = sett.write.dir_res;

d0       = [size(model.mu,1) size(model.mu,2) size(model.mu,3)];
z        = single(reshape(d./d0,[1 1 1 3]));
Mzoom    = oMmu\Mmu;
y        = reshape(reshape(Identity(d),[prod(d),3])*Mzoom(1:3,1:3)' + Mzoom(1:3,4)',[d 3]);
% --- Template ------------------------------------------------------------
model.mu = spm_diffeo('pullc',model.mu,y);
% --- Subspace ------------------------------------------------------------
if do_pca
    if isa(model.U, 'file_array')
        [path,fname,ext] = fileparts(model.U.fname);
        tmpname = fullfile(path, [fname '_tmp' ext]);
        ok = copyfile(model.U.fname, tmpname);
        if ~ok, error('Failed to copy the subspace file'); end
        nii = nifti(tmpname);
        nii.dat.fname = model.U.fname;
        nii.dat.dim(1:3) = d;
        nii.mat = Mmu;
        create(nii);
        if nii.dat.dim(4) == 1, nii.dat.dim(4) = []; end
        U = nii.dat;
    elseif prod(d)*size(model.U,4)*size(model.U,5)*4 > max_mem
        fname = fullfile(dir_res, 'subspace.nii');
        fa = file_array(fname, [d(1:3) 1 size(model.U,4) size(model.U,5)], 'float32');
        nii = nifti;
        nii.dat = fa;
        nii.mat = Mmu;
        nii.descrip = 'Subspace';
        create(nii);
        if nii.dat.dim(4) == 1, nii.dat.dim(4) = []; end
        U = nii.dat;
    else
        U = zeros([d(1:3) size(model.U,4) size(model.U,5)], 'single');
    end
    for k=1:size(model.U,4)
        for l=1:size(model.U,5)
            U(:,:,:,k,l) = spm_diffeo('pullc',single(model.U(:,:,:,k,l)),y);
        end
    end
    if isa(model.U, 'file_array')
        delete(tmpname);
    end
    model.U = U;
end
% --- Subject data --------------------------------------------------------
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
function model = SuffStatVelocities(dat,model,sett)
% Update velocity-related sufficient statistics for the model energy.
%   . ss.trLVV - trace(L*(V-WZ)*(V-WZ)')
%   . ss.trLSV - trace(L*Var[V])

% Parse function settings
v_settings = sett.var.v_settings;
do_pca     = sett.do.pca;

model.ss.trLVV = 0;
model.ss.trLSV = 0;
for n=1:numel(dat)
    v = spm_mb_io('GetData',dat(n).v);
    if do_pca
        v0 = 0;
        for l=1:size(model.U,5)
            v0 = v0 + spm_mb_io('GetData',model.U,5,l) * dat(n).z(l);
        end
        v = v - v0;
    end
    u0  = spm_diffeo('vel2mom', v, v_settings); % Initial momentum
    model.ss.trLVV = model.ss.trLVV + sum(u0(:).*v(:), 'double');
    model.ss.trLSV = model.ss.trLSV + dat(n).ss.trLSv;
end
end
%==========================================================================

%==========================================================================
% TemplateEnergy()
function model = SuffStatTemplate(model,sett)

% Parse function settings
do_updt_template = sett.do.updt_template;
mu_settings      = sett.var.mu_settings;

if do_updt_template
    g = spm_field('vel2mom', model.mu, mu_settings);
    model.ss.mLm = sum(model.mu(:).*g(:));
else
    model.ss.mLm = 0;
end
end
%==========================================================================

%==========================================================================
% ShapeEnergy()
function E = ShapeEnergy(model,sett)
% Compute energy of the PCA shape model from sufficient statistics.
%
%   . All sufficient statistics must have been updated prior to caling this
%     function.
%   . By energy we mean 'free energy' or 'negative evidence lower bound', 
%     in the variational Bayes context.
%   . When the PCA model is used (sett.do.pca == true), the velocity energy
%     is part of it. Else, it is stored in each subject (dat.E(3)).

E = [EnergyTemplate(model,sett) ...
     EnergyVelocity(model,sett) ...
     EnergyLatent(model,sett) ...
     EnergyLatentPrecision(model,sett) ...
     EnergyResidualPrecision(model,sett) ...
     EnergySubspace(model,sett)];

end
%==========================================================================


%==========================================================================
function E = EnergyTemplate(model, sett)
% Negative log-likelihood of a multivariate Gaussian distribution

E = 0.5*model.ss.mLm;

% TODO: we could learn the template regularisation
% (like we learn the residual precision in the PCA model)
% KL divergence between two multivariate Gaussian distributions.
% Log-determinants of the prior and posterior precision cannot be computed
% exactly so are dropped from the energy.

end
%==========================================================================


%==========================================================================
function E = EnergyVelocity(model, sett)
% KL divergence between two multivariate Gaussian distributions.
% Log-determinants of the prior and posterior precision cannot be computed
% exactly so are dropped from the energy.

if ~sett.do.pca
    % Mode estimate -> negative log-likelihood
    E = model.ss.trLVV;       % velocity prior (sum subj)
else
    % Variational posterior -> KL
    N      = size(model.Z,2);
    D      = size(model.U,1)*size(model.U,2)*size(model.U,3);
    F      = size(model.U,4);
    Su     = model.Su;        % subspace posterior covariance wrt PCs
    Sz     = model.Sz;        % latent posterior covariance (per subj)
    trLVV  = model.ss.trLVV;  % expected velocity prior (sum subj)
    trLSV  = model.ss.trLSV;  % velocity posterior covariance: Tr(L*Var[Z]) (sum subj)
    lam    = model.lam;       % residual noise precision (mean)
    nlam   = model.nlam;      % residual noise precision (deg. freedom)
    loglam = ELogGamma(lam,nlam);
    
    E = N*D*F*(loglam - 1) + lam*(trLVV + N*D*F*trace(Su*Sz) + trLSV);

end
E = 0.5*E;

end
%==========================================================================

%==========================================================================
function E = EnergyLatent(model, sett)
% KL divergence between two multivariate Gaussian distributions (xN subj)

if ~sett.do.pca, E = 0; return; end
    
N  = size(model.Z,2);  % Number of subjects
L  = size(model.U,5);  % Number of principal components
ZZ = model.ZZ;         % latent 2nd order moment     (sum subj)
Sz = model.Sz;         % latent posterior covariance (per subj)
A  = model.A;          % residual precision matrix (mean)
nA = model.nA;         % residual precision matrix (deg. freedom)
    
E =     trace(A*ZZ) ...
  - N * LogDetChol(Sz) ...
  - N * ELogDetWishart(A,nA) ...
  - N * L;
E = E/2;
end
%==========================================================================

%==========================================================================
function E = EnergySubspace(model, sett)
% KL divergence between two multivariate Gaussian distributions (xDF vox)
% The prior is set to be a standard Normal (zero-mean, identity covariance)

if ~sett.do.pca, E = 0; return; end

D   = size(model.U,1)*size(model.U,2)*size(model.U,3); % Number of voxels
F   = size(model.U,4);   % Should be 3
L   = size(model.U,5);   % Number of principal components
ULU = model.ULU;         % subspace 2nd order moment
Su  = model.Su;          % subspace posterior covariance wrt PCs

E =         trace(ULU) ...
    - D*F * LogDetChol(Su) ...
    - D*F * L;
E = E/2;
end
%==========================================================================

%==========================================================================
function E = EnergyLatentPrecision(model, sett)
% KL divergence between two Wishart distributions.

if ~sett.do.pca, E = 0; return; end

L = size(model.U,5);

n0 = sett.pca.latent_df;
A0 = sett.pca.latent_prior;
n  = model.nA;
A  = model.A;

if ~isfinite(n0) || n0 == 0  % Fixed value | Maximum likelihood
    E = 0;
else
    E =   n0 * (LogDetChol(A0) - LogDetChol(A))...
        + n0 * L * (log(n) - log(n0)) ...
        + n0 * trace(A0\A) ...
        - n * L ...
        + 2 * (LogGamma(n0/2,L) - LogGamma(n/2,L)) ...
        + (n-n0) * DiGamma(n/2,L);
    E = E/2;
end
end
%==========================================================================

%==========================================================================
function E = EnergyResidualPrecision(model, sett)
% KL divergence between two Gamma distributions.
%
%   The Gamma distributions are parameterised by lam and n, where
%   alpha = n*D/2, beta = n*D/(2*lam) and D is the number of elements in a
%   velocity field (D=nvox*3). This parameterisation allows lam to act as 
%   an expected value and n to act as a degrees of freedom.

if ~sett.do.pca, E = 0; return; end

D = size(model.U,1)*size(model.U,2)*size(model.U,3);
F = size(model.U,4);

n0   = sett.pca.res_df;
lam0 = sett.pca.res_prior;
n    = model.nlam;
lam  = model.lam;

if ~isfinite(n0) || n0 == 0  % Fixed value | Maximum likelihood
    E = 0;
else
    E =   D*F * n0 * log((lam0/n0)/(lam/n)) ...
        + D*F * n * ((lam/n)/(lam0/n0) - 1) ...
        +   2 * (LogGamma(n0*D*F/2) - LogGamma(n*D*F/2)) ...
        + D*F * (n-n0) * DiGamma(n*D*F/2);
    E = E/2;
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
if dims(1,3) == 1, B = AffineBases('SE(2)');
else,              B = AffineBases('SE(3)');
end

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
prct  = 0.05;            % percentage to remove (in each direction)
o     = -prct*(mx - mn); % offset -> make template a bit smaller (for using less memory!)
% o     = ones(3,1);
if dims(1,3) == 1, o(3) = 0; end
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
    clear Jdp
    
    u = spm_diffeo('pushc',u1,id+v/T);

    % v_t \gets L^g u_t
    v = spm_shoot_greens(u,kernel.F,kernel.v_settings); % Convolve with Greens function of L    
    
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
clear psi0  psi1

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
    clear Gm
end
clear J mu

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

[a,w] = Push1(f - Softmax(mu1,4),psi,d);
clear mu1 psi f

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
function datn = UpdateVelocitiesSub(datn,model,G,H0,sett)

% Parse function settings
B          = sett.registr.B;
d          = sett.var.d;
Mmu        = sett.var.Mmu;
s_settings = sett.shoot.s_settings;
scal       = sett.optim.scal_v;
v_settings = sett.var.v_settings;
do_pca     = sett.do.pca;

v  = spm_mb_io('GetData',datn.v);
v0 = 0;
if do_pca
    % Compute subject-specific mean velocity from PCA
    for l=1:size(model.U,5)
        v0 = v0 + spm_mb_io('GetData',model.U,5,l) * datn.z(l);
    end
    % Modulate precision
    v_settings(4:end) = v_settings(4:end) * model.lam;
    v = v - v0;
end
q         = datn.q;
Mn        = datn.Mat;
Mr        = spm_dexpm(q,B);
Mat       = Mmu\Mr*Mn;
df        = spm_mb_io('GetSize',datn.f);
psi       = Compose(spm_mb_io('GetData',datn.psi),Affine(df,Mat));
mu        = Pull1(model.mu,psi);
[f,datn]  = spm_mb_io('GetClasses',datn,mu,sett);
[a,w]     = Push1(f - Softmax(mu,4),psi,d);
clear psi f mu

g         = reshape(sum(a.*G,4),[d 3]);
H         = w.*H0;
clear a w

g         = g + spm_diffeo('vel2mom', v, v_settings);                             % Prior term
v         = v + v0 - scal*spm_diffeo('fmg', H, g, [v_settings s_settings]); % Gauss-Newton update

if do_pca
    % Uncertainty term for lower bound: tr(L*Sv)
    datn.ss.trLSv = double(spm_diffeo('trapprox', H, v_settings));
    datn.ss.trLSv = datn.ss.trLSv(1) / model.lam;
    % lam is removed from L so it can be updated after.
end

if d(3)==1, v(:,:,:,3) = 0; end % If 2D
if v_settings(4)==0             % Mean displacement should be 0
    avg = mean(mean(mean(v,1),2),3);
    v   = v - avg;
end
datn.v   = spm_mb_io('SetData',datn.v,v);
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
if accel>0, s  = Softmax(mu,4); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
H1 = zeros(d,'single');
H2 = H1;
H3 = H1;
H4 = H1;
H5 = H1;
H6 = H1;
for m1=1:M
    Gm11 = G(:,:,:,m1,1);
    Gm12 = G(:,:,:,m1,2);
    Gm13 = G(:,:,:,m1,3);
    if accel==0
        tmp = Ab(m1,m1);
    else
        sm1 = s(:,:,:,m1);
        tmp = (max(sm1.*(1-sm1),0))*accel + (1-accel)*Ab(m1,m1);
    end
    H1 = H1 + tmp.*Gm11.*Gm11;
    H2 = H2 + tmp.*Gm12.*Gm12;
    H3 = H3 + tmp.*Gm13.*Gm13;
    H4 = H4 + tmp.*Gm11.*Gm12;
    H5 = H5 + tmp.*Gm11.*Gm13;
    H6 = H6 + tmp.*Gm12.*Gm13;
    for m2=(m1+1):M
        if accel==0
            tmp = Ab(m1,m2);
        else
            sm2 = s(:,:,:,m2);
            tmp = (-sm1.*sm2)*accel + (1-accel)*Ab(m1,m2);
        end
        Gm21 = G(:,:,:,m2,1);
        Gm22 = G(:,:,:,m2,2);
        Gm23 = G(:,:,:,m2,3);
        H1 = H1 + 2*tmp.*Gm11.*Gm21;
        H2 = H2 + 2*tmp.*Gm12.*Gm22;
        H3 = H3 + 2*tmp.*Gm13.*Gm23;
        H4 = H4 + tmp.*(Gm11.*Gm22 + Gm21.*Gm12);
        H5 = H5 + tmp.*(Gm11.*Gm23 + Gm21.*Gm13);
        H6 = H6 + tmp.*(Gm12.*Gm23 + Gm22.*Gm13);
    end
end
clear Gm11 Gm12 Gm13 Gm21 Gm22 Gm23 sm1 sm2 tmp
H = cat(4, H1, H2, H3, H4, H5, H6);
end
%==========================================================================

%==========================================================================
function [T, iT] = JointDiag(ZZ, ULU)
% FORMAT [T, iT] = JointDiag(ZZ,ULU)
%
% ** Input **
% ezz - E[ZZ']  = E[Z]E[Z]'   + inv(Az)
% ww  - E[U'LU] = E[U]'LE[U]' + inv(Au)
% ** Output **
% T   - Orthogonalisation matrix,  s.t.   T * ZZ  * T' ~ diag
% iT  - Almost inverse of T,       s.t. iT' * ULU * iT ~ diag


[Vz, Dz2]  = svd(double(ZZ));
[Vw, Dw2]  = svd(double(ULU));
Dz         = diag(sqrt(diag(Dz2) + eps));
Dw         = diag(sqrt(diag(Dw2) + eps));
[U, D, V]  = svd(Dw * Vw' * Vz * Dz');
% Dz         = LoadDiag(Dz);
% Dw         = LoadDiag(Dw);
T          = D * V' * (Dz \ Vz');
iT         = Vw * (Dw \ U);
end
%==========================================================================
    
%==========================================================================
function [Q, iQ, q] = ScaleSubspace(ULU, ZZ, A0, n0, ndat, nvox, q0)
% Gauss-Newton optimisation of the scaling factor between the subspace (U)
% and the latent variables (Z)
%
% FORMAT [Q, iQ] = ScaleSubspace(ULU, ZZ, A0, n0, ndat, nvox, q0)
% ULU  - E[U'LU] (subspace suff stat)
%        > Must have been orthogonalised before.
% ZZ   - E[Z*Z'] (latent suff stat)
%        > Must have been orthogonalised before.
% A0   - Expected value of Wishart prior (n0*V0)
% n0   - Degrees of freedom of Wishart prior
% ndat - Number of subjects
% nvox - Number of voxels * 3
% q0   - Initial parameters [defaut: -0.5*log(nvox)]

M = size(ZZ, 1);

if nargin < 7 || isempty(q0)
    q0    = zeros(M,1)-0.5*log(nvox);
end

q = min(max(q0,-10),10);  % Heuristic to avoid bad starting estimate
Q = diag(exp(q));
A = UpdateWishart(Q*ZZ*Q, ndat, A0, n0);
E = 0.5*(trace(Q*ZZ*Q*A) + trace(ULU/(Q*Q))) ...
  + (nvox - ndat) * LogDetChol(Q);
% fprintf('[%3d %2d] %15.6g',0,0,E)

all_E0 = E;
for iter=1:100 % EM loop
    A   = UpdateWishart(Q*ZZ*Q, ndat, A0, n0);
    oE0 = E;

    all_E = E;
    for subit=1:10 % Gauss-Newton loop
        R  = A.*ZZ'+A'.*ZZ;
        g1 = Q*R*diag(Q);
        g2 =-2*(Q^2\diag(ULU));
        g  = g1 + g2 + 2 * (nvox - ndat);

        H1 = Q*R*Q + diag(g1);
        H2 = 4*(Q^2\ULU);
        H  = H1+H2;

        H  = LoadDiag(H);
        q  = q - H\g;
        q  = min(max(q,-10),10); % Heuristic to avoid overshoot
        Q  = diag(exp(q));

        oE = E;
        E  = 0.5*(trace(Q*ZZ*Q*A) + trace(ULU/(Q*Q))) ...
              + (nvox - ndat) * LogDetChol(Q);
        all_E = [all_E E];

        % fprintf('\n[%3d %2d] %15.6g (%8.1e)',iter,subit,E, ...
        %     abs(oE-E)/(max(all_E)-min(all_E)+eps))
        if abs(oE-E)/(max(all_E)-min(all_E)+eps) < 1e-8, break; end
    end
    all_E0 = [all_E0 E];
    % fprintf(' (%8.1e)', abs(oE0-E)/(max(all_E0)-min(all_E0)+eps))
    if abs(oE0-E)/(max(all_E0)-min(all_E0)+eps) < 1e-20, break; end
end
iQ = inv(Q);
% fprintf('\n');
end

% % Code for working out the gradients and Hessians
% q   = sym('q',[3,1],'real');
% Q   = diag(exp(q));
% % A   = diag(sym('a',[3,1],'real'));
% % ZZ  = diag(sym('z',[3,1],'real'));
% % WW   = diag(sym('w',[3,1],'real'));
% A   = sym('a',[3,3],'real');
% ZZ  = sym('z',[3,3],'real');
% WW   = sym('w',[3,3],'real');
% DmN = sym('DmN','real');
% %%
% E   = trace(Q*ZZ*Q*A) + trace(WW*inv(Q*Q)) + 2 * DmN * log(det(Q));
% %%
% pretty(simplify(diff(E,sym('q1')),1000))
% pretty(simplify(diff(diff(E,sym('q1')),sym('q2')),1000))
% pretty(simplify(diff(diff(E,sym('q1')),sym('q1')),1000))
% %%
% g1 =  Q*(A.*ZZ'+A'.*ZZ)*diag(Q);
% g2 = -Q^2\diag(WW)*2;
% g  =  g1+g2+2*DmN;
% H1 =  Q*(A'.*ZZ + A.*ZZ')*Q +diag(g1);
% H2 =  4*WW*Q^(-2);
% H  =  H1+H2;
% %%
% d1  = simplify(g(1)  -diff(E,sym('q1')),1000)
% d11 = simplify(H(1,1)-diff(diff(E,sym('q1')),sym('q1')),1000)
%==========================================================================


% === LogGamma ============================================================
function lg = LogGamma(a, p)
% Multivariate log Gamma function.
%
% FORMAT lg = LogGamma(a, p)
% a - value
% p - dimension [1]
if nargin < 2
    p = 1;
end
lg = (p*(p-1)/4)*log(pi);
for i=1:p
    lg = lg + gammaln(a + (1-p)/2);
end
end

% === DiGamma =============================================================
function dg = DiGamma(a, p, k)
% Multivariate di/tri/etc Gamma function.
%
% FORMAT lg = DiGamma(a, p, k)
% a - value
% p - dimension  [1]
% k - derivative [0=digamma], 1=trigamma, 2=tetragamma, etc.
if nargin < 3
    k = 0;
    if nargin < 2
        p = 1;
    end
end
dg = 0;
for i=1:p
    dg = dg + psi(k, a + (1-i)/2);
end
end

% === LogDetChol ==========================================================
function ld = LogDetChol(A)
% Log-determinant of a positive-definite matrix.
% Cholesky factorisation is used to compute a more stable log-determinant.
%
% FORMAT ld = LogDetChol(A)
% A  - A postive-definite square matrix
% ld - Logarithm of determinant of A

% Cholseki decomposition of A (A = C' * C, with C upper-triangular)
[C, p] = chol(double(A));

if p > 0
   % A should usually be positive definite, but check anyway.
   warning(['Attempting to compute log determinant of matrix ' ...
            'that is not positive definite (p=%d).'], p);
end

% Because C is triangular, |C| = prod(diag(C))
% Hence: log|C| = sum(log(diag(C)))
% And:   log|A| = log|C'*C| = log(|C|^2) = 2 * sum(log(diag(C)))
ld = 2 * sum(log(diag(C)));
end

% === LoadDiag ===========================================================
function A = LoadDiag(A)
% A  - A square matrix
%
% Load A's diagonal until it is well conditioned for inversion.

factor = 1e-7;
while rcond(A) < 1e-5
    A = A + factor * max([diag(A); eps]) * eye(size(A));
    factor = 10 * factor;
end
end

% === ELogGamma ===========================================================
function ld = ELogGamma(lam, n, k)
% Expected log of a Gamma random variable.
%   
%   /!\ It is parameterised by the expected value and degrees of freedom,
%       not the shape and scale.
%   If n == 0 or n == Inf, returns simple log.
%
% FORMAT ld = ELogGamma(A,n,[k])
% A  - Expected value of the Wishart distribution
% n  - Degrees of freedom
% k  - Dimension [1]
% ld - Expected logarithm of G(nk/2,nk/(2lam))
if nargin < 3, k = 1; end
ld = log(lam);
if isfinite(n) && n > 0
    ld = ld + DiGamma(k*n/2) - log(k*n/2);
end
end

% === ELogDetWishart ======================================================
function ld = ELogDetWishart(A, n)
% Expected log-determinant of a Wishart random variable.
%   
%   /!\ It is parameterised by the expected value, not the scale matrix.
%   If n == 0 or n == Inf, returns simple log-determinant.
%
% FORMAT ld = LogDetWishart(A,n)
% A  - Expected value of the Wishart distribution
% n  - Degrees of freedom
% ld - Expected log determinant of W(A/n,n)
ld = LogDetChol(A);
if isfinite(n) && n > 0
    k  = size(A, 1);
    ld = ld + DiGamma(n/2, k) - k*log(n/2);
end
end


% === UpdateWishart =======================================================
function [A,n] = UpdateWishart(ZZ, ndat, A0, n0)
% Posterior parameters of a Wishart distribution
%
% FORMAT [A,n] = UpdateWishart(ZZ, ndat, A0, n0)
% ZZ   - 2nd order moment of observed data (assumes zero mean)
% ndat - 0-th order moment of observed data
% A0   - Expected value of the prior
% n0   - Degrees of freedom of the prior
n = n0 + ndat;
if n0 > 0,  A = (n*A0)/(A0*ZZ +  n0*eye(size(ZZ)));
else,       A = n*inv(ZZ + 1E-3*eye(size(ZZ))); % slight regularisation
end
if nargout > 1 && n0 == 0
    n = Inf;
end
end