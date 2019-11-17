function varargout = spm_multireg_updt(varargin)
%__________________________________________________________________________
%
% Update functions for spm_multireg.
%
% FORMAT dat      = spm_multireg_updt('UpdateAffines',dat,mu,sett)
% FORMAT dat      = spm_multireg_updt('UpdateBiasField',dat,mu,sett)
% FORMAT dat      = spm_multireg_updt('UpdateGMM',dat,mu,sett)
% FORMAT dat      = spm_multireg_updt('UpdateIntensity',dat, sett)
% FORMAT [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett)
% FORMAT dat      = spm_multireg_updt('UpdateSimpleAffines',dat,mu,sett)
% FORMAT [mu,dat] = spm_multireg_updt('UpdateSimpleMean',dat, mu, sett)
% FORMAT dat      = spm_multireg_updt('UpdateVelocities',dat,mu,sett)
% FORMAT dat      = spm_multireg_updt('UpdateWarps',dat,sett)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_updt
    error('Not enough argument. Type ''help spm_multireg_updt'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'UpdateAffines'
        [varargout{1:nargout}] = UpdateAffines(varargin{:});
    case 'UpdateBiasField'
        [varargout{1:nargout}] = UpdateBiasField(varargin{:});           
    case 'UpdateGMM'  
        [varargout{1:nargout}] = UpdateGMM(varargin{:});        
    case 'UpdateIntensity'
        [varargout{1:nargout}] = UpdateIntensity(varargin{:});        
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
    otherwise
        help spm_multireg_updt
        error('Unknown function %s. Type ''help spm_multireg_updt'' for help.', id)
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
threads     = sett.gen.threads;

if ~do_updt_aff, return; end

% Update the affine parameters
spm_multireg_util('SetBoundCond');
if ~isempty(B)
    if threads>1 && numel(dat)>1
        % Memory = loads
        parfor(n=1:numel(dat),threads) % PARFOR
            spm_multireg_util('SetBoundCond');
            dat(n) = UpdateAffinesSub(dat(n),mu,sett);
        end
    else
        for n=1:numel(dat)
            dat(n) = UpdateAffinesSub(dat(n),mu,sett);
        end
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
% UpdateBiasField()
function dat = UpdateBiasField(dat,mu,sett)

% Parse function settings
do_bf_norm = sett.do.bf_norm;
do_updt_bf = sett.do.updt_aff;
threads    = sett.gen.threads;

if ~do_updt_bf, return; end

if threads>1 && numel(dat)>1    
    parfor(n=1:numel(dat),threads) % PARFOR
        spm_multireg_util('SetBoundCond');
        dat(n) = UpdateBiasFieldSub(dat(n),mu,sett);
    end
else
    for n=1:numel(dat)
        dat(n) = UpdateBiasFieldSub(dat(n),mu,sett);
    end
end    

if do_bf_norm
    % Zero-mean the field DC component
    dat = ZeroMeanDC(dat,mu,sett);
    
    % Update GMM after mean correction    
    dat = UpdateGMM(dat,mu,sett);
end
end
%==========================================================================

%==========================================================================
% UpdateGMM()
function dat = UpdateGMM(dat,mu,sett)

% Parse function settings
threads = sett.gen.threads;

if threads>1 && numel(dat)>1    
    parfor(n=1:numel(dat),threads) % PARFOR
        spm_multireg_util('SetBoundCond');
        dat(n) = UpdateGMMSub(dat(n),mu,sett);
    end
else
    for n=1:numel(dat)
        dat(n) = UpdateGMMSub(dat(n),mu,sett);
    end
end 
end
%==========================================================================

%==========================================================================
% UpdateIntensity()
function dat = UpdateIntensity(dat, sett)

% Parse function settings
fig_name = sett.show.figname_int;

if ~isfield(dat(1),'mog'), return; end

N  = numel(dat);
po = cell(1,N);
for n=1:N
    po{n}{1}{1} = dat(n).mog.po.m;
    po{n}{1}{2} = dat(n).mog.po.b;
    po{n}{2}{1} = dat(n).mog.po.V;
    po{n}{2}{2} = dat(n).mog.po.n;
end
pr = dat(1).mog.pr;
pr = {pr.m,pr.b,pr.V,pr.n};
pr = spm_gmm_lib('updatehyperpars',po,pr,...
                 'verbose', true, ...
                 'figname', fig_name);
for n=1:N
    dat(n).mog.pr.m = pr{1};
    dat(n).mog.pr.b = pr{2};
    dat(n).mog.pr.V = pr{3};
    dat(n).mog.pr.n = pr{4};
end             
end
%==========================================================================

%==========================================================================
% UpdateMean()
function [mu,dat] = UpdateMean(dat, mu, sett)

% Parse function settings
accel       = sett.gen.accel;
mu_settings = sett.var.mu_settings;
s_settings  = sett.shoot.s_settings;
threads     = sett.gen.threads;

spm_multireg_util('SetBoundCond');
g  = spm_field('vel2mom', mu, mu_settings);
M  = size(mu,4);
H  = zeros([sett.var.d M*(M+1)/2],'single');
H0 = spm_multireg_der('AppearanceHessian',mu,accel);
if threads>1 && numel(dat)>1
    % Memory = 4*prod(sett.var.d)*2*(M*(M+1)/2 + M) + max(prod(sett.var.d)*3+3*prod(dn)*2,(3+4*M)*prod(dn))); % Needs more work
    parfor(n=1:numel(dat),threads) % PARFOR
        spm_multireg_util('SetBoundCond');
        [gn,Hn,dat(n)] = UpdateMeanSub(dat(n),mu,H0,sett);
        g              = g + gn;
        H              = H + Hn;
    end
else
    for n=1:numel(dat)
        [gn,Hn,dat(n)] = UpdateMeanSub(dat(n),mu,H0,sett);
        g              = g + gn;
        H              = H + Hn;
    end
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
threads     = sett.gen.threads;

if ~do_updt_aff, return; end

% Update the affine parameters
G  = spm_diffeo('grad',mu);
H0 = spm_multireg_der('VelocityHessian',mu,G,accel);

if ~isempty(B)
    if threads>1 && numel(dat)>1
        % Memory = loads
        parfor(n=1:numel(dat),threads) % PARFOR
            dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
        end
    else
        for n=1:numel(dat)
            dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
        end
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
threads     = sett.gen.threads;

spm_multireg_util('SetBoundCond');
w  = zeros(sett.var.d,'single');
gf = zeros(size(mu),'single');
if threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),threads) % PARFOR
        [gn,wn,dat(n)] = UpdateSimpleMeanSub(dat(n),mu,sett);
        gf             = gf + gn;
        w              = w  + wn;
    end
else
    for n=1:numel(dat)
        [gn,wn,dat(n)] = UpdateSimpleMeanSub(dat(n),mu,sett);
        gf             = gf + gn;
        w              = w  + wn;
    end
end
for it=1:ceil(4+2*log2(numel(dat)))
    H  = w.*spm_multireg_der('AppearanceHessian',mu,accel);
    g  = w.*spm_multireg_util('softmax',mu,4) - gf;
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
threads     = sett.gen.threads;

if ~do_updt_vel, return; end

spm_multireg_util('SetBoundCond');
G  = spm_diffeo('grad',mu);
H0 = spm_multireg_der('VelocityHessian',mu,G,accel);
if size(G,3) == 1
    % Data is 2D
    H0(:,:,:,3) = H0(:,:,:,3) + mean(reshape(H0(:,:,:,[1 2]),[],1));
end
if threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),threads) % PARFOR
        spm_multireg_util('SetBoundCond');
        dat(n) = UpdateVelocitiesSub(dat(n),mu,G,H0,sett);
    end
else
    for n=1:numel(dat)
        dat(n) = UpdateVelocitiesSub(dat(n),mu,G,H0,sett);
    end
end
end
%==========================================================================

%==========================================================================
% UpdateWarps()
function dat = UpdateWarps(dat,sett)

% Parse function settings
groupwise  = sett.model.groupwise;
threads    = sett.gen.threads;
v_settings = sett.var.v_settings;

if groupwise
    % Total initial velocity should be zero (Khan & Beg)
    avg_v = single(0);
    for n=1:numel(dat)
        avg_v = avg_v + spm_multireg_io('GetData',dat(n).v); % For mean correcting initial velocities
    end
    avg_v = avg_v/numel(dat);
    d     = [size(avg_v,1) size(avg_v,2) size(avg_v,3)];
else
    avg_v = [];
    d     = spm_multireg_io('GetSize',dat(1).v);
end
kernel = spm_multireg_util('Shoot',d,v_settings);
if threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),threads) % PARFOR
        dat(n) = UpdateWarpsSub(dat(n),avg_v,sett,kernel);
    end
else
    for n=1:numel(dat)
        dat(n) = UpdateWarpsSub(dat(n),avg_v,sett,kernel);
    end
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% Mask()
function f = Mask(f,msk)
f(~isfinite(f)) = 0;
f = f.*msk;
end
%==========================================================================

%==========================================================================
function dat = ZeroMeanDC(dat,mu,sett)

% Parse function settings
threads = sett.gen.threads;

% Get correction factor
sm_dc_ln  = 0;
sm_dc_int = 0;
for n=1:numel(dat)
    sm_dc_ln  = sm_dc_ln  + dat(n).bf.dc.ln;
    sm_dc_int = sm_dc_int + dat(n).bf.dc.int;
end   
mn_dc_ln  = sm_dc_ln./numel(dat);
mn_dc_int = sm_dc_int./numel(dat);
C         = numel(mn_dc_ln); 

if 0     
    % Some verbose (should tend towards 1, for all channels)
    fprintf('mn_dc_int = [');
    for c=1:C - 1
        fprintf('%4.5f ',mn_dc_int(c));
    end
    fprintf('%4.5f',mn_dc_int(C))
    fprintf(']\n');
end

% Adjust bias field DC component
scl = (1./mn_dc_int)';
for n=1:numel(dat)
    for c=1:C
        dat(n).bf.chan(c).T(1,1,1) = dat(n).bf.chan(c).T(1,1,1) - mn_dc_ln(c);
    end
end
   
% Correct GMM posterior parameters and lower bound/objective
if threads>1 && numel(dat)>1    
    parfor(n=1:numel(dat),threads) % PARFOR
        spm_multireg_util('SetBoundCond');
        dat(n) = ZeroMeanDCSub(dat(n),mu,scl,sett);
    end
else
    for n=1:numel(dat)
        dat(n) = ZeroMeanDCSub(dat(n),mu,scl,sett);
    end
end    
end
%==========================================================================

%==========================================================================
function datn = ZeroMeanDCSub(datn,mu,scl,sett)

% Parse function settings
B   = sett.registr.B;
Mmu = sett.var.Mmu;

L     = datn.bf.L;
[d,C] = spm_multireg_io('GetSize',datn.f);

% GMM posterior
m = datn.mog.po.m;
b = datn.mog.po.b;
V = datn.mog.po.V;
n = datn.mog.po.n;

% GMM prior
m0 = datn.mog.pr.m;
b0 = datn.mog.pr.b;
V0 = datn.mog.pr.V;
n0 = datn.mog.pr.n;

% Get suffstats
lSS0 = datn.bf.lSS0;
lSS1 = datn.bf.lSS1;
lSS2 = datn.bf.lSS2;

% Rescale suffstats
K = size(m,2);
A = bsxfun(@times,V,reshape(n,[1 1 K]));
for l=2:numel(L)
    obs = spm_gmm_lib('code2bin', L(l), C);
    lSS1{l} = bsxfun(@times,lSS1{l},scl(obs));
    for k=1:size(lSS2{l},3)
        lSS2{l}(:,:,k) = (scl(obs)*scl(obs)').*lSS2{l}(:,:,k);
    end     
end

% Update GMM posteriors
[SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', lSS0, lSS1, lSS2, {m,A}, L);        
[m,~,b,V,n]  = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2, {m0,b0,V0,n0});    

% Save updated GMM posteriors
datn.mog.po.m = m;
datn.mog.po.b = b;
datn.mog.po.V = V;
datn.mog.po.n = n;

% Update lower bound/objective
%--------------------------------------------------------------------------

% Get subject-space template (softmaxed K + 1)
q    = double(datn.q);
Mn   = datn.Mat;
Mr   = spm_dexpm(q,B);
psi1 = spm_multireg_io('GetData',datn.psi);
psi0 = spm_multireg_util('Affine',d,Mmu\Mr*Mn);
psi  = spm_multireg_util('Compose',psi1,psi0);
clear psi0 psi1
mu   = spm_multireg_util('Pull1',mu,psi);
mu   = log(spm_multireg_util('softmaxmu',mu,4));
clear psi
mu   = reshape(mu,[prod(d(1:3)) size(mu,4)]);

% Get bias field
[bf,pr_bf] = spm_multireg_io('GetBiasField',datn.bf.chan,d);

% Get image(s)
fn = spm_multireg_io('GetData',datn.f);
fn = reshape(fn,[prod(d(1:3)) C]);

% Get responsibilities
fn   = spm_multireg_util('MaskF',fn);
code = spm_gmm_lib('obs2code', fn);
zn   = spm_multireg_io('ComputeResponsibilities',datn,bf.*fn,mu,code);
clear mu

lx  = spm_multireg_energ('LowerBound','X',bf.*fn,zn,code,{m,b},{V,n});
clear zn fn

lxb = spm_multireg_energ('LowerBound','XB',bf);
clear bf

datn.mog.lb.XB(end + 1) = lxb;
datn.mog.lb.X(end  + 1) = lx;
datn.mog.lb             = spm_multireg_energ('SumLowerBound',datn.mog.lb);

datn.E(1) = -datn.mog.lb.sum(end);
datn.E(3) = -sum(pr_bf);
end
%==========================================================================

%==========================================================================
% SclFromBiasFieldDC()
function scl = SclFromBiasFieldDC(chan)
C   = numel(chan);
scl = zeros(1,C);
for c=1:C
    b1 = chan(c).B1(1,1);
    b2 = chan(c).B2(1,1);
    b3 = chan(c).B3(1,1);
    t1 = chan(c).T(1,1,1);
    
    scl(c) = exp(b1*b2*b3*t1);
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

df   = spm_multireg_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,B);
dM   = zeros(12,size(B,3));
for m=1:size(B,3)
    tmp     = Mmu\dM3(:,:,m)*Mn;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi1 = spm_multireg_io('GetData',datn.psi);
psi0 = spm_multireg_util('Affine',df,Mmu\Mr*Mn);
J    = spm_diffeo('jacobian',psi1);
J    = reshape(spm_multireg_util('Pull1',reshape(J,[d 3*3]),psi0),[df 3 3]);
psi  = spm_multireg_util('Compose',psi1,psi0);
clear psi0 psi1

mu1  = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu1,sett);
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
    clear Gm;
end
clear J mu

msk       = all(isfinite(f),4);
a         = Mask(f - spm_multireg_util('softmax',mu1,4),msk);
[H,g]     = spm_multireg_der('AffineHessian',mu1,G,a,single(msk),accel);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-6 + 0.1);
q         = q + scal*(H\g);
datn.q    = q;
end
%==========================================================================

%==========================================================================
% UpdateBiasFieldSub()
function datn = UpdateBiasFieldSub(datn,mu,sett)

% Parse function settings
B          = sett.registr.B;
do_bf_norm = sett.do.bf_norm;
Mmu        = sett.var.Mmu;
nit_bf     = sett.nit.bf;
nit_ls     = sett.optim.nls_bf;

% Get subject-space template (softmaxed K + 1)
[d,C] = spm_multireg_io('GetSize',datn.f);
q     = double(datn.q);
Mn    = datn.Mat;
Mr    = spm_dexpm(q,B);
psi1  = spm_multireg_io('GetData',datn.psi);
psi0  = spm_multireg_util('Affine',d,Mmu\Mr*Mn);
psi   = spm_multireg_util('Compose',psi1,psi0);
clear psi0 psi1
mu    = spm_multireg_util('Pull1',mu,psi);
mu    = log(spm_multireg_util('softmaxmu',mu,4));
clear psi

% Get bias field
chan       = datn.bf.chan;
kron       = @(a,b) spm_krutil(a,b);
[bf,pr_bf] = spm_multireg_io('GetBiasField',chan,d);

% Get responsibilities
[zn,datn,code] = spm_multireg_io('GetClasses',datn,mu,sett,true);
L              = unique(code);
nL             = numel(L);
mu             = reshape(mu,[prod(d(1:3)) size(mu,4)]);

% Get image(s)
fn = spm_multireg_io('GetData',datn.f);
fn = reshape(fn,[prod(d(1:3)) C]);

% GMM posterior
m  = datn.mog.po.m;
b  = datn.mog.po.b;
V  = datn.mog.po.V;
n  = datn.mog.po.n;
K1 = size(m,2);

lx  = spm_multireg_energ('LowerBound','X',bf.*fn,zn,code,{m,b},{V,n});
lxb = spm_multireg_energ('LowerBound','XB',bf);
            
% -----------------------------------------------------------------
% Update bias field parameters
for it=1:nit_bf
    
    % -----------------------------------------------------------------
    % Update bias field parameters for each channel separately
    for c=1:C % Loop over channels

        % -----------------------------------------------------------------
        % Compute gradient and Hessian    
        gr_l = zeros(d(1:3));
        H_l  = zeros(d(1:3));

        % -----------------------------------------------------------------
        % For each combination of missing voxels
        for l=1:nL

            % -------------------------------------------------------------
            % Get mask of missing modalities (with this particular code)        
            observed_channels = spm_gmm_lib('code2bin', L(l), C);
            missing_channels  = ~observed_channels;
            if missing_channels(c), continue; end
            if isempty(code), selected_voxels = ones(size(code), 'logical');
            else,                   selected_voxels = (code == L(l));
            end
            nb_channels_missing  = sum(missing_channels);
            nb_voxels_coded      = sum(selected_voxels);
            if nb_voxels_coded == 0, continue; end

            % -------------------------------------------------------------
            % Convert channel indices to observed indices
            mapped_c     = 1:C;
            mapped_c     = mapped_c(observed_channels);
            mapped_c     = find(mapped_c == c);
            cc           = mapped_c; % short alias

            selected_obs = bf(selected_voxels,observed_channels).*fn(selected_voxels,observed_channels);
            gi = 0; % Gradient accumulated accross clusters
            Hi = 0; % Hessian accumulated accross clusters
            for k=1:K1

                % ---------------------------------------------------------
                % Compute expected precision (see GMM+missing data)
                Voo = V(observed_channels,observed_channels,k);
                Vom = V(observed_channels,missing_channels,k);
                Vmm = V(missing_channels,missing_channels,k);
                Vmo = V(missing_channels,observed_channels,k);
                Ao  = Voo - Vom*(Vmm\Vmo);
                Ao  = (n(k)-nb_channels_missing) * Ao;
                MUo = m(observed_channels,k);

                % ---------------------------------------------------------
                % Compute statistics
                gk = bsxfun(@minus, selected_obs, MUo.') * Ao(cc,:).';
                Hk = Ao(cc,cc);

                selected_resp = zn(selected_voxels,k);
                gk = bsxfun(@times, gk, selected_resp);
                Hk = bsxfun(@times, Hk, selected_resp);
                clear selected_resp

                % ---------------------------------------------------------
                % Accumulate across clusters
                gi = gi + gk;
                Hi = Hi + Hk;
                clear sk1x sk2x
            end

            % -------------------------------------------------------------
            % Multiply with bias corrected value (chain rule)
            gi = gi .* selected_obs(:,cc);
            Hi = Hi .* (selected_obs(:,cc).^2);
            clear selected_obs

            % -------------------------------------------------------------
            % Normalisation term
            gi = gi - 1;

            % -------------------------------------------------------------
            % Accumulate across missing codes
            gr_l(selected_voxels) = gr_l(selected_voxels) + gi;
            H_l(selected_voxels)  = H_l(selected_voxels) + Hi;
            clear selected_voxels
        end
        clear zn
        
        d3 = numel(chan(c).T); % Number of DCT parameters
        H  = zeros(d3,d3);     
        gr = zeros(d3,1);      
        for z=1:d(3)
            b3 = double(chan(c).B3(z,:)');
            gr = gr + kron(b3,spm_krutil(gr_l(:,:,z),double(chan(c).B1),double(chan(c).B2),0));
            H  = H  + kron(b3*b3',spm_krutil(H_l(:,:,z),double(chan(c).B1),double(chan(c).B2),1));
        end
        clear b3                       

        % -------------------------------------------------------------
        % Gauss-Newton update of bias field parameters
        Update = reshape((H + chan(c).C)\(gr + chan(c).C*chan(c).T(:)),size(chan(c).T));
        clear H gr

        % Line-search
        %------------------------------------------------------------------    
        armijo = 1;        

        % Old parameters
        oT     = chan(c).T;       
        obf    = bf;
        opr_bf = pr_bf;
        olxb   = lxb;   
        olx    = lx;

        for ls=1:nit_ls

            % Update bias-field parameters
            chan(c).T = chan(c).T - armijo*Update;

            % Compute new bias-field (only for channel c)
            [bf,pr_bf] = spm_multireg_io('GetBiasField',chan,d,obf,c,opr_bf);

            % Recompute responsibilities (with updated bias field)
            zn = spm_multireg_io('ComputeResponsibilities',datn,bf.*fn,mu,code);
            
            % Compute new lower bound
            lx  = spm_multireg_energ('LowerBound','X',bf.*fn,zn,code,{m,b},{V,n});
            lxb = spm_multireg_energ('LowerBound','XB',bf);
            
            % Check new lower bound
            if (lx + lxb + sum(pr_bf)) > (olx + olxb + sum(opr_bf))       
                datn.mog.lb.XB(end + 1) = lxb;
                datn.mog.lb.X(end  + 1) = lx;
                break;
            else                                
                armijo    = armijo*0.5;
                chan(c).T = oT;

                if ls == nit_ls                  
                    bf    = obf;         
                    pr_bf = opr_bf;
                    zn    = spm_multireg_io('ComputeResponsibilities',datn,bf.*fn,mu,code);
                end
            end
        end
        clear oT Update obf    
    end
end

% Set output (update datn)
datn.bf.chan = chan;
datn.mog.lb  = spm_multireg_energ('SumLowerBound',datn.mog.lb);
datn.E(1)    = -datn.mog.lb.sum(end);
datn.E(3)    = -sum(pr_bf);

if do_bf_norm
    [lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', bf.*fn, zn, 1, {code,L});   
    datn.bf.lSS0     = lSS0;
    datn.bf.lSS1     = lSS1;
    datn.bf.lSS2     = lSS2;
    datn.bf.L        = L;

    % Get DC component
    dc    = struct;
    dc.ln = zeros(1,C);
    for c=1:C
        dc.ln(c) = chan(c).T(1,1,1);
    end
    dc.int     = SclFromBiasFieldDC(chan);
    datn.bf.dc = dc;
end

if 0
    % Show stuff    
    z  = ceil(d(3).*0.5);
    fn = reshape(fn,[d(1:3) C]);
    bf = reshape(bf,[d(1:3) C]);
    
    figure(664);
    for c=1:C
        subplot(C,3,1 + (c - 1)*3); imagesc(fn(:,:,z,c)'); 
        title('fn');
        axis image xy off; drawnow

        subplot(C,3,2 + (c - 1)*3); imagesc(bf(:,:,z,c)'); 
        title('bf');
        axis image xy off; drawnow

        subplot(C,3,3 + (c - 1)*3); imagesc(bf(:,:,z,c)'.*fn(:,:,z,c)'); 
        title('bf.*fn');
        axis image xy off; drawnow
    end
end

end
%==========================================================================

%==========================================================================
% UpdateGMMSub()
function datn = UpdateGMMSub(datn,mu,sett)

% Parse function settings
B   = sett.registr.B;
Mmu = sett.var.Mmu;

% Get subject-space template (softmaxed K + 1)
d    = spm_multireg_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
Mr   = spm_dexpm(q,B);
psi1 = spm_multireg_io('GetData',datn.psi);
psi0 = spm_multireg_util('Affine',d,Mmu\Mr*Mn);
psi  = spm_multireg_util('Compose',psi1,psi0);
clear psi0 psi1
mu   = spm_multireg_util('Pull1',mu,psi);
clear psi

% Update GMM parameters
[~,datn] = spm_multireg_io('GetClasses',datn,mu,sett);
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

df  = spm_multireg_io('GetSize',datn.f);
q   = double(datn.q);
Mn  = datn.Mat;
psi = spm_multireg_util('Compose',spm_multireg_io('GetData',datn.psi), ...
                                  spm_multireg_util('Affine',df, Mmu\spm_dexpm(q,B)*Mn));
mu  = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu,sett);
if isempty(H0)
    g     = spm_multireg_util('Push1',spm_multireg_util('softmax',mu,4) - f,psi,d);
    H     = spm_multireg_util('Push1',spm_multireg_der('AppearanceHessian',mu,accel),psi,d);
else
    % Faster approximation - but might be unstable
    % If there are problems, then revert to the slow
    % way.
    [g,w] = spm_multireg_util('Push1',spm_multireg_util('softmax',mu,4) - f,psi,d);
    H     = w.*H0;
end
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

df   = spm_multireg_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,B);
dM   = zeros(12,size(B,3));
for m=1:size(B,3)
    tmp     = Mmu\dM3(:,:,m)*Mmu;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi      = spm_multireg_util('Affine',df,Mmu\Mr*Mn);
mu1      = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu1,sett);

[a,w]     = spm_multireg_util('Push1',f - spm_multireg_util('softmax',mu1,4),psi,d);
clear mu1 psi f

[H,g]     = spm_multireg_der('SimpleAffineHessian',mu,G,H0,a,w);
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

df    = spm_multireg_io('GetSize',datn.f);
q     = double(datn.q);
Mn    = datn.Mat;
psi   = spm_multireg_util('Compose',spm_multireg_io('Affine',datn.psi), ...
                                    spm_multireg_util('affine',df, Mmu\spm_dexpm(q,B)*Mn));
mu    = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu,sett);
[g,w] = spm_multireg_util('Push1',f,psi,d);
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

v         = spm_multireg_io('GetData',datn.v);
q         = datn.q;
Mn        = datn.Mat;
Mr        = spm_dexpm(q,B);
Mat       = Mmu\Mr*Mn;
df        = spm_multireg_io('GetSize',datn.f);
psi       = spm_multireg_util('Compose',spm_multireg_io('GetData',datn.psi), spm_multireg_util('Affine',df,Mat));
mu        = spm_multireg_util('Pull1',mu,psi);
[f,datn]  = spm_multireg_io('GetClasses',datn,mu,sett);
[a,w]     = spm_multireg_util('Push1',f - spm_multireg_util('softmax',mu,4),psi,d);
clear psi f mu

g         = reshape(sum(a.*G,4),[d 3]);
H         = w.*H0;
clear a w

u0        = spm_diffeo('vel2mom', v, v_settings);                          % Initial momentum
datn.E(2) = 0.5*sum(u0(:).*v(:));                                          % Prior term
v         = v - scal*spm_diffeo('fmg',H, g + u0, [v_settings s_settings]); % Gauss-Newton update

if d(3)==1, v(:,:,:,3) = 0; end % If 2D
if v_settings(1)==0             % Mean displacement should be 0
    avg = mean(mean(mean(v,1),2),3);
    v   = v - avg;
end
datn.v = spm_multireg_io('SetData',datn.v,v);
end
%==========================================================================

%==========================================================================
% UpdateWarpsSub()
function datn = UpdateWarpsSub(datn,avg_v,sett,kernel)
v        = spm_multireg_io('GetData',datn.v);
if ~isempty(avg_v)
    v    = v - avg_v;
end
datn.v   = spm_multireg_io('SetData',datn.v,v);
psi1     = spm_multireg_util('Shoot',v, kernel, sett.shoot.args); % Geodesic shooting
datn.psi = spm_multireg_io('SetData',datn.psi,psi1);
end
%==========================================================================