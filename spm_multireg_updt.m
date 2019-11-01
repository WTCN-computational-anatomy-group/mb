function varargout = spm_multireg_updt(varargin)
%__________________________________________________________________________
%
% Update functions for spm_multireg.
%
% FORMAT dat      = spm_multireg_updt('UpdateAffines',dat,mu,sett)
% FORMAT dat      = spm_multireg_updt('UpdateGMM',dat,F,mu,adjust,sett)
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
    case 'UpdateGMM'
        [varargout{1:nargout}] = UpdateGMM(varargin{:});        
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

if ~sett.do.updt_aff, return; end

% Update the affine parameters
spm_multireg_util('SetBoundCond');
if ~isempty(sett.registr.B)
    if sett.gen.threads>1 && numel(dat)>1
        % Memory = loads
        parfor(n=1:numel(dat),sett.gen.threads) % PARFOR
            spm_multireg_util('SetBoundCond');
            dat(n) = UpdateAffinesSub(dat(n),mu,sett);
        end
    else
        for n=1:numel(dat)
            dat(n) = UpdateAffinesSub(dat(n),mu,sett);
        end
    end

    if sett.model.groupwise
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
% UpdateGMM()
function datn = UpdateGMM(datn,Fn,mu,adjust,sett)
K   = size(mu,2);
mog = datn.mog; % Get GMM parameters
p   = zeros([numel(Fn),K + 1],'single');
if isempty(adjust)
    maxmu  = max(max(mu,[],2),0);
    adjust = sum(log(sum(exp(mu-maxmu),2) + exp(-maxmu))+maxmu);
end
for it=1:sett.nit.gmm
    
    for k=1:(K+1)
        p(:,k) = -((Fn-mog.mu(k)).^2)/(2*mog.sig2(k)+eps) - 0.5*log(2*pi*mog.sig2(k)+eps);
        if k<=K, p(:,k) = p(:,k) + mu(:,k); end
    end
    
    pmx = max(p,[],2);
    p   = p - pmx;
    p   = exp(p);
    sp  = sum(p,2);
    prevE     = datn.E(1);
    datn.E(1) = -sum(log(sp)+pmx,1) + adjust; % doesn't account for priors (so can increase)
   %fprintf(' %g', datn.E(1));
   
    p   = p./sp;
    for k=1:(K+1)
        sp          = sum(p(:,k))+eps;
        mog.mu(k)   = sum(Fn.*p(:,k))./sp;
        mog.sig2(k) = (sum((Fn-mog.mu(k)).^2.*p(:,k))+10000*80^2)./(sp+10000); % ad hoc "Wishart priors"
    end
    
   %disp([mog.mu; mog.sig2])
    if it>1 && abs(prevE-datn.E(1))/size(p,1) < 1e-4, break; end
end

datn.mog = mog; % Update GMM parameters

end
%==========================================================================

%==========================================================================
% UpdateMean()
function [mu,dat] = UpdateMean(dat, mu, sett)
spm_multireg_util('SetBoundCond');
g  = spm_field('vel2mom', mu, sett.var.mu_settings);
M  = size(mu,4);
H  = zeros([sett.var.d M*(M+1)/2],'single');
H0 = spm_multireg_der('AppearanceHessian',mu,sett.gen.accel);
if sett.gen.threads>1 && numel(dat)>1
    % Memory = 4*prod(sett.var.d)*2*(M*(M+1)/2 + M) + max(prod(sett.var.d)*3+3*prod(dn)*2,(3+4*M)*prod(dn))); % Needs more work
    parfor(n=1:numel(dat),sett.gen.threads) % PARFOR
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
mu = mu - spm_field(H, g, [sett.var.mu_settings sett.shoot.s_settings]);
end
%==========================================================================

%==========================================================================
% UpdateSimpleAffines()
function dat = UpdateSimpleAffines(dat,mu,sett)

if ~sett.do.updt_aff, return; end

% Update the affine parameters
G  = spm_diffeo('grad',mu);
H0 = spm_multireg_der('VelocityHessian',mu,G,sett.gen.accel);

if ~isempty(sett.registr.B)
    if sett.gen.threads>1 && numel(dat)>1
        % Memory = loads
        parfor(n=1:numel(dat),sett.gen.threads) % PARFOR
            dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
        end
    else
        for n=1:numel(dat)
            dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
        end
    end

    if sett.model.groupwise
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
spm_multireg_util('SetBoundCond');
w  = zeros(sett.var.d,'single');
gf = zeros(size(mu),'single');
if sett.gen.threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),sett.gen.threads) % PARFOR
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
    H  = w.*spm_multireg_der('AppearanceHessian',mu,sett.gen.accel);
    g  = w.*spm_multireg_util('softmax',mu) - gf;
    g  = g  + spm_field('vel2mom', mu, sett.var.mu_settings);
    mu = mu - spm_field(H, g, [sett.var.mu_settings sett.shoot.s_settings]);
end
end
%==========================================================================

%==========================================================================
% UpdateVelocities()
function dat = UpdateVelocities(dat,mu,sett)
spm_multireg_util('SetBoundCond');
G  = spm_diffeo('grad',mu);
H0 = spm_multireg_der('VelocityHessian',mu,G,sett.gen.accel);
if sett.gen.threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),sett.gen.threads) % PARFOR
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
if sett.model.groupwise
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
kernel = spm_multireg_util('Shoot',d,sett.var.v_settings);
if sett.gen.threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),sett.gen.threads) % PARFOR
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
% UpdateAffinesSub()
function datn = UpdateAffinesSub(datn,mu,sett)
% This could be made more efficient.
d    = spm_multireg_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,sett.registr.B);
dM   = zeros(12,size(sett.registr.B,3));
for m=1:size(sett.registr.B,3)
    tmp     = sett.var.Mmu\dM3(:,:,m)*Mn;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi1 = spm_multireg_io('GetData',datn.psi);
psi0 = spm_multireg_util('Affine',d,sett.var.Mmu\Mr*Mn);
J    = spm_diffeo('jacobian',psi1);
J    = reshape(spm_multireg_util('Pull1',reshape(J,[sett.var.d 3*3]),psi0),[d 3 3]);
psi  = spm_multireg_util('Compose',psi1,psi0);
clear psi0 psi1
mu1  = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu1,sett);
M    = size(mu,4);
G    = zeros([d M 3],'single');
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
clear J

msk       = all(isfinite(f),4);
a         = spm_multireg_util('Mask',f - spm_multireg_util('softmax',mu1),msk);
[H,g]     = spm_multireg_der('AffineHessian',mu1,G,a,single(msk),sett.gen.accel);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-6 + 0.1);
q         = q + sett.optim.scal*(H\g);
datn.q    = q;
end
%==========================================================================

%==========================================================================
% UpdateMeanSub()
function [g,H,datn] = UpdateMeanSub(datn,mu,H0,sett)
d    = spm_multireg_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
psi  = spm_multireg_util('Compose',spm_multireg_io('GetData',datn.psi), spm_multireg_util('Affine',d, sett.var.Mmu\spm_dexpm(q,sett.registr.B)*Mn));
mu1  = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu1,sett);
if isempty(H0)
    g     = spm_multireg_util('Push1',spm_multireg_util('softmax',mu1) - f,psi,sett.var.d);
    H     = spm_multireg_util('Push1',spm_multireg_der('AppearanceHessian',mu1,sett.gen.accel),psi,sett.var.d);
else
    % Faster approximation - but might be unstable
    % If there are problems, then revert to the slow
    % way.
    [g,w] = spm_multireg_util('Push1',spm_multireg_util('softmax',mu1)-f,psi,sett.var.d);
    H     = w.*H0;
end
end
%==========================================================================

%==========================================================================
% UpdateSimpleAffinesSub()
function datn = UpdateSimpleAffinesSub(datn,mu,G,H0,sett)
d    = spm_multireg_io('GetSize',datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,sett.registr.B);
dM   = zeros(12,size(sett.registr.B,3));
for m=1:size(sett.registr.B,3)
    tmp     = sett.var.Mmu\dM3(:,:,m)*sett.var.Mmu;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi      = spm_multireg_util('Affine',d,sett.var.Mmu\Mr*Mn);
mu1      = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu1,sett);

[a,w]     = spm_multireg_util('Push1',f - spm_multireg_util('softmax',mu1),psi,sett.var.d);
[H,g]     = spm_multireg_der('SimpleAffineHessian',mu,G,H0,a,w);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-6 + 0.1);
q         = q + sett.optim.scal*(H\g);
datn.q    = q;
end
%==========================================================================

%==========================================================================
% UpdateSimpleMeanSub()
function [g,w,datn] = UpdateSimpleMeanSub(datn,mu,sett)
d     = spm_multireg_io('GetSize',datn.f);
q     = double(datn.q);
Mn    = datn.Mat;
psi   = spm_multireg_util('Compose',spm_multireg_io('Affine',datn.psi), spm_multireg_util('affine',d, sett.var.Mmu\spm_dexpm(q,sett.registr.B)*Mn));
mu1   = spm_multireg_util('Pull1',mu,psi);
[f,datn] = spm_multireg_io('GetClasses',datn,mu1,sett);
[g,w] = spm_multireg_util('Push1',f,psi,sett.var.d);
end
%==========================================================================

%==========================================================================
% UpdateVelocitiesSub()
function datn = UpdateVelocitiesSub(datn,mu,G,H0,sett)
v         = spm_multireg_io('GetData',datn.v);
q         = datn.q;
Mn        = datn.Mat;
Mr        = spm_dexpm(q,sett.registr.B);
Mat       = sett.var.Mmu\Mr*Mn;
d         = spm_multireg_io('GetSize',datn.f);
psi       = spm_multireg_util('Compose',spm_multireg_io('GetData',datn.psi), spm_multireg_util('Affine',d,Mat));
mu1       = spm_multireg_util('Pull1',mu,psi);
[f,datn]  = spm_multireg_io('GetClasses',datn,mu1,sett);
[a,w]     = spm_multireg_util('Push1',f - spm_multireg_util('softmax',mu1),psi,sett.var.d);
g         = reshape(sum(a.*G,4),[sett.var.d 3]);
H         = w.*H0;
u0        = spm_diffeo('vel2mom', v, sett.var.v_settings);                                  % Initial momentum
datn.E(2) = 0.5*sum(u0(:).*v(:));                                                       % Prior term
v         = v - sett.optim.scal*spm_diffeo('fmg',H, g+u0, [sett.var.v_settings sett.shoot.s_settings]); % Gauss-Newton update

if sett.var.d(3)==1, v(:,:,:,3) = 0; end % If 2D
if sett.var.v_settings(1)==0             % Mean displacement should be 0
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