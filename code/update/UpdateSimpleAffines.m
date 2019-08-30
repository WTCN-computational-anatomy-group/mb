function dat = UpdateSimpleAffines(dat,mu,sett)
% Update the affine parameters

G  = spm_diffeo('grad',mu);
H0 = VelocityHessian(mu,G,sett.accel);

if ~isempty(sett.B)
    if sett.threads>1 && numel(dat)>1
        % Memory = loads
        parfor(n=1:numel(dat),sett.threads) % PARFOR
            dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
        end
    else
        for n=1:numel(dat)
            dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
        end
    end

    if sett.groupwise
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
function datn = UpdateSimpleAffinesSub(datn,mu,G,H0,sett)
d    = GetSize(datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,sett.B);
dM   = zeros(12,size(sett.B,3));
for m=1:size(sett.B,3)
    tmp     = sett.Mmu\dM3(:,:,m)*sett.Mmu;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi      = affine(d,sett.Mmu\Mr*Mn);
mu1      = Pull1(mu,psi);
[f,datn] = GetClasses(datn,mu1,sett);

[a,w]     = Push1(f-softmax(mu1),psi,sett.d);
[H,g]     = SimpleAffineHessian(mu,G,H0,a,w);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-6 + 0.1);
q         = q + sett.scal*(H\g);
datn.q    = q;
end
%==========================================================================