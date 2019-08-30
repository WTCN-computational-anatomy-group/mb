function dat = UpdateAffines(dat,mu,sett)
% Update the affine parameters
set_bound;
if ~isempty(sett.B)
    if sett.threads>1 && numel(dat)>1
        % Memory = loads
        parfor(n=1:numel(dat),sett.threads) % PARFOR
            set_bound;
            dat(n) = UpdateAffinesSub(dat(n),mu,sett);
        end
    else
        for n=1:numel(dat)
            dat(n) = UpdateAffinesSub(dat(n),mu,sett);
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
function datn = UpdateAffinesSub(datn,mu,sett)
% This could be made more efficient.
d    = GetSize(datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
[Mr,dM3] = spm_dexpm(q,sett.B);
dM   = zeros(12,size(sett.B,3));
for m=1:size(sett.B,3)
    tmp     = sett.Mmu\dM3(:,:,m)*Mn;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi1 = GetData(datn.psi);
psi0 = affine(d,sett.Mmu\Mr*Mn);
J    = spm_diffeo('jacobian',psi1);
J    = reshape(Pull1(reshape(J,[sett.d 3*3]),psi0),[d 3 3]);
psi  = compose(psi1,psi0);
clear psi0 psi1
mu1  = Pull1(mu,psi);
[f,datn] = GetClasses(datn,mu1,sett);
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
a         = mask(f - softmax(mu1),msk);
[H,g]     = AffineHessian(mu1,G,a,single(msk),sett.accel);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-6 + 0.1);
q         = q + sett.scal*(H\g);
datn.q    = q;
end
%==========================================================================