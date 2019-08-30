function [P,datn] = GetClasses(datn,mu1,~)
% Could be extended to include GMM stuff
f      = GetData(datn.f);
d      = [size(f) 1];
d      = d(1:3);
mog    = datn.mog;
K      = size(mu1,4);
%for k=1:K, subplot(3,3,k); imagesc(mu1(:,:,100,k)'); axis image xy off; end; drawnow
msk    = find(f~=0 & isfinite(f) & isfinite(mu1(:,:,:,1)));
mu     = reshape(mu1,[prod(d) K]);
mu     = mu(msk,:);
fm     = f(msk);
maxmu  = max(max(mu,[],2),0);
adjust = sum(log(sum(exp(mu-maxmu),2)+exp(-maxmu))+maxmu);
%adjust= log(sum(exp(mu),2)+1);
p      = zeros([numel(fm),K+1],'single');
for it=1:60
    for k=1:(K+1)
        p(:,k) = -((fm-mog.mu(k)).^2)/(2*mog.sig2(k)+eps) - 0.5*log(2*pi*mog.sig2(k)+eps);
        if k<=K, p(:,k) = p(:,k) + mu(:,k); end
    end
    pmx = max(p,[],2);
    p   = p - pmx;
    p   = exp(p);
    sp  = sum(p,2);
    prevE     = datn.E(1);
    datn.E(1) = -sum(log(sp)+pmx,1)+adjust; % doesn't account for priors (so can increase)
   %fprintf(' %g', datn.E(1));
    p   = p./sp;
    for k=1:(K+1)
        sp          = sum(p(:,k))+eps;
        mog.mu(k)   = sum(fm.*p(:,k))./sp;
        mog.sig2(k) = (sum((fm-mog.mu(k)).^2.*p(:,k))+10000*80^2)./(sp+10000); % ad hoc "Wishart priors"
    end
   %disp([mog.mu; mog.sig2])
    if it>1 && abs(prevE-datn.E(1))/size(p,1) < 1e-4, break; end
end
datn.mog = mog;

for k=1:(K+1)
    p(:,k) = -((fm-mog.mu(k)).^2)/(2*mog.sig2(k)+eps) - 0.5*log(2*pi*mog.sig2(k)+eps);
    if k<=K, p(:,k) = p(:,k) + mu(:,k); end
end
pmx = max(p,[],2);
p   = p - pmx;
p   = exp(p);
sp  = sum(p,2);
datn.E(1) = -sum(log(sp)+pmx,1)+adjust; % doesn't account for priors (so can increase)
%fprintf(' %g\n', datn.E(1));
p   = p./sp;
P   = zeros(size(mu1),'single')+NaN;
for k=1:K
    Pk         = P(:,:,:,k);
    Pk(msk)    = p(:,k);
    P(:,:,:,k) = Pk;
%   subplot(3,3,k); imagesc(Pk(:,:,100)'); axis image xy off; drawnow
end
end
%==========================================================================