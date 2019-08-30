function [mu,dat] = UpdateMean(dat, mu, sett)
set_bound;
g  = spm_field('vel2mom', mu, sett.mu_settings);
M  = size(mu,4);
H  = zeros([sett.d M*(M+1)/2],'single');
H0 = AppearanceHessian(mu,sett.accel);
if sett.threads>1 && numel(dat)>1
    % Memory = 4*prod(sett.d)*2*(M*(M+1)/2 + M) + max(prod(sett.d)*3+3*prod(dn)*2,(3+4*M)*prod(dn))); % Needs more work
    parfor(n=1:numel(dat),sett.threads) % PARFOR
        set_bound;
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
mu = mu - spm_field(H, g, [sett.mu_settings sett.s_settings]);
end
%==========================================================================

%==========================================================================
function [g,H,datn] = UpdateMeanSub(datn,mu,H0,sett)
d    = GetSize(datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
psi  = compose(GetData(datn.psi), affine(d, sett.Mmu\spm_dexpm(q,sett.B)*Mn));
mu1  = Pull1(mu,psi);
[f,datn] = GetClasses(datn,mu1,sett);
if isempty(H0)
    g     = Push1(softmax(mu1)-f,psi,sett.d);
    H     = Push1(AppearanceHessian(mu1,sett.accel),psi,sett.d);
else
    % Faster approximation - but might be unstable
    % If there are problems, then revert to the slow
    % way.
    [g,w] = Push1(softmax(mu1)-f,psi,sett.d);
    H     = w.*H0;
end
end
%==========================================================================