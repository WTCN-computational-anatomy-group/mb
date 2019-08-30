function [mu,dat] = UpdateSimpleMean(dat, mu, sett)
set_bound;
w  = zeros(sett.d,'single');
gf = zeros(size(mu),'single');
if sett.threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),sett.threads) % PARFOR
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
    H  = w.*AppearanceHessian(mu,sett.accel);
    g  = w.*softmax(mu) - gf;
    g  = g  + spm_field('vel2mom', mu, sett.mu_settings);
    mu = mu - spm_field(H, g, [sett.mu_settings sett.s_settings]);
end
end
%==========================================================================


%==========================================================================
function [g,w,datn] = UpdateSimpleMeanSub(datn,mu,sett)
d     = GetSize(datn.f);
q     = double(datn.q);
Mn    = datn.Mat;
psi   = compose(GetData(datn.psi), affine(d, sett.Mmu\spm_dexpm(q,sett.B)*Mn));
mu1   = Pull1(mu,psi);
[f,datn] = GetClasses(datn,mu1,sett);
[g,w] = Push1(f,psi,sett.d);
end
%==========================================================================