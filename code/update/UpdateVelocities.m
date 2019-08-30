function dat = UpdateVelocities(dat,mu,sett)
set_bound;
G  = spm_diffeo('grad',mu);
H0 = VelocityHessian(mu,G,sett.accel);
if sett.threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),sett.threads) % PARFOR
        set_bound;
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
function datn = UpdateVelocitiesSub(datn,mu,G,H0,sett)
v         = GetData(datn.v);
q         = datn.q;
Mn        = datn.Mat;
Mr        = spm_dexpm(q,sett.B);
Mat       = sett.Mmu\Mr*Mn;
d         = GetSize(datn.f);
psi       = compose(GetData(datn.psi), affine(d,Mat));
mu1       = Pull1(mu,psi);
[f,datn]  = GetClasses(datn,mu1,sett);
[a,w]     = Push1(f-softmax(mu1),psi,sett.d);
g         = reshape(sum(a.*G,4),[sett.d 3]);
H         = w.*H0;
u0        = spm_diffeo('vel2mom', v, sett.v_settings);                                  % Initial momentum
datn.E(2) = 0.5*sum(u0(:).*v(:));                                                       % Prior term
v         = v - sett.scal*spm_diffeo('fmg',H, g+u0, [sett.v_settings sett.s_settings]); % Gauss-Newton update

if sett.d(3)==1, v(:,:,:,3) = 0; end % If 2D
if sett.v_settings(1)==0             % Mean displacement should be 0
    avg = mean(mean(mean(v,1),2),3);
    v   = v - avg;
end
datn.v = SetData(datn.v,v);
end
%==========================================================================