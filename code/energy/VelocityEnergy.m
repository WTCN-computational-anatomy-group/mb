function dat = VelocityEnergy(dat,sett)
set_bound;
v_settings = sett.v_settings;
if sett.threads>1 && numel(dat)>1
    % Memory  = 4*2*3*prod(GetSize(dat(1).v));
    % nthread = min(sett.threads,floor(sett.maxmem/Memory));
    parfor(n=1:numel(dat),sett.threads) % PARFOR
        set_bound;
        v           = GetData(dat(n).v);
        u0          = spm_diffeo('vel2mom', v, v_settings); % Initial momentum
        dat(n).E(2) = 0.5*sum(u0(:).*v(:));                 % Prior term
    end
else
    for n=1:numel(dat)
        v           = GetData(dat(n).v);
        u0          = spm_diffeo('vel2mom', v, v_settings); % Initial momentum
        dat(n).E(2) = 0.5*sum(u0(:).*v(:));                 % Prior term
    end
end
end
%==========================================================================