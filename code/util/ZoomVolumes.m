function [dat,mu] = ZoomVolumes(dat,mu,sett,oMmu)
set_bound;
d     = sett.d;
d0    = [size(mu,1) size(mu,2) size(mu,3)];
Mmu   = sett.Mmu;
z     = single(reshape(d./d0,[1 1 1 3]));
Mzoom = oMmu\Mmu;
y     = reshape(reshape(identity(d),[prod(d),3])*Mzoom(1:3,1:3)'+Mzoom(1:3,4)',[d 3]);
mu    = spm_diffeo('pullc',mu,y);

if sett.threads>1 && numel(dat)>1
    % Should attempt to change number of threads accrding to how much memory
    % each one requires
    % Memory = 4*(prod(d)*(3+3)+prod(d0)*(3+3));
    parfor(n=1:numel(dat),sett.threads) % PARFOR
        v          = GetData(dat(n).v);
        v          = spm_diffeo('pullc',v.*z,y);
        dat(n).v   = ResizeFile(dat(n).v  ,d,Mmu);
        dat(n).psi = ResizeFile(dat(n).psi,d,Mmu);
        dat(n).v   = SetData(dat(n).v,v);
    end
else
    for n=1:numel(dat)
        v          = GetData(dat(n).v);
        v          = spm_diffeo('pullc',v,y).*z;
        dat(n).v   = ResizeFile(dat(n).v  ,d,Mmu);
        dat(n).psi = ResizeFile(dat(n).psi,d,Mmu);
        dat(n).v   = SetData(dat(n).v,v);
    end
end

end
%==========================================================================
