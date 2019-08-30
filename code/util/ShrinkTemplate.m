function mu1 = ShrinkTemplate(mu,oMmu,sett)
d       = [size(mu,1) size(mu,2) size(mu,3)];
Mmu     = sett.Mmu;
Mzoom   = Mmu\oMmu;
if norm(Mzoom-eye(4))<1e-4 && all(d==sett.d)
    mu1 = mu;
else
    y       = reshape(reshape(identity(d),[prod(d),3])*Mzoom(1:3,1:3)'+Mzoom(1:3,4)',[d 3]);
    [mu1,c] = Push1(mu,y,sett.d);
    mu1     = mu1./(c+eps);
end
end
%==========================================================================
