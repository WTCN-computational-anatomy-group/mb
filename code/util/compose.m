function psi = compose(psi1,psi0)
set_bound;
psi = spm_diffeo('comp',psi1,psi0);
end
%==========================================================================