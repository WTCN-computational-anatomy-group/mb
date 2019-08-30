function E = TemplateEnergy(mu,sett)
set_bound;
g    = spm_field('vel2mom', mu, sett.mu_settings);
E    = 0.5*mu(:)'*g(:);
end
%==========================================================================