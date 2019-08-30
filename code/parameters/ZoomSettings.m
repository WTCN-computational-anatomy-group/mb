function sz = ZoomSettings(d, Mmu, v_settings, mu_settings, n)
[dz{1:n}] = deal(d);
sz        = struct('Mmu',Mmu,'d',dz,...
                   'v_settings', v_settings,...
                   'mu_settings',mu_settings);

% I'm still not entirely sure how best to deal with regularisation
% when dealing with different voxel sizes.
scale     = 1/abs(det(Mmu(1:3,1:3)));
for i=1:n
    sz(i).d           = ceil(d/(2^(i-1)));
    z                 = d./sz(i).d;
    sz(i).Mmu         = Mmu*[diag(z), (1-z(:))*0.5; 0 0 0 1];
    vx                = sqrt(sum(sz(i).Mmu(1:3,1:3).^2));
    sz(i).v_settings  = [vx v_settings *(scale*abs(det(sz(i).Mmu(1:3,1:3))))];
    sz(i).mu_settings = [vx mu_settings*(scale*abs(det(sz(i).Mmu(1:3,1:3))))];
end
end
%==========================================================================
