function dat = SaveImages(dat,mu,sett)
for n=1:numel(dat)
    dat(n).psi = ResavePsiSub(dat(n),sett);
end

if ~isempty(mu)
    fa       = file_array('logTPM.nii',size(mu),'float32',0);
    Nmu      = nifti;
    Nmu.dat  = fa;
    Nmu.mat  = sett.Mmu;
    Nmu.mat0 = sett.Mmu;
    Nmu.descrip = 'Mean parameters (log)';
    create(Nmu);
    Nmu.dat(:,:,:,:) = mu;
    
    fa       = file_array('softmaxTPM.nii',size(mu),'float32',0);
    Nmu      = nifti;
    Nmu.dat  = fa;
    Nmu.mat  = sett.Mmu;
    Nmu.mat0 = sett.Mmu;
    Nmu.descrip = 'Mean parameters (softmax)';
    create(Nmu);
    Nmu.dat(:,:,:,:) = softmax(mu);
end
end
%==========================================================================