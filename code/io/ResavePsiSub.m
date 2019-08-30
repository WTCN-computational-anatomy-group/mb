function psi = ResavePsiSub(datn,sett)
d    = GetSize(datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
psi1 = GetData(datn.psi);
psi  = compose(psi1,affine(d,sett.Mmu\spm_dexpm(q,sett.B)*Mn));
psi  = reshape(reshape(psi,[prod(d) 3])*sett.Mmu(1:3,1:3)' + sett.Mmu(1:3,4)',[d 1 3]);
if isa(datn.psi(1),'nifti')
    to_delete = datn.psi(1).dat.fname;
    [~,nam,~] = fileparts(datn.f(1).dat.fname);
    datn.psi(1).dat.fname = fullfile('.',['y_' nam '.nii']);
    datn.psi(1).dat.dim = [d 1 3];
    datn.psi(1).mat = datn.f(1).mat0; % For working with "imported" images;
   %datn.psi(1).mat = Mn;
    datn.psi(1).descrip = 'Deformation';
    create(datn.psi(1));
    datn.psi(1).dat(:,:,:,:,:) = psi;
    delete(to_delete);
end
end
%==========================================================================