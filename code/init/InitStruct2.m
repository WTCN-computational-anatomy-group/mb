function dat = InitStruct2(dat,sett)
v    = zeros([sett.d,3],'single');
psi1 = identity(sett.d);
for n=1:numel(dat)
    dat(n).q = zeros(size(sett.B,3),1);
    if isnumeric(dat(n).f)
        dat(n).v   = v;
        dat(n).psi = psi1;
    else
        if isa(dat(n).f,'nifti')
            [~,nam,~] = fileparts(dat(n).f(1).dat.fname);
            vname    = fullfile('.',['v_' nam '.nii']);
            pname    = fullfile('.',['psi_' nam '.nii']);
            fa       = file_array(vname,[sett.d(1:3) 1 3],'float32',0);
            nii      = nifti;
            nii.dat  = fa;
            nii.mat  = sett.Mmu;
            nii.mat0 = sett.Mmu;
            nii.descrip = 'Velocity';
            create(nii);
            nii.dat(:,:,:,:) = v;
            dat(n).v    = nii;

            nii.dat.fname = pname;
            nii.descrip = 'Deformation (WIP)';
            create(nii);
            nii.dat(:,:,:,:) = psi1;
            dat(n).psi  = nii;
        end
    end
end
end
%==========================================================================