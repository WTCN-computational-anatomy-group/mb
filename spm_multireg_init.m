function varargout = spm_multireg_init(varargin)
%__________________________________________________________________________
%
% Initialisation functions for spm_multireg.
%
% FORMAT dat = spm_multireg_init('InitDat',F,K,sett)
% FORMAT dat = spm_multireg_init('InitDef',dat,sett)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_init
    error('Not enough argument. Type ''help spm_multireg_init'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'InitDat'
        [varargout{1:nargout}] = InitDat(varargin{:});   
    case 'InitDef'
        [varargout{1:nargout}] = InitDef(varargin{:});           
    otherwise
        help spm_multireg_init
        error('Unknown function %s. Type ''help spm_multireg_init'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% InitDat()
function dat = InitDat(F,K,sett)

if sett.do.gmm
    % Need a better way of initialising these parameters
    K1    = K + 1;
    mog   = struct('mu',(1:-(1/K1):(1/K1))*2000,'sig2',ones(1,K1)/3*sqrt(1000)); % Random
    %mog  = struct('mu',ones(1,K1)*500,'sig2',ones(1,K1)*500^2); % Same (for existing TPM)
end

M0 = eye(4);
if sett.do.gmm
    dat = struct('f',F,'M',M0, 'q',zeros(6,1), 'v',[], 'psi',[], 'E',[0 0],'mog',mog);
else
    dat = struct('f',F,'M',M0, 'q',zeros(6,1), 'v',[], 'psi',[], 'E',[0 0]);
end

for n=1:numel(dat)
    if isnumeric(dat(n).f)
        dat(n).Mat = eye(4); % Should really do this better
    else
        if isa(dat(n).f,'char')
            dat(n).f   = nifti(dat(n).f);
        end
        if isa(dat(n).f,'nifti')
            dat(n).Mat = dat(n).f(1).mat;
        end
    end
end
end
%==========================================================================

%==========================================================================
% InitDef()
function dat = InitDef(dat,sett)
v    = zeros([sett.var.d,3],'single');
psi1 = spm_multireg_util('Identity',sett.var.d);
for n=1:numel(dat)
    dat(n).q = zeros(size(sett.registr.B,3),1);
    if isnumeric(dat(n).f)
        dat(n).v   = v;
        dat(n).psi = psi1;
    else
        if isa(dat(n).f,'nifti')
            [~,nam,~] = fileparts(dat(n).f(1).dat.fname);
            vname    = fullfile('.',['v_' nam '.nii']);
            pname    = fullfile('.',['psi_' nam '.nii']);
            fa       = file_array(vname,[sett.var.d(1:3) 1 3],'float32',0);
            nii      = nifti;
            nii.dat  = fa;
            nii.mat  = sett.var.Mmu;
            nii.mat0 = sett.var.Mmu;
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