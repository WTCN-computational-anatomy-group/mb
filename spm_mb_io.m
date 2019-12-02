function varargout = spm_mb_io(varargin)
%__________________________________________________________________________
%
% Functions for I/O related.
%
% FORMAT to         = spm_mb_io('CopyFields',from,to)
% FORMAT [P,datn]   = spm_mb_io('GetClasses',datn,mu,sett)
% FORMAT out        = spm_mb_io('GetData',in)
% FORMAT Mat        = spm_mb_io('GetMat',fin)
% FORMAT [d,M]      = spm_mb_io('GetSize',fin)
% FORMAT dat        = spm_mb_io('InitDat',F,sett)
% FORMAT [psi,fpth] = spm_mb_io('SavePsiSub',datn,sett) 
% FORMAT dat        = spm_mb_io('SaveTemplate',dat,mu,sett)
% FORMAT              spm_mb_io('SetBoundCond')
% FORMAT fout       = spm_mb_io('SetData',fin,f) 
% FORMAT              spm_mb_io('SetPath')
% FORMAT              spm_mb_io('WriteNii',f,img,Mmu,descrip);
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_io
    error('Not enough argument. Type ''help spm_mb_io'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'CopyFields'
        [varargout{1:nargout}] = CopyFields(varargin{:});       
    case 'GetClasses'
        [varargout{1:nargout}] = GetClasses(varargin{:});    
    case 'GetData'
        [varargout{1:nargout}] = GetData(varargin{:});            
    case 'GetMat'
        [varargout{1:nargout}] = GetMat(varargin{:});             
    case 'GetSize'
        [varargout{1:nargout}] = GetSize(varargin{:});    
    case 'InitDat' 
        [varargout{1:nargout}] = InitDat(varargin{:});            
    case 'SavePsiSub' 
        [varargout{1:nargout}] = SavePsiSub(varargin{:});           
    case 'SaveTemplate'
        [varargout{1:nargout}] = SaveTemplate(varargin{:});    
    case 'SetBoundCond'
        [varargout{1:nargout}] = SetBoundCond(varargin{:});        
    case 'SetData'
        [varargout{1:nargout}] = SetData(varargin{:});                
    case 'SetPath'
        [varargout{1:nargout}] = SetPath(varargin{:});        
    case 'WriteNii'
        [varargout{1:nargout}] = WriteNii(varargin{:});            
    otherwise
        help spm_mb_io
        error('Unknown function %s. Type ''help spm_mb_io'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% CopyFields()
function to = CopyFields(from,to)
fn = fieldnames(from);
for i=1:numel(fn)
    to.(fn{i}) = from.(fn{i});
end
end
%==========================================================================

%==========================================================================
% GetClasses()
function [P,datn] = GetClasses(datn,mu,sett)

if ~isfield(datn,'mog')
    P = GetData(datn.f);

    if nargout > 1
        % Make mask
        msk = sum(P,4) > 0; % Removes voxels that sums to zero across classes..
        msk = msk & sum(~isfinite(P),4) == 0; % ..and voxels that are not finite in segmentation..
        msk = msk & sum(~isfinite(mu),4) == 0; % ..and voxels that are not finite in template

        % Compute subject-specific categorical cross-entropy loss between
        % segmentation and template
        tmp       = sum(P.*mu,4) - spm_mb_shape('LSE',mu,4);  
        datn.E(1) = -sum(tmp(msk));
    end
else
    [P,datn] = spm_mb_appearance('Update',datn,mu,sett);
end

if 0
    % Show stuff
    d  = size(mu);
    d  = [d 1 1];
    K  = d(4);        
    nr = floor(sqrt(K));
    nc = ceil(K/nr);  
        
    for k=1:K    
        % Show template    
        figure(664); subplot(nr,nc,k); imagesc(mu(:,:,ceil(d(3).*0.55),k)'); 
        title('mu');
        axis image xy off; drawnow
        
        % Show segmentation
        figure(665); subplot(nr,nc,k); imagesc(P(:,:,ceil(d(3).*0.55),k)'); 
        title('Seg')
        axis image xy off; drawnow
    end
end
end
%==========================================================================

%==========================================================================
% GetData()
function out = GetData(in)
if isnumeric(in)
    out = single(in);
    return
end
if isa(in,'char')
    in = nifti(in);
end
if isa(in,'nifti')
    M = numel(in);
    d = size(in(1).dat,[1 2 3 4 5]);
    if M>1
        d(4) = M;
        out = zeros(d,'single');
        for m=1:M
            out(:,:,:,m) = single(in(m).dat(:,:,:,:,:));
        end
    else
        out = single(in.dat(:,:,:,:,:));
        if numel(d)>4 && d(4)==1
            out = reshape(out,[d(1:3) d(5)]);
        end
    end
    return
end
error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
% GetMat()
function Mat = GetMat(fin)
if isnumeric(fin)
    Mat = eye(4);
    return;
end
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    Mat  = fin(1).mat;
    return
end
error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
% GetSize()
function [d,M] = GetSize(fin)
d = [GetDimensions(fin) 1 1 1];
M = d(4);
d = d(1:3);
end
%==========================================================================

%==========================================================================
% InitDat()
function dat = InitDat(F,sett)

% Parse function settings
run2d = sett.gen.run2d;

M0 = eye(4);
for n=1:numel(F)
    
    % Init datn.f
    if iscell(F(n)) && isnumeric(F{n})
        % Input F is numeric -> store as numeric
        
        if run2d            
            % Get 2D slice from 3D data
            dat(n).f = GetSlice(F{n},run2d);
        else
            dat(n).f = single(F{n});
        end
    elseif isa(F(n),'nifti') || (iscell(F(n)) && (isa(F{n},'char') || isa(F{n},'nifti')))
        % Input F is nifti (path or object) -> store as nifti
                       
        if isa(F(n),'nifti')
            Nii      = F(n);
            dat(n).f = Nii;        
        elseif iscell(F(n)) 
            if isa(F{n},'char')
                Nii      = nifti(F{n});
                dat(n).f = Nii;        
            elseif isa(F{n},'nifti')                
                Nii      = F{n};
                C        = numel(Nii);
                dat(n).f = nifti;
                for c=1:C
                    dat(n).f(c) = Nii(c);
                end
            end
        end
        
        if run2d
            % Get 2D slice from 3D data
            fn       = spm_mb_io('GetData',dat(n).f);
            dat(n).f = GetSlice(fn,run2d);
        end
    end
    
    % Other parameters
    dat(n).M   = M0;    
    dat(n).q   = zeros(6,1);    
    dat(n).v   = [];    
    dat(n).psi = [];    
    dat(n).E   = [0 0 0]; % Px Pv Pbf
    dat(n).bf  = [];
                     
    % Orientation matrix (image voxel-to-world)    
    dat(n).Mat = eye(4); % Should really do this better           
    if isa(F(n),'nifti') || (iscell(F(n)) && (isa(F{n},'char') || isa(F{n},'nifti')))
        Mat = Nii(1).mat;
        vx  = sqrt(sum(Mat(1:3,1:3).^2));    
        if run2d
            dat(n).Mat = diag([vx 1]);
        else
            dat(n).Mat = Nii(1).mat;        
        end
    end
end
end
%==========================================================================

%==========================================================================
% SavePsiSub()
function [psi,fpth] = SavePsiSub(datn,sett)

% Parse function settings
B       = sett.registr.B;
dir_res = sett.write.dir_res;
Mmu     = sett.var.Mmu;

d    = GetSize(datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
psi1 = GetData(datn.psi);
psi  = spm_mb_shape('Compose',psi1,spm_mb_shape('Affine',d,Mmu\spm_dexpm(q,B)*Mn));
psi  = reshape(reshape(psi,[prod(d) 3])*Mmu(1:3,1:3)' + Mmu(1:3,4)',[d 1 3]);
fpth = '';
if isa(datn.psi(1),'nifti')
    to_delete = datn.psi(1).dat.fname;
    [~,nam,~] = fileparts(datn.f(1).dat.fname);
    fpth = fullfile(dir_res,['y_' nam '.nii']);
    datn.psi(1).dat.fname = fpth;
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

%==========================================================================
% SaveTemplate()
function dat = SaveTemplate(dat,mu,sett)

% Parse function settings
dir_res = sett.write.dir_res;
Mmu     = sett.var.Mmu;

if ~isempty(mu)
    % Save mu (log)
    f        = fullfile(dir_res ,'mu_log.nii');
    fa       = file_array(f,size(mu),'float32',0);
    Nmu      = nifti;
    Nmu.dat  = fa;
    Nmu.mat  = Mmu;
    Nmu.mat0 = Mmu;
    Nmu.descrip = 'Mean parameters (log)';
    create(Nmu);
    Nmu.dat(:,:,:,:) = mu;
    
    % Save mu (softmax)    
    d        = size(mu);
    d        = [d 1];
    mu       = cat(4,mu,zeros(d(1:3),'single'));
    mu       = spm_mb_shape('Softmax',mu,4);
    f        = fullfile(dir_res ,'mu_softmax.nii');        
    fa       = file_array(f,size(mu),'float32',0);
    Nmu      = nifti;
    Nmu.dat  = fa;
    Nmu.mat  = Mmu;
    Nmu.mat0 = Mmu;
    Nmu.descrip = 'Mean parameters (softmax)';
    create(Nmu);
    Nmu.dat(:,:,:,:) = mu;
end
end
%==========================================================================

%==========================================================================
% SetBoundCond()
function SetBoundCond
spm_diffeo('bound',0);
spm_field('bound',1);
end
%==========================================================================

%==========================================================================
% SetData()
function fout = SetData(fin,f)
fout = fin;
if isnumeric(fin)
    fout = f;
    return
end
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    M    = numel(fin);
    if M>1
        for m=1:M
            fout(m).dat(:,:,:,1,:) = f(:,:,:,m,:);
        end
    else
        fout(1).dat(:,:,:,:,:) = reshape(f,size(fout(1).dat));
    end
    return
end
end
%==========================================================================

%==========================================================================
% SetPath()
function SetPath
pth=fileparts(which('spm'));
addpath(pth);
addpath(fullfile(pth,'toolbox','Longitudinal'));
addpath(fullfile(pth,'toolbox','Shoot'));
end
%==========================================================================

%==========================================================================
% WriteNii()
function WriteNii(f,img,M,descrip)
fa       = file_array(f,size(img),'float32',0);
Nii      = nifti;
Nii.dat  = fa;
Nii.mat  = M;
Nii.mat0 = M;
Nii.descrip = descrip;
create(Nii);
d = size(img);
Nii.dat(:,:,:,:,:,:) = img;
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% GetDimensions()
function d = GetDimensions(fin)
if isnumeric(fin)
    d = size(fin);
    d = [d 1 1];
    return
end
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    M    = numel(fin);
    d    = size(fin(1).dat,[1 2 3 4 5]);
    if M>1
        d(4) = M;
    else
        if numel(d)>4 && d(4)==1
            d = [d(1:3) d(5)];
        end
    end
    return
end
end
%==========================================================================

%==========================================================================
% GetSlice()
function fn = GetSlice(fn,direction)
d  = size(fn);
d  = [d 1];
ix = round(d(1:3)*0.5);

if d(3) == 1, return; end

if direction == 1
    fn = single(fn(ix(1),:,:,:));
elseif direction == 2
    fn = single(fn(:,ix(2),:,:));
elseif direction == 3
    fn = single(fn(:,:,ix(3),:));
end

% Reshape
C  = d(4);
ix = 1:3;
d  = d(1:3);
fn = reshape(fn, [d(ix ~= direction) 1 C]);
end
%==========================================================================   