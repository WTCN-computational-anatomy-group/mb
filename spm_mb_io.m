function varargout = spm_mb_io(varargin)
%__________________________________________________________________________
%
% Functions for I/O related.
%
% FORMAT to         = spm_mb_io('CopyFields',from,to)
% FORMAT [P,datn]   = spm_mb_io('GetClasses',datn,mu,sett)
% FORMAT [out,M]    = spm_mb_io('GetData',in)
% FORMAT Mat        = spm_mb_io('GetMat',fin)
% FORMAT [d,M]      = spm_mb_io('GetSize',fin)
% FORMAT s          = spm_mb_io('GetScale',fin,sett);
% FORMAT dat        = spm_mb_io('InitDat',data,sett)
% FORMAT model      = spm_mb_io('MakeModel',dat,model,sett)
% FORMAT [psi,fpth] = spm_mb_io('SavePsiSub',datn,sett) 
% FORMAT              spm_mb_io('SaveTemplate',dat,mu,sett)
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
    case 'GetScale'
        [varargout{1:nargout}] = GetScale(varargin{:});
    case 'InitDat' 
        [varargout{1:nargout}] = InitDat(varargin{:});
    case 'MakeModel' 
        [varargout{1:nargout}] = MakeModel(varargin{:});        
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
        % Compute subject-specific categorical cross-entropy loss between
        % segmentation and template
        msk       = all(isfinite(P),4);
        tmp       = sum(P.*mu,4) - spm_mb_shape('LSE',mu,4);
        datn.E(1) = -sum(tmp(msk(:)));
    end
else
    % Update appearance model
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
function out = GetScale(in,sett)
% Return a scale for adding random numbers

% Parse function settings
run2d = sett.gen.run2d;

if isnumeric(in)
    if isa(in,'integer') || run2d
        out = ones([1,1,1,size(in,4)]);
    else
        out = zeros([1,1,1,size(in,4)]);
    end
    return;
end
if isa(in,'char')
    in = nifti(in);
end
if isa(in,'nifti')
    C   = numel(in);
    d   = size(in(1).dat,[1 2 3 4 5]);
    if d(4)>1 && C>1, error('Don''t know what to do with this image data.'); end
    d(4) = max(d(4),C);
    out  = zeros([1,1,1,d(4)]);
    if C>1
        for c=1:C
            dt1 = in(c).dat.dtype(1);
            if dt1=='I' || dt1=='U'
                out(c) = in(c).dat.scl_slope(1);
            end
        end
    else
        dt1 = in(1).dat.dtype(1);
        if dt1=='I' || dt1=='U'
            out(:) = in(1).dat.scl_slope(1);
        end
    end
    return
end
error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
% GetData()
function [out,Mn] = GetData(in)
Mn = eye(4);
if isnumeric(in)
    out = single(in);
    return
end
if isa(in,'char')
    in = nifti(in);
end
if isa(in,'nifti')
    C  = numel(in);
    d  = size(in(1).dat,[1 2 3 4 5]);
    Mn = in(1).mat;
    if C>1
        d(4) = C;
        out = zeros(d,'single');
        for m=1:C
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
function dat = InitDat(data,sett)

% Parse function settings
do_gmm = sett.do.gmm;
run2d  = sett.gen.run2d;

% Initialise for each subject
N  = numel(data);
M0 = eye(4);
for n=1:N

    if isstruct(data(n)) && isfield(data(n),'F'), F = data(n).F;
    else,                                         F = data(n);
    end

    % Init datn.f
    if iscell(F) && isnumeric(F{1})
        % Input F is numeric -> store as numeric
        
        if run2d            
            % Get 2D slice from 3D data
            dat(n).f = GetSlice(F{1},run2d);
        else
            dat(n).f = single(F{1});
        end
    elseif isa(F,'nifti') || (iscell(F) && (isa(F{1},'char') || isa(F{1},'nifti')))
        % Input F is nifti (path or object) -> store as nifti
                       
        if isa(F,'nifti')
            Nii      = F;
            dat(n).f = Nii;        
        elseif iscell(F) 
            if isa(F{1},'char')
                Nii      = nifti(F{1});
                dat(n).f = Nii;        
            elseif isa(F{1},'nifti')                
                Nii      = F{1};
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
    
    % Get number of channels
    [~,C] = spm_mb_io('GetSize',dat(n).f);
    
    % Parameters
    dat(n).M     = M0;    
    dat(n).q     = zeros(6,1);    
    dat(n).v     = [];    
    dat(n).psi   = [];    
    dat(n).E     = [0 0 0]; % Px Pv Pbf
    dat(n).bf    = [];   
       
    if do_gmm
        dat(n).mog = [];
        dat(n).bf  = [];                
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subject-level 'extras'
    %------------------------------------------------------------
    
    % Do bias field (per channel)
    if isstruct(data(n)) && isfield(data(n),'do_bf') && ~isempty(data(n).do_bf)
        dat(n).do_bf = data(n).do_bf;
    else
        dat(n).do_bf = true;        
    end
    if numel(dat(n).do_bf) < C 
        dat(n).do_bf = repmat(dat(n).do_bf,[1 C]); 
    end
    
    % Intensity prior index
    if isstruct(data(n)) && isfield(data(n),'ix_pop') && ~isempty(data(n).ix_pop)
        dat(n).ix_pop = data(n).ix_pop;
    else
        dat(n).ix_pop = 1;
    end
    
    % Is CT data
    if isstruct(data(n)) && isfield(data(n),'is_ct') && ~isempty(data(n).is_ct)
        dat(n).is_ct = data(n).is_ct;
    else
        dat(n).is_ct = false;
    end
    if numel(dat(n).is_ct) < C 
        dat(n).is_ct = repmat(dat(n).is_ct,[1 C]); 
    end
    
    % Labels in a cell array as {nifti,cm_map}
    if isstruct(data(n)) && isfield(data(n),'labels') && (~isempty(data(n).labels{1}) && ~isempty(data(n).labels{2}))
        dat(n).labels = data(n).labels;
        
        if run2d
            % Get 2D slice from 3D data
            labels           = spm_mb_io('GetData',dat(n).labels{1});
            dat(n).labels{1} = GetSlice(labels,run2d);
        end
    else
        dat(n).labels = [];        
    end
            
    % Orientation matrix (image voxel-to-world)    
    dat(n).Mat = eye(4);
    if isa(F,'nifti') || (iscell(F) && (isa(F{1},'char') || isa(F{1},'nifti')))
        dat(n).Mat = Nii(1).mat;        
        if run2d
            vx         = sqrt(sum(dat(n).Mat(1:3,1:3).^2));
            dat(n).Mat = [diag(vx) zeros(3,1); 0 0 0 1];
        end
    end
end
end
%==========================================================================

%==========================================================================
% MakeModel()
function model = MakeModel(dat,model,sett)

% Parse function settings
do_gmm           = sett.do.gmm;
do_updt_int      = sett.do.updt_int;
do_updt_template = sett.do.updt_template;
dir_res          = sett.write.dir_res;
write_model      = sett.write.model;

if isempty(dir_res) 
    pth     = fileparts(dat(1).f(1).dat.fname);
    dir_res = pth; 
end

if do_updt_template
    % Shape related
    f                    = fullfile(dir_res ,'mu_log_spm_mb.nii');
    model.shape.template = f;
end

if do_updt_int && do_gmm
    % Appearance related
    p_ix = spm_mb_appearance('GetPopulationIdx',dat);
    Npop = numel(p_ix);

    model.appear = containers.Map;
    for p=1:Npop
        n      = p_ix{p}(1);
        datn   = dat(n);
        ix_pop = datn.ix_pop;
        pr     = datn.mog.pr;
        model.appear(num2str(ix_pop)) = pr;
    end
end

if write_model && (do_updt_template || do_updt_int)
    % Save model
    save(fullfile(dir_res,'model_spm_mb.mat'),'model')
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
    to_delete   = datn.psi(1).dat.fname;
    [pth,nam,~] = fileparts(datn.f(1).dat.fname);
    if isempty(dir_res), dir_res = pth; end
    fpth        = fullfile(dir_res,['y_' nam '.nii']);
    
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
function SaveTemplate(dat,mu,sett)

% Parse function settings
dir_res          = sett.write.dir_res;
do_updt_template = sett.do.updt_template;
Mmu              = sett.var.Mmu;
write_mu         = sett.write.mu;

if isempty(dir_res) 
    pth     = fileparts(dat(1).f(1).dat.fname);
    dir_res = pth; 
end

if ~isempty(mu) && do_updt_template
    if write_mu(1)
        % Log
        f        = fullfile(dir_res ,'mu_log_spm_mb.nii');
        fa       = file_array(f,size(mu),'float32',0);
        Nmu      = nifti;
        Nmu.dat  = fa;
        Nmu.mat  = Mmu;
        Nmu.mat0 = Mmu;
        Nmu.descrip = 'Mean parameters (log)';
        create(Nmu);
        Nmu.dat(:,:,:,:) = mu;
    end
    
    if write_mu(2)
        % Softmax   
        mu = spm_mb_shape('TemplateK1',mu,4);
        mu = exp(mu);    

        f        = fullfile(dir_res ,'mu_softmax_spm_mb.nii');        
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

nslices = 0; % 1 + 2*nslices

if direction == 1
    fn = single(fn(ix(1) - nslices:ix(1) + nslices,:,:,:));
elseif direction == 2
    fn = single(fn(:,ix(2) - nslices:ix(2) + nslices,:,:));
elseif direction == 3
    fn = single(fn(:,:,ix(3) - nslices:ix(3) + nslices,:));
end

% Reshape
C  = d(4);
ix = 1:3;
d  = d(1:3);
fn = reshape(fn, [d(ix ~= direction) 1 + 2*nslices C]);
end
%==========================================================================   

