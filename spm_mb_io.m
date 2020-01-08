function varargout = spm_mb_io(varargin)
%__________________________________________________________________________
%
% Functions for I/O related.
%
% FORMAT to         = spm_mb_io('CopyFields',from,to)
% FORMAT pth        = spm_mb_io('CropLearnedTemplate',pth,centre,bb,do_crop)
% FORMAT [P,datn]   = spm_mb_io('GetClasses',datn,mu,sett)
% FORMAT [out,M]    = spm_mb_io('GetData',in,[sl,idx])
% FORMAT [out,M]    = spm_mb_io('MemMapData',in)
% FORMAT Mat        = spm_mb_io('GetMat',fin)
% FORMAT [d,K,L]    = spm_mb_io('GetSize',fin)
% FORMAT s          = spm_mb_io('GetScale',fin,sett);
% FORMAT dat        = spm_mb_io('InitDat',data,sett)
% FORMAT model      = spm_mb_io('MakeModel',dat,model,sett)
% FORMAT [psi,fpth] = spm_mb_io('SavePsiSub',datn,sett) 
% FORMAT              spm_mb_io('SaveSubspace',U,sett)
% FORMAT              spm_mb_io('SaveTemplate',mu,sett)
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
    case 'CropLearnedTemplate'
        [varargout{1:nargout}] = CropLearnedTemplate(varargin{:});        
    case 'GetClasses'
        [varargout{1:nargout}] = GetClasses(varargin{:});
    case 'GetData'
        [varargout{1:nargout}] = GetData(varargin{:});
    case 'MemMapData'
        [varargout{1:nargout}] = MemMapData(varargin{:});
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
    case 'SaveSubspace'
        [varargout{1:nargout}] = SaveSubspace(varargin{:});
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
% Copy fields from one structure to another
fn = fieldnames(from);
for i=1:numel(fn)
    to.(fn{i}) = from.(fn{i});
end
end
%==========================================================================

%==========================================================================
% CropLearnedTemplate()
function fout = CropLearnedTemplate(fin,c,mrg,do)
if nargin < 2, c   = []; end
if nargin < 3, mrg = 0; end
if nargin < 4, do  = true; end

if isscalar(mrg), mrg = mrg*ones(1,3); end

if do
    [pth,nam,ext] = fileparts(fin);
    fout          = fullfile(pth,['cpy_' nam ext]);
    copyfile(fin,fout);
    
    c = [192 189 201];
    bb = [c - [100 120 180]; c + [100 140 110]]; % [l b d], [r f u]
    
    V  = spm_vol(fout);
    for k=1:numel(V)
        VO = SubVol(V(k),bb);
    end
else
%     Nii       = nifti(pth);
%     dim       = Nii.dat.dim;    
%     tis_class = 1;
%     centre_ix = [1,2; 3,2; 1,3];
%     
%     figure(666)
%     for i=1:3
%         ci       = c(centre_ix(i,:));
%         slice_ix = c(round(setdiff(1:3,centre_ix(i,:))));
%         
%         subplot(1,3,i)
%         if     i == 1, imagesc(Nii.dat(:,:,slice_ix,tis_class)); 
%         elseif i == 2, imagesc(squeeze(Nii.dat(:,slice_ix,:,tis_class))); 
%         else,          imagesc(squeeze(Nii.dat(slice_ix,:,:,tis_class))); 
%         end
%         axis image xy               
%         
%         hold on
%         plot(ci(1),ci(2),'rx'); 
%         bb1 = [ci(1) + mrg(i), ci(1) - mrg(i), ci(1) - mrg(i), ci(1) + mrg(i), ci(1) + mrg(i)];
%         bb2 = [ci(2) + mrg(i), ci(2) + mrg(i), ci(2) - mrg(i), ci(2) - mrg(i), ci(2) + mrg(i)];
%         plot(bb1, bb2, 'r-', 'LineWidth', 1);
%         hold off
%     end
end
end
%==========================================================================
% GetClasses()
function [P,datn] = GetClasses(datn,mu,sett)
% Update GMM parameters of an image and return responsibilities.
% If the image contains pre-computed responsibilities, just return it.
%
% FORMAT [P,datn] = GetClasses(datn,mu,sett)
% datn - Data structure of one subject
% mu   - Log-template
% sett - Option structure
% P    - Responsibilities [dx dy dz K]
%
% TODO (YB): I find it strange that this is in io even though it implies an
% update of the GMM parameters. Having it if appearance would make more
% sense. (I understand 'io' as a rather static operation)

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
%
% FORMAT out = GetScale(in,sett)
% in   - array or file_array
% sett - settings
% out  - scale for each channel

% Parse function settings
run2d = sett.gen.run2d;

if isnumeric(in)
    if isa(in,'integer') || run2d
        out = ones([1 1 1 size(in,4)]);
    else
        out = zeros([1 1 1 size(in,4)]);
    end
    return;
end
if isa(in,'file_array')
    d = size(in);
    d((end+1):4) = 1;
    out = zeros([1 1 1 d(4)]);
    dtypes = {in.dtype};
    slopes = {in.scl_slope};
    dims   = {in.dim};
    j = 1;
    for i=1:numel(dtypes)
        dt = dtypes{i}(1);
        C  = dims{i};
        C((end+1):4) = 1;
        C = C(4);
        if dt(1)=='I' || dt(1)=='U'
            out(j:(j+C-1)) = slopes{i};
        end
        j = j + C;
    end
    return
end
error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
% GetData()
function [out,Mn] = GetData(in, slice, idx)
% FORMAT [out,Mn] = GetData(in, [dim, idx])
% in    - array or file_array or nifti or filename
% slice - Dimension along which to slice [none]
% idx   - Index to read in slice [all]

if nargin < 2, slice = 0; end
Mn = eye(4);
if isnumeric(in) || isa(in, 'file_array')
    if slice > 0
        [subs{1:numel(size(in))}] = deal(':');
        subs{slice} = idx;
        subs = struct('type', '()', 'subs', {subs});
        out = single(subsref(in, subs));
    else
        out = single(in());
    end
    return
end
if isa(in,'char')
    in = nifti(in);
end
if isa(in,'nifti')
    C  = numel(in);
    d  = size(in(1).dat);
    Mn = in(1).mat;
    if C>1
        d(end+1:4) = 1;
        d(4) = C;
        if slice > 0
            d(end+1:slice) = 1;
            d(slice) = 1;
        end
        out = zeros(d,'single');
        [subsin{1:numel(d)}] = deal(':');
        subsout = subsin;
        if slice > 0
            subsin{slice} = idx;
        end
        subsin  = struct('type', '()', 'subs', {subsin});
        subsout = struct('type', '()', 'subs', {subsout});
        for m=1:C
            subsout.subs{4} = m;
            out = subsasgn(out, subsout, subsref(in(m).dat, subsin));
        end
    else
        if slice > 0
            [subs{1:numel(d)}] = deal(':');
            subs{slice} = idx;
            subs = struct('type', '()', 'subs', {subs});
            out = single(subsref(in.dat, subs));
        else
            out = single(in.dat());
        end
        if numel(d)>4 && d(4)==1
            out = reshape(out,[d(1:3) d(5:end)]);
        end
    end
    return
end
error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
% MemMapData()
function [out,Mn] = MemMapData(in,dimcat)
% Memory map a nifti file (or a series of nifti files).
%
%   . This function outputs a file_array (eventually a virtual 
%     concatenation of file_array objects, built using `cat`).
%   . If `in` contains a series of files, they must all have the same
%     dimensions. If the concatenation dimension is not specified, it is  
%     the 4-th by default.
%   . If `in` contains a series of files, they are assumed to all share th
%     same orientation matrix.
%   . The Nifti standard reserves the 4-th dimension for time (which does
%     not exist in our model). This time dimension is removed from the
%     mapped data (so the 4-th mapped dimension may correspond to the 5-th 
%     file dimension).
%
% FORMAT [out,Mn] = MemMapData(in,dimcat)
% in     - [cell of] nifti, filename, array or file_array
% dimcat - dimension along which to concatenate [4]
% out    - file_array
% Mn     - Orientation matrix (of the first volume)

if nargin < 2
    dimcat = 4;
end
if ~iscell(in), in = {in}; end
if isnumeric(in{1}),         out = in{1}; return; end
if isa(in{1}, 'file_array'), out = cat(dimcat, in{:}); return; end
if isa(in{1},'char'),        in  = {nifti(in{:})}; end
if isa(in{1}, 'nifti'),      in  = cat(2,in{:}); end
if isa(in,'nifti')
    C    = numel(in);
    Mn   = in(1).mat;
    dats = cell(1,C);
    for c=1:C
        dats{c} = in(c).dat;
        if numel(size(in(c).dat)) > 4 && size(in(c).dat,4) == 1
            dats{c}.dat(4) = [];
        end
    end
    out = cat(dimcat, dats{:});
    return
end
error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
% GetMat()
function Mat = GetMat(fin)
% Returns the orientation matrix of an array, file_array or [series of]
% nifti files.
if isnumeric(fin)
    Mat = eye(4);
    return;
end
if isa(fin,'char')
    fin = nifti(fin);
elseif isa(fin, 'file_array')
    fin = nifti({fin.fname});
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
function [d,M,N] = GetSize(fin)
% Returns the size of the 3D lattice of a volume, and its 4-th and 5-th
% dimensions (usually, the 4-th dimension maps to tissue classes or
% deformation direction; the 5-th dimension maps to repeats of these
% features, e.g., principal components).
d = [GetDimensions(fin) 1 1 1];
M = d(4);
N = d(5);
d = d(1:3);
end
%==========================================================================

%==========================================================================
% InitDat()
function dat = InitDat(data,sett)
% Convert the input data structure (that usually contains filenames) into a
% conveniant structure that cen be used during processing.
%
% The input data structure must be a:
% . a cell array (of arrays, filenames, nifti objects or file_array objects)
% . a char array (of filenames)
% . a string array (of filenames)
% . a struct array with field 'F'
%
% Furthermore, if it is a struct array, it can have the fields:
% . do_bf:  do/don't update the bias field for this subject
% . is_ct:  is/isn't a CT image
% . ix_pop: index of the population it belongs to
% . labels: an image of ground truth labels
% . Mat:    orientation matrix
%
% FORMAT odat = InitDat(idat,sett)
% idat - Input data strucutre
% odat - Output data structure
% sett - Option structure

% Parse function settings
do_gmm = sett.do.gmm;
do_pca = sett.do.pca;
run2d  = sett.gen.run2d;  % Extract 2D slices from 3D volumes

% Initialise for each subject
N   = numel(data);
dat = struct;
dat(N).f = [];
for n=1:N

    % Init datn.f
    if isstruct(data(n)) && isfield(data(n),'F'), F = data(n).F;
    else,                                         F = data(n);
    end
    if ~iscell(F), F = {F}; end
    dat(n).f = MemMapData(F);
    if run2d, dat(n).f = GetSlice(dat(n).f,run2d); end
    
    % Get number of channels
    [~,C] = spm_mb_io('GetSize',dat(n).f);
    
    % Parameters
    dat(n).M     = GetMat(dat(n).f);  % Orientation matrix
    dat(n).q     = zeros(6,1);        % Affine parameters
    dat(n).v     = [];                % Velocity field
    dat(n).psi   = [];                % Warp
    dat(n).E     = [0 0 0];           % Energy: (gmm) (velocity) (bias field) 
    dat(n).ss    = struct;            % Sufficent statistics for energy
    
    if do_pca
        dat(n).z   = [];              % PCA latent variable
    end
        
    if do_gmm
        dat(n).mog = [];              % GMM parameters
        dat(n).bf  = [];              % Bias field  
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
    if isstruct(data(n)) && isfield(data(n),'labels') && ~isempty(data(n).labels{1})
        dat(n).labels    = data(n).labels;
        dat(n).labels{1} = MemMapData(dat(n).labels{1});
        if run2d, dat(n).labels{1} = GetSlice(dat(n).labels{1},run2d); end
    else
        dat(n).labels = [];        
    end
            
    % Orientation matrix (image voxel-to-world)
    dat(n).Mat = GetMat(dat(n).f);    
    if run2d
        vx         = sqrt(sum(dat(n).Mat(1:3,1:3).^2));
        dat(n).Mat = [diag(vx) zeros(3,1); 0 0 0 1];
    end
end
end
%==========================================================================

%==========================================================================
% MakeModel()
function model = MakeModel(dat,shape,model,sett)

% Parse function settings
do_gmm               = sett.do.gmm;
do_pca               = sett.do.pca;
do_updt_int          = sett.do.updt_int;
do_updt_template     = sett.do.updt_template;
do_updt_subspace     = sett.do.updt_subspace;
do_updt_latent_prior = sett.do.updt_latent_prior;
do_updt_res_prior    = sett.do.updt_res_prior;
dir_res              = sett.write.dir_res;
mg_ix                = sett.model.mg_ix;
write_model          = sett.write.model;

if isempty(dir_res) 
    pth     = fileparts(dat(1).f(1).dat.fname);
    dir_res = pth; 
end

% -----------
% Shape model

% TODO: Maybe better to save relative path (i.e., filename) rather than
% absolute path, so that people can copy/paste their model files without
% breaking anything. the only requirement would be for the model nifti
% files (template, subspace) to be located next to the model mat file.

if do_updt_template
    model.shape.template = fullfile(dir_res ,'mu_log_spm_mb.nii');
end
if do_pca
    if do_updt_subspace
        model.shape.subspace = fullfile(dir_res ,'subspace_spm_mb.nii');
        model.shape.subspace_cov = shape.Su;
        model.shape.v_settings = sett.var.v_settings;
        if numel(model.shape.v_settings) == 8
            model.shape.v_settings = model.shape.v_settings(4:8);
        end
    end
    if do_updt_latent_prior
        model.shape.latent_prior = shape.A;
        model.shape.latent_df    = shape.nA;
    end
    if do_updt_res_prior
        model.shape.res_prior = shape.lam;
        model.shape.res_df    = shape.nlam;
    end
end

if do_updt_int && do_gmm
    % Appearance related
    p_ix = spm_mb_appearance('GetPopulationIdx',dat);
    Npop = numel(p_ix);

    model.appear       = struct('pr',[],'mg_ix',[]);
    model.appear.pr    = containers.Map;    
    model.appear.mg_ix = mg_ix;
    for p=1:Npop
        n      = p_ix{p}(1);
        datn   = dat(n);
        ix_pop = datn.ix_pop;
        pr     = datn.mog.pr;
        
        model.appear.pr(num2str(ix_pop)) = pr;
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
function SaveTemplate(mu,sett)

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
% SaveSubspace()
function SaveSubspace(U,sett)

% Parse function settings
dir_res          = sett.write.dir_res;
do_pca           = sett.do.pca;
do_updt_subspace = sett.do.updt_subspace;
Mmu              = sett.var.Mmu;
write_subspace   = sett.write.subspace;

if ~do_pca, return; end

if isempty(dir_res) 
    pth     = fileparts(dat(1).f(1).dat.fname);
    dir_res = pth; 
end

if ~isempty(U) && do_updt_subspace
    if write_subspace
        fname       = fullfile(dir_res ,'subspace_spm_mb.nii');
        fa          = file_array(fname,size(U),'float32');
        Nii         = nifti;
        Nii.dat     = fa;
        Nii.mat     = Mmu;
        Nii.mat0    = Mmu;
        Nii.descrip = 'Principal subspace';
        create(Nii);
        for l=1:size(U,5)
            Nii.dat(:,:,:,:,l) = single(U(:,:,:,:,l));
        end
        
        % 4-th dimension reserved for time
        dim = size(U);
        dim = [dim(1:3) 1 dim(4:end)];
        Nii.dat.dim = dim;
        create(Nii);
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
% Returns the dimensions of an array, file_array or [series of] nifti
% files. The function ensures that the output has at least 5 elements by
% one-padding.
if isnumeric(fin) || isa(fin, 'file_array')
    d = size(fin);
    d(end+1:5) = 1;
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
    d(end+1:5) = 1;
    return
end
end
%==========================================================================

%==========================================================================
% GetSlice()
function fn = GetSlice(fn,direction)
% Exract the central slice from a 3D Volume in a given direction.
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

%==========================================================================
function VO = SubVol(V,bb,prefix,deg)
% Extract a subvolume
% FORMAT VO = subvol(V,bb,prefix)
% V      - SPM volume object
% bb     - bounding box
% prefix - file prefix (if empty -> overwrites)
% VO     - resized image
%
% Example:
%     V = spm_vol(spm_select(1,'image'));
%     subvol(V,[32 64 ; 1 64 ; 1 48]');
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin < 3, prefix = 'sv'; end
if nargin < 4, deg    = 0;     end

bb      = round(bb);
bb      = sort(bb);
bb(1,:) = max(bb(1,:),[1 1 1]);
bb(2,:) = min(bb(2,:),V.dim(1:3));

[~,~,oext] = fileparts(prefix);

VO            = V;
[pth,nam,ext] = fileparts(V.fname);
if ~isempty(oext)
    VO.fname  = prefix;
else
    VO.fname  = fullfile(pth,[prefix nam ext]);
end
VO.dim(1:3)   = diff(bb)+1;
VO.mat        = V.mat*spm_matrix((bb(1,:)-1));

VO = spm_create_vol(VO);
for z=1:VO.dim(3)
    M   = V.mat\VO.mat*spm_matrix([0 0 z]);
    img = spm_slice_vol(V,M,VO.dim(1:2),deg);
    VO  = spm_write_plane(VO,img,z);
end
end
%==========================================================================