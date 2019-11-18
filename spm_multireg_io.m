function varargout = spm_multireg_io(varargin)
%__________________________________________________________________________
%
% I/O functions for spm_multireg.
%
% FORMAT [zn,lx,lz]    = spm_multireg_io('ComputeResponsibilities',datn,fn,mu,code)
% FORMAT to            = spm_multireg_io('CopyFields',from,to)
% FORMAT [bfn,lln]     = spm_multireg_io('GetBiasField',chan,d,varargin)
% FORMAT chan          = spm_multireg_io('GetBiasFieldStruct',C,d,Mat,reg,fwhm,scl,T)
% FORMAT [P,datn,code] = spm_multireg_io('GetClasses',datn,mu,sett,get_k1)
% FORMAT out           = spm_multireg_io('GetData',in)
% FORMAT Mat           = spm_multireg_io('GetMat',fin)
% FORMAT [d,M]         = spm_multireg_io('GetSize',fin)
% FORMAT psi           = spm_multireg_io('ResavePsiSub',datn,sett)
% FORMAT dat           = spm_multireg_io('SaveImages',dat,mu,sett)
% FORMAT fout          = spm_multireg_io('SetData',fin,f)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_io
    error('Not enough argument. Type ''help spm_multireg_io'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'ComputeResponsibilities'
        [varargout{1:nargout}] = ComputeResponsibilities(varargin{:});   
    case 'CopyFields'
        [varargout{1:nargout}] = CopyFields(varargin{:});                   
    case 'GetBiasField'
        [varargout{1:nargout}] = GetBiasField(varargin{:});  
    case 'GetBiasFieldStruct'
        [varargout{1:nargout}] = GetBiasFieldStruct(varargin{:});  
    case 'GetClasses'
        [varargout{1:nargout}] = GetClasses(varargin{:});    
    case 'GetData'
        [varargout{1:nargout}] = GetData(varargin{:});            
    case 'GetMat'
        [varargout{1:nargout}] = GetMat(varargin{:});             
    case 'GetSize'
        [varargout{1:nargout}] = GetSize(varargin{:});    
    case 'ResavePsiSub'
        [varargout{1:nargout}] = ResavePsiSub(varargin{:});    
    case 'SaveImages'
        [varargout{1:nargout}] = SaveImages(varargin{:});    
    case 'SetData'
        [varargout{1:nargout}] = SetData(varargin{:});        
    otherwise
        help spm_multireg_io
        error('Unknown function %s. Type ''help spm_multireg_io'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% ComputeResponsibilities()
function [zn,lx,lz] = ComputeResponsibilities(datn,fn,mu,code)

% Is there missing data?
L       = unique(code);
do_miss = numel(L) > 1;

% Posterior
m = datn.mog.po.m;
b = datn.mog.po.b;
V = datn.mog.po.V;
n = datn.mog.po.n;

if do_miss, const = spm_gmm_lib('Const', {m,b}, {V,n}, L);
else,       const = spm_gmm_lib('Const', {m,b}, {V,n});
end

fn = spm_gmm_lib('Marginal', fn, {m,V,n}, const, {code,L});
zn = spm_gmm_lib('Responsibility', fn, mu);

if nargout > 1
    lx = spm_multireg_energ('LowerBound','X',fn,zn,code,{m,b},{V,n}); 
    lz = spm_gmm_lib('KL', 'Categorical', zn, 1, mu);
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
% GetBiasField()
function [bfn,lln] = GetBiasField(chan,d,varargin)
C  = numel(chan);
I  = prod(d);
Iz = prod(d(1:2));
nz = d(3);   
if numel(varargin) == 0
    % Compute full bias field (for all channels)
    bfn = zeros([I C],'single');
    lln = zeros(1,C);
    for c=1:C           
        lln(c) = double(-0.5*chan(c).T(:)'*chan(c).C*chan(c).T(:));    

        for z=1:nz
            ix        = IndexSlice2Vol(z,Iz);
            bf_c      = transf(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);
            bf_c      = bf_c(:);      
            bfn(ix,c) = single(exp(bf_c));        
        end
    end
else
    % Compute just for one channel
    bfn = varargin{1};
    c   = varargin{2};
    lln = varargin{3};

    lln(c) = double(-0.5*chan(c).T(:)'*chan(c).C*chan(c).T(:));    

    for z=1:nz
        ix        = IndexSlice2Vol(z,Iz);
        bf_c      = transf(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);
        bf_c      = bf_c(:);      
        bfn(ix,c) = single(exp(bf_c));        
    end
end
end
%==========================================================================

%==========================================================================
% GetBiasFieldStruct()
function chan = GetBiasFieldStruct(C,d,Mat,reg,fwhm,scl,T,samp)
if nargin < 7, T    = {}; end
if nargin < 8, samp = 1; end

cl   = cell(C,1);    
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
chan = struct(args{:});

vx = sqrt(sum(Mat(1:3,1:3).^2));                     
sd = vx(1)*d(1)/fwhm; d3(1) = ceil(sd*2);
sd = vx(2)*d(2)/fwhm; d3(2) = ceil(sd*2);
sd = vx(3)*d(3)/fwhm; d3(3) = ceil(sd*2);

% Precision (inverse covariance) of Gaussian prior on bias field parameters
ICO = spm_bias_lib('regulariser','bending',d,d3,vx);
ICO = single(ICO*reg);

samp      = max([1 1 1],round(samp*[1 1 1]./vx));
[x0,y0,~] = ndgrid(single(1:samp(1):d(1)),single(1:samp(2):d(2)),1);
z0        = single(1:samp(3):d(3));

for c=1:C
    % GAUSSIAN REGULARISATION for bias correction                        
    chan(c).C = ICO;

    % Basis functions for bias correction
    chan(c).B3 = single(spm_dctmtx(d(3),d3(3),z0));
    chan(c).B2 = single(spm_dctmtx(d(2),d3(2),y0(1,:)'));
    chan(c).B1 = single(spm_dctmtx(d(1),d3(1),x0(:,1)));

    if isempty(T)
        % Initial parameterisation of bias field
        chan(c).T = zeros(d3,'single');
    else
        % Parameterisation given
        chan(c).T = T{c};
    end

    if ~isempty(scl)
        % Change DC component of bias field to make intensities more
        % simillar between MR images.
        b1               = chan(c).B1(1,1);
        b2               = chan(c).B2(1,1);
        b3               = chan(c).B3(1,1);
        chan(c).T(1,1,1) = 1/(b1*b2*b3)*log(scl(c));
    end
end   
end
%========================================================================== 

%==========================================================================
% GetClasses()
function [P,datn,code] = GetClasses(datn,mu,sett,get_k1)
if nargin < 4, get_k1 = false; end

code = [];
if ~isfield(datn,'mog')
    P = GetData(datn.f);

    if nargout > 1
        % Make mask
        msk = sum(P,4) > 0; % Removes voxels that sums to zero across classes..
        msk = msk & sum(~isfinite(P),4) == 0; % ..and voxels that are not finite in segmentation..
        msk = msk & sum(~isfinite(mu),4) == 0; % ..and voxels that are not finite in template

        % Compute subject-specific categorical cross-entropy loss between
        % segmentation and template
        tmp       = sum(P.*mu,4) - spm_multireg_util('lse',mu,4);  
        datn.E(1) = -sum(tmp(msk));
    end
else
    if nargout > 1
        [P,datn,code] = GetClassesFromGMM(datn,mu,sett,get_k1);
    else
        P = GetClassesFromGMM(datn,mu,sett,get_k1);
    end
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
% ResavePsiSub()
function psi = ResavePsiSub(datn,sett)

% Parse function settings
B       = sett.registr.B;
dir_res = sett.write.dir_res;
Mmu     = sett.var.Mmu;

d    = GetSize(datn.f);
q    = double(datn.q);
Mn   = datn.Mat;
psi1 = GetData(datn.psi);
psi  = spm_multireg_util('Compose',psi1,spm_multireg_util('Affine',d,Mmu\spm_dexpm(q,B)*Mn));
psi  = reshape(reshape(psi,[prod(d) 3])*Mmu(1:3,1:3)' + Mmu(1:3,4)',[d 1 3]);
if isa(datn.psi(1),'nifti')
    to_delete = datn.psi(1).dat.fname;
    [~,nam,~] = fileparts(datn.f(1).dat.fname);
    datn.psi(1).dat.fname = fullfile(dir_res,['y_' nam '.nii']);
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
% SaveImages()
function dat = SaveImages(dat,mu,sett)

% Parse function settings
dir_res = sett.write.dir_res;
Mmu     = sett.var.Mmu;

for n=1:numel(dat)
    dat(n).psi = ResavePsiSub(dat(n),sett);
end

if ~isempty(mu)
    % Save mu (log)
    fa       = file_array(fullfile(dir_res ,'mu_log.nii'),size(mu),'float32',0);
    Nmu      = nifti;
    Nmu.dat  = fa;
    Nmu.mat  = Mmu;
    Nmu.mat0 = Mmu;
    Nmu.descrip = 'Mean parameters (log)';
    create(Nmu);
    Nmu.dat(:,:,:,:) = mu;
    
    % Save mu (softmax)    
    mu       = spm_multireg_util('softmaxmu',mu,4);    
    fa       = file_array(fullfile(dir_res ,'mu_softmax.nii'),size(mu),'float32',0);
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
%
% Utility functions
%
%==========================================================================

%==========================================================================
% GetClassesFromGMM()
function [zn,datn,code] = GetClassesFromGMM(datn,mu,sett,get_k1)
if nargin < 4, get_k1 = false; end

% Parse function settings
nit_gmm      = sett.nit.gmm;
nit_gmm_miss = sett.nit.gmm_miss;
samp         = sett.gen.samp_gmm;
updt_bf      = sett.do.updt_bf;

fn     = GetData(datn.f);
[d0,C] = GetSize(datn.f);
fn     = reshape(fn,[prod(d0(1:3)) C]);
Mat    = datn.Mat;
W      = 1;

if updt_bf, bf = spm_multireg_io('GetBiasField',datn.bf.chan,d0);
else,       bf = ones(1,C);
end

% Missing data stuff
fn      = spm_multireg_util('MaskF',fn);
code    = spm_gmm_lib('obs2code', fn);
fn      = bf.*fn;
L       = unique(code);
do_miss = numel(L) > 1;

% GMM posterior
m  = datn.mog.po.m;
b  = datn.mog.po.b;
V  = datn.mog.po.V;
n  = datn.mog.po.n;

% GMM prior
m0 = datn.mog.pr.m;
b0 = datn.mog.pr.b;
V0 = datn.mog.pr.V;
n0 = datn.mog.pr.n;

% Lower bound
lb = datn.mog.lb;

% Make softmaxed K + 1 template
K1 = numel(b);
K  = K1 - 1;
if size(mu,4) < K1
    mu = log(spm_multireg_util('softmaxmu',mu,4));
end
mu = reshape(mu,[prod(d0(1:3)) K1]);

if nargout > 1    
    % Update GMM and get responsibilities
               
    if samp > 1
        % Subsample (runs faster, lower bound is corrected by scalar W)              
        [code0,code,fn0,fn,mu0,mu,W] = spm_multireg_util('SubSample',samp,Mat,d0,code,fn,mu);
    end        
    
    [zn,mog,~,lb] = spm_gmm_loop({fn,W},{{m,b},{V,n}},{'LogProp', mu}, ...
                                 'GaussPrior',   {m0,b0,V0,n0}, ...
                                 'Missing',      do_miss, ...
                                 'LowerBound',   lb, ...
                                 'MissingCode',  {code,L}, ...
                                 'IterMax',      nit_gmm, ...
                                 'Tolerance',    1e-4, ...
                                 'SubIterMax',   nit_gmm_miss, ...
                                 'SubTolerance', 1e-4, ...
                                 'Verbose',      0);
    clear mu fn

    % Update datn
    datn.mog.po.m = mog.MU; % GMM posteriors
    datn.mog.po.b = mog.b;
    datn.mog.po.V = mog.V;
    datn.mog.po.n = mog.n;        
    
    if samp > 1
        % Compute responsibilities on original data         
        zn = ComputeResponsibilities(datn,fn0,mu0,code0);
    end    
        
    datn.mog.lb = lb; % Lower bound            
    datn.E(1)   = -sum(lb.sum(end)); % objective function
else
    % Just compute responsibilities    
    zn = ComputeResponsibilities(datn,fn,mu,code);
    
end

if ~get_k1
    % Get 4D versions of K1 - 1 classes
    zn = reshape(zn(:,1:K),[d0(1:3) K]);
end
end
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
% IndexSlice2Vol()
function ix = IndexSlice2Vol(z,Iz)
ix = ((z - 1)*Iz + 1):z*Iz;
end
%==========================================================================

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
end
%==========================================================================