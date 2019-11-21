function varargout = spm_multireg_init(varargin)
%__________________________________________________________________________
%
% Initialisation functions for spm_multireg.
%
% FORMAT dat = spm_multireg_init('InitBiasField',dat,sett)
% FORMAT dat = spm_multireg_init('InitDat',F,sett)
% FORMAT dat = spm_multireg_init('InitDef',dat,sett)
% FORMAT dat = spm_multireg_init('InitGMM',dat,K,sett)
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
    case 'InitBiasField'
        [varargout{1:nargout}] = InitBiasField(varargin{:});  
    case 'InitDat'
        [varargout{1:nargout}] = InitDat(varargin{:});   
    case 'InitDef'
        [varargout{1:nargout}] = InitDef(varargin{:});           
    case 'InitGMM'
        [varargout{1:nargout}] = InitGMM(varargin{:});                
    otherwise
        help spm_multireg_init
        error('Unknown function %s. Type ''help spm_multireg_init'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% InitBiasField()
function dat = InitBiasField(dat,sett)

% Parse function settings
do_bf_norm = sett.do.bf_norm;
fwhm       = sett.bf.fwhm;
reg        = sett.bf.reg;
updt_bf    = sett.do.updt_bf;

if ~updt_bf, return; end

for n=1:numel(dat)    
    [d,C]          = spm_multireg_io('GetSize',dat(n).f);
    if do_bf_norm
        val = 1e3;            
        fn  = spm_multireg_io('GetData',dat(n).f);                    
        fn  = reshape(fn,[prod(d(1:3)) C]);                      
        fn  = spm_multireg_util('MaskF',fn);
        scl = double(val./nanmean(fn,1)); 
    else        
        scl = ones(1,C);
    end
    
    % Get bias field parameterisation struct
    chan        = spm_multireg_io('GetBiasFieldStruct',C,d,dat(n).Mat,reg,fwhm,scl);
    dat(n).bf.T = {chan(:).T};
    
    % struct used for rescaling images using DC component of bias fields
    dc            = struct;
    dc.int        = zeros(1,C);
    dc.ln         = zeros(1,C);
    dat(n).bf.dc  = dc;
end
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
            dat(n).f = get_slice(F{n},run2d);
        else
            dat(n).f = single(F{n});
        end
    elseif isa(F(n),'nifti') || (iscell(F(n)) && (isa(F{n},'char') || isa(F{n},'nifti')))
        % Input F is nifti (path or object) -> store as nifti
                       
        if isa(F(n),'nifti')
            dat(n).f = F(n);        
        elseif iscell(F(n)) 
            if isa(F{n},'char')
                dat(n).f = nifti(F{n});        
            elseif isa(F{n},'nifti')
                dat(n).f = nifti;
                C        = numel(F{n});
                for c=1:C
                    dat(n).f(c) = F{n}(c);
                end
            end
        end
        
        if run2d
            % Get 2D slice from 3D data
            fn       = spm_multireg_io('GetData',dat(n).f);
            dat(n).f = get_slice(fn,run2d);
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
    if isa(dat(n).f,'nifti') && ~run2d
        dat(n).Mat = dat(n).f(1).mat;        
    end
end
end
%==========================================================================

%==========================================================================
% InitDef()
function dat = InitDef(dat,sett)

% Parse function settings
B       = sett.registr.B;
d       = sett.var.d;
dir_res = sett.write.dir_res;
Mmu     = sett.var.Mmu;

v    = zeros([d,3],'single');
psi1 = spm_multireg_util('Identity',d);
for n=1:numel(dat)
    dat(n).q = zeros(size(B,3),1);
    if isnumeric(dat(n).f)
        dat(n).v   = v;
        dat(n).psi = psi1;
    else
        if isa(dat(n).f,'nifti')
            [~,nam,~] = fileparts(dat(n).f(1).dat.fname);
            vname    = fullfile(dir_res,['v_' nam '.nii']);
            pname    = fullfile(dir_res,['psi_' nam '.nii']);
            
            fa       = file_array(vname,[d(1:3) 1 3],'float32',0);
            nii      = nifti;
            nii.dat  = fa;
            nii.mat  = Mmu;
            nii.mat0 = Mmu;
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

%==========================================================================
% InitGMM()
function dat = InitGMM(dat,K,sett)

% Parse function settings
do_gmm  = sett.do.gmm;
fwhm    = sett.bf.fwhm;
reg     = sett.bf.reg;
updt_bf = sett.do.updt_bf;

if ~do_gmm, return; end
    
K1 = K + 1;
lb = struct('sum', NaN, 'X', [], 'XB', [], ...
            'Z', [], 'P', [], 'MU', [], 'A', []);

[~,C] = spm_multireg_io('GetSize',dat(1).f);
mx    = zeros(C,numel(dat));
mn    = zeros(C,numel(dat));
vr    = zeros(C,numel(dat));
for n=1:numel(dat)            
    % GMM        
    [df,C] = spm_multireg_io('GetSize',dat(n).f);
    fn     = spm_multireg_io('GetData',dat(n).f);                    
    fn     = reshape(fn,[prod(df(1:3)) C]);                      
    fn     = spm_multireg_util('MaskF',fn);
    
    % Modulate with bias field
    if updt_bf
        chan = spm_multireg_io('GetBiasFieldStruct',C,df,dat(n).Mat,reg,fwhm,[],dat(n).bf.T);
        bf   = spm_multireg_io('GetBiasField',chan,df);        
    else
        bf = ones(1,C);
    end
    fn = bf.*fn;

    % Init GMM posterior
    [po,mx(:,n),mn(:,n),vr(:,n)] = init_gmm_po(fn,K1);        
    
    mog.po = po;
    
    % Init lower bound struct
    mog.lb = lb;
    
    dat(n).mog = mog;
end

% Init GMM empirical prior
pr = init_gmm_pr(mx,mn,vr,K1);
for n=1:numel(dat)                
    dat(n).mog.pr = pr;
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% get_slice()
function fn = get_slice(fn,direction)
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

%==========================================================================    
% init_gmm_po()
function [po,mx,mn,vr] = init_gmm_po(fn,K)
C  = size(fn,2);
mx = zeros(1,C);
mn = zeros(1,C);
vr = zeros(1,C);
mu = zeros(C,K);
A  = zeros(C,C,K);  
n  = 3;
for c=1:C
    mx(c)    = nanmax(fn(:,c));                         
    mn(c)    = nanmean(fn(:,c));
    vr(c)    = nanvar(fn(:,c));
    mu(c,:)  = (0:(K - 1))'*mx(c)/(1.5*K);
    A(c,c,:) = mx(c)/(1.5*K);      
%     vrc      = vr(c)/(K + 1);
%     mnc      = 1.0*mn(c);
%     sd       = sqrt(vrc);
%     mu(c,:)  = abs(linspace(mnc - n*sd,mnc + n*sd,K));    
%     A(c,c,:) = vrc;
    A(c,c,:) = 1/A(c,c,:);
end   

po   = struct('m',[],'b',[],'n',[],'V',[]);
po.m = mu;
po.b = ones(1,K);
po.n = C*ones(1,K);
po.V = bsxfun(@times, A, reshape(po.n, [1 1 K])); % Expected precision
end
%==========================================================================      

%==========================================================================    
% init_gmm_pr()
function pr = init_gmm_pr(mx,mn,vr,K)
C   = size(mx,1);
mmx = max(mx,[],2);
mmn = mean(mn,2);
mvr = mean(vr,2);

% pr   = struct('m',[],'b',[],'n',[],'V',[]);
% pr.m = zeros(C,K);
% pr.b = ones(1,K);
% pr.n = C*ones(1,K);
% pr.V = bsxfun(@times, repmat(eye(C),[1 1 K]), reshape(pr.n, [1 1 K]));
mu = zeros(C,K);
A  = zeros(C,C,K);        
n  = 3;
for c=1:C         
    mu(c,:)  = (0:(K - 1))'*mmx(c)/(1.5*K);
    A(c,c,:) = mmx(c)/(1.5*K);    
%     vrc      = mvr(c)/(K + 1);
%     mnc      = 1.0*mn(c);
%     sd       = sqrt(vrc);
%     mu(c,:)  = abs(linspace(mnc - n*sd,mnc + n*sd,K));    
%     A(c,c,:) = vrc;        
    A(c,c,:) = 1/A(c,c,:);
end   

pr   = struct('m',[],'b',[],'n',[],'V',[]);
pr.m = mu;
pr.b = ones(1,K);
pr.n = C*ones(1,K);
pr.V = bsxfun(@times, A, reshape(pr.n, [1 1 K])); % Expected precision
end
%==========================================================================   