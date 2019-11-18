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
    cl             = cell(C,1);    
    args           = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
    dat(n).bf.chan = struct(args{:});
    
    vx = sqrt(sum(dat(n).Mat(1:3,1:3).^2));                     
    sd = vx(1)*d(1)/fwhm; d3(1) = ceil(sd*2);
    sd = vx(2)*d(2)/fwhm; d3(2) = ceil(sd*2);
    sd = vx(3)*d(3)/fwhm; d3(3) = ceil(sd*2);
    
    % Precision (inverse covariance) of Gaussian prior on bias field parameters
    ICO = spm_bias_lib('regulariser','bending',d,d3,vx);
    ICO = single(ICO*reg);
    
    [x0,y0,~] = ndgrid(single(1:d(1)),single(1:d(2)),1);
    z0        = single(1:d(3));

    for c=1:C
        % GAUSSIAN REGULARISATION for bias correction                        
        dat(n).bf.chan(c).C = ICO;

        % Basis functions for bias correction
        dat(n).bf.chan(c).B3  = single(spm_dctmtx(d(3),d3(3),z0));
        dat(n).bf.chan(c).B2  = single(spm_dctmtx(d(2),d3(2),y0(1,:)'));
        dat(n).bf.chan(c).B1  = single(spm_dctmtx(d(1),d3(1),x0(:,1)));

        % Initial parameterisation of bias field
        dat(n).bf.chan(c).T = zeros(d3,'single');

        % Change DC component of bias field to make intensities more
        % simillar between MR images.
        b1 = dat(n).bf.chan(c).B1(1,1);
        b2 = dat(n).bf.chan(c).B2(1,1);
        b3 = dat(n).bf.chan(c).B3(1,1);

        dat(n).bf.chan(c).T(1,1,1) = 1/(b1*b2*b3)*log(scl(c));
    end    
    
    % struct used for rescaling images using DC component of bias fields
    dc           = struct;
    dc.int       = zeros(1,C);
    dc.ln        = zeros(1,C);
    dat(n).bf.dc = dc;
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
updt_bf = sett.do.updt_bf;

if ~do_gmm, return; end
    
K1 = K + 1;
lb = struct('sum', NaN, 'X', [], 'XB', [], ...
            'Z', [], 'P', [], 'MU', [], 'A', []);

for n=1:numel(dat)            
    % GMM        
    [d,C] = spm_multireg_io('GetSize',dat(n).f);
    fn    = spm_multireg_io('GetData',dat(n).f);                    
    fn    = reshape(fn,[prod(d(1:3)) C]);                      
    fn    = spm_multireg_util('MaskF',fn);
    
    % Modulate with bias field
    if updt_bf
        bf    = spm_multireg_io('GetBiasField',dat(n).bf.chan,d);
        lb.XB = spm_multireg_energ('LowerBound','XB',bf);
    else
        bf = ones(1,C);
    end
    fn = bf.*fn;

    % Initial means and precisions from image channel max
    mog    = init_gmm(fn,K1);        
    mog.lb = lb;
    
    dat(n).mog = mog;
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
% init_gmm()
function mog = init_gmm(fn,K)
C = size(fn,2);

% Posterior
mu = zeros(C,K);
A  = zeros(C,C,K);        
for c=1:C
    mx       = nanmax(fn(:,c));                            
    mu(c,:)  = (0:(K - 1))'*mx/(1.5*K);
    A(c,c,:) = mx/K;        
    A(c,c,:) = 1/A(c,c,:);
end   

mog.po.m = mu;
mog.po.b = ones(1,K);
mog.po.n = C*ones(1,K);
mog.po.V = bsxfun(@times, A, reshape(mog.po.n, [1 1 K])); % Expected precision

% Prior (uninformative)
mog.pr.m = zeros(C,K);
mog.pr.b = ones(1,K);
mog.pr.n = C*ones(1,K);
mog.pr.V = bsxfun(@times, repmat(eye(C),[1 1 K]), reshape(mog.pr.n, [1 1 K]));
end
%==========================================================================            