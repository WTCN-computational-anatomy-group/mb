function varargout = spm_mb_appearance(varargin)
%__________________________________________________________________________
%
% Functions for appearance model related.
%
% FORMAT [bfn,lln]  = spm_mb_appearance('BiasField',chan,d,varargin)
% FORMAT chan       = spm_mb_appearance('BiasFieldStruct',datn,C,df,reg,fwhm,dc,T,samp)
% FORMAT labels     = spm_mb_appearance('GetLabels',datn,sett,do_samp)
% FORMAT [nvx_obs,msk_allmiss] = spm_mb_appearance('GetNumObsVox',f,is_ct,ax)
% FORMAT p_ix       = spm_mb_appearance('GetPopulationIdx',dat)
% FORMAT [dat,sett] = spm_mb_appearance('Init',dat,model,K,sett)
% FORMAT fn         = spm_mb_appearance('Mask',fn,is_ct)
% FORMAT zn         = spm_mb_appearance('Responsibility',m,b,W,n,fn,mu,msk_chn)
% FORMAT [zn,datn]  = spm_mb_appearance('Update',datn,mun0,sett)
% FORMAT dat        = spm_mb_appearance('UpdatePrior',dat,mu,sett,add_po_observation)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_appearance
    error('Not enough argument. Type ''help spm_mb_appearance'' for help.');
end
id       = varargin{1};
varargin = varargin(2:end);
switch id
    case 'BiasField'
        [varargout{1:nargout}] = BiasField(varargin{:});      
    case 'BiasFieldStruct'
        [varargout{1:nargout}] = BiasFieldStruct(varargin{:});         
    case 'GetLabels'
        [varargout{1:nargout}] = GetLabels(varargin{:});            
    case 'GetNumVoxObserved'
        [varargout{1:nargout}] = GetNumObsVox(varargin{:});           
    case 'GetPopulationIdx'
        [varargout{1:nargout}] = GetPopulationIdx(varargin{:});           
    case 'Init'
        [varargout{1:nargout}] = Init(varargin{:});              
    case 'Mask'
        [varargout{1:nargout}] = Mask(varargin{:});                     
    case 'Responsibility'
        [varargout{1:nargout}] = Responsibility(varargin{:});             
    case 'Update'
        [varargout{1:nargout}] = Update(varargin{:});             
    case 'UpdatePrior'
        [varargout{1:nargout}] = UpdatePrior(varargin{:});            
    otherwise
        help spm_mb_appearance
        error('Unknown function %s. Type ''help spm_mb_appearance'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% BiasField()
function [bfn,lln] = BiasField(chan,d,varargin)
I  = prod(d);
Iz = prod(d(1:2));
nz = d(3);
if numel(varargin) == 0
    % Compute full bias field (for all channels)
    C   = numel(chan);
    bfn = zeros([I C],'single');
    cr  = 1:C;
    lln = zeros(1,C);
else
    % Compute just for one channel
    bfn = varargin{1};
    cr  = varargin{2};
    lln = varargin{3};
end
for c=cr
    lln(c) = double(-0.5*chan(c).T(:)'*chan(c).ICO*chan(c).T(:));

    for z=1:nz
        ix        = IndexSlice2Vol(z,Iz);
        bf_c      = TransformBF(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);
        bf_c      = bf_c(:);
        bfn(ix,c) = single(exp(bf_c));
    end
end
end
%==========================================================================

%==========================================================================
% BiasFieldStruct()
function chan = BiasFieldStruct(datn,C,df,reg,fwhm,dc,T,samp)
if nargin < 7, T    = {}; end
if nargin < 8, samp = 1; end

cl   = cell(C,1);
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
chan = struct(args{:});

do_bf = datn.do_bf;
Mn    = datn.Mat;

vx  = sqrt(sum(Mn(1:3,1:3).^2));
sd1 = vx(1)*df(1)/fwhm; d3(1) = ceil(sd1*2);
sd2 = vx(2)*df(2)/fwhm; d3(2) = ceil(sd2*2);
sd3 = vx(3)*df(3)/fwhm; d3(3) = ceil(sd3*2);

% Precision (inverse covariance) of Gaussian prior on bias field parameters
if 1    
    % Bending energy
    ICO = spm_bias_lib('regulariser','bending',df,d3,vx);
    ICO = ICO*reg;
else
    % spm_preproc8
    kron  = @(a,b) spm_krutil(a,b);
    krn_x = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
    krn_y = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
    krn_z = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
    reg   = 1e-3;
    ICO   = kron(krn_z,kron(krn_y,krn_x)).^(-2)*reg;
    ICO   = sparse(1:length(ICO),1:length(ICO),ICO,length(ICO),length(ICO));
end
    
samp      = max([1 1 1],round(samp*[1 1 1]./vx));
[x0,y0,~] = ndgrid(single(1:samp(1):df(1)),single(1:samp(2):df(2)),1);
z0        = single(1:samp(3):df(3));

for c=1:C
    % GAUSSIAN REGULARISATION for bias correction
    chan(c).ICO = ICO;

    % Basis functions for bias correction
    chan(c).B3 = spm_dctmtx(df(3),d3(3),z0);
    chan(c).B2 = spm_dctmtx(df(2),d3(2),y0(1,:)');
    chan(c).B1 = spm_dctmtx(df(1),d3(1),x0(:,1));

    if isempty(T) || ~do_bf(c)
        % Initial parameterisation of bias field
        chan(c).T = zeros(d3);
    else
        % Parameterisation given
        chan(c).T = T{c};
    end

    if ~isempty(dc) && do_bf(c)
        % Change DC component of bias field to make intensities more
        % simillar between MR images.
        b1               = chan(c).B1(1,1);
        b2               = chan(c).B2(1,1);
        b3               = chan(c).B3(1,1);
        chan(c).T(1,1,1) = 1/(b1*b2*b3)*dc(c);
    end
end   
end
%========================================================================== 

%==========================================================================
% GetLabels()
function labels = GetLabels(datn,sett,do_samp)
if nargin < 3, do_samp = false; end

% Parse function settings
K          = sett.model.K;
samp       = sett.gen.samp;
use_labels = sett.labels.use;

if ~use_labels || isempty(datn.labels) || isempty(datn.labels{1}) || isempty(datn.labels{2})
    % Do not use labels
    K1     = K + 1;
    labels = zeros(1,K1,'single');
    return
end

cm_map     = datn.labels{2}; % cell array that defines the confusion matrix (cm)
ix_nonlabl = numel(cm_map);  % label of nonlabelled data

% There might be more labels in the label image than what we want to use,
% here we get the number of labels we want to use from cm_map
labels2use = [];
for l=1:numel(cm_map) - 1
    if ~isempty(cm_map{l}), labels2use = [labels2use l]; end
end

% Load labels
labels = spm_mb_io('GetData',datn.labels{1});
if do_samp && samp > 1
    % Subsample labels
    Mn     = datn.M;
    labels = SubSample(labels,Mn,samp);    
end

% Use labels2use to keep only labels of interest
labels       = labels(:);
msk          = ismember(labels,labels2use);
labels(~msk) = ix_nonlabl;

% Get confusion matrix that maps from label value to probability value
cm = GetLabelConfMatrix(cm_map,sett);

% Build NxK1 label image using confusion matrix
labels = cm(labels,:);

% Make log probability
labels = log(labels);
end
%==========================================================================

%==========================================================================
% GetNumObsVox()
function [nvx_obs,msk_allmiss] = GetNumObsVox(fn,is_ct,ax)
fn          = ApplyMask(fn,is_ct);
msk_allmiss = all(isnan(fn),ax);
nvx_obs     = sum(msk_allmiss(:) == 0);
end
%==========================================================================

%==========================================================================
% GetPopulationIdx()
function p_ix = GetPopulationIdx(dat)
ix_pop = [dat.ix_pop];
un     = unique(ix_pop);
p_ix   = cell(1,numel(un));
cnt    = 1;
for i=un
    p_ix{cnt} = find(ix_pop == i);
    cnt       = cnt + 1;
end
end
%==========================================================================

%==========================================================================
% Init()
function [dat,sett] = Init(dat,model,K,sett)

% Parse function settings
do_gmm = sett.do.gmm;
mg_ix  = sett.model.mg_ix;
fwhm   = sett.bf.fwhm;
reg    = sett.bf.reg;

if ~do_gmm, return; end

% Get parameters
N   = numel(dat);
Kmg = numel(mg_ix);
K1  = K + 1;

% What model parameters were given?
[template_given,appear_given] = spm_mb_param('SetFit',model,sett);
if appear_given
    % Set mg_ix to mg_ix that was used when learning intensity model
    if isfield(model.appear,'mg_ix'), mg_ix = model.appear.mg_ix;
    else,                             mg_ix = 1:K + 1;
    end
    sett.model.mg_ix = mg_ix;
end

% Lower bound struct
lb  = struct('sum', NaN, 'X', [], 'XB', [], 'Z', [], 'P', [], 'MU', [], ...
             'A', [], 'pr_v', [], 'pr_bf',[]);
             
if template_given
    % Build a weigthed mean (w_mu) that can be used to initi bias field DC scaling
    mu = spm_mb_io('GetData',model.shape.template);                    
        
    w_mu = spm_mb_shape('TemplateK1',mu,4);
    w_mu = exp(w_mu);
    w_mu = sum(sum(sum(w_mu,1),2),3);
    w_mu = w_mu./sum(w_mu);
    w_mu = reshape(w_mu,[1 numel(w_mu)]);
    w_mu = w_mu(mg_ix)./arrayfun(@(x) sum(x == mg_ix), mg_ix); % incorporate multiple Gaussians per tissue
else
    mu   = [];
    w_mu = 1;
end

if appear_given && template_given    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % TODO
    %------------------------------------------------------------     
    
    % TODO
    pr = model.appear.pr;
    dc = w_mu.*pr.m;    
    dc = sum(dc,2); 
    
    for n=1:N
        [df,C] = spm_mb_io('GetSize',dat(n).f);
    
        % Load image data
        fn = spm_mb_io('GetData',dat(n).f);
        fn = reshape(fn,[prod(df(1:3)) C]);
        fn = spm_mb_appearance('Mask',fn,dat(n).is_ct);
        
        for c=1:C
            msk   = isfinite(fn(:,c));
            dc(c) = dc(c)./mean(fn(msk,c));
        end
        dc = log(dc);
        
        if any(dat(n).do_bf == true)
            % Get bias field parameterisation struct
            chan        = spm_mb_appearance('BiasFieldStruct',dat(n),C,df,reg,fwhm,dc);
            dat(n).bf.T = {chan(:).T};
            
            % Get bias field
            bf = spm_mb_appearance('BiasField',chan,df);
        else
            bf = ones([1 C],'single');
        end
        
        % Modulate with bias field
        fn = bf.*fn;
        
        po = InitSimplePosteriorGMM(dat(n),fn,mu,pr,K1,mg_ix,sett);
         
        mog.po   = po;
        mog.pr   = pr;
        mog.lb   = lb;
        mog.mg_w = ones([1 K1]);
   
        dat(n).mog = mog;                        
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % TODO
    %------------------------------------------------------------

    % TODO
    [ix_ct,ix_mri1,ix_mri2] = spm_mb_io('GetCTandMRI',dat,sett);

    if ~isempty(ix_ct) && (isempty(ix_mri1) && isempty(ix_mri2))
        % TODO
        
        dat(ix_ct) = InitGMM(dat(ix_ct),sett);
    else
        % TODO
        
        if ~isempty(ix_ct)
            % Init CT subjects    
            Nct = numel(ix_ct);
            mx  = 3000;
            mu  = zeros([1 K1]);
            Sig = diag((mx/K1).^2).*ones([1,1,K1]);

            % Compute precision
            W                    = zeros(size(Sig));
            for k=1:K1, W(:,:,k) = inv(Sig(:,:,k)); end

            b  = zeros(1,K1) + 0.01;
            nu = ones(1,K1);

            po       = struct('m',mu,'b',b,'W',W,'n',nu);      
            mog.po   = po;
            mog.pr   = po;
            mog.lb   = lb;
            mog.mg_w = ones([1 K1]);

            for n=1:Nct   
                n1          = ix_ct(n);        
                dat(n1).mog = mog;
            end
        end

        % Run on all MRI subjects
        dat(ix_mri1) = InitGMM(dat(ix_mri1),sett);

        if ~isempty(ix_mri2)
            % TODO
            
            Nmri = numel(ix_mri2);

            b  = zeros(1,K1) + 0.01;
            nu = ones(1,K1);

            for n=1:Nmri   

                n1    = ix_mri2(n);       
                [df,C] = spm_mb_io('GetSize',dat(n1).f);

                po       = struct('m',ones(C,K1),'b',b,'W',repmat(eye(C),[1 1 K1])/C,'n',C*nu);      
                mog.po   = po;
                mog.pr   = po;
                mog.lb   = lb;
                mog.mg_w = ones([1 K1]);

                dat(n1).mog = mog;        

                fn = spm_mb_io('GetData',dat(n1).f);
                fn = reshape(fn,[prod(df(1:3)) C]);
                fn = spm_mb_appearance('Mask',fn,dat(n1).is_ct);

                % Make mean close to val
                val = 1e3*ones(1,C);                        
                dc  = zeros(1,C);
                for c=1:C
                    msk   = isfinite(fn(:,c));
                    dc(c) = val(c)./mean(fn(msk,c));
                end
                dc = log(dc);

                % Get bias field parameterisation struct
                chan         = spm_mb_appearance('BiasFieldStruct',dat(n1),C,df,reg,fwhm,dc);
                dat(n1).bf.T = {chan(:).T};
            end    
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Modify posteriors and priors, if using multiple Gaussians per tissue
%------------------------------------------------------------
if K1 < Kmg            
    for n=1:N              
        is_ct = dat(n).is_ct;
        
        % Posterior
        po            = dat(n).mog.po;
        gmm           = spm_gmm_lib('extras', 'more_gmms', {po.m,po.b,po.W,po.n}, mg_ix);        
        if ~is_ct, gmm{1} = abs(gmm{1}); end % make sure non-negative
        po.m          = gmm{1};
        po.b          = gmm{2};
        po.W          = gmm{3};
        po.n          = gmm{4};
        dat(n).mog.po = po;
        
        % Prior
        pr            = dat(n).mog.pr;
        gmm           = spm_gmm_lib('extras', 'more_gmms', {pr.m,pr.b,pr.W,pr.n}, mg_ix);        
        if ~is_ct, gmm{1} = abs(gmm{1}); end % make sure non-negative
        pr.m          = gmm{1};
        pr.b          = gmm{2};
        pr.W          = gmm{3};
        pr.n          = gmm{4};
        dat(n).mog.pr = pr;
        
        dat(n).mog.mg_w   = ones(1,Kmg)./arrayfun(@(x) sum(x == mg_ix), mg_ix);
    end
end
end
%==========================================================================

%==========================================================================
% Mask()
function fn = Mask(fn,is_ct)
C = size(fn,2);
for c=1:C
    fn(:,c) = ApplyMask(fn(:,c),is_ct(c));
end
end
%==========================================================================

%==========================================================================
% Responsibility()
function zn = Responsibility(m,b,W,n,fn,mu,msk_chn)
% Compute responsibilities.
%
% FORMAT zn = Responsibility(m,b,W,n,fn,mu,L,code)
% m       - GMM Means
% b       - GMM Mean d.f.
% W       - GMM Scale matrices
% n       - GMM Scale d.f.
% fn      - Bias-corrected observed image in matrix form [nbvox nbchannel]
% mu      - Deformed and exponentiated template
% msk_chn - Mask of observed channels per code
% zn      - Image of responsibilities [nbvox K]

const = spm_gmm_lib('Normalisation', {m,b}, {W,n}, msk_chn);
fn    = spm_gmm_lib('Marginal', fn, {m,W,n}, const, msk_chn);
zn    = spm_gmm_lib('Responsibility', fn, mu);
end
%==========================================================================

%==========================================================================
% Update()
function [datn,zn] = Update(datn,mun0,sett)
% Update appearance model for a single subject (GMM & bias field)
%
% FORMAT [datn,zn] = Update(datn,mun0,sett)
% datn - Structure holding data for a single subject
% mun0 - Log template
% sett - Structure of settings

% Parse function settings
do_updt_bf   = sett.do.updt_bf;
fwhm         = sett.bf.fwhm;
mg_ix        = sett.model.mg_ix;
nit_bf       = sett.nit.bf;
nit_gmm      = sett.nit.gmm;
nit_gmm_miss = sett.nit.gmm_miss;
nit_appear   = sett.nit.appear;
nit_lsbf     = sett.optim.nls_bf;
reg          = sett.bf.reg;
samp         = sett.gen.samp;
tol_gmm      = sett.appear.tol_gmm;
scal_bf      = sett.optim.scal_bf;

% Parameters
[df,C]     = spm_mb_io('GetSize',datn.f);
K          = size(mun0,4);
K1         = K + 1;
Kmg        = numel(mg_ix);
Mn         = datn.Mat;
scl_samp   = 1; % sampling scaling (defined as W = prod(d0(1:3))/prod(d(1:3)), when samp > 1)
do_bf      = datn.do_bf;
is_ct      = datn.is_ct;
mg_w         = datn.mog.mg_w;

% Get image data
fn = spm_mb_io('GetData',datn.f);

% % Store template voxels for where there are no observations in the image
% % data. These values will be used at the end of this function to fill in
% % responsibilities with NaNs.
% [~,msk_allmiss] = GetNumVoxObserved(fn,is_ct);
% bg_mun          = zeros([nnz(msk_allmiss) K],'single');
% for k=1:K
%     kbg_mun     = mun0(:,:,:,k);
%     kbg_mun     = kbg_mun(msk_allmiss);
%     bg_mun(:,k) = kbg_mun;
% end
% clear kbg_mun
% bg_mun = spm_mb_shape('TemplateK1',bg_mun,2);
% bg_mun = exp(bg_mun);
% bg_mun = bg_mun(:,1:end - 1);

% Get amount to jitter by
jitter = spm_mb_io('GetScale',datn.f,sett);
jitter = reshape(jitter,[1 C]);

% GMM posterior
m  = datn.mog.po.m;
b  = datn.mog.po.b;
W  = datn.mog.po.W;
n  = datn.mog.po.n;

% GMM prior
m0 = datn.mog.pr.m;
b0 = datn.mog.pr.b;
W0 = datn.mog.pr.W;
n0 = datn.mog.pr.n;

% Lower bound
lb = datn.mog.lb;

if samp > 1
    % Subsample (runs faster, lower bound is corrected by scl_samp)
    
    % Image data        
    nobs0    = GetNumObsVox(fn,is_ct,4);
    [fn,d]   = SubSample(fn,Mn,samp);            
    nobs1    = GetNumObsVox(fn,is_ct,4);
    scl_samp = nobs0/nobs1; % get voxel ratio between original and subsamped image(s)
    
    % Template
    mun = SubSample(mun0,Mn,samp);      
    mun = reshape(mun,[prod(d(1:3)) K]);        
else
    d    = df;
    mun  = reshape(mun0,[prod(d(1:3)) K]);
    clear mun0
end
fn = reshape(fn,[prod(d(1:3)) C]);

% If labels are provided, use these
labels = GetLabels(datn,sett,true); % true -> subsample labels (if samp > 1)

% Missing data stuff
fn     = Mask(fn,is_ct);
msk_vx = ~isnan(fn);

if any(jitter~=0)
    % Data is an integer type, so to prevent aliasing in the histogram, small
    % random values are added.
    rng('default'); rng(1);
    fn = fn + bsxfun(@times,rand(size(fn)) - 1/2,jitter);
end

% Make K + 1 template
mun = spm_mb_shape('TemplateK1',mun,2);

% Add labels and template
mun    = mun + labels;
clear labels

% Expand, if using multiple Gaussians per tissue
mun = mun(:,mg_ix);

% Bias field related
if any(do_bf == true) 
    chan       = BiasFieldStruct(datn,C,df,reg,fwhm,[],datn.bf.T,samp);
    [bf,pr_bf] = BiasField(chan,d);
    bffn       = bf.*fn;
else
    bffn       = fn;
end

% Format for spm_gmm
[bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', bffn);
mun                       = spm_gmm_lib('obs2cell', mun, code_image, false);
code_list                 = unique(code_image);
code_list                 = code_list(code_list ~= 0);

for it_appear=1:nit_appear

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GMM and get responsibilities (zn)
    %------------------------------------------------------------
    
    [zn,mog,~,lb,mg_w] = spm_gmm_lib('loop',bffn,scl_samp,{{m,b},{W,n}},{'LogProp', mun}, ...
                                   'GaussPrior',   {m0,b0,W0,n0}, ...
                                   'Missing',      msk_chn, ...
                                   'LowerBound',   lb, ...
                                   'IterMax',      nit_gmm, ...
                                   'Tolerance',    tol_gmm, ...
                                   'SubIterMax',   nit_gmm_miss, ...
                                   'SubTolerance', tol_gmm, ...
                                   'Verbose',      [0 0], ...
                                   'MultGaussPi',  {mg_ix,mg_w});
    m = mog.MU;
    b = mog.b;
    W = mog.V;
    n = mog.n;           

    if it_appear > 1 && ((lb.sum(end) - lb.sum(end - 1))/abs(lb.sum(end)) > -eps('single')*10 || it_appear == nit_appear)
        % Finished
        break
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update bias field parameters
    % This computes the derivatives of the negative logarithm of the
    % joint probability distribution
    %------------------------------------------------------------

    if do_updt_bf && any(do_bf == true)        

        % Make sure to use the latest responsibilties         
        zn = Responsibility(m,b,W,n,bffn,ReWeightMu(mun,log(mg_w)),msk_chn);

        % Recompute parts of objective function that depends on bf
        lx  = LowerBound('ln(P(X|Z))',bffn,zn,msk_chn,{m,b},{W,n},scl_samp);
        lxb = scl_samp*LowerBound('ln(|bf|)',bf,msk_vx);                     

        for it_bf=1:nit_bf

            % Update bias field parameters for each channel separately
            for c=1:C % Loop over channels

                if ~datn.do_bf(c)
                    continue; 
                end

                % Compute gradient and Hessian (in image space)
                gr_im = zeros(d(1:3),'single');
                H_im  = zeros(d(1:3),'single');                    
                for l=1:size(msk_chn,1) % loop over combinations of missing voxels

                    % Get mask of missing modalities (with this particular code)    
                    ixo = msk_chn(l,:);         % Observed channels                    
                    ixm = ~ixo;                 % Missing channels
                    nm  = sum(ixm);              % Number of missing channels
                    if ~ixo(c), continue; end

                    % Convert channel indices to observed indices
                    ixc = 1:C; % mapped_c
                    ixc = ixc(ixo);
                    ixc = find(ixc == c);
                    
                    go = 0; % Gradient accumulated accross clusters
                    Ho = 0; % Hessian accumulated accross clusters
                    for k=1:Kmg
                        % Compute expected precision (see GMM + missing data)
                        Woo = W(ixo,ixo,k);
                        Wom = W(ixo,ixm,k);
                        Wmm = W(ixm,ixm,k);
                        Wmo = W(ixm,ixo,k);
                        Ao  = Woo - Wom*(Wmm\Wmo);
                        Ao  = (n(k) - nm) * Ao;
                        mo  = m(ixo,k);

                        % Compute statistics
                        gk = bsxfun(@minus, bffn{l}, mo.') * Ao(ixc,:).';
                        Hk = Ao(ixc,ixc);

                        gk   = bsxfun(@times, gk, zn{l}(:,k));
                        Hk   = bsxfun(@times, Hk, zn{l}(:,k));

                        % Accumulate across clusters
                        go = go + gk;
                        Ho = Ho + Hk;
                    end

                    % Multiply with bias corrected value (chain rule)
                    obffn = bffn{l}(:,ixc);
                    go    = go .* obffn;
                    Ho    = Ho .* (obffn.^2);
                    clear obffn
                    
                    % Add terms related to the normalisation (log(b))
                    go = go - 1;
                    Ho = Ho + 1; % Comes from subs(H,g,0)

                    % Accumulate across missing codes
                    ixvx        = (code_image == code_list(l));
                    gr_im(ixvx) = gr_im(ixvx) + go;
                    H_im(ixvx)  = H_im(ixvx)  + Ho;
                    clear ixvx
                end
                clear zn

                 % Compute gradient and Hessian (transform from image space to parameter space)
                d3 = numel(chan(c).T); % Number of DCT parameters
                H  = zeros(d3,d3);     
                gr = zeros(d3,1);      
                for z=1:d(3)
                    b3 = double(chan(c).B3(z,:)');
                    gr = gr + kron(b3,spm_krutil(double(gr_im(:,:,z)),double(chan(c).B1),double(chan(c).B2),0));
                    H  = H  + kron(b3*b3',spm_krutil(double(H_im(:,:,z)),double(chan(c).B1),double(chan(c).B2),1));
                end
                clear b3

                % Gauss-Newton update of bias field parameters
                Update = reshape((H + chan(c).ICO)\(gr + chan(c).ICO*chan(c).T(:)),size(chan(c).T));
                clear H gr

                % Line-search
                armijo = scal_bf;        
                oT     = chan(c).T;  
                opr_bf = pr_bf;
                olxb   = lxb;   
                olx    = lx;

                for ls=1:nit_lsbf

                    % Update bias-field parameters
                    chan(c).T = chan(c).T - armijo*Update;

                    % Compute new bias-field (only for channel c)
                    [bf,pr_bf] = BiasField(chan,d,bf,c,opr_bf);                        
                    bffn       = bf.*fn;
                    bffn       = spm_gmm_lib('obs2cell', bffn, code_image, true);
                    
                    % Recompute responsibilities (with updated bias field)
                    zn = Responsibility(m,b,W,n,bffn,ReWeightMu(mun,log(mg_w)),msk_chn);

                    % Compute new lower bound
                    lx  = LowerBound('ln(P(X|Z))',bffn,zn,msk_chn,{m,b},{W,n},scl_samp);            
                    lxb = scl_samp*LowerBound('ln(|bf|)',bf,msk_vx);

                    % Check new lower bound
                    if  ((lx + lxb + sum(pr_bf)) - (olx + olxb + sum(opr_bf)))/abs(lx + lxb + sum(pr_bf)) > -eps('single')*10
                        % Converged
%                         fprintf('it2=%i\tc=%i\tls=%i\tarmijo=%0.7f\tnl=%0.7f\tgain=%0.7f :o)\n',it_bf,c,ls,armijo,lx + lxb + sum(pr_bf),lx + lxb + sum(pr_bf) - olx + olxb + sum(opr_bf));
                        lb.XB(end + 1) = lxb;
                        lb.X(end  + 1) = lx;
                        break;
                    else                                
                        armijo    = armijo*0.5;
                        chan(c).T = oT;
                        if ls == nit_lsbf   
                            % Did not converge -> reset
%                             fprintf('it2=%i\tc=%i\tls=%i :o(\n',it_bf,c,ls);
                            lx    = olx;
                            lxb   = olxb;
                            bf    = BiasField(chan,d,bf,c,opr_bf);    
                            bffn  = bf.*fn;
                            bffn  = spm_gmm_lib('obs2cell', bffn, code_image, true);
                            pr_bf = opr_bf;
                            if nit_appear == 1 || nit_bf > 1
                                % Recompute responsibilities
                                zn = Responsibility(m,b,W,n,bffn,ReWeightMu(mun,log(mg_w)),msk_chn);
                            end
                        end
                    end
                end
                clear oT Update
            end
        end   

        % Update datn     
        datn.bf.T = {chan(:).T};
        datn.E(3) = -sum(pr_bf); % global objective function is negative log-likelihood..
        
        % Ensure correct lower bound
        lb.pr_bf(end + 1) = -datn.E(3);
    end
end
clear fn bf mun

% Update datn     
datn.E(1)     = -lb.sum(end);
datn.mog.po.m = m;
datn.mog.po.b = b;
datn.mog.po.W = W;
datn.mog.po.n = n;          
datn.mog.lb   = lb;  
datn.mog.mg_w = mg_w;

if nargout > 1
    % Compute full size resps
    
    if samp > 1
        % Compute responsibilities on original data
        fn = spm_mb_io('GetData',datn.f);
        fn = reshape(fn,[prod(df(1:3)) C]);
        fn = Mask(fn,is_ct);

        % Bias field
        if any(do_bf == true)
            % Get full-sized bias field
            chan = BiasFieldStruct(datn,C,df,reg,fwhm,[],datn.bf.T);
            bf   = BiasField(chan,df);
            bffn = bf.*fn;
            clear bf
        else
            bffn = fn;
        end
        clear fn
        [bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', bffn);    

        % Template
        mun0 = reshape(mun0,[prod(df(1:3)) K]);
        mun0 = spm_mb_shape('TemplateK1',mun0,2);

        % Integrate labels and multiple Gaussians per tissue
        labels = GetLabels(datn,sett);    
        mun0   = mun0 + labels;
        clear labels
        mun0   = mun0(:,mg_ix);
        mun0   = mun0 + log(mg_w);    

        % Compute full-size resps
        mun0 = spm_gmm_lib('obs2cell', mun0, code_image, false);            
        zn   = Responsibility(m,b,W,n,bffn,mun0,msk_chn); 
        clear mun0
    end       
    zn = spm_gmm_lib('cell2obs', zn, code_image, msk_chn);

    % Get K1 - 1 classes resp
    zn = zn(:,mg_ix <= K);

    % If using multiple Gaussians per tissue, collapse so that zn is of
    % size K
    if size(zn,2) > K
        for k=1:K, zn(:,k) = sum(zn(:,mg_ix==k),2); end
        zn(:,K + 1:end)    = [];
    end

    % % Fill in resps with no observations using template
    % for k=1:K, zn(msk_allmiss,k) = bg_mun(:,k); end
    % clear bg_mun msk_zn

    % Make 4D
    zn = reshape(zn,[df(1:3) K]);
end
end
%==========================================================================

%==========================================================================
% UpdateAllGMMs()
function dat = UpdateAllGMMs(dat, mu0, sett)

% Parse function settings
B     = sett.registr.B;
Mmu   = sett.var.Mmu;
mu_bg = sett.model.mu_bg;

for n=1:numel(dat)
    
    % Warp template to subject
    df  = spm_mb_io('GetSize',dat(n).f);
    q   = double(dat(n).q);
    Mn  = dat(n).Mat;
    psi = spm_mb_shape('Compose',spm_mb_io('GetData',dat(n).psi),spm_mb_shape('Affine',df, Mmu\spm_dexpm(q,B)*Mn));
    mu  = spm_mb_shape('Pull1',mu0,psi,mu_bg);
    clear psi
    
    % Update GMM parameters
    dat(n) = Update(dat(n),mu,sett);
end
end
%==========================================================================

%==========================================================================
% UpdatePrior()
function dat = UpdatePrior(dat, mu, sett, add_po_observation)
if nargin < 4, add_po_observation = false; end

% Parse function settings
mg_ix = sett.model.mg_ix;

% Make sure we use the latest GMM parameters
dat = UpdateAllGMMs(dat, mu, sett);

% Get population indices
p_ix = GetPopulationIdx(dat);

for p=1:numel(p_ix) % Loop over populations

    N = numel(p_ix{p}); % Number of subjects in population

    % Get old prior
    pr = dat(p_ix{p}(1)).mog.pr;
    pr = {pr.m,pr.b,pr.W,pr.n};
    C  = size(pr{1},1);

    % Get all posteriors
    K     = size(pr{1},2);
    K1    = numel(mg_ix);
    po    = cell(1,N);
    for n=1:N
        n1          = p_ix{p}(n);
        po{n}{1}{1} = dat(n1).mog.po.m;
        po{n}{1}{2} = dat(n1).mog.po.b;
        po{n}{2}{1} = dat(n1).mog.po.W;
        po{n}{2}{2} = dat(n1).mog.po.n;
    end

if add_po_observation
    % Add one artificial observation (increases numerical stability)
    
    % Get overall mean and variance for regularising
    avgmn = 0;
    sum_b = 0;
    avgvr = 0;
    sum_n = 0;
    for n=1:N
        pon = dat(p_ix{p}(n)).mog.po;
        for k=1:K1
            avgmn = avgmn + pon.m(:,k)*pon.b(k);
            avgvr = avgvr + inv(pon.W(:,:,k));
            sum_n = sum_n + pon.n(k);
            sum_b = sum_b + pon.b(k);
        end
    end
    avgvr = avgvr/sum_n;
    avgmn = avgmn/sum_b;
    avgpr = diag(1./diag(avgvr));

    % Add the artificial observation
    po1{1}{1} = repmat(avgmn,[1 K1]);     % m
    po1{1}{2} = zeros(1,K1) + 0.01;       % b
    po1{2}{1} = repmat(avgpr/C,[1 1 K1]); % W
    po1{2}{2} = C*ones(1,K1);             % n
    po{end+1} = po1;
end
 
    % Update prior
    pr = spm_gmm_lib('updatehyperpars',po,pr);

if false
    sum_m = 0;
    sum_b = 0;
    sum_P = 0;
    sum_n = 0;
    for k=1:K1
        sum_m = sum_m + pr{1}(:,k)*pr{2}(k);
        sum_P = sum_P + inv(pr{3}(:,:,k));
        sum_n = sum_n + pr{4}(k);
        sum_b = sum_b + pr{2}(k);
    end
    b_extra = 0;
    m_extra = sum_m/sum_b;
    n_extra = C;
    W_extra = C*diag(1./diag(sum_P/sum_n));              % Double check
    W_new   = W_extra;
    P_extra = inv(W_extra);
    for k=1:K1
        W_new(:,:,k) = inv(P_extra + inv(pr{3}(:,:,k))); % Double check
    end

    % Assign new prior
    for n=p_ix{p}
        dat(n).mog.pr.m = (pr{1} + m_extra)./reshape(pr{2} + b_extra,[1 K1]);
        dat(n).mog.pr.b =  pr{2} + b_extra;
        dat(n).mog.pr.W =  W_new;
        dat(n).mog.pr.n =  pr{4} + n_extra;
    end
end

    % Assign new prior
    for n=p_ix{p}
        dat(n).mog.pr.m = pr{1};
        dat(n).mog.pr.b = pr{2};
        dat(n).mog.pr.W = pr{3};
        dat(n).mog.pr.n = pr{4};
    end
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% ApplyMask()
function f = ApplyMask(f,is_ct)
if is_ct, f(~isfinite(f) | f == 0 | f < - 1020 | f > 3000) = NaN;
else,     f(~isfinite(f) | f == 0)                         = NaN;
end
end
%==========================================================================

%==========================================================================
% GetLabelConfMatrix()
function cm = GetLabelConfMatrix(cm_map,sett)
% FORMAT CM = get_label_cm(cm_map,opt)
% cm_map - Defines the confusion matrix
% sett   - Options structure
% cm     - confusion matrix
%
% Build Rater confusion matrix for one subject.
% This matrix maps template classes to manually segmented classes.
% Manual labels often do not follow the same convention as the Template, 
% and not all regions may be labelled. Therefore, a manual label may 
% correspond to several Template classes and, conversely, one Template
% class may correspond to several manual labels.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Here, we assume that all subjects from the same population (e.g.,
% a publicily available dataset) have the same labelling protocole and 
% confusion matrix.
% We allow the rater's sensitivity to change every few acquistion. We would
% typically start with a high sensitivity, to weight the labels strongly,
% and then decrease this value to allow the model to correct the rater's
% mistakes (especially near boundaries).

% Parse function settings
K = sett.model.K;
w = sett.labels.w;

K1 = K + 1;
L  = numel(cm_map); % Number of labels
cm = zeros([L K1],'single'); % Allocate confusion matrix
for l=1:L % Loop over labels    
    ix            = false(1,K1);
    ix(cm_map{l}) = true;   
    
    cm(l,ix)  = w/nnz(ix); 
    cm(l,~ix) = (1 - w)/nnz(~ix);
end
cm = bsxfun(@rdivide,cm,sum(cm,2));
end
%==========================================================================

%==========================================================================
% InitPopulation()
function dat = InitPopulation(dat,mu,pr,K,use_initgmm,w_mu,mg_ix,sett)

% Parse function settings
fwhm = sett.bf.fwhm;
reg  = sett.bf.reg;

% Parameters
N     = numel(dat);
K1    = K + 1;
Kmg   = numel(mg_ix);
[~,C] = spm_mb_io('GetSize',dat(1).f);

% Lower bound struct
lb  = struct('sum', NaN, 'X', [], 'XB', [], 'Z', [], 'P', [], 'MU', [], ...
             'A', [], 'pr_v', [], 'pr_bf',[]);
     
if ~use_initgmm
    % A simple scheme based on image statistics is used to init the GMM
    % parameters
    mx  = zeros(C,numel(dat));
    mn  = zeros(C,numel(dat));
    avg = zeros(C,numel(dat));
    vr  = zeros(C,numel(dat));
else
    % The InitGMM function is used to init GMM posterior and prior, as well
    % as the bias field DC scaling
    [po,pr,dc_all] = InitGMM(dat,sett);
end 

% Loop over subjects in population(s)
for n=1:N
    [df,C] = spm_mb_io('GetSize',dat(n).f);
    
    if ~use_initgmm                
        % Load image data
        fn = spm_mb_io('GetData',dat(n).f);
        fn = reshape(fn,[prod(df(1:3)) C]);
        fn = spm_mb_appearance('Mask',fn,dat(n).is_ct);
        
        % Set bias field DC component based on making images in
        % population closs to a mean value given by val
        if any(dat(n).do_bf == true)            
            if ~isempty(pr)        
                % Weighted mean based on template and prior mean                
                val = w_mu.*pr.m;    
                val = sum(val,2); 
            else                
                val = 1e3*ones(1,C);
            end
            
            % Make mean close to val
            dc = zeros(1,C);
            for c=1:C
                msk   = isfinite(fn(:,c));
                dc(c) = val(c)./mean(fn(msk,c));
            end
            dc = log(dc);
        end
    else
        dc = dc_all(:,n);
    end

    if any(dat(n).do_bf == true)
        % Get bias field parameterisation struct
        chan        = spm_mb_appearance('BiasFieldStruct',dat(n),C,df,reg,fwhm,dc);
        dat(n).bf.T = {chan(:).T};

        % Get bias field
        bf = spm_mb_appearance('BiasField',chan,df);
    else
        bf = ones([1 C],'single');
    end

    if ~use_initgmm
        % Modulate with bias field
        fn = bf.*fn;
    
        % Init GMM
        [po,mx(:,n),mn(:,n),avg(:,n),vr(:,n)] = InitSimplePosteriorGMM(dat(n),fn,mu,pr,K1,mg_ix,sett);
    end
    
    % Assign GMM
    mog.po     = po;
    mog.lb     = lb;
    mog.mg_w   = ones(1,Kmg)./arrayfun(@(x) sum(x == mg_ix), mg_ix);
    dat(n).mog = mog;
end

if isempty(pr)
    % Init GMM empirical prior
    pr = InitSimplePriorGMM(mx,mn,avg,vr,mu,K1);
end

% Assign prior
for n=1:numel(dat)                
    dat(n).mog.pr = pr; 
end

if K1 < Kmg && numel(pr.n) ~= Kmg
    % Modify posteriors and priors for when using multiple Gaussians per
    % tissue
    for n=1:N        
        
        is_ct = dat(n).is_ct;
        
        % Posterior
        po            = dat(n).mog.po;
        gmm           = spm_gmm_lib('extras', 'more_gmms', {po.m,po.b,po.W,po.n}, mg_ix);        
        if ~is_ct, gmm{1} = abs(gmm{1}); end % make sure non-negative
        po.m          = gmm{1};
        po.b          = gmm{2};
        po.W          = gmm{3};
        po.n          = gmm{4};
        dat(n).mog.po = po;
        
        % Prior
        pr            = dat(n).mog.pr;
        gmm           = spm_gmm_lib('extras', 'more_gmms', {pr.m,pr.b,pr.W,pr.n}, mg_ix);        
        if ~is_ct, gmm{1} = abs(gmm{1}); end % make sure non-negative
        pr.m          = gmm{1};
        pr.b          = gmm{2};
        pr.W          = gmm{3};
        pr.n          = gmm{4};
        dat(n).mog.pr = pr;
    end
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
% LowerBound()
function lb = LowerBound(type,varargin)
% Compute parts of the lower bound
%
% FORMAT lb = LowerBound('ln(|bf|)',bf,obs_msk)
% bf      - Exponentiated bias field [one channel]
% msk_vx  - Mask of observed values [one channel]
% lb      - Sum of the log bias field
%   >> This is part of the normalisation term in the GMM
%
% FORMAT lb = LowerBound('ln(P(X|Z))',fn,zn,code,mean,prec,W,L)
% fn       - Bias corrected observed image in matrix form [nbvox nbchannel]
% zn       - Responsibilities in matrix form [nbvox nbclass]
% msk_chn  - Mask of observed channels per code
% mean     - Mean parameters {m b}
% prec     - Precision parameters {W n}
% scl_samp - Image of weights [nbvox 1]
% lb       - Marginal log-likelihood
%   >> Marginal log-likelihood of the observed data, without the
%      bias-related normalisation.
if strcmpi(type,'ln(|bf|)')    
    bf     = varargin{1};
    msk_vx = varargin{2};
        
    lb = sum(log(bf(msk_vx)),'double');
elseif strcmpi(type,'ln(P(X|Z))')    
    fn       = varargin{1};
    zn       = varargin{2};
    msk_chn  = varargin{3};
    mean     = varargin{4};
    prec     = varargin{5};
    scl_samp = varargin{6};
    
    [lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', fn, zn, scl_samp, msk_chn);
    lb               = spm_gmm_lib('MarginalSum', lSS0, lSS1, lSS2, mean, prec, msk_chn);
else
    error('Undefined type!');
end
end
%==========================================================================

%==========================================================================
% InitGMM()
function dat = InitGMM(dat,sett)
% Code for initialising a GMM over lots of images of different channels and
% that can be scaled differently.
%
% FORMAT [po,pr,dc] = InitGMM(dat,sett)
% dat  - Structure holding data of N subjects
% sett - Structure of settings
%
% Adds dat.mog and dat.bf to the dat struct, for each subject. The mog
% field holds GMM parameters and the bf field the parameterisation of the
% bias field model. These parameters are estimated by this function.

% Parse function settings
dmu     = sett.dmu;
K       = sett.model.K;
Mmu     = sett.Mmu;     
fwhm    = sett.bf.fwhm;
reg     = sett.bf.reg;
fignam  = sett.show.figname_gmm;
figs    = sett.show.figs;
use_lab = sett.labels.use_initgmm;

% Parameters
tol     = 1e-4;  % Convergence tolerance
nit     = 256;   % Max number of iterations
nit_sub = 64;    % Max number of sub-iterations to update GMM mu and Sigma
do_dc   = true;
wp_reg  = 0.01;  % Regularises the GMM proportion by a percentage of the sampled number of voxels
verbose = any(strcmp(figs,'InitGMM'));

Ndat = numel(dat);   % Total number of subjects
K1   = K + 1;        % Number of Gaussians

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get number of populations, number of subjects of each population, and
% number of channels
%------------------------------------------------------------

pop_id  = []; % vector of populations ids [1,2,3,...]
pop_chn = []; % vector with number of image channels in each population
pop_cnt = []; % vector with count of number of subjects in each population
pop_ix  = {}; % cell array, where each element is a vector with the dat struct indices of a population
for n=1:Ndat
    if ~any(pop_id == dat(n).ix_pop) 
        pop_id = [pop_id dat(n).ix_pop];
        
        cnt1 = 0;
        nix  = [];
        for n1=1:Ndat
            if dat(n1).ix_pop == dat(n).ix_pop 
                cnt1 = cnt1 + 1;
                nix  = [nix n1];
            end
        end
        pop_ix{end + 1} = nix;
        pop_cnt         = [pop_cnt cnt1];
        
        [~,C]   = spm_mb_io('GetSize',dat(n).f); 
        pop_chn = [pop_chn C];
    else
        continue
    end
end

% cell array. where each element is a vector that maps that populations
% channels to...
pop_cnh_ix = cell([1 numel(pop_chn)]);
for i=1:numel(pop_chn)
    if i == 1, pop_cnh_ix{i} = 1:pop_chn(i);
    else,      pop_cnh_ix{i} = sum(pop_chn(1:i - 1)) + 1:sum(pop_chn(1:i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init model
%------------------------------------------------------------

% Size of model data
N = max(pop_cnt); % Number of population observations
C = sum(pop_chn); % Number of channels

% We sample a subset of voxels for each image, these points are chosen in
% template space, and then projected to each image's space, so that voxel
% locations corresponds. 'sampmu' decides on the distance between sampled 
% voxels in template space (smaller -> uses more memory, but better estimates)
sampmu       = 3;
vxmu         = sqrt(sum(Mmu(1:3,1:3).^2));
sk           = max([1 1 1],round(sampmu*[1 1 1]./vxmu));
sk(dmu == 1) = 1; % defines sampling grid

% Number of voxels to sample
Nvx        = min(round(prod(dmu(1:3))/N), round(prod(dmu(1:3)./sk))); 
r          = randperm(prod(dmu(1:3)),Nvx);
[x1,x2,x3] = ind2sub(dmu(1:3),r(:)); % voxel indicies in template space
ymu        = cat(2,single(x1),single(x2),single(x3),ones([Nvx 1],'single')); 
clear r x1 x2 x3

% Get sample of image data (and labels, if provided)
F = cell(1,N); % holds imaging data
L = cell(1,N); % holds possible label data
for n=1:N, L{n} = zeros([K1 1],'single'); end

labels_given = false;

for n=1:N % Loop over subjects
    
    fn    = NaN([C Nvx],'single');    
    cnt_c = 1; % channel count
    cnt_l = single(0); % label count
    for c=1:numel(pop_ix)
        if n > numel(pop_ix{c})
            % Population has no more images
            n1      = pop_ix{c}(1);
            [~,C1]  = spm_mb_io('GetSize',dat(n1).f);            
            cnt_c = cnt_c + C1;
            continue; 
        end
        
        % Parameters
        n1      = pop_ix{c}(n);
        [df,C1] = spm_mb_io('GetSize',dat(n1).f);
        is_ct   = dat(n1).is_ct;  
        
        % Get image data
        f1 = spm_mb_io('GetData',dat(n1).f);        
        if any(is_ct == true), f1(f1 < -1020 | f1 > 3000) = 0; end
        
        % Move template space sample points to subject space        
        Mn = dat(n1).Mat;  
        M  = Mn\Mmu;
        yf = ymu*M';       
        yf = reshape(yf(:,1:3),[Nvx 1 1 3]);        
        if df(3) == 1, yf(:,:,:,3) = 1; end
        
        % Get sample of image(s)        
        f1 = spm_diffeo('pull',f1,yf);        
        f1 = reshape(f1,[Nvx C1]);          
        
        if ~is_ct
            % Set leq to zero as NaN
            f1(f1 <= 0) = NaN;
        end
    
        for c1=1:C1
            fn(cnt_c,:) = f1(:,c1);  
            cnt_c       = cnt_c + 1;
        end
        
        % Deal with (possible) labels
        l1 = spm_mb_appearance('GetLabels',dat(n1),sett);            
        if use_lab && size(l1,1) > 1
            labels_given = true;
            
            l1      = reshape(l1,[df K1]);
            l1      = spm_diffeo('pull',l1,yf);
            l1      = reshape(l1,[Nvx K1]);      
            msk     = isnan(l1);
            l1(msk) = 0;             % we don't want NaN in the transformed labels
            L{n}    = L{n} + l1';
            cnt_l   = cnt_l + ~msk'; % for normalising labels across populations
            clear msk
        else                
            L{n} = L{n} + l1';
        end
    end
        
    if numel(cnt_l) > 1
        % Makes sure that labels sum to one (if more than one population of labels)
        cnt_l(cnt_l == 0) = 1;
        L{n}              = L{n}./cnt_l;
        clear cnt_l
    end 
        
    % Add to F array
    F{n} = fn;     
end
clear fn d1 mask l1 yf ymu

% Init bias field DC component
dc = zeros(C,N);
mn = zeros(C,1);
mx = zeros(C,1);
for n=1:N % Loop over subjects
    fn         = F{n};
    mn         = min(mn,min(fn,[],2,'omitnan'));    
    mx         = max(mx,max(fn,[],2,'omitnan'));
    fn(mn<0,:) = NaN; % For CT so to not optimise dc
    dc(:,n)    = -log(mean(fn,2,'omitnan'));
end
            
% Make DC component zero mean, across N
for c=1:C
    dcn     = dc(c,:);
    msk_dcn = isfinite(dcn);
    if all(msk_dcn == 0) || sum(msk_dcn) == 1, continue; end
    
    mn_dcn  = mean(dcn(msk_dcn));
    dc(c,:) = dc(c,:) - mn_dcn;
end

% Init GMM parameters
gam = ones(1,K1)./K1;

if labels_given
    mu = rand(C,K1).*(mx + mn); 
else
    mu = zeros(C,K1);
    for c=1:C
        rng     = linspace(mn(c),mx(c),K1);
        rng     = -sum(rng<0):sum(rng>=0) - 1;
        mu(c,:) = rng'*mx(c)/(1.0*K1);
    end    
end

Sig = diag(((mx + mn)/(2*K1)).^2).*ones([1,1,K1]);

% Compute precision
prec = zeros(size(Sig));
for k=1:K1, prec(:,:,k) = inv(Sig(:,:,k)); end
    
% Regularise the GMM proportion (as in spm_preproc8)
wp_reg = wp_reg*Nvx;  

% Get combinations of missing data, and the number of such combinations (Cm)
fn             = [];
for n=1:N, fn  = [fn; F{n}']; end
[~,~,msk_chnm] = spm_gmm_lib('obs2cell',fn);
clear fn
Cm             = size(msk_chnm,1);

% Cast images and labels into spm_gmm_lib format
for n=1:N
    fn                      = F{n}'; 
    [fn,code_image,msk_chn] = spm_gmm_lib('obs2cell',fn);
    F{n}                    = {fn,msk_chn};
    if verbose > 0, F{n}{end + 1} = code_image; end % store code_image so that we can use the image dat for final visualisation
    ln = L{n};
    ln = ln';
    if size(ln,1) > 1            
        L{n} = spm_gmm_lib('obs2cell', ln, code_image, false);
    end 
end
clear fn code_image msk_chn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit model
%------------------------------------------------------------

if ~do_dc, dc = zeros(size(dc)); end

ll = -Inf;
for it=1:nit
    oll = ll(it);
    
    if verbose > 1 && (mod(it,20) == 0 || it == 1)
        % Visualise
        figure(665); clf
        nr  = floor(sqrt(2*C));
        nc  = ceil(2*C/nr);  
        for c=1:C        
            x = linspace(mn(c),mx(c),1000);
            p = 0;            
            subplot(nr,nc,c);        
            hold on
            for k=1:K1
                pk = exp(log(gam(k)) - 0.5*log(Sig(c,c,k)) - 0.5*log(2*pi) - 0.5*(x-mu(c,k)).^2/(Sig(c,c,k)));
                p  = p + pk;
                plot(x,pk,'b-','LineWidth',1);
            end
            plot(x,p,'r-','LineWidth',3);
            hold off

            subplot(nr,nc,C + c);          
            bar(exp(dc(c,:)));
        end
        drawnow;       
    end
           
    % For summing suffstats across subjects
    SS0m = cell(1,Cm);
    SS1m = cell(1,Cm);
    SS2m = cell(1,Cm);    
    for c=1:Cm        
        Co      = sum(msk_chnm(c,:));                                        
        SS0m{c} = zeros(1,K1);
        SS1m{c} = zeros(Co,K1);
        SS2m{c} = zeros(Co,Co,K1);
    end
    
    %---------------------
    % Update GMM parameters (mu,Sig)
    %---------------------
    
    sll = 0;
    for n=1:N % Loop over subjects
        
        % Get image data
        fn      = F{n}{1};
        msk_chn = F{n}{2};        
        Cmn     = size(msk_chn,1);
        
        % Modulate with scaling parameter
        dcn                 = dc(:,n);
        dcn(~isfinite(dcn)) = 0;   
        for c=1:Cmn
            fn{c} = fn{c}.*exp(dcn(msk_chn(c,:)))';
        end 
        
        % Get labels and add bias normilisation term
        ln = L{n};
        if iscell(ln)
            for c=1:Cmn
                ln{c} = ln{c}  + sum(dcn);
            end 
        else
            ln = ln + sum(dcn);
            ln = ln';
        end              
                
        % Compute responsibilities
        norm_term        = spm_gmm_lib('normalisation',mu, prec, msk_chn);
        logpX            = spm_gmm_lib('marginal',fn, [{mu} prec], norm_term, msk_chn);
        zn               = spm_gmm_lib('responsibility',logpX, log(gam(:)'), ln);
        [SS0n,SS1n,SS2n] = spm_gmm_lib('suffstat','base',fn, zn, 1, msk_chn); 
        clear norm_term logpX fn
        
        % Lower bound
        sll = sll + spm_gmm_lib('marginalsum',SS0n, SS1n, SS2n, mu, prec, msk_chn);                
        sll = sll + spm_gmm_lib('kl','categorical',zn, 1, log(gam(:)'), ln);
        clear zn ln
        
        % Add nth suffstats to total suffstats
        for c=1:Cmn    
            Co        = msk_chn(c,:);
            [~,index] = ismember(Co,msk_chnm,'rows');
                                      
            for k=1:K1
                SS0m{index}(k)     = SS0m{index}(k)     + SS0n{c}(k);
                SS1m{index}(:,k)   = SS1m{index}(:,k)   + SS1n{c}(:,k);
                SS2m{index}(:,:,k) = SS2m{index}(:,:,k) + SS2n{c}(:,:,k);
            end
        end
    end
    ll = [ll sll];
    
    for it1=1:nit_sub               
                        
        % Infer missing suffstat
        [SS0i,SS1i,SS2i] = spm_gmm_lib('suffstat','infer',SS0m, SS1m, SS2m, {mu,prec}, msk_chnm);        

        % Update GMM
        SS0i = SS0i';        
        mu   = SS1i./SS0i';
        Sig0 = (sum(SS2i,3) - (SS0i'.*mu)*mu')/sum(SS0i);
        for k=1:K1
            Sig(:,:,k)  = (SS2i(:,:,k) - SS0i(k)*mu(:,k)*mu(:,k)' + C*Sig0)./(SS0i(k) + C);
        end
                    
        % Update precision
        oprec = prec;
        prec  = zeros(size(Sig));
        for k=1:K1, prec(:,:,k) = inv(Sig(:,:,k)); end
        
        for k=1:K1
            [~,cholp] = chol(prec(:,:,k));
            if cholp ~= 0
                prec(:,:,k) = oprec(:,:,k);
            end            
        end
    end
    clear SS0m SS1m SS2m
    
    % Update proportion
    gam = (SS0i + wp_reg)./(sum(SS0i) + wp_reg*K1);
    
    if do_dc
        %---------------------
        % Update rescaling (s)
        %---------------------

        sll = 0;
        for n=1:N % Loop over subjects

            % Get image data
            fn      = F{n}{1};
            msk_chn = F{n}{2};        
            Cmn     = size(msk_chn,1);

            % Modulate with scaling parameter
            dcn                 = dc(:,n);
            dcn(~isfinite(dcn)) = 0;   
            for c=1:Cmn
                fn{c} = fn{c}.*exp(dcn(msk_chn(c,:)))';
            end 

            % Get labels and add bias normilisation term
            ln = L{n};
            if iscell(ln)
                for c=1:Cmn
                    ln{c} = ln{c}  + sum(dcn);
                end 
            else
                ln = ln + sum(dcn);
                ln = ln';
            end                

            % Compute responsibilities
            norm_term        = spm_gmm_lib('normalisation',mu, prec, msk_chn);
            logpX            = spm_gmm_lib('marginal',fn, [{mu} prec], norm_term, msk_chn);
            zn               = spm_gmm_lib('responsibility',logpX, log(gam(:)'), ln);   
            [SS0n,SS1n,SS2n] = spm_gmm_lib('suffstat','base',fn, zn, 1, msk_chn); 
            clear norm_term logpX

            % Lower bound
            sll = sll + spm_gmm_lib('marginalsum',SS0n, SS1n, SS2n, mu, prec, msk_chn);                
            sll = sll + spm_gmm_lib('kl','categorical',zn, 1, log(gam(:)'), ln);

            % Compute gradient and Hessian
            msk_dc = isfinite(dc(:,n))';
            g      = zeros(C,1);
            H      = zeros(C,C);        
            for l=1:Cmn % loop over combinations of missing channels                        
                ixo = msk_chn(l,:) & msk_dc;
                if all(ixo == 0), continue; end

                mskfc                    = ixo;
                mskfc(msk_chn(l,:) == 0) = [];          

                fnc = fn{l}(:,mskfc);            
                fnc = fnc';                       
                znc = zn{l};
                Ns  = size(fnc,2);

                go = 0;
                Ho = 0;
                for k=1:K1
                    rk = znc(:,k);
                    rk = rk';

                    muk = mu(ixo,k);
                    Sk  = Sig(ixo,ixo,k);

                    go = go + Sk\(rk.*(fnc - muk));
                    Ho = Ho + ((rk.*fnc)*fnc').*inv(Sk);
                end
                go = sum(fnc.*go,2) - Ns;
                Ho = Ho + Ns*eye(nnz(ixo));

                g(ixo)     = g(ixo)     + go;
                H(ixo,ixo) = H(ixo,ixo) + Ho;
            end        
            clear zn fn

            % Gauss-Newton update     
            dcn(dcn == 0) = NaN;
            dcn(msk_dc)   = dcn(msk_dc) - ((N - 1)/N)*(H(msk_dc,msk_dc)\g(msk_dc));            
            dc(:,n)       = dcn;
        end
        ll = [ll sll];

        % Make DC component zero mean, across N
        for c=1:C
            dcn     = dc(c,:);
            msk_dcn = isfinite(dcn);
            if all(msk_dcn == 0) || sum(msk_dcn) == 1, continue; end

            mn_dcn  = mean(dcn(msk_dcn));
            dc(c,:) = dc(c,:) - mn_dcn;
        end
    end
    
    if verbose > 1 && mod(it,20) == 0      
        % Visualise
        figure(666);
        subplot(121)
        plot(1:numel(ll),ll,'b-','LineWidth',2)
        titll = sprintf('it=%i, ll - oll=%0.6f',it,(ll(end) - oll)/abs(ll(end) + oll));
        title(titll)
        subplot(122)
        bar(gam)
        drawnow;       
    end    

    if (ll(end) - oll)/abs(ll(end) + oll) < tol
        % Finished
        break; 
    end
end

if verbose > 0
    % Visualise
    f  = findobj('Type', 'Figure', 'Name', fignam);
    if isempty(f), f = figure('Name', fignam, 'NumberTitle', 'off'); end
    set(0, 'CurrentFigure', f);   

    nr = floor(sqrt((C + 2)));
    nc = ceil((C + 2)/nr);  
    
    % GMM fit
    for c=1:C
        subplot(nr,nc,c)                
        for n=1:N            
            fn  = spm_gmm_lib('cell2obs', F{n}{1}, F{n}{3}, F{n}{2}); 
            fn  = fn';
            fnc = fn(c,:);
            fnc = fnc(isfinite(fnc));
            if isempty(fnc), continue; end
            
            dcn                 = dc(c,n);
            dcn(~isfinite(dcn)) = 0;   
        
            xn  = fnc.*exp(dcn);
            x   = min(xn):20:max(xn);
            h   = hist(xn,x);
            h   = h/sum(h)/20;
            plot(x,h,'k.','MarkerSize',1);
            hold on
        end
        x = linspace(min(xn),max(xn),1000);
        p = 0;
        for k=1:K1
            pk = exp(log(gam(k)) - 0.5*log(Sig(c,c,k)) - 0.5*log(2*pi) - 0.5*(x-mu(c,k)).^2/(Sig(c,c,k)));
            p  = p + pk;
            plot(x,pk,'b-','LineWidth',1);
        end
        plot(x,p,'r-','LineWidth',3);
        hold off
        title(['C=' num2str(c)])
    end
    
    % Log-likelihood
    subplot(nr,nc,c + 1)
    plot(1:numel(ll),ll,'b-','LineWidth',2)
    title('ll');
    
    % GMM prop
    subplot(nr,nc,c + 2)
    bar(gam)
    title('prop')
    
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make function output
%------------------------------------------------------------

% Lower bound struct
lb  = struct('sum', NaN, 'X', [], 'XB', [], 'Z', [], 'P', [], 'MU', [], ...
             'A', [], 'pr_v', [], 'pr_bf',[]);
         
% Set GMM and bias field parameters for all subjects
for n=1:N
    for c=1:numel(pop_ix)
        if n > numel(pop_ix{c}), continue; end
        
        n1      = pop_ix{c}(n);
        [df,C1] = spm_mb_io('GetSize',dat(n1).f);                        
        p       = dat(n1).ix_pop;    
        chn     = pop_cnh_ix{p};
        
        m     = mu(chn,:);
        b     = zeros(1,K1) + 0.01;
        ico   = zeros(C1,C1,K1);
        mnico = zeros(C1,C1);
        for k=1:K1
            for c1=1:C1
                ico(c1,c1,k) = inv(Sig(chn(c1),chn(c1),k));
                mnico(c1,c1) = mnico(c1,c1) + Sig(chn(c1),chn(c1),k);
            end
        end
        mnico = inv(mnico/K);
        
        W  = ico/C1;
        nu = C1*ones(1,K1);
        po = struct('m',m,'b',b,'W',W,'n',nu);    

        mog.po      = po;
        mog.pr      = po;
%         mog.pr      = struct('m',repmat(mean(m,2),[1 K1]),'b',b,'W',repmat(mnico/C1,[1 1 K1]),'n',nu); % uninformative prior
        mog.lb      = lb;
        mog.mg_w    = ones([1 K1]);
        dat(n1).mog = mog;
        
        if any(dat(n1).do_bf == true)
            chan         = spm_mb_appearance('BiasFieldStruct',dat(n1),C1,df,reg,fwhm,dc(chn,n));
            dat(n1).bf.T = {chan(:).T};
        end
    end
end
end
%==========================================================================

%==========================================================================    
% InitSimplePriorGMM()
function pr = InitSimplePriorGMM(mx,mn,avg,vr,mu0,K1)
% Initialise the prior parameters of a Gaussian mixture model. 
%
%   This function is only used to initialise the GMM parameters at the
%   beginning of the algorithm.
% 
% FORMAT pr = InitSimplePriorGMM(mx,mn,vr,mu0,K)
% mx   - Maximum observed value [per channel]
% mn   - Minimum observed value [per channel]
% avg  - Mean observed value [per channel]
% vr   - Variance of observed value [per channel]
% mu0  - Log template
% K1   - Number of classes
% pr   - Structure holding prior GMM parameters (m, b, W, n)

mvr  = mean(vr,2); % mean variance across all subjects in population
mmx  = mean(mx,2); % mean mean across all subjects in population
mmn  = mean(mn,2); % mean mean across all subjects in population
mavg = mean(avg,2); % mean mean across all subjects in population

C   = size(mx,1);
m   = zeros(C,K1);
ico = zeros(C,C,K1);   
for c=1:C      
    rng    = linspace(mmn(c),mmx(c),K1);
    rng    = -sum(rng<0):sum(rng>=0) - 1;
    m(c,:) = rng'*mmx(c)/(1.0*K1);
    
    ico(c,c,:) = (mmx(c)/(1.0*K1))^2;
    ico(c,c,:) = 1/ico(c,c,:); % precision
end

if ~isempty(mu0)
    % If template is given, make sure that pr.m is the same for all classes
    m = repmat(mean(m,2),[1 K1]);
end

% Define prior
pr   = struct('m',[],'b',[],'n',[],'W',[]);
pr.m = m;
pr.b = zeros(1,K1) + 0.01;
pr.n = C*ones(1,K1);
pr.W = ico/C;

if 0
    fig_name = 'Prior';
    spm_gmm_lib('plot','gaussprior',{pr.m,pr.b,pr.W,pr.n},[],fig_name);
end

end
%==========================================================================   

%==========================================================================    
% InitSimplePosteriorGMM()
function [po,mx,mn,avg,vr] = InitSimplePosteriorGMM(datn,fn,mu,pr,K1,mg_ix,sett)
% Initialise the posterior parameters of a Gaussian mixture model. 
%
%   This function is only used to initialise the GMM parameters at the
%   beginning of the algorithm.
% 
% FORMAT [po,mx,mn,vr] = InitSimplePosteriorGMM(datn,fn,mu,pr,K,sett)
% datn - Structure holding data of a single subject
% fn   - Bias corrected observed image, in matrix form [nbvox nbchannel]
% mu   - Log template
% pr   - Structure holding prior GMM parameters (m, b, W, n)
% K1   - Number of classes
% sett - Structure of settings
% po   - Structure of posterior parameters (m, b, W, n)
% mx   - Maximum observed value [per channel]
% mn   - Minimum observed value [per channel]
% avg  - Mean observed value [per channel]
% vr   - Variance of observed value [per channel]

% Parse function settings
B   = sett.registr.B;
Mmu = sett.Mmu;

C   = size(fn,2);
mx  = zeros(C,1);
mn  = zeros(C,1);
avg = zeros(C,1);
vr  = zeros(C,1);
m   = zeros(C,K1);
ico = zeros(C,C,K1);
for c=1:C
    msk    = isfinite(fn(:,c));
    mx(c)  = max(fn(msk,c));
    mn(c)  = min(fn(msk,c));
    avg(c) = mean(fn(msk,c));
    vr(c)  = var(fn(msk,c));
    
    rng    = linspace(mn(c),mx(c),K1);
    rng    = -sum(rng<0):sum(rng>=0) - 1;
    m(c,:) = rng'*mx(c)/(1.0*K1);
    
    ico(c,c,:) = (mx(c)/(1.0*K1))^2;
    ico(c,c,:) = 1/ico(c,c,:); % precision
end

% Initialise posterior
if isempty(pr)
    % No prior given, initialise from image statistics (mean, min, max, var)
    po   = struct('m',[],'b',[],'n',[],'W',[]);
    po.m = m;
    po.b = zeros(1,K1) + 0.01;
    po.n = C*ones(1,K1);
    po.W = ico/C;
else
    % Use given prior
    po.m = pr.m;
    po.b = pr.b;
    po.n = pr.n;
    po.W = pr.W;
end

if (isempty(pr) && ~isempty(mu))
    % Template is given, but not prior, so make sure that m parameter is the
    % same for all classes
    po.m = repmat(mean(po.m,2),[1 K1]);
end
   
if ~isempty(mu)
    
    mu_bg = sett.model.mu_bg;
    
    % Use template as resposibilities to compute values for GMM
    % posterior from suffstats
    df = spm_mb_io('GetSize',datn.f);          
    q  = double(datn.q);
    Mr = spm_dexpm(q,B);
    Mn = datn.Mat;    

    % Warp template
    mu = spm_mb_shape('Pull1',mu,spm_mb_shape('Affine',df,Mmu\Mr*Mn),mu_bg);                

    % Add class, then softmax -> can now be used as resps
    mu = spm_mb_shape('TemplateK1',mu,4);
    mu = exp(mu);
    mu = reshape(mu,[prod(df(1:3)) K1]);        
        
    if K1 < numel(pr.m)
        % Prior has more Gaussians then tissue classes, incorporate this
        % into template model
        mu = mu(:,mg_ix);
    end
    
    % Compute posterior estimates
    [SS0,SS1,SS2] = spm_gmm_lib('SuffStat',fn,mu,1);
    [MU,~,b,W,n]  = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2, {po.m,po.b,po.W,po.n});
    
    po.m = MU;
    po.b = b;%zeros(1,K1) + 0.01;
%     A    = po.W.*reshape(po.n,[1 1 K1]);
    po.n = n;%C*ones(1,K1);    
    po.W = W;%A/C;
end

if 0
    fig_name = 'Posterior';
    spm_gmm_lib('plot','gaussprior',{po.m,po.b,po.W,po.n},[],fig_name);
end
end
%==========================================================================  

%==========================================================================
% ReWeightMu()
function mun = ReWeightMu(mun,logmg_w)
if sum(logmg_w) == 0, return; end
for i=1:numel(mun)
    mun{i} = mun{i} + logmg_w; 
end
end
%==========================================================================

%==========================================================================
% SubSample()
function [ofn,d,w] = SubSample(fn,Mn,samp)
% Subsample a multichannel volume.
%
% FORMAT [ofn,d,scl_samp] = SubSample(fn,Mn,samp,deg,bc)
% fn   - Original volume
% Mn   - Original orientation matrix
% samp - Sampling distance in mm
% ofn  - Resampled volume
% d    - Output dimensions
% w    - Proportion of sampled voxels

% Input image properties
vx = sqrt(sum(Mn(1:3,1:3).^2));
df = size(fn);
df = [df 1];

% Output image properties
samp = max([1 1 1],round(samp*[1 1 1]./vx));
D    = diag([1./samp 1]);
d    = ceil(D(1:3,1:3)*df(1:3)')'; % New dimensions

ofn = fn(1:samp(1):end,1:samp(2):end,1:samp(3):end,:); 

% For weighting data parts of lowerbound with factor based on amount of
% downsampling  
w = prod(df(1:3))/prod(d(1:3));
end
%==========================================================================

%==========================================================================
% TransformBF()
function t = TransformBF(B1,B2,B3,T)
% Create an image-space log bias field from its basis function encoding.
%
% FORMAT t = TransformBF(B1,B2,B3,T)
% B1 - x-dim DCT basis [nx kx]
% B2 - y-dim DCT basis [ny ky]
% B3 - z-dim DCT basis [nz kz]
% T  - DCT encoding of the log bias field [kx ky kz]
% t  - Reconstructed log bias field [nx ny nz]
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
end
%==========================================================================