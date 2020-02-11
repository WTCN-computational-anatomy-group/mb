function varargout = spm_mb_appearance(varargin)
%__________________________________________________________________________
%
% Functions for appearance model related.
%
% FORMAT chan       = spm_mb_appearance('BiasBasis',T,df,Mat,reg,samp)
% FORMAT [bfn,lln]  = spm_mb_appearance('BiasField',T,chan,d,varargin)
% FORMAT labels     = spm_mb_appearance('GetLabels',datn,sett,do_samp)
% FORMAT [nvx_obs,msk_allmiss] = spm_mb_appearance('GetNumObsVox',f,is_ct,ax)
% FORMAT p_ix       = spm_mb_appearance('GetPopulationIdx',dat)
% FORMAT [dat,sett] = spm_mb_appearance('IntroduceMG',dat,sett)
% FORMAT zn         = spm_mb_appearance('Responsibility',m,b,W,n,fn,mu,msk_chn)
% FORMAT [zn,datn]  = spm_mb_appearance('Update',datn,mun0,sett)
% FORMAT dat        = spm_mb_appearance('UpdatePrior',dat,sett)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_appearance
    error('Not enough argument. Type ''help spm_mb_appearance'' for help.');
end
id       = varargin{1};
varargin = varargin(2:end);
switch id
    case 'BiasBasis'
        [varargout{1:nargout}] = BiasBasis(varargin{:});
    case 'BiasField'
        [varargout{1:nargout}] = BiasField(varargin{:});
    case 'GetLabels'
        [varargout{1:nargout}] = GetLabels(varargin{:});
    case 'GetNumVoxObserved'
        [varargout{1:nargout}] = GetNumObsVox(varargin{:});
    case 'GetPopulationIdx'
        [varargout{1:nargout}] = GetPopulationIdx(varargin{:});
    case 'IntroduceMG'
        [varargout{1:nargout}] = IntroduceMG(varargin{:});
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
function [bfn,lln] = BiasField(T,chan)
d  = [size(chan(1).B1,1) size(chan(1).B2,1) size(chan(1).B3,1)];
nz = d(3);
C  = numel(T);
cr = 1:C;

% Compute full bias field (for all channels)
bfn = zeros([d C],'single');
lln = zeros(1,C);

for c=cr
    lln(c) = double(-0.5*T{c}(:)'*chan(c).ICO*T{c}(:));
    for z=1:nz
        bf_c         = TransformBF(chan(c).B1,chan(c).B2,chan(c).B3(z,:),T{c});
        bfn(:,:,z,c) = single(exp(bf_c));
    end
end
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

%==========================================================================
function chan = BiasBasis(T,df,Mat,reg,samp)
if nargin<5, samp = 0; end
cl   = cell(1, numel(T));
chan = struct('ICO', cl, 'B1',cl, 'B2',cl, 'B3',cl);
ind  = SampleInd(df,Mat,samp);
vx   = sqrt(sum(Mat(1:3,1:3).^2,1));
for c=1:numel(T)
    d3 = [size(T{c}) 1];
    d3 = d3(1:3);

    % GAUSSIAN REGULARISATION for bias correction
    chan(c).ICO = reg*spm_bias_lib('regulariser','bending',df,d3,vx);

    % Basis functions for bias correction
    chan(c).B1 = spm_dctmtx(df(1),d3(1),ind{1});
    chan(c).B2 = spm_dctmtx(df(2),d3(2),ind{2});
    chan(c).B3 = spm_dctmtx(df(3),d3(3),ind{3});
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
    Mn     = datn.Mat;
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
% GetNumObsVox()
function [nvx_obs,msk_allmiss] = GetNumObsVox(fn,ax)
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
% IntroduceMG()
function [dat,sett] = IntroduceMG(dat,sett)
% Modify posteriors and priors, if using multiple Gaussians per tissue

% Parse function settings
mg_ix  = sett.model.mg_ix_intro;

K1 = numel(dat(1).mog.po.n);
K  = K1 - 1;
if isscalar(mg_ix)
    mg_ix = repelem(1:K + 1,mg_ix);
end
sett.model.mg_ix = mg_ix;

% Get parameters
N   = numel(dat);
Kmg = numel(mg_ix);
K1  = K + 1;

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
Kmg        = numel(mg_ix);
Mn         = datn.Mat;
scl_samp   = 1; % sampling scaling (defined as W = prod(d0(1:3))/prod(d(1:3)), when samp > 1)
do_bf      = datn.do_bf;
mg_w       = datn.mog.mg_w;

% Get image data
fn0 = spm_mb_io('GetImage',datn);

% For visualising results (spm_mb_show)
spm_mb_io('Write2Visualise',datn,mun0,'mu',sett);

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
jitter = reshape(jitter,[1 1 1 C]);

% GMM posterior
m  = datn.mog.po.m;
b  = datn.mog.po.b;
W  = datn.mog.po.W;
n  = datn.mog.po.n;

% Lower bound
lb = datn.mog.lb;

if samp > 0 
    % Subsample (runs faster, lower bound is corrected by scl_samp)

    % Image data
    nobs0    = GetNumObsVox(fn0,4);
    [fn,d]   = SubSample(fn0,Mn,samp);
    nobs1    = GetNumObsVox(fn,4);
    scl_samp = nobs0/nobs1; % get voxel ratio between original and subsamped image(s)

    % Template
    mun = SubSample(mun0,Mn,samp);
    mun = reshape(mun,[prod(d(1:3)) K]);
else
    d    = df;
    mun  = reshape(mun0,[prod(d(1:3)) K]);
    clear mun0
    fn   = fn0;
end

% If labels are provided, use these
labels = GetLabels(datn,sett,true); % true -> subsample labels (if samp > 1)

% Missing data stuff
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
mun = mun + labels;
clear labels

% Expand, if using multiple Gaussians per tissue
mun = mun(:,mg_ix);

% Bias field related
if any(do_bf == true)
    T          = datn.T;
    chan       = BiasBasis(T,df,Mn,reg,samp);
    [bf,llpbf] = BiasField(T,chan);
    bffn       = bf.*fn;    
else
    bffn = fn;
end

% For visualising results (spm_mb_show)    
if samp == 1
    spm_mb_io('Write2Visualise',datn,bffn,'bff',sett);
    if any(do_bf == true)
        spm_mb_io('Write2Visualise',datn,fn,'f',sett);
        spm_mb_io('Write2Visualise',datn,bf,'bf',sett);
    end
end

if isempty(lb.XB)
    % Make sure bias field part of lower bound is correct
    lxb            = scl_samp*LowerBound('ln(|bf|)',vol2vec(bf),msk_vx) + sum(llpbf);
    lb.XB(end + 1) = lxb;
end

% Format for spm_gmm
[bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', vol2vec(bffn));
mun                       = spm_gmm_lib('obs2cell', mun, code_image, false);
code_list                 = unique(code_image);
code_list                 = code_list(code_list ~= 0);

for it_appear=1:nit_appear

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GMM and get responsibilities (zn)
    %------------------------------------------------------------
    pr = datn.mog.pr;
    [zn,mog,~,lb,mg_w] = spm_gmm_lib('loop',bffn,scl_samp,{{m,b},{W,n}},{'LogProp', mun}, ...
                                     'GaussPrior',   {pr.m,pr.b, pr.W,pr.n}, ...
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

    if any(do_bf == true)

        % Update bias field parameters for each channel separately
        for c=1:C % Loop over channels

            if isempty(T{c}), continue; end

            % Compute gradient and Hessian (in image space)
            gr_im = zeros(d(1:3),'single');
            H_im  = zeros(d(1:3),'single');
            for l=1:size(msk_chn,1) % loop over combinations of missing voxels

                % Get mask of missing modalities (with this particular code)
                ixo = msk_chn(l,:);         % Observed channels
                ixm = ~ixo;                 % Missing channels
                nm  = sum(ixm);             % Number of missing channels
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
            d3 = numel(T{c}); % Number of DCT parameters
            H  = zeros(d3,d3);
            gr = zeros(d3,1);
            for z=1:d(3)
                b3 = double(chan(c).B3(z,:)');
                gr = gr + kron(b3,spm_krutil(double(gr_im(:,:,z)),double(chan(c).B1),double(chan(c).B2),0));
                H  = H  + kron(b3*b3',spm_krutil(double(H_im(:,:,z)),double(chan(c).B1),double(chan(c).B2),1));
            end
            clear b3

            % Gauss-Newton update of bias field parameters
            Update = reshape((H + chan(c).ICO)\(gr + chan(c).ICO*T{c}(:)),size(T{c}));
            clear H gr

            % Line-search
            armijo = scal_bf;
            oT     = T{c};
            ollpbf = llpbf;
            olx    = lb.X(end);
            olxb   = lb.XB(end);

            for ls=1:nit_lsbf

                % Update bias-field parameters
                T{c} = T{c} - armijo*Update;

                % Compute new bias-field (only for channel c)
                [bf(:,:,:,c),llpbf(c)] = BiasField(T(c),chan(c));
                bffn       = spm_gmm_lib('obs2cell', vol2vec(bf.*fn), code_image, true);

                % Recompute responsibilities (with updated bias field)
                zn = Responsibility(m,b,W,n,bffn,ReWeightMu(mun,log(mg_w)),msk_chn);

                % Compute new lower bound
                lx  = LowerBound('ln(P(X|Z))',bffn,zn,msk_chn,{m,b},{W,n},scl_samp);
                lxb = scl_samp*LowerBound('ln(|bf|)',bf,msk_vx) + sum(llpbf);

                % Check new lower bound
                if  ((lx + lxb) - (olx + olxb))/abs(lx + lxb) > -eps('single')*10
                    % Converged
                    %fprintf('it2=%i\tc=%i\tls=%i\tarmijo=%0.7f\tnl=%0.7f\tgain=%0.7f :o)\n',it_bf,c,ls,armijo,lx + lxb,lx + lxb - olx + olxb);
                    lb.X(end  + 1) = lx;
                    lb.XB(end + 1) = lxb; % bias field part of lower bound
                    break;
                else
                    armijo = armijo*0.5;
                    T{c}   = oT;
                    if ls == nit_lsbf
                        % Did not converge -> reset
                        %fprintf('it2=%i\tc=%i\tls=%i :o(\n',it_bf,c,ls);
                        llpbf       = ollpbf;
                        bf(:,:,:,c) = BiasField(T(c),chan(c));
                        bffn        = spm_gmm_lib('obs2cell', vol2vec(bf.*fn), code_image, true);
                        if nit_appear == 1 || nit_bf > 1
                            % Recompute responsibilities
                            zn = Responsibility(m,b,W,n,bffn,ReWeightMu(mun,log(mg_w)),msk_chn);
                        end
                    end
                end
            end
            clear oT Update
        end

        % Update datn
        datn.T = T;
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
        fn = fn0;

        % Bias field
        if any(do_bf == true)
            % Get full-sized bias field
            chan = BiasBasis(T,df,Mn,reg);
            bf   = BiasField(T,chan);
            bffn = bf.*fn;            
        else
            bffn = fn;
        end
                
        % For visualising results (spm_mb_show)
        spm_mb_io('Write2Visualise',datn,bffn,'bff',sett);
        if any(do_bf == true)
            spm_mb_io('Write2Visualise',datn,fn,'f',sett);
            spm_mb_io('Write2Visualise',datn,bf,'bf',sett);
            clear bf
        end
        clear fn 
    
        [bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', vol2vec(bffn));

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
    
    % For visualising results (spm_mb_show)    
    spm_mb_io('Write2Visualise',datn,zn,'z',sett);
end
end
%==========================================================================

%==========================================================================
function X2d = vol2vec(X4d)
d   = [size(X4d) 1 1];
X2d = reshape(X4d,[prod(d(1:3)) d(4)]);
end
%==========================================================================

%==========================================================================
% UpdatePrior()
function dat = UpdatePrior(dat, settn)
if ~isfield(dat,'mog'), return; end

% Parse function settings
mg_ix = sett.model.mg_ix;

% Get population indices
p_ix = GetPopulationIdx(dat);

for p=1:numel(p_ix) % Loop over populations

    N = numel(p_ix{p}); % Number of subjects in population

    % Get old prior
    pr = dat(p_ix{p}(1)).mog.pr;
    pr = {pr.m,pr.b,pr.W,pr.n};

    % Get all posteriors
    K1    = numel(mg_ix);
    po    = cell(1,N);
    for n=1:N
        n1          = p_ix{p}(n);
        po{n}{1}{1} = dat(n1).mog.po.m;
        po{n}{1}{2} = dat(n1).mog.po.b;
        po{n}{2}{1} = dat(n1).mog.po.W;
        po{n}{2}{2} = dat(n1).mog.po.n;
    end

    % Update prior
    pr = spm_gmm_lib('updatehyperpars',po,pr);

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
% IndexSlice2Vol() % unused
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
    bf       = varargin{1};
    msk_vx   = varargin{2};
    lb       = sum(log(bf(msk_vx)),'double');
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
df   = size(fn);
df   = [df 1];
df   = df(1:3);

ind  = SampleInd(df,Mn,samp);
d    = cellfun(@length,ind);  % New dimensions
ofn  = fn(ind{:},:,:);

% For weighting data parts of lowerbound with factor based on amount of
% downsampling
w    = prod(df)/prod(d);
end
%==========================================================================

%==========================================================================
function ind = SampleInd(df,Mn,samp)
vx   = sqrt(sum(Mn(1:3,1:3).^2));
df   = [df 1];
df   = df(1:3);
sk   = max([1 1 1],round(samp*[1 1 1]./vx));
ind  = {1:sk(1):df(1), 1:sk(2):df(2), 1:sk(3):df(3)};
end
%==========================================================================

%==========================================================================
% SumLowerBound()
function lb = SumLowerBound(lb)

if ~isfinite(lb.sum(1)), return; end

fields = fieldnames(lb);
lb.sum(end+1) = 0;
for i=1:numel(fields)
    field = fields{i};
    if ~any(strcmpi(field, {'sum' 'last'})) && ~isempty(lb.(field)) && ~isnan(lb.(field)(end))
        lb.sum(end) = lb.sum(end) + sum(lb.(field)(:,end));
    end
end
end
%==========================================================================
