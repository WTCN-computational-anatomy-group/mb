function varargout = spm_mb_appearance(varargin)
%__________________________________________________________________________
%
% Functions for appearance model related.
%
% FORMAT [bfn,lln] = spm_mb_appearance('BiasField',chan,d,varargin)
% FORMAT chan      = spm_mb_appearance('BiasFieldStruct',datn,C,df,reg,fwhm,scl,T,samp)
% FORMAT labels    = spm_mb_appearance('GetLabels',datn,sett,do_samp)
% FORMAT p_ix      = spm_mb_appearance('GetPopulationIdx',dat)
% FORMAT dat       = spm_mb_appearance('Init',dat,model,K,sett)
% FORMAT fn        = spm_mb_appearance('Mask',fn,is_ct)
% FORMAT zn        = spm_mb_appearance('Responsibility',m,b,W,n,fn,mu,msk_chn)
% FORMAT [zn,datn] = spm_mb_appearance('Update',datn,mun0,sett)
% FORMAT dat       = spm_mb_appearance('UpdatePrior',dat,mu,sett)
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
function chan = BiasFieldStruct(datn,C,df,reg,fwhm,scl,T,samp)
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

    if ~isempty(scl) && do_bf(c)
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
    labels = zeros(1,K1);
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
function dat = Init(dat,model,K,sett)
p_ix = spm_mb_appearance('GetPopulationIdx',dat);
Np   = numel(p_ix);

appear_given   = isfield(model,'appear');
template_given = (isfield(model,'shape') && isfield(model.shape,'template'));

% Get template
if template_given, mu = spm_mb_io('GetData',model.shape.template);
else,              mu = [];
end

for p=1:Np
    datn = dat(p_ix{p});
    if appear_given, pr = model.appear(num2str(datn(1).ix_pop));
    else,            pr = [];
    end
    dat(p_ix{p}) = InitPopulation(datn,mu,pr,K,sett);
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
function [zn,datn] = Update(datn,mun0,sett)
% Update appearance model for a single subject (GMM & bias field)
%
% FORMAT [zn,datn] = Update(datn,mun0,sett)
% datn - Structure holding data for a single subject
% mun0 - Log template
% sett - Structure of settings

% Parse function settings
do_updt_bf   = sett.do.updt_bf;
fwhm         = sett.bf.fwhm;
nit_bf       = sett.nit.bf;
nit_gmm      = sett.nit.gmm;
nit_gmm_miss = sett.nit.gmm_miss;
nit_appear   = sett.nit.appear;
nit_lsbf     = sett.optim.nls_bf;
reg          = sett.bf.reg;
samp         = sett.gen.samp;
tol_gmm      = sett.appear.tol_gmm;

% Parameters
[df,C]     = spm_mb_io('GetSize',datn.f);
K          = size(mun0,4);
K1         = K + 1;
Mn         = datn.Mat;
scl_samp   = 1; % sampling scaling (defined as W = prod(d0(1:3))/prod(d(1:3)), when samp > 1)
do_bf      = datn.do_bf;
is_ct      = datn.is_ct;

% Get image data
fn = spm_mb_io('GetData',datn.f);

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
    nvx_full = GetNumVoxObserved(fn,is_ct);  % get number of obseved voxels in input image(s)
    [fn,d]   = SubSample(fn,Mn,samp);        
    nvx_samp = GetNumVoxObserved(fn,is_ct);  % get number of obseved voxels in subsampled image(s)  
    scl_samp = nvx_full/nvx_samp;            % get voxel ratio between original and subsamped image(s)
    % Template
    mun = SubSample(mun0,Mn,samp);      
    mun = reshape(mun,[prod(d(1:3)) K]);        
else
    d    = df;
    mun  = reshape(mun0,[prod(d(1:3)) K]);
    mun0 = [];
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
labels = [];

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

ol = lb.sum(end);
for it_appear=1:nit_appear

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GMM and get responsibilities (zn)
    %------------------------------------------------------------

    [zn,mog,~,lb] = spm_gmm_lib('loop',bffn,scl_samp,{{m,b},{W,n}},{'LogProp', mun}, ...
                                 'GaussPrior',   {m0,b0,W0,n0}, ...
                                 'Missing',      msk_chn, ...
                                 'LowerBound',   lb, ...
                                 'IterMax',      nit_gmm, ...
                                 'Tolerance',    tol_gmm, ...
                                 'SubIterMax',   nit_gmm_miss, ...
                                 'SubTolerance', tol_gmm, ...
                                 'Verbose',      [0 0]);
    m = mog.MU;
    b = mog.b;
    W = mog.V;
    n = mog.n;           

    nl = lb.sum(end);        
%     fprintf('it1=%i\tnl=%0.7f\tgain=%0.7f\n',it_likel,nl,nl - ol);
    if it_appear > 1 && ((nl - ol)/abs(nl) > -eps('single')*10 || it_appear == nit_appear)
        % Finished
        break
    end
    ol = nl;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update bias field parameters
    % This computes the derivatives of the negative logarithm of the
    % joint probability distribution
    %------------------------------------------------------------

    if do_updt_bf && any(do_bf == true)        

        % Make sure to use the latest responsibilties
        zn = Responsibility(m,b,W,n,bffn,mun,msk_chn);

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
                    for k=1:K1
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
                    obffn = [];
                    
                    % Add terms related to the normalisation (log(b))
                    go = go - 1;
                    Ho = Ho + 1; % Comes from subs(H,g,0)

                    % Accumulate across missing codes
                    ixvx        = (code_image == code_list(l));
                    gr_im(ixvx) = gr_im(ixvx) + go;
                    H_im(ixvx)  = H_im(ixvx)  + Ho;
                    ixvx        = [];
                end
                zn = [];

                 % Compute gradient and Hessian (transform from image space to parameter space)
                d3 = numel(chan(c).T); % Number of DCT parameters
                H  = zeros(d3,d3);     
                gr = zeros(d3,1);      
                for z=1:d(3)
                    b3 = double(chan(c).B3(z,:)');
                    gr = gr + kron(b3,spm_krutil(double(gr_im(:,:,z)),double(chan(c).B1),double(chan(c).B2),0));
                    H  = H  + kron(b3*b3',spm_krutil(double(H_im(:,:,z)),double(chan(c).B1),double(chan(c).B2),1));
                end
                b3 = [];                    

                % Gauss-Newton update of bias field parameters
                Update = reshape((H + chan(c).ICO)\(gr + chan(c).ICO*chan(c).T(:)),size(chan(c).T));
                H      = []; gr = [];

                % Line-search
                armijo = 1;        
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
                    zn = Responsibility(m,b,W,n,bffn,mun,msk_chn);

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
                            zn    = Responsibility(m,b,W,n,bffn,mun,msk_chn);
                        end
                    end
                end
                oT = []; Update = [];
            end
        end   

        % Update datn     
        datn.bf.T = {chan(:).T};
        datn.E(3) = -sum(pr_bf); % global objective function is negative log-likelihood..
    end
end
fn = []; bf = []; mun = [];

if samp > 1
    % Compute responsibilities on original data
    fn = spm_mb_io('GetData',datn.f);
    fn = reshape(fn,[prod(df(1:3)) C]);
    fn = Mask(fn,is_ct);
    
    mun0 = reshape(mun0,[prod(df(1:3)) K]);
    if any(do_bf == true)
        % Get full-sized bias field
        chan = BiasFieldStruct(datn,C,df,reg,fwhm,[],datn.bf.T);
        bf   = BiasField(chan,df);
        bffn = bf.*fn;
        bf   = [];
    else
        bffn = fn;
    end
    fn = [];
       
    [bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', bffn);    

    mun0 = spm_mb_shape('TemplateK1',mun0,2);
    
    labels = GetLabels(datn,sett);    
    mun0   = mun0 + labels;
    labels = [];
    
    mun0 = spm_gmm_lib('obs2cell', mun0, code_image, false);            
    zn   = Responsibility(m,b,W,n,bffn,mun0,msk_chn);
end       
zn = spm_gmm_lib('cell2obs', zn, code_image, msk_chn);

% Get 4D versions of K1 - 1 classes
zn = reshape(zn(:,1:K),[df(1:3) K]);

% Update datn     
datn.E(1)     = -lb.sum(end);
datn.mog.po.m = m;
datn.mog.po.b = b;
datn.mog.po.W = W;
datn.mog.po.n = n;          
datn.mog.lb   = lb;  
end
%==========================================================================

%==========================================================================
% UpdatePrior()
function dat = UpdatePrior(dat, mu, sett)

if ~sett.do.updt_int,      return; end
if ~isfield(dat(1),'mog'), return; end

% Get population indices
p_ix = GetPopulationIdx(dat);

for p=1:numel(p_ix) % Loop over populations

    N = numel(p_ix{p}); % Number of subjects in population

    % Get old prior
    pr = dat(p_ix{p}(1)).mog.pr;
    pr = {pr.m,pr.b,pr.W,pr.n};
    C  = size(pr{1},1);

    % Get all posteriors
    K     = size(mu,4);
    K1    = K + 1;
    po    = cell(1,N);
    for n=1:N
        n1          = p_ix{p}(n);
        po{n}{1}{1} = dat(n1).mog.po.m;
        po{n}{1}{2} = dat(n1).mog.po.b;
        po{n}{2}{1} = dat(n1).mog.po.W;
        po{n}{2}{2} = dat(n1).mog.po.n;
    end

if false
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

    % Add one artificial observation
    po1{1}{1} = repmat(avgmn,[1 K1]);     % m
    po1{1}{2} = zeros(1,K1) + 0.01;       % b
    po1{2}{1} = repmat(avgpr/C,[1 1 K1]); % W
    po1{2}{2} = C*ones(1,K1);             % n
%   po{N+1}   = po1; %%% DISABLED
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
if is_ct, f(~isfinite(f) | f == 0 | f <= - 1020) = NaN;
else,     f(~isfinite(f) | f == 0)               = NaN;
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
% GetMask()
function nvx = GetNumVoxObserved(f,is_ct)
f   = ApplyMask(f,is_ct);
msk = ~isnan(f);
msk = sum(msk,4);
nvx = sum(msk(:) > 0);
end
%==========================================================================

%==========================================================================
% InitPopulation()
function dat = InitPopulation(dat,mu,pr,K,sett)

% Parse function settings
do_gmm = sett.do.gmm;
fwhm   = sett.bf.fwhm;
reg    = sett.bf.reg;

if ~do_gmm, return; end

N  = numel(dat);
K1 = K + 1;
lb = struct('sum', NaN, 'X', [], 'XB', [], ...
            'Z', [], 'P', [], 'MU', [], 'A', []);

[~,C] = spm_mb_io('GetSize',dat(1).f);
mx    = zeros(C,numel(dat));
mn    = zeros(C,numel(dat));
avg   = zeros(C,numel(dat));
vr    = zeros(C,numel(dat));

for n=1:N
    [df,C] = spm_mb_io('GetSize',dat(n).f);
    fn     = spm_mb_io('GetData',dat(n).f);
    fn     = reshape(fn,[prod(df(1:3)) C]);
    fn     = spm_mb_appearance('Mask',fn,dat(n).is_ct);
    if any(dat(n).do_bf == true)
        val = 1e3;
        scl = ones(1,C);
        for c=1:C
            msk    = isfinite(fn(:,c));
            scl(c) = val./mean(fn(msk,c));
        end
    else
        scl = ones(1,C);
    end
 
    if any(dat(n).do_bf == true)
        % Get bias field parameterisation struct
        chan        = spm_mb_appearance('BiasFieldStruct',dat(n),C,df,reg,fwhm,scl);
        dat(n).bf.T = {chan(:).T};

        % Get bias field
        bf = spm_mb_appearance('BiasField',chan,df);
    else
        bf = ones([1 C],'single');
    end
 
    % Modulate with bias field
    fn = bf.*fn;
 
    % Init GMM
    [po,mx(:,n),mn(:,n),avg(:,n),vr(:,n)] = InitPosteriorGMM(dat(n),fn,mu,pr,K1,sett);
    mog.po     = po;
    mog.lb     = lb;
    dat(n).mog = mog;
end

if isempty(pr)
    % Init GMM empirical prior
    pr = InitPriorGMM(mx,mn,avg,vr,mu,K1);
end
for n=1:numel(dat)                
    dat(n).mog.pr = pr; 
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
% InitPriorGMM()
function pr = InitPriorGMM(mx,mn,avg,vr,mu0,K)
% Initialise the prior parameters of a Gaussian mixture model. 
%
%   This function is only used to initialise the GMM parameters at the
%   beginning of the algorithm.
% 
% FORMAT pr = InitPriorGMM(mx,mn,vr,mu0,K)
% mx   - Maximum observed value [per channel]
% mn   - Minimum observed value [per channel]
% avg  - Mean observed value [per channel]
% vr   - Variance of observed value [per channel]
% mu0  - Log template
% K    - Number of classes
% pr   - Structure holding prior GMM parameters (m, b, W, n)

mvr  = mean(vr,2); % mean variance across all subjects in population
mmx  = mean(mx,2); % mean mean across all subjects in population
mmn  = mean(mn,2); % mean mean across all subjects in population
mavg = mean(avg,2); % mean mean across all subjects in population

C   = size(mx,1);
m   = zeros(C,K);
ico = zeros(C,C,K);   
for c=1:C      
    rng    = linspace(mmn(c),mmx(c),K);
    rng    = -sum(rng<0):sum(rng>=0) - 1;
    m(c,:) = rng'*mmx(c)/(1.5*K);
    
    ico(c,c,:) = mmx(c)/(1.5*K);
    ico(c,c,:) = 1/ico(c,c,:); % precision
end

if ~isempty(mu0)
    % If template is given, make sure that pr.m is the same for all classes
    m = repmat(mean(m,2),[1 K]);
end

% Define prior
pr   = struct('m',[],'b',[],'n',[],'W',[]);
pr.m = m;
pr.b = zeros(1,K) + 0.01;
pr.n = C*ones(1,K);
pr.W = ico/C;

if 0
    fig_name = 'Prior';
    spm_gmm_lib('plot','gaussprior',{pr.m,pr.b,pr.W,pr.n},[],fig_name);
end

end
%==========================================================================   

%==========================================================================    
% InitPosteriorGMM()
function [po,mx,mn,avg,vr] = InitPosteriorGMM(datn,fn,mu,pr,K,sett)
% Initialise the posterior parameters of a Gaussian mixture model. 
%
%   This function is only used to initialise the GMM parameters at the
%   beginning of the algorithm.
% 
% FORMAT [po,mx,mn,vr] = InitPosteriorGMM(datn,fn,mu,pr,K,sett)
% datn - Structure holding data of a single subject
% fn   - Bias corrected observed image, in matrix form [nbvox nbchannel]
% mu   - Log template
% pr   - Structure holding prior GMM parameters (m, b, W, n)
% K    - Number of classes
% sett - Structure of settings
% po   - Structure of posterior parameters (m, b, W, n)
% mx   - Maximum observed value [per channel]
% mn   - Minimum observed value [per channel]
% avg  - Mean observed value [per channel]
% vr   - Variance of observed value [per channel]

% Parse function settings
B   = sett.registr.B;
Mmu = sett.var.Mmu;

C   = size(fn,2);
mx  = zeros(C,1);
mn  = zeros(C,1);
avg = zeros(C,1);
vr  = zeros(C,1);
m   = zeros(C,K);
ico = zeros(C,C,K);
for c=1:C
    msk    = isfinite(fn(:,c));
    mx(c)  = max(fn(msk,c));
    mn(c)  = min(fn(msk,c));
    avg(c) = mean(fn(msk,c));
    vr(c)  = var(fn(msk,c));
    
    rng    = linspace(mn(c),mx(c),K);
    rng    = -sum(rng<0):sum(rng>=0) - 1;
    m(c,:) = rng'*mx(c)/(1.5*K);
    
    ico(c,c,:) = mx(c)/(1.5*K);
    ico(c,c,:) = 1/ico(c,c,:); % precision
end

% Initialise posterior
if isempty(pr)
    % No prior given, initialise from image statistics (mean, min, max, var)
    po   = struct('m',[],'b',[],'n',[],'W',[]);
    po.m = m;
    po.b = zeros(1,K) + 0.01;
    po.n = C*ones(1,K);
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
    po.m = repmat(mean(po.m,2),[1 K]);
end
   
if ~isempty(mu)
    % Use template as resposibilities to compute values for GMM
    % posterior from suffstats
    df = spm_mb_io('GetSize',datn.f);          
    q  = double(datn.q);
    Mr = spm_dexpm(q,B);
    Mn = datn.Mat;    
    
    % Warp template
    psi1 = spm_mb_io('GetData',datn.psi);
    psi  = spm_mb_shape('Compose',psi1,spm_mb_shape('Affine',df,Mmu\Mr*Mn));
    mu   = spm_mb_shape('Pull1',mu,psi);            
    psi1 = [];

    % Add class, then softmax -> can now be used as resps
    mu = cat(4,mu,zeros(df(1:3),'single'));
    mu = spm_mb_shape('Softmax',mu,4);
    mu = reshape(mu,[prod(df(1:3)) K]);        
        
    % Compute posterior estimates
    [SS0,SS1,SS2] = spm_gmm_lib('SuffStat',fn,mu,1);
    [MU,~,b,W,n]  = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2, {po.m,po.b,po.W,po.n});
    
    po.m = MU;
    po.b = b;
    po.n = n;
    po.W = W;
end

if 0
    fig_name = 'Posterior';
    spm_gmm_lib('plot','gaussprior',{po.m,po.b,po.W,po.n},[],fig_name);
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
