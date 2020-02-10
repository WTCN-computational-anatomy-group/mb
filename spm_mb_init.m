function varargout = spm_mb_init(varargin)
%__________________________________________________________________________
%
% Functions for initialising appearance model.
%
% FORMAT datn       = spm_mb_init('InitBias',datn,fwhm,dc)
% FORMAT [dat,sett] = spm_mb_init('Init',dat,model,K,sett)
% FORMAT [dat,mu]   = PropagateTemplate(dat,mu,sz,sett)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_appearance
    error('Not enough argument. Type ''help spm_mb_appearance'' for help.');
end
id       = varargin{1};
varargin = varargin(2:end);
switch id
    case 'InitBias'
        [varargout{1:nargout}] = InitBias(varargin{:});
    case 'Init'
        [varargout{1:nargout}] = Init(varargin{:});
    case 'PropagateTemplate'
        [varargout{1:nargout}] = PropagateTemplate(varargin{:});
    otherwise
        help spm_mb_appearance
        error('Unknown function %s. Type ''help spm_mb_appearance'' for help.', id)
end
end
%==========================================================================

%==========================================================================
function datn = InitBias(datn,fwhm,dc)
[df,C] = spm_mb_io('GetSize',datn.f);
T      = cell(1,C);
vx     = sqrt(sum(datn.Mat(1:3,1:3).^2));
for c=1:C
    if datn.do_bf(c)
        d3 = ceil(2*vx.*df/fwhm);
    else
        d3 = [1 1 1];
    end
    T{c} = zeros(d3);

    if ~isempty(dc)
        % Change DC component of bias field to make intensities more
        % simillar between MR images.
        bbb         = spm_dctmtx(df(1),1,1)*spm_dctmtx(df(2),1,1)*spm_dctmtx(df(3),1,1);
        T{c}(1,1,1) = dc(c)/bbb;
    end
end
datn.T = T;
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

if isscalar(mg_ix)
    mg_ix            = repelem(1:K + 1,mg_ix);
    sett.model.mg_ix = mg_ix;
end

% Get parameters
N   = numel(dat);
K1  = K + 1;

% What model parameters were given?
[template_given,appear_given] = spm_mb_param('SetFit',model,sett);
if appear_given
    % Set mg_ix to mg_ix that was used when learning intensity model
    if isfield(model.appear,'mg_ix'), mg_ix = model.appear.mg_ix;
    else,                             mg_ix = 1:K + 1;
    end
    sett.model.mg_ix       = mg_ix;
    sett.model.mg_ix_intro = mg_ix;

    Kmg  = numel(mg_ix);
    mg_w = ones(1,Kmg)./arrayfun(@(x) sum(x == mg_ix), mg_ix);
else
    mg_w = ones([1 K1]);
end

% Lower bound struct
lb  = struct('sum', NaN, 'X', [], 'XB', [], 'Z', [], 'P', [], 'MU', [], 'A', []);

mu   = [];
w_mu = 1;
if template_given
    % Get template
    mu = spm_mb_io('GetData',model.shape.template);

    % Build a weigthed mean (w_mu) that can be used to initi bias field DC scaling
    w_mu = spm_mb_shape('TemplateK1',mu,4);
    w_mu = exp(w_mu);
    w_mu = sum(sum(sum(w_mu,1),2),3);
    w_mu = w_mu./sum(w_mu);
    w_mu = reshape(w_mu,[1 numel(w_mu)]);
    w_mu = w_mu(mg_ix)./arrayfun(@(x) sum(x == mg_ix), mg_ix); % incorporate multiple Gaussians per tissue
end

if appear_given && template_given
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO
    %------------------------------------------------------------

    % TODO
    pr  = model.appear.pr;
    dc0 = w_mu.*pr.m;
    dc0 = sum(dc0,2);

    for n=1:N
        [df,C] = spm_mb_io('GetSize',dat(n).f);
        do_dc  = dat(n).do_dc;

        % Load image data
        fn = spm_mb_io('GetImage',dat(n));        
        fn = reshape(fn,[prod(df(1:3)) C]);
        
        dc = zeros([1 C]);
        if do_dc
            for c=1:C
                msk   = isfinite(fn(:,c));
                dc(c) = dc0(c)./mean(mean(mean(fn(msk,c))));
            end
            dc = log(dc);
        end

        dat(n) = InitBias(dat(n),fwhm,dc);
        if any(dat(n).do_bf == true)
            % Modulate with bias field
            chan = spm_mb_appearance('BiasBasis',dat(n).T,df,dat(n).Mat,reg);
            bf   = spm_mb_appearance('BiasField',dat(n).T,chan);
            bf   = reshape(bf,[prod(df(1:3)) C]);
            fn   = bf.*fn;
        end

        % TODO        
        po = InitSimplePosteriorGMM(dat(n),fn,mu,pr,K1,mg_ix,sett);

        mog.po   = po;
        mog.pr   = pr;
        mog.lb   = lb;
        mog.mg_w = mg_w;

        dat(n).mog = mog;
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO
    %------------------------------------------------------------

    % TODO
    [ix_ct,ix_mri1,ix_mri2] = spm_mb_io('GetCTandMRI',dat,sett);

    b  = zeros(1,K1) + 0.01;
    nu = ones(1,K1);

    if ~isempty(ix_ct) && (isempty(ix_mri1) && isempty(ix_mri2))
        % TODO

        dat(ix_ct) = InitGMM(dat(ix_ct),sett);
    else
        % TODO

        % Run on 'ix_mri1' subjects
        dat(ix_mri1) = InitGMM(dat(ix_mri1),sett);        

        if ~isempty(ix_mri2)
            % TODO

            Nmri2 = numel(ix_mri2);
            for n=1:Nmri2

                n1     = ix_mri2(n);
                [df,C] = spm_mb_io('GetSize',dat(n1).f);
                do_dc  = dat(n1).do_dc;

                po       = struct('m',ones(C,K1),'b',b,'W',repmat(eye(C),[1 1 K1])/C,'n',C*nu);
                mog.po   = po;
                mog.pr   = po;
                mog.lb   = lb;
                mog.mg_w = mg_w;

                dat(n1).mog = mog;
                
                fn = spm_mb_io('GetImage',dat(n1));
                fn = reshape(fn,[prod(df(1:3)) C]);
                
                % Set bias field dc scaling
                if do_dc
                    dc = ones(1,C);
                    for c=1:C
                        msk   = isfinite(fn(:,c));
                        dc(c) = dc(c)./mean(fn(msk,c));
                    end
                    dc = -log(dc);
                else
                    dc = zeros(1,C);
                end
                dat(n1) = InitBias(dat(n1),fwhm,dc);
            end
        end

        if ~isempty(ix_ct)
            % Init CT subjects

            Nct = numel(ix_ct);
            mn  = -1020;
            mx  = 3000;
            mu  = 50*ones([1 K1]);
            Sig = diag(((mx + mn)/K1).^2).*ones([1,1,K1]);

            % Compute precision
            W                    = zeros(size(Sig));
            for k=1:K1, W(:,:,k) = inv(Sig(:,:,k)); end

            po       = struct('m',mu,'b',b,'W',W,'n',nu);
            mog.po   = po;
            mog.pr   = po;
            mog.lb   = lb;
            mog.mg_w = mg_w;

            for n=1:Nct
                n1          = ix_ct(n);
                dat(n1).mog = mog;
                if any(dat(n1).do_bf == true)
                    dat(n1) = InitBias(dat(n1),fwhm,dc);
                end
            end
        end
    end
end
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
        fn = spm_mb_io('GetImage',dat(n));

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

    dat(n) = InitBias(dat(n),fwhm,dc);
    if any(dat(n).do_bf == true)
        % Modulate with bias field
        chan = spm_mb_appearance('BiasBasis',dat(n).T,df,dat(n).Mat,reg,samp);
        bf   = spm_mb_appearance('BiasField',dat(n).T,chan);
        fn   = bf.*fn;
    end

    if ~use_initgmm
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

% Number of voxels to sample, and what indices (defined in template space)
crp_prct = 0.3; % to skip sampling to many air voxels, which are usually in the outer parts of the images
dmu_crp  = ceil((1 - 2*crp_prct)*dmu);
diff_mu  = ceil((dmu - dmu_crp)/2);
Nvx      = min(round(prod(dmu(1:3))/N), round(prod(dmu(1:3)./sk)));
r        = randperm(prod(dmu_crp(1:3)),Nvx);

% Voxel indicies in template space
[x1,x2,x3] = ind2sub(dmu(1:3),r(:)); % use full template dimensions
x1         = x1 + diff_mu(1);        % then 'crop' according to 'crp_prct'
x2         = x2 + diff_mu(2);
x3         = x3 + diff_mu(3);
ymu        = cat(2,single(x1),single(x2),single(x3),ones([Nvx 1],'single'));
clear r x1 x2 x3

% Get sample of image data (and labels, if provided)
F = cell(1,N); % holds imaging data
L = cell(1,N); % holds possible label data
for n=1:N, L{n} = zeros([K1 1],'single'); end

labels_given = false;

for n=1:N % Loop over subjects

    fn    = NaN([C Nvx],'single');
    cnt_c = 1;         % channel count
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
        f1 = spm_mb_io('GetImage',dat(n1));

        % Move template space sample points to subject space
        Mn = dat(n1).Mat;
        M  = Mmu\Mn;
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
        mu(c,:) = rng'*mx(c)/K1;
    end
end

Sig = diag(((mx + mn)/K1).^2).*ones([1,1,K1]);

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

% Set NaNs to zero
dc(~isfinite(dc)) = 0;

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
            x   = min(xn):max(xn);
            h   = hist(xn,x);
            h   = h/sum(h);
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
lb  = struct('sum', NaN, 'X', [], 'XB', [], 'Z', [], 'P', [], 'MU', [], 'A', []);

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

        dat(n1) = InitBias(dat(n1),fwhm,dc);
    end
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
    m(c,:) = rng'*mx(c)/K1;

    ico(c,c,:) = (mx(c)/K1)^2;
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
% PropagateTemplate()
function [dat,mu] = PropagateTemplate(dat,mu,sz,sett)

% Parse function settings
nit_init_mu = sett.nit.init_mu;

% Partion CT and MR images
[ix_ct,ix_mri1,ix_mri2] = spm_mb_io('GetCTandMRI',dat,sett);

if ~isempty(ix_mri1) && (~isempty(ix_ct) || ~isempty(ix_mri2))

    sett.gen.samp = numel(sz); % coarse-to-fine sampling of observed data

    for it=1:nit_init_mu
        [mu,dat(ix_mri1)] = spm_mb_shape('UpdateMean',dat(ix_mri1), mu, sett);
        dat(ix_mri1)      = spm_mb_appearance('UpdatePrior',dat(ix_mri1), sett);
    end

    if ~isempty(ix_mri2)
        for it=1:nit_init_mu
            [mu,dat([ix_mri1 ix_mri2])] = spm_mb_shape('UpdateMean',dat([ix_mri1 ix_mri2]), mu, sett);
            dat([ix_mri1 ix_mri2])      = spm_mb_appearance('UpdatePrior',dat([ix_mri1 ix_mri2]), sett);
        end
    end
end
end
%==========================================================================

