function varargout = spm_mb_appearance(varargin)
%__________________________________________________________________________
%
% Functions for appearance model related.
%
% FORMAT [bfn,lln] = spm_mb_appearance('BiasField',chan,d,varargin)
% FORMAT chan      = spm_mb_appearance('BiasFieldStruct',datn,C,d,reg,fwhm,scl,T,samp)
% FORMAT p_ix      = spm_mb_appearance('GetPopulationIdx',dat)
% FORMAT dat       = spm_mb_appearance('Init',dat,K,sett)
% FORMAT fn        = spm_mb_appearance('Mask',fn)
% FORMAT zn        = spm_mb_appearance('Responsibility',m,b,V,n,fn,mu,L,code)
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
            bf_c      = TransformBF(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);
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
        bf_c      = TransformBF(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);
        bf_c      = bf_c(:);      
        bfn(ix,c) = single(exp(bf_c));        
    end
end
end
%==========================================================================

%==========================================================================
% BiasFieldStruct()
function chan = BiasFieldStruct(datn,C,d,reg,fwhm,scl,T,samp)
if nargin < 7, T    = {}; end
if nargin < 8, samp = 1; end

cl   = cell(C,1);    
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
chan = struct(args{:});

do_bf = datn.do_bf;
Mn    = datn.Mat;

vx = sqrt(sum(Mn(1:3,1:3).^2));                     
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

    if isempty(T) || ~do_bf(c)
        % Initial parameterisation of bias field
        chan(c).T = zeros(d3,'single');
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
% GetPopulationIdx()
function p_ix = GetPopulationIdx(dat)
int_pr_ix = [dat.int_pr_ix];
un        = unique(int_pr_ix);
p_ix      = cell(1,numel(un));
for i=un
    p_ix{i} = find(int_pr_ix == i);
end
end
%==========================================================================

%==========================================================================
% Init()
function dat = Init(dat,K,sett)

% Parse function settings
do_bf_norm = sett.do.bf_norm;
do_gmm     = sett.do.gmm;
fwhm       = sett.bf.fwhm;
reg        = sett.bf.reg;
do_updt_bf = sett.do.updt_bf;

if (~do_gmm && ~do_updt_bf), return; end

N  = numel(dat);
K1 = K + 1;
lb = struct('sum', NaN, 'X', [], 'XB', [], ...
            'Z', [], 'P', [], 'MU', [], 'A', []);

[~,C] = spm_mb_io('GetSize',dat(1).f);
mx    = zeros(C,numel(dat));
mn    = zeros(C,numel(dat));
vr    = zeros(C,numel(dat));

for n=1:N
    [df,C] = spm_mb_io('GetSize',dat(n).f);
    fn     = spm_mb_io('GetData',dat(n).f);                    
    fn     = reshape(fn,[prod(df(1:3)) C]);                      
    fn     = spm_mb_appearance('Mask',fn);
    if do_updt_bf
        val = 1e3;                    
        scl = double(val./nanmean(fn,1)); 
    else        
        scl = ones(1,C);
    end
    
    if do_updt_bf && any(dat(n).do_bf == true)
        % Get bias field parameterisation struct
        chan        = spm_mb_appearance('BiasFieldStruct',dat(n),C,df,reg,fwhm,scl);
        dat(n).bf.T = {chan(:).T};

        % struct used for rescaling images using DC component of bias fields
        dc            = struct;
        dc.int        = zeros(1,C);
        dc.ln         = zeros(1,C);
        dat(n).bf.dc  = dc;    
        
        % Get bias field
        bf = spm_mb_appearance('BiasField',chan,df);     
    else
        bf = ones([1 C],'single');
    end
    
    % Modulate with bias field    
    fn = bf.*fn;
    
    % Init GMM
    [po,mx(:,n),mn(:,n),vr(:,n)] = PosteriorGMM(fn,K1,sett);            
    mog.po     = po;    
    mog.lb     = lb;    
    dat(n).mog = mog;
end

% Init GMM empirical prior
pr = PriorGMM(mx,mn,vr,K1,sett);
for n=1:numel(dat)                
    dat(n).mog.pr = pr;
end
end
%==========================================================================

%==========================================================================
% Mask()
function fn = Mask(fn)
C = size(fn,2);
for c=1:C
    fn(:,c) = ApplyMask(fn(:,c));
end
end
%==========================================================================

%==========================================================================
% Responsibility()
function zn = Responsibility(m,b,V,n,fn,mu,L,code)

% Is there missing data?
do_miss = numel(L) > 1;

if do_miss, const = spm_gmm_lib('Const', {m,b}, {V,n}, L);
else,       const = spm_gmm_lib('Const', {m,b}, {V,n});
end

fn = spm_gmm_lib('Marginal', fn, {m,V,n}, const, {code,L});
zn = spm_gmm_lib('Responsibility', fn, mu);
end
%==========================================================================

%==========================================================================
% Update()
function [zn,datn] = Update(datn,mun0,sett)

% Parse function settings
do_bf_norm   = sett.do.bf_norm;
do_updt_bf   = sett.do.updt_bf;
fwhm         = sett.bf.fwhm;
nit_bf       = sett.nit.bf;
nit_gmm      = sett.nit.gmm;
nit_gmm_miss = sett.nit.gmm_miss;
nit_appear   = sett.nit.appear;
nit_lsbf     = sett.optim.nls_bf;
reg          = sett.bf.reg;
samp         = sett.gen.samp;
tol_bf       = sett.bf.tol;
tol_appear   = sett.appear.tol;

% Parameters
[df,C]     = spm_mb_io('GetSize',datn.f);
K          = size(mun0,4);
K1         = K + 1;
Mn         = datn.Mat;
W          = 1;

% Get image data
fn = spm_mb_io('GetData',datn.f);

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

if nargout > 1
    % Lower bound
    lb = datn.mog.lb;

    if samp > 1
        % Subsample (runs faster, lower bound is corrected by scalar W)              
        [fn,mun,W,d] = SubSample(samp,Mn,fn,mun0);
        mun          = reshape(mun,[prod(d(1:3)) K]);
    else
        d    = df;    
        mun  = reshape(mun0,[prod(d(1:3)) K]);
        mun0 = [];
    end        
    nm = prod(d(1:3));
    fn = reshape(fn,[prod(d(1:3)) C]);

    % Missing data stuff
    fn      = Mask(fn);
    code    = spm_gmm_lib('obs2code', fn);
    L       = unique(code);
    nL      = numel(L);
    do_miss = numel(L) > 1;

    % Make K + 1 template
    mun = cat(2,mun,zeros([prod(d(1:3)) 1],'single'));

    % Bias field related
    if do_updt_bf    
        chan       = BiasFieldStruct(datn,C,df,reg,fwhm,[],datn.bf.T,samp);
        [bf,pr_bf] = BiasField(chan,d);
        bffn       = bf.*fn;
    else
        bffn       = fn;
    end

    ol = lb.sum(end);
    for it_appear=1:nit_appear
        % Update GMM and get responsibilities (zn)
        [zn,mog,~,lb] = spm_gmm_loop({bffn,W},{{m,b},{V,n}},{'LogProp', mun}, ...
                                     'GaussPrior',   {m0,b0,V0,n0}, ...
                                     'Missing',      do_miss, ...
                                     'LowerBound',   lb, ...
                                     'MissingCode',  {code,L}, ...
                                     'IterMax',      nit_gmm, ...
                                     'Tolerance',    1e-4, ...
                                     'SubIterMax',   nit_gmm_miss, ...
                                     'SubTolerance', 1e-4, ...
                                     'Verbose',      [0 0]);
        m = mog.MU;
        b = mog.b;
        V = mog.V;
        n = mog.n;    

        nl = lb.sum(end);        
    %     fprintf('it1=%i\tnl=%0.7f\tgain=%0.7f\n',it_likel,nl,nl - ol);
        if it_appear > 1 && ((nl - ol) < 2*nm*tol_appear || it_appear == nit_appear)
            % Finished
            break
        end
        ol = nl;

        if do_updt_bf
            % Update bias field parameters                    
            lx   = W*LowerBound('X',bffn,zn,code,{m,b},{V,n});
            lxb  = W*LowerBound('XB',bf);       
            done = false(1,C);
            for it_bf=1:nit_bf

                % Update bias field parameters for each channel separately
                for c=1:C % Loop over channels

                    if done(c) || ~datn.do_bf(c)
                        % Channel c finished
                        %fprintf('Done! c=%i, it=%i\n',c,it);
                        continue; 
                    end

                    % Compute gradient and Hessian    
                    gr_l = zeros(d(1:3),'single');
                    H_l  = zeros(d(1:3),'single');

                    % For each combination of missing voxels
                    for l=1:nL

                        % Get mask of missing modalities (with this particular code)        
                        observed_channels = spm_gmm_lib('code2bin', L(l), C);
                        missing_channels  = ~observed_channels;
                        if missing_channels(c), continue; end
                        if isempty(code), selected_voxels = ones(size(code), 'logical');
                        else,                   selected_voxels = (code == L(l));
                        end
                        nb_channels_missing  = sum(missing_channels);
                        nb_voxels_coded      = sum(selected_voxels);
                        if nb_voxels_coded == 0, continue; end

                        % Convert channel indices to observed indices
                        mapped_c     = 1:C;
                        mapped_c     = mapped_c(observed_channels);
                        mapped_c     = find(mapped_c == c);
                        cc           = mapped_c; % short alias

                        selected_obs = bffn(selected_voxels,observed_channels);
                        gi = 0; % Gradient accumulated accross clusters
                        Hi = 0; % Hessian accumulated accross clusters
                        for k=1:K1

                            % Compute expected precision (see GMM+missing data)
                            Voo = V(observed_channels,observed_channels,k);
                            Vom = V(observed_channels,missing_channels,k);
                            Vmm = V(missing_channels,missing_channels,k);
                            Vmo = V(missing_channels,observed_channels,k);
                            Ao  = Voo - Vom*(Vmm\Vmo);
                            Ao  = (n(k) - nb_channels_missing) * Ao;
                            MUo = m(observed_channels,k);

                            % Compute statistics
                            gk = bsxfun(@minus, selected_obs, MUo.') * Ao(cc,:).';
                            Hk = Ao(cc,cc);

                            selected_resp = zn(selected_voxels,k);
                            gk = bsxfun(@times, gk, selected_resp);
                            Hk = bsxfun(@times, Hk, selected_resp);
                            selected_resp = [];

                            % Accumulate across clusters
                            gi = gi + gk;
                            Hi = Hi + Hk;
                        end

                        % Multiply with bias corrected value (chain rule)
                        gi = gi .* selected_obs(:,cc);
                        Hi = Hi .* (selected_obs(:,cc).^2);
                        selected_obs = [];

                        % Normalisation term
                        gi = gi - 1;

                        % Accumulate across missing codes
                        gr_l(selected_voxels) = gr_l(selected_voxels) + gi;
                        H_l(selected_voxels)  = H_l(selected_voxels) + Hi;
                        selected_voxels = [];
                    end
                    zn = [];

                    d3 = numel(chan(c).T); % Number of DCT parameters
                    H  = zeros(d3,d3);     
                    gr = zeros(d3,1);      
                    for z=1:d(3)
                        b3 = double(chan(c).B3(z,:)');
                        gr = gr + kron(b3,spm_krutil(double(gr_l(:,:,z)),double(chan(c).B1),double(chan(c).B2),0));
                        H  = H  + kron(b3*b3',spm_krutil(double(H_l(:,:,z)),double(chan(c).B1),double(chan(c).B2),1));
                    end
                    b3 = [];                    

                    % Gauss-Newton update of bias field parameters
                    Update = reshape((H + chan(c).C)\(gr + chan(c).C*chan(c).T(:)),size(chan(c).T));
                    H = []; gr = [];

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

                        % Recompute responsibilities (with updated bias field)
                        zn = Responsibility(m,b,V,n,bffn,mun,L,code);

                        % Compute new lower bound
                        lx  = W*LowerBound('X',bffn,zn,code,{m,b},{V,n});            
                        lxb = W*LowerBound('XB',bf);

                        % Check new lower bound
                        if (lx + lxb + sum(pr_bf)) > (olx + olxb + sum(opr_bf))                                                                          
                            lb.XB(end + 1) = lxb;
                            lb.X(end  + 1) = lx;

                            nl = lx + lxb + sum(pr_bf);
                            ol = olx + olxb + sum(opr_bf);  
    %                         fprintf('it2=%i\tc=%i\tls=%i\tnl=%0.7f\tgain=%0.7f :o)\n',it_bf,c,ls,nl,nl - ol);
                            if (nl - ol) < nm*tol_bf
                                % Finished for channel c
                                done(c) = true;
                            end
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
                                pr_bf = opr_bf;
                                zn    = Responsibility(m,b,V,n,bffn,mun,L,code);
                            end
                        end
                    end
                    oT = []; Update = [];      
                end
            end   

            % Update datn     
            datn.bf.T = {chan(:).T};
            datn.E(3) = -sum(pr_bf);
        end
    end
    fn = []; bf = []; mun = [];

    if do_bf_norm
        [lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', bffn, zn, 1, {code,L});   
        datn.bf.lSS0     = lSS0;
        datn.bf.lSS1     = lSS1;
        datn.bf.lSS2     = lSS2;
        datn.bf.L        = L;

        % Get DC component
        dc    = struct;
        dc.ln = zeros(1,C);
        for c=1:C
            dc.ln(c) = chan(c).T(1,1,1);
        end
        dc.int     = SclFromBiasFieldDC(chan);
        datn.bf.dc = dc;
    end
end

if samp > 1 || nargout == 1
    % Compute responsibilities on original data
    fn   = spm_mb_io('GetData',datn.f);
    fn   = reshape(fn,[prod(df(1:3)) C]);
    fn   = Mask(fn);
    code = spm_gmm_lib('obs2code', fn);
    mun0 = reshape(mun0,[prod(df(1:3)) K]);
    if do_updt_bf
        % Get full-sized bias field
        chan = BiasFieldStruct(datn,C,df,reg,fwhm,[],datn.bf.T);
        bf   = BiasField(chan,df);
        bffn = bf.*fn;
        bf   = [];
    else
        bffn = fn;
    end
    fn   = [];
    L    = unique(code);
    mun0 = cat(2,mun0,zeros([prod(df(1:3)) 1],'single'));
    zn   = Responsibility(m,b,V,n,bffn,mun0,L,code);
end       

if nargout > 1
    % Get 4D versions of K1 - 1 classes
    zn = reshape(zn(:,1:K),[df(1:3) K]);

    % Update datn     
    datn.E(1)     = -lb.sum(end);
    datn.mog.po.m = m;
    datn.mog.po.b = b;
    datn.mog.po.V = V;
    datn.mog.po.n = n;          
    datn.mog.lb   = lb;     
end
end
%==========================================================================

%==========================================================================
% UpdatePrior()
function dat = UpdatePrior(dat, mu, sett)

% Parse function settings
do_bf_norm  = sett.do.bf_norm;
do_updt_int = sett.do.updt_int;

if ~isfield(dat(1),'mog'), return; end

% Get population indices
p_ix = GetPopulationIdx(dat);

for p=1:numel(p_ix) % Loop over populations
    if do_bf_norm
        % Zero-mean the bias field DC component
        dat(p_ix{p}) = ZeroMeanDC(dat(p_ix{p}),mu,sett);
    end

    if do_updt_int
        N = numel(p_ix{p}); % Number of subjects in population
        
        % Get all posteriors
        po = cell(1,N);
        for n=p_ix{p}
            po{n}{1}{1} = dat(n).mog.po.m;
            po{n}{1}{2} = dat(n).mog.po.b;
            po{n}{2}{1} = dat(n).mog.po.V;
            po{n}{2}{2} = dat(n).mog.po.n;
        end
        
        % Get old prior
        pr = dat(p_ix{p}(1)).mog.pr;
        pr = {pr.m,pr.b,pr.V,pr.n};
        
        % Update prior
        pr = DoIntensityPriorUpdate(pr,po);
                     
        % Assign new prior
        for n=p_ix{p}
            dat(n).mog.pr.m = pr{1};
            dat(n).mog.pr.b = pr{2};
            dat(n).mog.pr.V = pr{3};
            dat(n).mog.pr.n = pr{4};
        end             
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
function f = ApplyMask(f)
f(~isfinite(f) | f == 0 | f == min(f(:))) = NaN;
end
%==========================================================================

%==========================================================================
% DoIntensityPriorUpdate()
function pr = DoIntensityPriorUpdate(pr,po)
m0 = pr{1};
b0 = pr{2};
W0 = pr{3};
n0 = pr{4};

N = numel(po);
D = size(m0,1);
K = size(m0,2);

for k=1:K % loop over classes
    
    % Compute m_0
    g = zeros(D,1);
    H = zeros(D,D);
    for i=1:N    
        [m,b,W,n] = GetPosterior(po,i);

        g = g + b0(k)*n(k)*W(:,:,k)*(m(:,k)-m0(:,k));
        H = H + b0(k)*n(k)*W(:,:,k);
    end
    m0(:,k) = m0(:,k) + H\g;



    % Compute \beta_0
    g_const = 0;
    for i=1:N
        [m,b,W,n] = GetPosterior(po,i);

        g_const = g_const - 0.5*(D/b(k) + n(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)));
    end
    for subit=1:100
        % Diff w.r.t. b0
        g  = 0.5*N*D/b0(k) + g_const; % Gradient
        H  = 0.5*N*D/b0(k)^2;         % Hessian
        b0(k) = max(b0(k) + H\g,1e-5);
        if norm(g)==0, break; end
    end

    % Set up some constants
    nW_const = zeros(D);
    for i=1:N
        [m,b,W,n] = GetPosterior(po,i);

        nW_const = nW_const + n(k)*W(:,:,k);
    end

    ElogLam = N*D*log(2);
    for i=1:N
        [m,b,W,n] = GetPosterior(po,i);

        ElogLam = ElogLam + 2*sum(log(diag(chol(W(:,:,k)))));
        for j=1:D
            ElogLam = ElogLam + psi((n(k)+1-j)/2);
        end
    end

    % convergence = [];
    E = -realmax;

    for it=1:1000    
        % Compute objective function (Equation 10.74 of Bishop)

        oE = E;
        logB = -n0(k)*sum(log(diag(chol(W0(:,:,k))))) - n0(k)*D/2*log(2) - D*(D-1)/4*log(pi);
        for j=1:D
            logB = logB - gammaln((n0(k)+1-j)/2); 
        end
        E = (0.5*D*log(b0(k)/(2*pi)) + logB)*N + 0.5*(n0(k)-D-1)*ElogLam;
        for i=1:N
            [m,b,W,n] = GetPosterior(po,i);

            e = 0.5*(-D*b0(k)/b(k) - b0(k)*n(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)))...
              - 0.5*n(k)*trace(W0(:,:,k)\W(:,:,k));
            E = E + e;
        end
        %if E-oE<abs(E)*eps*D^2, break; end
        if E-oE==0, break; end

    %     convergence = [convergence E];
    %     plot(convergence,'.-'); drawnow;


        % Compute \nu_0

        % Objective function terms containing n0:
        % NlogB = -n0*N*(sum(log(diag(chol(W0)))) + D/2*log(2));
        % for j=1:D, NlogB = NlogB - N*gammaln((n0+1-j)/2); end
        % E = NlogB + n0*0.5*ElogLam

        g = (sum(log(diag(chol(W0(:,:,k))))) + D/2*log(2))*N - 0.5*ElogLam;
        H = 0;
        for j=1:D
            g = g + N*psi(  (n0(k)+1-j)/2)/2;
            H = H + N*psi(1,(n0(k)+1-j)/2)/4;
        end
        n0(k) = max(n0(k) - H\g,D-0.99999);

        % Compute W_0

        % Objective function terms containing W0:
        % E = -n0*N*sum(log(diag(chol(W0))));
        % for i=1:N
        %    E = E - 0.5*n(i)*trace(W0\W(:,:,i));
        % end

        C = inv(chol(W0(:,:,k)));

        % Objective function terms containing W0, after
        % re-expressing using C = inv(chol(W0)):
        % E = n0*N*sum(log(diag(C)));
        % for i=1:N
        %    E = E - 0.5*n(i)*trace(C'*W(:,:,i)*C);
        % end

        G  = -n0(k)*N*diag(1./diag(C)) + nW_const*C;
        for d=1:D
            c        = C(1:d,d);
            g        = G(1:d,d);
            H        = nW_const(1:d,1:d);
            H(d,d)   = H(d,d) + n0(k)*N/c(d)^2;
            C(1:d,d) = c - H\g;
        end
        C         = inv(C);
        W0(:,:,k) = C'*C;
    end
end

% Assign new prior
pr{1} = m0;
pr{2} = b0;
pr{3} = W0;
pr{4} = n0;
end
%==========================================================================

%==========================================================================
% GetPosterior()
function [m,b,W,n] = GetPosterior(po,i)
m = po{i}{1}{1};
b = po{i}{1}{2};
W = po{i}{2}{1};
n = po{i}{2}{2};
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
if strcmpi(type,'XB')    
    bf = varargin{1};
    
    bf(isnan(bf)) = 1;        
    bf            = log(prod(bf,2)); 
    lb            = sum(double(bf));
elseif strcmpi(type,'X')    
    fn   = varargin{1};
    zn   = varargin{2};
    code = varargin{3};
    mean = varargin{4};
    prec = varargin{5};
    L    = unique(code);
    
    [lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', fn, zn, 1, {code,L});
    lb               = spm_gmm_lib('MarginalSum', lSS0, lSS1, lSS2, mean, prec, L);
else
    error('Undefined type!');
end
end
%==========================================================================

%==========================================================================    
% PriorGMM()
function pr = PriorGMM(mx,mn,vr,K,sett)

% Parse function settings
fig_name = sett.show.figname_int;

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
%     mu(c,:)  = (0:(K - 1))'*mmx(c)/(1.5*K);
%     A(c,c,:) = mmx(c)/(1.5*K);    
    vrc      = mvr(c)/(K + 1);
    mnc      = mn(c);
    sd       = sqrt(vrc);
    mu(c,:)  = abs(linspace(mnc - n*sd,mnc + n*sd,K));    
    A(c,c,:) = vrc;        
    A(c,c,:) = 1/A(c,c,:);
end   

pr   = struct('m',[],'b',[],'n',[],'V',[]);
pr.m = mu;
pr.b = ones(1,K);
pr.n = C*ones(1,K);
pr.V = bsxfun(@times, A, reshape(pr.n, [1 1 K])); % Expected precision

if 0
    spm_gmm_lib('plot','gaussprior',{pr.m,pr.b,pr.V,pr.n},[],fig_name);
end

end
%==========================================================================   

%==========================================================================    
% PosteriorGMM()
function [po,mx,mn,vr] = PosteriorGMM(fn,K,sett)

% Parse function settings
fig_name = sett.show.figname_int;

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

if 0
    spm_gmm_lib('plot','gaussprior',{po.m,po.b,po.V,po.n},[],fig_name);
end

end
%==========================================================================  

%==========================================================================
% SclFromBiasFieldDC()
function scl = SclFromBiasFieldDC(chan)
C   = numel(chan);
scl = zeros(1,C);
for c=1:C
    b1 = chan(c).B1(1,1);
    b2 = chan(c).B2(1,1);
    b3 = chan(c).B3(1,1);
    t1 = chan(c).T(1,1,1);
    
    scl(c) = exp(b1*b2*b3*t1);
end
end
%==========================================================================

%==========================================================================
% SubSample()
function varargout = SubSample(samp,Mat,varargin)
vx        = sqrt(sum(Mat(1:3,1:3).^2));
samp      = max([1 1 1],round(samp*[1 1 1]./vx));
N         = numel(varargin);
varargout = cell(1,N + 2);
for n=1:N
    f  = varargin{n};    
    d0 = size(f);
    d0 = [d0 1];
    
    f = f(1:samp:end,1:samp:end,1:samp:end,:);    
    d = size(f); 
    d = [d 1];    
    
    varargout{n} = f; 
end

% For weighting data parts of lowerbound with factor based on amount of
% downsampling  
W                = prod(d0(1:3))/prod(d(1:3));
varargout{N + 1} = W;
varargout{N + 2} = d(1:3);
end
%==========================================================================

%==========================================================================
% SumLowerBound()
function lb = SumLowerBound(lb)
fields          = fieldnames(lb);
lb.sum(end + 1) = 0;
for i=1:numel(fields)
    field = fields{i};
    if ~any(strcmpi(field, {'sum' 'last'})) && ~isempty(lb.(field)) && ~isnan(lb.(field)(end))
        lb.sum(end) = lb.sum(end) + sum(lb.(field)(:,end));
    end
end
end
%==========================================================================

%==========================================================================
% TransformBF()
function t = TransformBF(B1,B2,B3,T)
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
% ZeroMeanDC()
function dat = ZeroMeanDC(dat,mu,sett)

% Get correction factor
sm_dc_ln  = 0;
sm_dc_int = 0;
for n=1:numel(dat)
    sm_dc_ln  = sm_dc_ln  + dat(n).bf.dc.ln;
    sm_dc_int = sm_dc_int + dat(n).bf.dc.int;
end   
mn_dc_ln  = sm_dc_ln./numel(dat);
mn_dc_int = sm_dc_int./numel(dat);
C         = numel(mn_dc_ln); 

if 0     
    % Some verbose (should tend towards 1, for all channels)
    fprintf('mn_dc_int = [');
    for c=1:C - 1
        fprintf('%4.5f ',mn_dc_int(c));
    end
    fprintf('%4.5f',mn_dc_int(C))
    fprintf(']\n');
end

% Adjust bias field DC component
scl = (1./mn_dc_int)';
for n=1:numel(dat)
    for c=1:C
        dat(n).bf.T{c}(1,1,1) = dat(n).bf.T{c}(1,1,1) - mn_dc_ln(c);
    end
end
   
% Correct GMM posterior parameters and lower bound/objective
for n=1:numel(dat)
    if all(dat(n).do_bf == false), continue; end
    
    dat(n) = ZeroMeanDCSub(dat(n),mu,scl,sett);
end
end
%==========================================================================

%==========================================================================
% ZeroMeanDCSub()
function datn = ZeroMeanDCSub(datn,mu,scl,sett)

% Parse function settings
B    = sett.registr.B;
fwhm = sett.bf.fwhm;
Mmu  = sett.var.Mmu;
reg  = sett.bf.reg;
samp = sett.gen.samp;

[df,C] = spm_mb_io('GetSize',datn.f);
Mat    = datn.Mat;
q      = double(datn.q);
Mn     = datn.Mat;
W      = 1;
Mr     = spm_dexpm(q,B);

% Get image(s)
fn = spm_mb_io('GetData',datn.f);
if samp > 1      
    [fn,W,d] = SubSample(samp,Mat,fn);
else
    d = df;
end
fn = reshape(fn,[prod(d(1:3)) C]);

% GMM posterior
m = datn.mog.po.m;
b = datn.mog.po.b;
V = datn.mog.po.V;
n = datn.mog.po.n;

% GMM prior
m0 = datn.mog.pr.m;
b0 = datn.mog.pr.b;
V0 = datn.mog.pr.V;
n0 = datn.mog.pr.n;

% Get suffstats
lSS0 = datn.bf.lSS0;
lSS1 = datn.bf.lSS1;
lSS2 = datn.bf.lSS2;
L    = datn.bf.L;

% Rescale suffstats
K1 = size(m,2);
K  = K1 - 1;
A  = bsxfun(@times,V,reshape(n,[1 1 K1]));
for l=2:numel(L)
    obs = spm_gmm_lib('code2bin', L(l), C);
    lSS1{l} = bsxfun(@times,lSS1{l},scl(obs));
    for k=1:size(lSS2{l},3)
        lSS2{l}(:,:,k) = (scl(obs)*scl(obs)').*lSS2{l}(:,:,k);
    end     
end

% Update GMM posteriors
[SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', lSS0, lSS1, lSS2, {m,A}, L);        
[m,~,b,V,n]   = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2, {m0,b0,V0,n0});    

% Save updated GMM posteriors
datn.mog.po.m = m;
datn.mog.po.b = b;
datn.mog.po.V = V;
datn.mog.po.n = n;

% Get subject-space template (softmaxed K + 1)
psi1 = spm_mb_io('GetData',datn.psi);
psi0 = spm_mb_shape('Affine',df,Mmu\Mr*Mn);
psi  = spm_mb_shape('Compose',psi1,psi0);
psi0 = []; psi1 = [];  
mu   = spm_mb_shape('Pull1',mu,psi);
psi  = [];

if samp > 1      
    mu = SubSample(samp,Mat,mu);
end
mu = reshape(mu,[prod(d(1:3)) K]);
mu = cat(2,mu,zeros([prod(d(1:3)) 1],'single'));

% Get bias field
chan       = BiasFieldStruct(datn,C,df,reg,fwhm,[],datn.bf.T,samp);
[bf,pr_bf] = BiasField(chan,d);

% Get responsibilities
fn   = Mask(fn);
code = spm_gmm_lib('obs2code', fn);
zn   = Responsibility(m,b,V,n,bf.*fn,mu,L,code);
mu   = [];

% Update lower bound/objective
lx  = W*LowerBound('X',bf.*fn,zn,code,{m,b},{V,n});
lxb = W*LowerBound('XB',bf);

datn.mog.lb.XB(end + 1) = lxb;
datn.mog.lb.X(end  + 1) = lx;
datn.mog.lb             = SumLowerBound(datn.mog.lb);

datn.E(1) = -datn.mog.lb.sum(end);
datn.E(3) = -sum(pr_bf);
end
%==========================================================================