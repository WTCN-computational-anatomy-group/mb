function varargout = spm_mb_show(varargin)
%__________________________________________________________________________
%
% Functions for visualisation related.
%
% FORMAT spm_mb_show('Clear',sett)
% FORMAT spm_mb_show('All',dat,mu,Objective,N,sett)
% FORMAT spm_mb_show('IntensityPrior',dat,sett,p)
% FORMAT spm_mb_show('Model',mu,Objective,N,sett)
% FORMAT spm_mb_show('Subjects',dat,mu,sett,p,show_extras)
% FORMAT spm_mb_show('Tissues',im,do_softmax,num_montage,fig_nam)
% FORMAT spm_mb_show('Speak',nam,sett,varargin)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_mb_show
    error('Not enough argument. Type ''help spm_mb_show'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id    
    case 'All'
        [varargout{1:nargout}] = All(varargin{:});
    case 'Clear'
        [varargout{1:nargout}] = Clear(varargin{:});        
    case 'IntensityPrior'
        [varargout{1:nargout}] = IntensityPrior(varargin{:});         
    case 'Model'
        [varargout{1:nargout}] = Model(varargin{:});          
    case 'Subjects'
        [varargout{1:nargout}] = Subjects(varargin{:});         
    case 'Tissues'
        [varargout{1:nargout}] = Tissues(varargin{:});        
    case 'Speak'
        [varargout{1:nargout}] = Speak(varargin{:});                
    otherwise
        help spm_mb_show
        error('Unknown function %s. Type ''help spm_mb_show'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% All()
function All(dat,mu,Objective,N,sett)
if sett.show.level >= 2
    % Template and negative loglikelihood
    Model(mu,Objective,N,sett);
end
if sett.show.level >= 3
    % Template space images
    Subjects(dat,mu,sett,false,false); 
end
if sett.show.level == 4 || sett.show.level == 5
    % Segmentations and warped template
    Subjects(dat,mu,sett,false);    
end
if sett.show.level >= 5
    % Intensity prior fit   
    IntensityPrior(dat,sett);
end
if sett.show.level >= 6
    % Subject parameters
    Subjects(dat,mu,sett,true);        
end
end
%==========================================================================

%==========================================================================
% Clear()
function Clear(sett)
fn = {sett.show.figname_bf, sett.show.figname_int, sett.show.figname_model, ...
      sett.show.figname_subjects, sett.show.figname_parameters,sett.show.figname_imtemplatepace};
for i=1:numel(fn)
    f = findobj('Type', 'Figure', 'Name', fn{i});
    if ~isempty(f), clf(f); drawnow; end
    p = 1;
    while true
        fnp = [fn{i} ' (p=' num2str(p) ')'];
        f   = findobj('Type', 'Figure', 'Name', fnp);
        if ~isempty(f) 
            clf(f); 
            drawnow; 
        else 
            break
        end    
        p = p + 1;
    end
end
end
%==========================================================================

%==========================================================================
% IntensityPrior()
function IntensityPrior(dat,sett)
    p_ix = spm_mb_appearance('GetPopulationIdx',dat);
    Np   = numel(p_ix);
    for p=1:Np
        if Np == 1, pp = []; 
        else,       pp = p;
        end
        ShowIntensityPrior(dat(p_ix{p}),sett,pp);
    end
end
%==========================================================================

%==========================================================================
% Model()
function Model(mu,Objective,N,sett)

% Parse function settings
fig_name = sett.show.figname_model;

mu  = spm_mb_shape('TemplateK1',mu,4);
mu  = exp(mu);   
nam = ['K1=' num2str(size(mu,4)) ', N=' num2str(N) ' (softmaxed)'];
if size(mu,3) > 1
    % 3D
    ShowCat(mu,1,2,3,1,fig_name);
    ShowCat(mu,2,2,3,2,fig_name);
    title(nam);
    ShowCat(mu,3,2,3,3,fig_name,true); % true -> show colorbar    
    subplot(2,1,2);     
else
    % 2D
    ShowCat(mu,3,1,2,1,fig_name,true); % true -> show colorbar
    title(nam);
    subplot(1,2,2); 
end
plot(Objective,'.-');        
title('Negative log-likelihood')
drawnow
end
%==========================================================================

%==========================================================================
% Speak()
function Speak(nam,sett,varargin)

% Parse function settings
do_updt_int      = sett.do.updt_int;
do_updt_template = sett.do.updt_template;
nit              = sett.nit.init;
show_level       = sett.show.level;

if show_level < 1, return; end

switch nam
    case 'Finished'
        t = varargin{1}; 
        
        fprintf('Algorithm finished in %.1f seconds.\n', t);
    case 'InitAff'
        fprintf('Optimising parameters at largest zoom level (nit = %i)\n',nit)
    case 'Iter'
        nzm = varargin{1}; 
        
        fprintf('Optimising parameters at decreasing zoom levels (nzm = %i)\n',nzm)        
    case 'Start'
        N = varargin{1};
        K = varargin{2};        

        fprintf('=========================================================================\n')
        fprintf('| Algorithm starting (N = %i, K = %i, updt_intprior = %i, updt_template = %i)\n',N,K,do_updt_int,do_updt_template)
        fprintf('=========================================================================\n\n')
    otherwise
        error('Unknown input!')
end
end
%==========================================================================

%==========================================================================
% ShowIntensityPrior()
function ShowIntensityPrior(dat,sett,p)
if nargin < 3, p = []; end

% Parse function settings
fig_name = sett.show.figname_int;
if ~isempty(p), fig_name = [fig_name ' (p=' num2str(p) ')']; end
mg_ix    = sett.model.mg_ix;

if ~isfield(dat(1),'mog'), return; end

n  = 1;
m0 = dat(n).mog.pr.m;
b0 = dat(n).mog.pr.b;
W0 = dat(n).mog.pr.W;
n0 = dat(n).mog.pr.n;

spm_gmm_lib('plot','gaussprior',{m0,b0,W0,n0},mg_ix,fig_name);
end
%==========================================================================

%==========================================================================
% ShowNativeSubjects()
function ShowNativeSubjects(dat,mu0,sett,p,show_extras)
if nargin < 4, p           = [];   end
if nargin < 5, show_extras = true; end

% Parse function settings
axis_3d       = sett.show.axis_3d;
B             = sett.registr.B;
c             = sett.show.channel;
fig_name_bf   = sett.show.figname_bf;
fig_name_par  = sett.show.figname_parameters;
fig_name_tiss = sett.show.figname_subjects;
fwhm          = sett.bf.fwhm;
mg_ix         = sett.model.mg_ix;
Mmu           = sett.var.Mmu;
mx_subj       = sett.show.mx_subjects;
reg           = sett.bf.reg;

if ~isempty(p), fig_name_bf   = [fig_name_bf ' (p=' num2str(p) ')']; end
if ~isempty(p), fig_name_par  = [fig_name_par ' (p=' num2str(p) ')']; end
if ~isempty(p), fig_name_tiss = [fig_name_tiss ' (p=' num2str(p) ')']; end

if isfield(dat,'mog'), nr_tiss = 3;
else,                  nr_tiss = 2;
end
nr_bf  = 3;
nr_par = 4;
if ~isfield(dat,'mog'), nr_par = nr_par - 2; end

[df,~] = spm_mb_io('GetSize',dat(1).f);
if df(3) == 1, ax = 3;
else,          ax = axis_3d;   
end
    
clr = {'r','g','b','y','m','c',};
K   = size(mu0,4);
K1  = K + 1;
Kmg = numel(mg_ix);
nd  = min(numel(dat),mx_subj);
for n=1:nd
    % Parameters
    [df,C] = spm_mb_io('GetSize',dat(n).f);          
    c      = min(c,C);    
    q      = double(dat(n).q);
    Mr     = spm_dexpm(q,B);
    Mn     = dat(n).Mat;
    do_bf  = dat(n).do_bf;
    is_ct  = dat(n).is_ct;    
    
    % Warp template
    psi1 = spm_mb_io('GetData',dat(n).psi);
    psi  = spm_mb_shape('Compose',psi1,spm_mb_shape('Affine',df,Mmu\Mr*Mn));
    mun  = spm_mb_shape('Pull1',mu0,psi);            
    psi1 = [];
    
    % Get bias field    
    if isfield(dat(n),'mog') && any(do_bf == true)
        chan = spm_mb_appearance('BiasFieldStruct',dat(n),C,df,reg,fwhm,[],dat(n).bf.T);
        bf   = spm_mb_appearance('BiasField',chan,df);        
    else
        bf = ones([1 C]);
    end  
    
    % Get template (K + 1)
    mun = spm_mb_shape('TemplateK1',mun,4);
    
    % Get segmentation
    if isfield(dat,'mog')
        fn   = spm_mb_io('GetData',dat(n).f);         
        fn   = reshape(fn,[prod(df(1:3)) C]);
        fn   = spm_mb_appearance('Mask',fn,is_ct);
        
        % Integrate labels and multiple Gaussians per tissue
        labels = spm_mb_appearance('GetLabels',dat(n),sett);
        mg_w   = dat(n).mog.mg_w;
        
        [bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', bf.*fn);
        mun                       = reshape(mun,[prod(df(1:3)) K + 1]);
        mun1                      = mun + labels;
        labels                    = [];
        mun1                      = mun1(:,mg_ix);
        mun1                      = mun1 + log(mg_w);        
        mun1                      = spm_gmm_lib('obs2cell', mun1, code_image, false);
                    
        % Get responsibility
        zn   = spm_mb_appearance('Responsibility',dat(n).mog.po.m,dat(n).mog.po.b, ...
                                    dat(n).mog.po.W,dat(n).mog.po.n,bffn,mun1,msk_chn);           
        zn   = spm_gmm_lib('cell2obs', zn, code_image, msk_chn);                
        mun1 = [];
        
        % Just to insert NaNs..
        mun     = spm_gmm_lib('obs2cell', mun, code_image, false);
        mun     = spm_gmm_lib('cell2obs', mun, code_image, msk_chn);                
        msk_chn = []; 
        
        % Reshape back
        zn  = reshape(zn,[df(1:3) Kmg]);
        fn  = reshape(fn,[df(1:3) C]);
        mun = reshape(mun,[df(1:3) K + 1]);
        
        % If using multiple Gaussians per tissue, collapse so that zn is of
        % size K1
        if Kmg > K1
            for k=1:K1, zn(:,:,:,k) = sum(zn(:,:,:,mg_ix==k),4); end
            zn(:,:,:,K1 + 1:end)    = [];
        end
    else
        zn = spm_mb_io('GetData',dat(n).f); 
        zn = cat(4,zn,1 - sum(zn,4));
    end    
    
    % Softmax template
    mun = exp(mun);

    if isfield(dat(n),'mog') && any(do_bf == true)
        bf = reshape(bf,[df C]);
    else
        bf = reshape(bf,[1 1 1 C]);
    end  
    
    % Show template, segmentation
    ShowCat(mun,ax,nr_tiss,nd,n,fig_name_tiss);
    ShowCat(zn,ax,nr_tiss,nd,n + nd,fig_name_tiss);        
    if isfield(dat,'mog')
        % and image (if using GMM)     
        ShowIm(bf(:,:,:,c).*fn(:,:,:,c),ax,nr_tiss,nd,n + 2*nd,fig_name_tiss,true,is_ct);        
    end
    zn = [];
    
    if show_extras
        % Show bias field fit, velocities, affine parameters, GMM fit,
        % lower bound
        if isfield(dat,'mog') && any(do_bf == true)
            % Show bias field
            ShowIm(fn(:,:,:,c),ax,nr_bf,nd,n,fig_name_bf,false,is_ct);
            ShowIm(bf(:,:,:,c),ax,nr_bf,nd,n + nd,fig_name_bf,false,is_ct);            
            ShowIm(bf(:,:,:,c).*fn(:,:,:,c),ax,nr_bf,nd,n + 2*nd,fig_name_bf,false,is_ct);
        end

        % Now show some other stuff
        fg = findobj('Type', 'Figure', 'Name', fig_name_par);
        if isempty(fg), fg = figure('Name', fig_name_par, 'NumberTitle', 'off'); end
        set(0, 'CurrentFigure', fg);   

        % Affine parameters    
        q    = spm_imatrix(Mr);            
        q    = q([1 2 6]);
        q(3) = 180/pi*q(3);
        q    = abs(q);
        sp = subplot(nr_par,nd,n);
        cla(sp); % clear subplot
        hold on
        for k=1:numel(q)
            bar(k, q(k), clr{k});
        end
        box on
        hold off            
        set(gca,'xtick',[])  

        % Velocities (x)
        v = spm_mb_io('GetData',dat(n).v);
        ShowIm(v(:,:,:,1),ax,nr_par,nd,n + 1*nd,fig_name_par)

        % Intensity histogram w GMM fit
        if isfield(dat,'mog')                
            % Here we get approximate class proportions from the (softmaxed K + 1)
            % tissue template
            mun       = reshape(mun,[prod(df(1:3)) size(mun,4)]);
            msk       = isfinite(mun);
            mun(~msk) = 0;
            mun       = sum(mun,1);
            mun       = mun./sum(mun);
            mun       = mun(mg_ix).*mg_w;
            
            % Plot GMM fit        
            ShowGMMFit(bf(:,:,:,c).*fn(:,:,:,c),mun,dat(n).mog,nr_par,nd,n + 2*nd,c,mg_ix);

            % Lower bound
            subplot(nr_par,nd,n + 3*nd) 
            plot(dat(n).mog.lb.sum,'-');
            axis off
        end    
    end
end
drawnow
end
%==========================================================================

%==========================================================================
% ShowTemplateSubjects()
function ShowTemplateSubjects(dat,mu,sett,p)

% Parse function settings
B        = sett.registr.B;
c        = sett.show.channel;
dmu      = sett.var.d;
fig_name = sett.show.figname_imtemplatepace;
fwhm     = sett.bf.fwhm;
Mmu      = sett.var.Mmu;
mx_subj  = sett.show.mx_subjects;
reg      = sett.bf.reg;

if ~isempty(p), fig_name   = [fig_name ' (p=' num2str(p) ')']; end

if size(mu,3) > 1, nr = 3;
else,              nr = 1;
end

nd = min(numel(dat),mx_subj);
for n=1:nd
    % Parameters
    [df,C] = spm_mb_io('GetSize',dat(n).f);          
    c      = min(c,C);    
    q      = double(dat(n).q);
    Mr     = spm_dexpm(q,B);
    Mn     = dat(n).Mat;        
    do_bf  = dat(n).do_bf;
    is_ct  = dat(n).is_ct;
    
    % Get forward deformation
    psi0 = spm_mb_io('GetData',dat(n).psi);
    psi  = spm_mb_shape('Compose',psi0,spm_mb_shape('Affine',df,Mmu\Mr*Mn));  
    psi0 = [];    
    if df(3) == 1, psi(:,:,:,3) = 1; end % 2D
    
    % Bias field
    if isfield(dat(n),'mog') && any(do_bf == true)
        chan = spm_mb_appearance('BiasFieldStruct',dat(n),C,df,reg,fwhm,[],dat(n).bf.T);
        bf   = spm_mb_appearance('BiasField',chan,df);  
        bf   = reshape(bf,[df(1:3) C]);
        bf   = bf(:,:,:,c);
    else
        bf = 1;
    end  
    
    % Warp observed data to template space
    fn = spm_mb_io('GetData',dat(n).f);
    sd = spm_mb_shape('SampDens',Mmu,Mn);
    if isfield(dat(n),'mog')
        fn       = bf.*fn(:,:,:,c);
        bf       = [];
        [fn,cnt] = spm_mb_shape('Push1',fn,psi,dmu,sd);
    else
        [fn,cnt] = spm_mb_shape('Push1',fn,psi,dmu,sd);
    end
    fn = fn./(cnt + eps('single'));
    
    % Show
    if isfield(dat(n),'mog')
        ShowIm(fn,3,nr,nd,n,fig_name,true,is_ct);         
        if size(mu,3) > 1
            ShowIm(fn,2,nr,nd,n + nd,fig_name,true,is_ct);   
            ShowIm(fn,1,nr,nd,n + 2*nd,fig_name,true,is_ct);        
        end
    else
        fn = cat(4,fn,1 - sum(fn,4));
        ShowCat(fn,3,nr,nd,n,fig_name);         
        if size(mu,3) > 1
            ShowCat(fn,2,nr,nd,n + nd,fig_name);   
            ShowCat(fn,1,nr,nd,n + 2*nd,fig_name);        
        end
    end
end
drawnow
end
%==========================================================================

%==========================================================================
% Subjects()
function Subjects(dat,mu,sett,show_extras,show_native)
if nargin < 4, show_extras = true; end
if nargin < 5, show_native = true; end

p_ix = spm_mb_appearance('GetPopulationIdx',dat);
Np   = numel(p_ix);
for p=1:Np
    if Np == 1, pp = []; 
    else,       pp = p;
    end
    if show_native
        % Subject space information
        ShowNativeSubjects(dat(p_ix{p}),mu,sett,pp,show_extras);
    else
        % Template space information
        ShowTemplateSubjects(dat(p_ix{p}),mu,sett,pp);        
    end
end
end
%==========================================================================

%==========================================================================
% Tissues()
function Tissues(im,do_softmax,num_montage,fig_nam)
if nargin < 2, do_softmax  = true; end
if nargin < 3, num_montage = 20; end
if nargin < 4, fig_nam     = '(spm_mb) Tissues'; end

f  = findobj('Type', 'Figure', 'Name', fig_nam);
if isempty(f)
    f = figure('Name', fig_nam, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f); 

if ischar(im),      im = nifti(im); end
if isa(im,'nifti'), im = im.dat();  end

dm = size(im);
K  = dm(4);

if do_softmax
    im = cat(4,im,zeros(dm(1:3),'single'));
    im = spm_mb_shape('Softmax',im,4); 
    K        = K + 1;
end

nr = floor(sqrt(K));
nc = ceil(K/nr);  
num_montage = 1:num_montage:dm(3);

for k=1:K
    subplot(nr,nc,k)
    montage(im(:,:,num_montage,k))
    axis off image xy;   
    colormap(gray)
end
drawnow
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% ShowCat()
function ShowCat(in,ax,nr,nc,np,fn,show_colorbar)
if nargin < 7, show_colorbar = false; end

f = findobj('Type', 'Figure', 'Name', fn);
if isempty(f)
    f = figure('Name', fn, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f); 

subplot(nr,nc,np);

dm = size(in);
if numel(dm) ~= 4
    dm = [dm 1 1];
end
K   = dm(4);
pal = hsv(K);
mid = ceil(dm(1:3).*0.5);

if ax == 1
    in = permute(in(mid(1),:,:,:),[3 2 1 4]);
elseif ax == 2
    in = permute(in(:,mid(2),:,:),[1 3 2 4]);
else
    in = in(:,:,mid(3),:);
end

tri = false;
if numel(size(in)) == 4 && size(in, 3) == 1
    tri = true;
    dm  = [size(in) 1 1];
    in   = reshape(in, [dm(1:2) dm(4)]);
end
if isa(pal, 'function_handle')
    pal = pal(size(in,3));
end

dm = [size(in) 1 1];
c   = zeros([dm(1:2) 3]); % output RGB image
s   = zeros(dm(1:2));     % normalising term

for k=1:dm(3)
    s = s + in(:,:,k);
    color = reshape(pal(k,:), [1 1 3]);
    c = c + bsxfun(@times, in(:,:,k), color);
end
if dm(3) == 1
    c = c / max(1, max(s(:)));
else
    c = bsxfun(@rdivide, c, s);
end

if tri
    c = reshape(c, [size(c, 1) size(c, 2) 1 size(c, 3)]);
end

c = squeeze(c(:,:,:,:));
imagesc(c); axis off image xy;   
colormap(pal);

if show_colorbar    
    cb = colorbar;
    set(gca, 'clim', [0.5 K+0.5]);
    set(cb, 'ticks', 1:K, 'ticklabels', 1:K); 
end
end
%==========================================================================

%==========================================================================
% ShowGMMFit()
function ShowGMMFit(f,PI,mog,nr,nd,n,c,mg_ix)
K      = size(mog.po.m,2);
colors = hsv(max(mg_ix));

sp = subplot(nr,nd,n);
cla(sp); % clear subplot
hold on

% ---------
% Histogram

% Input is a list of observations
[H, edges] = histcounts(f(f > 0 & isfinite(f)), 64, 'Normalization', 'pdf');
centres = (edges(1:end-1) + edges(2:end))/2;
bar(centres, H, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);

% -----------
% GMM Density
A     = bsxfun(@times, mog.po.W, reshape(mog.po.n, [1 1 K])); % Expected precision
xlims = [inf -inf];
for k=1:K
    MU   = mog.po.m(c,k);    
    sig2 = inv(A(c,c,k));
    
    x = linspace(MU - 3*sqrt(sig2), MU + 3*sqrt(sig2),100);
    y = PI(k)*spm_Npdf(x, MU, sig2);
    plot(x, y, 'Color', colors(mg_ix(k),:), 'LineWidth', 1)
    xlims = [min([xlims(1) x]) max([xlims(2) x])];
end

ymax = max(H);
xlim(xlims);
if ymax == 0
    ylim([0 1.1]);
else
    ylim([0 1.1*ymax]);
end

box on
hold off
set(gca,'ytick',[])  
end
%==========================================================================
        
%==========================================================================
% ShowIm()
function ShowIm(in,ax,nr,nc,np,fn,use_gray,is_ct)
if nargin < 7, use_gray = true; end
if nargin < 8, is_ct    = false; end

f  = findobj('Type', 'Figure', 'Name', fn);
if isempty(f)
    f = figure('Name', fn, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   

subplot(nr,nc,np);

dm  = size(in);
dm  = [dm 1];
mid = ceil(dm(1:3).*0.5);

if ax == 1
    in = permute(in(mid(1),:,:,:),[3 2 1 4]);
elseif ax == 2
    in = permute(in(:,mid(2),:,:),[1 3 2 4]);
else
    in = in(:,:,mid(3),:);
end

in = squeeze(in(:,:,:,:));
if is_ct
    % Assume CT..
    imagesc(in,[0 100]);
else
    imagesc(in);
end
axis off image xy;
if use_gray
    colormap(gray);
end
end 
%==========================================================================