function varargout = spm_multireg_show(varargin)
%__________________________________________________________________________
%
% Visualisation functions for spm_multireg.
%
% FORMAT spm_multireg_show('Clear',sett)
% FORMAT spm_multireg_show('ShowAll',dat,mu,Objective,N,sett)
% FORMAT spm_multireg_show('ShowBiasField',dat,sett)
% FORMAT spm_multireg_show('ShowIntensityModel',dat,sett)
% FORMAT spm_multireg_show('ShowModel',mu,Objective,N,sett)
% FORMAT spm_multireg_show('ShowParameters',dat,mu,sett)
% FORMAT spm_multireg_show('ShowSubjects',dat,mu,sett)
% FORMAT spm_multireg_show('Speak',nam,varargin)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_show
    error('Not enough argument. Type ''help spm_multireg_show'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'Clear'
        [varargout{1:nargout}] = Clear(varargin{:});
    case 'ShowAll'
        [varargout{1:nargout}] = ShowAll(varargin{:});
    case 'ShowBiasField'
        [varargout{1:nargout}] = ShowBiasField(varargin{:}); 
    case 'ShowIntensityModel'
        [varargout{1:nargout}] = ShowIntensityModel(varargin{:});         
    case 'ShowModel'
        [varargout{1:nargout}] = ShowModel(varargin{:});
    case 'ShowParameters'
        [varargout{1:nargout}] = ShowParameters(varargin{:});            
    case 'ShowSubjects'
        [varargout{1:nargout}] = ShowSubjects(varargin{:});        
    case 'Speak'
        [varargout{1:nargout}] = Speak(varargin{:});                
    otherwise
        help spm_multireg_show
        error('Unknown function %s. Type ''help spm_multireg_show'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% Clear()
function Clear(sett)
fn = {sett.show.figname_bf, sett.show.figname_int, sett.show.figname_model, ...
      sett.show.figname_subjects, sett.show.figname_parameters};
for i=1:numel(fn)
    f  = findobj('Type', 'Figure', 'Name', fn{i});
    if ~isempty(f), clf(f); drawnow; end    
end
end
%==========================================================================

%==========================================================================
% ShowAll()
function ShowAll(dat,mu,Objective,N,sett)
if sett.show.level >= 2
    ShowModel(mu,Objective,N,sett);
end
if sett.show.level >= 3
    ShowSubjects(dat,mu,sett);
    ShowParameters(dat,mu,sett);
    ShowBiasField(dat,sett);
    ShowIntensityModel(dat,sett);
end
end
%==========================================================================

%==========================================================================
% ShowBiasField()
function ShowBiasField(dat,sett)

% Parse function settings
axis_3d  = sett.show.axis_3d;
c        = sett.show.channel;
fig_name = sett.show.figname_bf;
fwhm     = sett.bf.fwhm;
mx_subj  = sett.show.mx_subjects;
reg      = sett.bf.reg;
updt_bf  = sett.do.updt_bf;

if ~updt_bf, return; end

[df,C] = spm_multireg_io('GetSize',dat(1).f);
nr = 3;
if df(3) == 1, ax = 3;
else,         ax = axis_3d;   
end
    
nd = min(numel(dat),mx_subj);
for n=1:nd
    df   = spm_multireg_io('GetSize',dat(n).f);
    fn   = spm_multireg_io('GetData',dat(n).f);   
    chan = spm_multireg_io('GetBiasFieldStruct',C,df,dat(n).Mat,reg,fwhm,[],dat(n).bf.T);
    bf   = spm_multireg_io('GetBiasField',chan,df);
    bf   = reshape(bf,[df C]);
    c    = min(c,C);
    
    % Show stuff
    ShowIm(fn(:,:,:,c),ax,nr,nd,n,fig_name,false);
    ShowIm(bf(:,:,:,c),ax,nr,nd,n + nd,fig_name,false);            
    ShowIm(bf(:,:,:,c).*fn(:,:,:,c),ax,nr,nd,n + 2*nd,fig_name,false);
end
drawnow
end
%==========================================================================

%==========================================================================
% ShowIntensityModel()
function ShowIntensityModel(dat,sett)

% Parse function settings
do_updt_int = sett.do.updt_int;
fig_name    = sett.show.figname_int;

if ~do_updt_int || ~isfield(dat(1),'mog'), return; end

n  = 1;
m0 = dat(n).mog.pr.m;
b0 = dat(n).mog.pr.b;
V0 = dat(n).mog.pr.V;
n0 = dat(n).mog.pr.n;

spm_gmm_lib('plot','gaussprior',{m0,b0,V0,n0},[],fig_name);
end
%==========================================================================

%==========================================================================
% ShowModel()
function ShowModel(mu,Objective,N,sett)

% Parse function settings
fig_name = sett.show.figname_model;

d  = size(mu);
mu = cat(4,mu,zeros(d(1:3),'single'));
mu = spm_multireg_util('softmax',mu,4);
% mu  = spm_multireg_util('softmaxmu',mu,4); % Gives K + 1 classes
nam = ['K=' num2str(size(mu,4)) ', N=' num2str(N) ' (softmaxed)'];
if size(mu,3) > 1
    ShowCat(mu,1,2,3,1,fig_name);
    ShowCat(mu,2,2,3,2,fig_name);
    ShowCat(mu,3,2,3,3,fig_name);
    title(nam);
    subplot(2,1,2); 
    plot(Objective,'.-');        
    title('Negative log-likelihood')
else
    ShowCat(mu,3,1,2,1,fig_name);
    title(nam);
    subplot(1,2,2); 
    plot(Objective,'.-');    
    title('Negative log-likelihood')
end
drawnow
end
%==========================================================================

%==========================================================================
% ShowParameters()
function ShowParameters(dat,mu0,sett)

% Parse function settings
B        = sett.registr.B;
c        = sett.show.channel;
fig_name = sett.show.figname_parameters;
fwhm     = sett.bf.fwhm;
Mmu      = sett.var.Mmu;
mx_subj  = sett.show.mx_subjects;
reg      = sett.bf.reg;
updt_bf  = sett.do.updt_bf;

fg = findobj('Type', 'Figure', 'Name', fig_name);
if isempty(fg)
    fg = figure('Name', fig_name, 'NumberTitle', 'off');
else
    clf(fg);
end
set(0, 'CurrentFigure', fg);   

nd = min(numel(dat),mx_subj);
nr = 3;
if ~isfield(dat,'mog')
    nr = nr - 1;
end
for n=1:nd              
    % Affine parameters
    q    = double(dat(n).q);
    M    = spm_dexpm(q,B);
    q    = spm_imatrix(M);            
    q    = q([1 2 6]);
    q(3) = 180/pi*q(3);
    q    = abs(q);
    subplot(nr,nd,n) 
    hold on
    for k=1:numel(q)
        bar(k, q(k));
    end
    box on
    hold off            
    set(gca,'xtick',[])  
    
    % Velocities (x)
    v = spm_multireg_io('GetData',dat(n).v);
    ShowIm(v(:,:,:,1),3,nr,nd,n + 1*nd,fig_name)
    
    if isfield(dat,'mog')
        % Intensity histogram
        
        % Here we get approximate class proportions from the (softmaxed K + 1)
        % tissue template
        [df,C] = spm_multireg_io('GetSize',dat(n).f);
        q      = double(dat(n).q);
        Mn     = dat(n).Mat;
        psi1   = spm_multireg_io('GetData',dat(n).psi);
        psi    = spm_multireg_util('Compose',psi1,spm_multireg_util('Affine',df,Mmu\spm_dexpm(q,B)*Mn));
        mu     = spm_multireg_util('Pull1',mu0,psi);               
        mu     = cat(4,mu,zeros(df(1:3),'single'));
        mu     = spm_multireg_util('softmax',mu,4);
        mu     = reshape(mu,[prod(df(1:3)) size(mu,4)]);
        mu     = sum(mu,1);
        mu     = mu./sum(mu);
        
        % Plot GMM fit        
        fn = spm_multireg_io('GetData',dat(n).f);        
        c  = min(c,size(fn,4));
        fn = fn(:,:,:,c);
        if updt_bf
            chan = spm_multireg_io('GetBiasFieldStruct',C,df,dat(n).Mat,reg,fwhm,[],dat(n).bf.T);
            bf   = spm_multireg_io('GetBiasField',chan,df);
            bf   = reshape(bf,[df C]);
            bf   = bf(:,:,:,c);
            fn   = bf.*fn;
        end
        ShowGMMFit(fn,mu,dat(n).mog,nr,nd,n + 2*nd,c);
    end
end
drawnow
end
%==========================================================================

%==========================================================================
% ShowSubjects()
function ShowSubjects(dat,mu0,sett)

% Parse function settings
axis_3d  = sett.show.axis_3d;
B        = sett.registr.B;
c        = sett.show.channel;
fig_name = sett.show.figname_subjects;
fwhm     = sett.bf.fwhm;
Mmu      = sett.var.Mmu;
mx_subj  = sett.show.mx_subjects;
reg      = sett.bf.reg;
updt_bf  = sett.do.updt_bf;

if isfield(dat,'mog'), nr = 3;
else,                  nr = 2;
end
if size(mu0,3) == 1, ax = 3;
else,                ax = axis_3d;   
end
    
nd = min(numel(dat),mx_subj);
for n=1:nd
    [df,C] = spm_multireg_io('GetSize',dat(n).f);
    q      = double(dat(n).q);
    Mn     = dat(n).Mat;
    psi1   = spm_multireg_io('GetData',dat(n).psi);
    psi    = spm_multireg_util('Compose',psi1,spm_multireg_util('Affine',df,Mmu\spm_dexpm(q,B)*Mn));
    mu     = spm_multireg_util('Pull1',mu0,psi);            
    
    fn = spm_multireg_io('GetClasses',dat(n),mu,sett);   
    fn = cat(4,fn,1 - sum(fn,4)); % Gives K + 1 classes
    
    mu = cat(4,mu,zeros(df(1:3),'single'));
    mu = spm_multireg_util('softmax',mu,4);
    
    % Show template, segmentation
    ShowCat(mu,ax,nr,nd,n,fig_name);
    ShowCat(fn,ax,nr,nd,n + nd,fig_name);        
    if isfield(dat,'mog')
        % and image (if using GMM)
        fn = spm_multireg_io('GetData',dat(n).f);
        c  = min(c,C);
        fn = fn(:,:,:,c);
        if updt_bf
            chan = spm_multireg_io('GetBiasFieldStruct',C,df,dat(n).Mat,reg,fwhm,[],dat(n).bf.T);
            bf   = spm_multireg_io('GetBiasField',chan,df);
            bf   = reshape(bf,[df C]);
            bf   = bf(:,:,:,c);
            fn   = bf.*fn;
        end
        ShowIm(fn,ax,nr,nd,n + 2*nd,fig_name);
    end
end
drawnow
end
%==========================================================================

%==========================================================================
% Speak()
function Speak(nam,varargin)
switch nam
    case 'Finished'
        t = varargin{1}; 
        
        fprintf('Algorithm finished in %.1f seconds.\n', t);
    case 'Init'
        nit = varargin{1}; 
        
        fprintf('Optimising parameters at largest zoom level (nit = %i)\n',nit)
    case 'Iter'
        nzm = varargin{1}; 
        
        fprintf('Optimising parameters at decreasing zoom levels (nzm = %i)\n',nzm)        
    case {'Register','Groupwise'}
        N = varargin{1};
        K = varargin{2};
        
        fprintf('------------------------------------\n')
        fprintf(' Begin %s (N = %i, K = %i)\n',nam,N,K)
        fprintf('------------------------------------\n\n')
    otherwise
        error('Unknown input!')
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% ShowCat()
function ShowCat(in,ax,nr,nc,np,fn)
f  = findobj('Type', 'Figure', 'Name', fn);
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

end
%==========================================================================

%==========================================================================
% ShowGMMFit()
function ShowGMMFit(f,PI,mog,nr,nd,n,c)
K      = size(mog.po.m,2);
colors = hsv(K);

subplot(nr,nd,n)
hold on

% ---------
% Histogram

% Input is a list of observations
[H, edges] = histcounts(f(f > 0 & isfinite(f)), 64, 'Normalization', 'pdf');
centres = (edges(1:end-1) + edges(2:end))/2;
bar(centres, H, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);

% -----------
% GMM Density
A     = bsxfun(@times, mog.po.V, reshape(mog.po.n, [1 1 K])); % Expected precision
xlims = [inf -inf];
for k=1:K
    MU   = mog.po.m(c,k);    
    sig2 = inv(A(c,c,k));
    
    x = linspace(MU - 3*sqrt(sig2), MU + 3*sqrt(sig2),100);
    y = PI(k)*normpdf(x, MU, sqrt(sig2));
    plot(x, y, 'Color', colors(k,:), 'LineWidth', 1)
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
function ShowIm(in,ax,nr,nc,np,fn,use_gray)
if nargin < 7, use_gray = true; end

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
imagesc(in); axis off image xy;
if use_gray
    colormap(gray);
end
end 
%==========================================================================