function varargout = spm_multireg_show(varargin)
%__________________________________________________________________________
%
% Visualisation functions for spm_multireg.
%
% FORMAT spm_multireg_show('ShowModel',mu,Objective,sett,N)
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
% ShowModel()
function ShowModel(mu,Objective,sett,N)
mu = spm_multireg_util('softmaxmu',mu,4); % Gives K + 1 classes    
if size(mu,3) > 1
    ShowCat(mu,1,2,3,1,sett.show.figname_model);
    ShowCat(mu,2,2,3,2,sett.show.figname_model);
    ShowCat(mu,3,2,3,3,sett.show.figname_model);
    subplot(2,1,2); plot(Objective,'.-');
    title(['K=' num2str(size(mu,4)) ', N=' num2str(N)]);
else
    ShowCat(mu,3,1,2,1,sett.show.figname_model);
    title(['K=' num2str(size(mu,4)) ', N=' num2str(N)]);
    subplot(1,2,2); plot(Objective,'.-');    
    title('Negative log-likelihood')
end
drawnow
end
%==========================================================================

%==========================================================================
% ShowParameters()
function ShowParameters(dat,mu0,sett)
fg = findobj('Type', 'Figure', 'Name', sett.show.figname_parameters);
if isempty(fg)
    fg = figure('Name', sett.show.figname_parameters, 'NumberTitle', 'off');
else
    clf(fg);
end
set(0, 'CurrentFigure', fg);  

nd = min(numel(dat),sett.show.mx_subjects);
nr = 3;
if ~isfield(dat,'mog')
    nr = nr - 1;
end
for n=1:nd              
    % Affine parameters
    q    = double(dat(n).q);
    M    = spm_dexpm(q,sett.registr.B);
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
    ShowIm(v(:,:,:,1),3,nr,nd,n + 1*nd,sett.show.figname_parameters)
    
    if isfield(dat,'mog')
        % Intensity histogram
        
        % Here we get approximate class proportions from the (softmaxed K + 1)
        % tissue template
        d    = spm_multireg_io('GetSize',dat(n).f);
        q    = double(dat(n).q);
        Mn   = dat(n).Mat;
        psi1 = spm_multireg_io('GetData',dat(n).psi);
        psi  = spm_multireg_util('Compose',psi1,spm_multireg_util('Affine',d,sett.var.Mmu\spm_dexpm(q,sett.registr.B)*Mn));
        mu   = spm_multireg_util('Pull1',mu0,psi);        
        mu   = spm_multireg_util('softmaxmu',mu,4);
        mu   = reshape(mu,[prod(d(1:3)) size(mu,4)]);
        mu   = sum(mu,1);
        mu   = mu./sum(mu);
        
        % Plot GMM fit
        f = spm_multireg_io('GetData',dat(n).f);
        ShowGMMFit(f, mu, dat(n).mog, nr, nd, n + 2*nd);
    end
end
drawnow
end
%==========================================================================

%==========================================================================
% ShowSubjects()
function ShowSubjects(dat,mu0,sett)
fg = findobj('Type', 'Figure', 'Name', sett.show.figname_subjects);
if isempty(fg)
    fg = figure('Name', sett.show.figname_subjects, 'NumberTitle', 'off');
else
    clf(fg);
end
set(0, 'CurrentFigure', fg);   

if isfield(dat,'mog'), nr = 3;
else,                  nr = 2;
end
if size(mu0,3) == 1, ax = 3;
else,                ax = sett.show.axis_3d;   
end
    
nd = min(numel(dat),sett.show.mx_subjects);
for n=1:nd
    d    = spm_multireg_io('GetSize',dat(n).f);
    q    = double(dat(n).q);
    Mn   = dat(n).Mat;
    psi1 = spm_multireg_io('GetData',dat(n).psi);
    psi  = spm_multireg_util('Compose',psi1,spm_multireg_util('Affine',d,sett.var.Mmu\spm_dexpm(q,sett.registr.B)*Mn));
    mu   = spm_multireg_util('Pull1',mu0,psi);            
    
    fn = spm_multireg_io('GetClasses',dat(n),mu,sett);   
    fn = cat(4,fn,1 - sum(fn,4)); % Gives K + 1 classes
    
    mu = spm_multireg_util('softmaxmu',mu,4); % Gives K + 1 classes
    
    if isfield(dat,'mog')
        f = spm_multireg_io('GetData',dat(n).f);
    end
    
    % Show template, segmentation
    ShowCat(mu,ax,nr,nd,n,sett.show.figname_subjects);
    ShowCat(fn,ax,nr,nd,n + nd,sett.show.figname_subjects);        
    if isfield(dat,'mog')
        % and image (if using GMM)
        ShowIm(f,ax,nr,nd,n + 2*nd,sett.show.figname_subjects);
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
mid = ceil(dm(1:3).*0.55);

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
function ShowGMMFit(f,PI,mog,nr,nd,n)
K      = numel(mog.mu);
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
xlims = [inf -inf];
for k=1:K
    x = linspace(mog.mu(k) - 3*sqrt(mog.sig2(k)), mog.mu(k) + 3*sqrt(mog.sig2(k)),100);
    y = PI(k)*normpdf(x, mog.mu(k), sqrt(mog.sig2(k)));
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
function ShowIm(in,ax,nr,nc,np,fn)
f  = findobj('Type', 'Figure', 'Name', fn);
if isempty(f)
    f = figure('Name', fn, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   

subplot(nr,nc,np);

dm  = size(in);
dm  = [dm 1];
mid = ceil(dm(1:3).*0.55);

if ax == 1
    in = permute(in(mid(1),:,:,:),[3 2 1 4]);
elseif ax == 2
    in = permute(in(:,mid(2),:,:),[1 3 2 4]);
else
    in = in(:,:,mid(3),:);
end

in = squeeze(in(:,:,:,:));
imagesc(in); axis off image xy;   
colormap(gray);
end 
%==========================================================================