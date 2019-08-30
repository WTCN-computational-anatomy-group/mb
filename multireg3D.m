function varargout = multireg3D(varargin)

pth=fileparts(which('spm'));
addpath(pth);
addpath(fullfile(pth,'toolbox','Longitudinal'));
addpath(fullfile(pth,'toolbox','Shoot'));

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'code')));

set_bound;

if nargin==1 && isa(varargin{1},'cell')
    [varargout{1:nargout}] = Groupwise3D(varargin{:});
elseif nargin==2
    [varargout{1:nargout}] = Register3D(varargin{:});
else
    error('Incorrect usage.');
end
end
%==========================================================================

%==========================================================================
function dat = Register3D(F,Mu)
% Rigid + nonlinear registration
sett           = Settings;
sett.groupwise = false;
sett.threads   = 1;

mu0  = GetData(Mu);
Mmu  = GetMat(Mu);
d    = [size(mu0,1) size(mu0,2) size(mu0,3)];
nz   = max(ceil(log2(min(d))-log2(8)),1);
sz   = ZoomSettings(d,Mmu,sett.v_settings,sett.mu_settings,nz);
sett = CopyFields(sz(end), sett);
sett.K   = size(mu0,4);
dat      = InitStruct1(F,sett);
dat(1).f = GetData(dat(1).f);
dat      = InitStruct2(dat,sett);

Objective = [];

zm  = numel(sz);
mu    = ShrinkTemplate(mu0,Mmu,sett);
E     = Inf;
prevt = Inf;
d     = GetSize(dat.f);
for iter=1:10
    Eold   = E; tic;
    dat    = UpdateSimpleAffines(dat,mu,sett);
    E      = sum(sum(cat(2,dat.E),2),1); % Cost function after diffeo update
    t      = toc;
    fprintf('%d 1 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
    prevt  = t;
    if (Eold-E)/prod(d)<1e-4, break; end

    ShowSubjects(dat,mu,sett);
    Objective = [Objective; E];
    ShowModel(mu,Objective,sett);   
end

for zm=numel(sz):-1:1
    mu    = ShrinkTemplate(mu0,Mmu,sett);
    niter = sett.nits;
    niter = niter + (zm-1);
    done  = false;
    dat   = VelocityEnergy(dat,sett);
    for iter=1:niter

        if iter==2 %iter<niter && zm~=numel(sz) && iter~=1
            Eold   = E; tic;
            dat    = UpdateAffines(dat,mu,sett);
            E      = sum(sum(cat(2,dat.E),2),1); % Cost function after diffeo update
            t      = toc;
            fprintf('%d 1 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
            prevt = t;
            Objective = [Objective; E];
        end
        Eold    = E; tic;
        dat     = VelocityEnergy(dat,sett);
        dat     = UpdateVelocities(dat,mu,sett);
        E       = sum(sum(cat(2,dat.E),2),1); % Cost function after rigid update

        if (iter==niter || done) && zm>1
            oMmu      = sett.Mmu;
            sett      = CopyFields(sz(zm-1), sett);
            [dat,mu]  = ZoomVolumes(dat,mu,sett,oMmu);
        end
        dat     = VelocityEnergy(dat,sett);
        dat     = UpdateWarps(dat,sett);
        t       = toc;
        fprintf('%d 2 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
        prevt   = t;

        save RegResults.mat dat mu sett

        ShowSubjects(dat,mu,sett);
        Objective = [Objective; E];
        ShowModel(mu,Objective,sett);  

        if done, break; end
    end
end
end
%==========================================================================

%==========================================================================
function [dat,mu] = Groupwise3D(F)

is3d = size(F{1},3) > 1;

sett           = Settings(is3d);
sett.K         = sett.model.K;
sett.groupwise = true;

N         = numel(F);
dat       = InitStruct1(F,sett);
%[~,M]     = GetSize(dat(1).f);
M         = numel(dat(1).mog.mu)-1;
[Mmu, d]  = SpecifyMean(dat,sett.model.vx);
nz        = max(ceil(log2(min(d(d~=1)))-log2(8)),1);
sz        = ZoomSettings(d,Mmu,sett.v_settings,sett.mu_settings,nz);
sett      = CopyFields(sz(end), sett);
dat       = InitStruct2(dat,sett);
%mu        = randn([sett.d M],'single')*10000;
mu        = zeros([sett.d M],'single');
Objective = [];

% zm    = numel(sz);
%mu    = ShrinkTemplate(mu,Mmu,sett);
E     = Inf;
prevt = Inf;
te    = TemplateEnergy(mu,sett);

for iter=1:10

    for subit=1:3
        Eold     = E; tic;
        [mu,dat] = UpdateMean(dat, mu, sett);
        E        = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after rigid update
        te       = TemplateEnergy(mu,sett);
        t        = toc;
        fprintf('%d 0 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
        prevt    = t;
        Objective = [Objective; E];
    end
   %if (Eold-E)/(numel(dat)*100^3)<1e-4, break; end

    Eold   = E; tic;
    dat    = UpdateSimpleAffines(dat,mu,sett);
    E      = sum(sum(cat(2,dat.E),2),1)+te;   % Cost function after mean update
    t      = toc;
    fprintf('%d 1 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
    prevt = t;

    ShowSubjects(dat,mu,sett);
    Objective = [Objective; E];
    ShowModel(mu,Objective,sett,N);        
end

tic
for zm=numel(sz):-1:1
   %if zm~=numel(sz), [mu,dat] = UpdateSimpleMean(dat, mu, sett); end
    if zm~=numel(sz),
        for i=1:4
            [mu,dat] = UpdateMean(dat, mu, sett);
        end
    end
    te       = TemplateEnergy(mu,sett);

    niter = sett.nits;
    niter = niter + (zm-1);
   %if zm==numel(sz), niter=niter+4; end
    E4    = Inf;
   %done  = false;
    for iter=1:niter

        % Might be an idea to run this multiple times
        [mu,dat] = UpdateMean(dat, mu, sett);
        E1       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after diffeo update
        te       = TemplateEnergy(mu,sett);

        % Might be an idea to run this less often (currently slow)
        dat      = UpdateAffines(dat,mu,sett);
        E2       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after mean update

        % Might be an idea to run this multiple times
        [mu,dat] = UpdateMean(dat, mu, sett);     % An extra mean iteration
        te       = TemplateEnergy(mu,sett);
        [mu,dat] = UpdateMean(dat, mu, sett);
        E3       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after rigid update
        te       = TemplateEnergy(mu,sett);

        dat      = VelocityEnergy(dat,sett);
        dat      = UpdateVelocities(dat,mu,sett);
        E4old    = E4;
        E4       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after mean update
        dat      = VelocityEnergy(dat,sett);

%       if (E4old-E4)/E4 < 3.5e-4, done = true; end
       %if (iter==niter || done) && zm>1
        if (iter==niter) && zm>1
            oMmu      = sett.Mmu;
            sett      = CopyFields(sz(zm-1), sett);
            [dat,mu]  = ZoomVolumes(dat,mu,sett,oMmu);
        end

        dat      = UpdateWarps(dat,sett);

        fprintf('%d\t%g\t%g\t%g\t%g\n', iter, E1, E2, E3, E4);
        save MultiRegResults.mat dat mu sett
        
        ShowSubjects(dat,mu,sett);
        Objective = [Objective; E4];        
        ShowModel(mu,Objective,sett,N);        

       %if done, break; end
    end
    fprintf('%g seconds\n\n', toc); tic;

end
[mu,dat] = UpdateMean(dat, mu, sett);

dat = SaveImages(dat,mu,sett);
end
%==========================================================================