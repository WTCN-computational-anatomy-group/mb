function varargout = spm_multireg(varargin)
%__________________________________________________________________________
%
% Subject- or population-level image segmentation + normalisation.
%
% FORMAT [dat,mu] = spm_multireg('Groupwise',F,sett)
% FORMAT dat      = spm_multireg('Register',F,mu,sett)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Set paths
spm_multireg_util('SetPath');

% Set boundary conditions
spm_multireg_util('SetBoundCond');

% Do population- or subject-level
if nargin==1 && isa(varargin{1},'cell')
    [varargout{1:nargout}] = Groupwise(varargin{:});
elseif nargin==2
    [varargout{1:nargout}] = Register(varargin{:});
else
    error('Incorrect usage.');
end
end
%==========================================================================

%==========================================================================
function [dat,mu] = Groupwise(F,sett)
if nargin < 2, sett = struct; end

[~,is3d] = spm_multireg_io('GetData',F{1});
N        = numel(F); % Number of subjects

%------------------
% Get algorithm settings
%------------------

sett                 = spm_multireg_par('Settings',sett,is3d);
sett.model.groupwise = true;

dat       = spm_multireg_init('InitDat',F,sett);
clear F

%[~,M]     = spm_multireg_io('GetSize',dat(1).f);
M         = numel(dat(1).mog.mu) - 1;
[Mmu, d]  = spm_multireg_util('SpecifyMean',dat,sett.model.vx);
nz        = max(ceil(log2(min(d(d~=1)))-log2(8)),1);
sz        = spm_multireg_par('ZoomSettings',d,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var  = spm_multireg_io('CopyFields',sz(end), sett.var);
dat       = spm_multireg_init('InitDef',dat,sett);
%mu        = randn([sett.var.d M],'single')*10000;
mu        = zeros([sett.var.d M],'single');
%mu    = spm_multireg_util('ShrinkTemplate',mu,Mmu,sett);

%------------------
% Start algorithm
%------------------

spm_multireg_show('Speak','Groupwise',N,sett.model.K);

Objective = [];
E         = Inf;
prevt     = Inf;
te        = spm_multireg_energ('TemplateEnergy',mu,sett);
for iter=1:sett.nit.init

    %------------------
    % Updates template, affine and GMM parameters (at largest template resolution)    
    %------------------
    
    for subit=1:sett.nit.init_mu
        Eold     = E; tic;
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
        E        = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after rigid update
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);
        t        = toc;
        fprintf('%d 0 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
        prevt    = t;
        Objective = [Objective; E];
    end
   %if (Eold-E)/(numel(dat)*100^3)<1e-4, break; end

    Eold   = E; tic;
    dat    = spm_multireg_updt('UpdateSimpleAffines',dat,mu,sett);
    E      = sum(sum(cat(2,dat.E),2),1)+te;   % Cost function after mean update
    t      = toc;
    fprintf('%d 1 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
    prevt = t;

    spm_multireg_show('ShowSubjects',dat,mu,sett);
    Objective = [Objective; E];
    spm_multireg_show('ShowModel',mu,Objective,sett,N);        
end

%------------------
% Iteratively decrease the template resolution
%------------------

tic
for zm=numel(sz):-1:1
   %if zm~=numel(sz), [mu,dat] = spm_multireg_updt('UpdateSimpleMean',dat, mu, sett); end
    if zm~=numel(sz),
        for i=1:4
            [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
        end
    end
    te       = spm_multireg_energ('TemplateEnergy',mu,sett);

    niter = sett.nit.zm;
    niter = niter + (zm-1);
   %if zm==numel(sz), niter=niter+4; end
    E4    = Inf;
   %done  = false;
    for iter=1:niter

        % Might be an idea to run this multiple times
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
        E1       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after diffeo update
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);

        % Might be an idea to run this less often (currently slow)
        dat      = spm_multireg_updt('UpdateAffines',dat,mu,sett);
        E2       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after mean update

        % Might be an idea to run this multiple times
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);     % An extra mean iteration
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
        E3       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after rigid update
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);

        dat      = spm_multireg_energ('VelocityEnergy',dat,sett);
        dat      = spm_multireg_updt('UpdateVelocities',dat,mu,sett);
        E4old    = E4;
        E4       = sum(sum(cat(2,dat.E),2),1)+te; % Cost function after mean update
        dat      = spm_multireg_energ('VelocityEnergy',dat,sett);

%       if (E4old-E4)/E4 < 3.5e-4, done = true; end
       %if (iter==niter || done) && zm>1
        if (iter==niter) && zm>1
            oMmu      = sett.var.Mmu;
            sett.var  = spm_multireg_io('CopyFields',sz(zm-1), sett.var);
            [dat,mu]  = spm_multireg_util('ZoomVolumes',dat,mu,sett,oMmu);
        end

        dat      = spm_multireg_updt('UpdateWarps',dat,sett);

        fprintf('%d\t%g\t%g\t%g\t%g\n', iter, E1, E2, E3, E4);
        save(fullfile(sett.write.dir_res,'results_Groupwise.mat'),'dat','mu','sett')
        
        spm_multireg_show('ShowSubjects',dat,mu,sett);
        Objective = [Objective; E4];        
        spm_multireg_show('ShowModel',mu,Objective,sett,N);        

       %if done, break; end
    end
    fprintf('%g seconds\n\n', toc); tic;

end
[mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);

dat = spm_multireg_io('SaveImages',dat,mu,sett);
end
%==========================================================================

%==========================================================================
function dat = Register(F,mu,sett)
if nargin < 3, sett = struct; end

[~,is3d] = spm_multireg_io('GetData',F{1});
N        = numel(F); % Number of subjects

%------------------
% Get algorithm settings
%------------------

sett                 = spm_multireg_par('Settings',sett,is3d);
sett.model.groupwise = false;
sett.gen.threads     = 1;

%------------------
% Get template
%------------------

mu0          = spm_multireg_io('GetData',mu);
Mmu          = spm_multireg_io('GetMat',mu);
sett.model.K = size(mu0,4);

%------------------
% Get zoom (multi-scale) settings
%------------------

d        = [size(mu0,1) size(mu0,2) size(mu0,3)];
nz       = max(ceil(log2(min(d))-log2(8)),1);
sz       = spm_multireg_par('ZoomSettings',d,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_multireg_io('CopyFields',sz(end), sett.var);

%------------------
% Init dat (f, M, q, v, psi, E, mog, Mat)
%------------------

dat = spm_multireg_init('InitDat',F,sett);
clear F
dat = spm_multireg_init('InitDef',dat,sett);

%------------------
% Shrink template
%------------------

mu = spm_multireg_util('ShrinkTemplate',mu0,Mmu,sett);

%------------------
% Start algorithm
%------------------

spm_multireg_show('Speak','Register',N,sett.model.K);

Objective = [];
E         = Inf;
prevt     = Inf;
d         = spm_multireg_io('GetSize',dat(1).f); % dim of one input image
for iter=1:sett.nit.init
    
    %------------------
    % Updates affine and GMM parameters (at largest template resolution)    
    %------------------

    Eold      = E; tic;
    dat       = spm_multireg_updt('UpdateSimpleAffines',dat,mu,sett); % GMM parameters are updated in here (GetClasses())
    E         = sum(sum(cat(2,dat.E),2),1); % Cost function after affine update
    t         = toc;
    Objective = [Objective; E];
    
    fprintf('%d 1 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);        
    prevt = t;
    
    % Finished?
    if (Eold-E)/prod(d) < 1e-4, break; end        

    % Show stuff
    spm_multireg_show('ShowSubjects',dat,mu,sett);    
    spm_multireg_show('ShowModel',mu,Objective,sett,N);   
end

%------------------
% Iteratively decrease the template resolution
%------------------

for zm=numel(sz):-1:1
    
    %------------------    
    % Updates affine, velocity and GMM parameters
    %------------------
    
    mu    = spm_multireg_util('ShrinkTemplate',mu0,Mmu,sett);
    niter = sett.nit.zm;
    niter = niter + (zm-1);
    done  = false;
    dat   = spm_multireg_energ('VelocityEnergy',dat,sett); % Update velocity energy
    
    for iter=1:niter

        if iter==2 %iter<niter && zm~=numel(sz) && iter~=1
            Eold   = E; tic;
            dat    = spm_multireg_updt('UpdateAffines',dat,mu,sett);
            E      = sum(sum(cat(2,dat.E),2),1); % Cost function after diffeo update
            t      = toc;
            
            fprintf('%d 1 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
            prevt = t;
            
            Objective = [Objective; E];
        end
        
        Eold    = E; tic;
        dat     = spm_multireg_energ('VelocityEnergy',dat,sett);
        dat     = spm_multireg_updt('UpdateVelocities',dat,mu,sett);
        E       = sum(sum(cat(2,dat.E),2),1); % Cost function after rigid update

        if (iter == niter || done) && zm>1
            oMmu      = sett.var.Mmu;
            sett.var  = spm_multireg_io('CopyFields',sz(zm-1), sett.var);
            [dat,mu]  = spm_multireg_util('ZoomVolumes',dat,mu,sett,oMmu);
        end
        
        dat = spm_multireg_energ('VelocityEnergy',dat,sett);
        dat = spm_multireg_updt('UpdateWarps',dat,sett);
        t   = toc;
        
        fprintf('%d 2 \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
        prevt     = t;
        
        Objective = [Objective; E];

        % Save stuff
        save(fullfile(sett.write.dir_res,'results_Register.mat'),'dat','mu','sett')        

        % Show stuff
        spm_multireg_show('ShowSubjects',dat,mu,sett);        
        spm_multireg_show('ShowModel',mu,Objective,sett,N);  

        % Finished?
        if done, break; end
    end
end
end
%==========================================================================