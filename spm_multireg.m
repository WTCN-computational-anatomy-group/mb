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
if nargin == 0
    help spm_multireg
    error('Not enough argument. Type ''help spm_multireg'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'Groupwise'
        [varargout{1:nargout}] = Groupwise(varargin{:});   
    case 'Register'
        [varargout{1:nargout}] = Register(varargin{:});           
    otherwise
        help spm_multireg
        error('Unknown function %s. Type ''help spm_multireg'' for help.', id)
end
end
%==========================================================================

%==========================================================================
function [dat,mu] = Groupwise(F,sett)
if nargin < 2, sett = struct; end

t0 = tic;

%------------------
% Get algorithm settings
%------------------

sett = spm_multireg_par('Settings',sett);

dir_res     = sett.write.dir_res;
do_gmm      = sett.do.gmm;
do_updt_aff = sett.do.updt_aff;
do_zoom     = sett.do.zoom;
K           = sett.model.K; 
nit_init    = sett.nit.init;
nit_init_mu = sett.nit.init_mu;
nit_zm      = sett.nit.zm;
vx          = sett.model.vx;

sett.model.groupwise = true;

%------------------
% Init dat
%------------------

dat = spm_multireg_init('InitDat',F,sett); clear F
N   = numel(dat); % Number of subjects

% Get number of template classes (if not using GMM)
if ~do_gmm, [~,K] = spm_multireg_io('GetSize',dat(1).f); end

%------------------
% Get template size and orientation
%------------------

[Mmu, d] = spm_multireg_util('SpecifyMean',dat,vx);

%------------------
% Get zoom (multi-scale) settings
%------------------

nz       = max(ceil(log2(min(d(d~=1))) - log2(8)),1);
if ~do_zoom, nz = 1; end
sz       = spm_multireg_par('ZoomSettings',d,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_multireg_io('CopyFields',sz(end), sett.var);

%------------------
% Init deformation, bias field and GMM
%------------------

dat = spm_multireg_init('InitDef',dat,sett);
dat = spm_multireg_init('InitBiasField',dat,sett);
dat = spm_multireg_init('InitGMM',dat,K,sett);

%------------------
% Start algorithm (Groupwise)
%------------------

spm_multireg_show('Speak','Groupwise',N,K);

Objective = [];
E         = Inf;
prevt     = Inf;
mu        = zeros([sett.var.d K],'single'); % Initial template (uniform)

if do_updt_aff
    spm_multireg_show('Speak','Init',sett.nit.init);
    for iter=1:nit_init

        %------------------
        % Updates template, affine and GMM parameters (at largest template resolution)    
        %------------------

        for subit=1:nit_init_mu
            % Update template, bias field and intensity model
            Eold     = E; tic;
            [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
            dat      = spm_multireg_updt('UpdateBiasField',dat,mu,sett);
            dat      = spm_multireg_updt('UpdateIntensity',dat, sett);
            te       = spm_multireg_energ('TemplateEnergy',mu,sett);
            E        = sum(sum(cat(2,dat.E),2),1) + te;
            t        = toc;
                               
            % Print stuff
            fprintf('it=%i mu \t%g\t%g\t%g\n', iter, E, t, (Eold - E)/prevt);
            prevt     = t;
            Objective = [Objective; E];
            
            % Show stuff
            spm_multireg_show('ShowAll',dat,mu,Objective,N,sett);
        end

        % Update affine
        Eold = E; tic;
        dat  = spm_multireg_updt('UpdateSimpleAffines',dat,mu,sett);
        E    = sum(sum(cat(2,dat.E),2),1) + te;
        t    = toc;
        
        % Print stuff
        fprintf('it=%i q  \t%g\t%g\t%g\n', iter, E, t, (Eold - E)/prevt);
        prevt = t;
        Objective = [Objective; E];

        % Show stuff
        spm_multireg_show('ShowAll',dat,mu,Objective,N,sett);
    end
end

%------------------
% Iteratively decrease the template resolution (Groupwise)
%------------------

spm_multireg_show('Speak','Iter',numel(sz)); tic;
for zm=numel(sz):-1:1 % loop over zoom levels
    
    E0 = 0;
    if zm ~= numel(sz) || zm == 1
        % Runs only at finest resolution
        for i=1:nit_init_mu
            % Update template, bias field and intensity model
            [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
            dat      = spm_multireg_updt('UpdateBiasField',dat,mu,sett);
            dat      = spm_multireg_updt('UpdateIntensity',dat, sett);
            
            % Show stuff
            spm_multireg_show('ShowAll',dat,mu,Objective,N,sett);
        end
        te = spm_multireg_energ('TemplateEnergy',mu,sett);
        E0 = sum(sum(cat(2,dat.E),2),1) + te;
    end    
        
    E4 = Inf;
    for iter=1:nit_zm

        % Update template, bias field and intensity model
        % Might be an idea to run this multiple times
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);
        dat      = spm_multireg_updt('UpdateBiasField',dat,mu,sett);
        dat      = spm_multireg_updt('UpdateIntensity',dat, sett);
        E1       = sum(sum(cat(2,dat.E),2),1) + te;        
                           
        % Update affine
        % (Might be an idea to run this less often - currently slow)
        dat      = spm_multireg_updt('UpdateAffines',dat,mu,sett);
        E2       = sum(sum(cat(2,dat.E),2),1) + te;

        % Update template, bias field and intensity model
        % (Might be an idea to run this multiple times)
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett); % An extra mean iteration
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);
        dat      = spm_multireg_updt('UpdateBiasField',dat,mu,sett);
        dat      = spm_multireg_updt('UpdateIntensity',dat, sett);        
                    
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);
        dat      = spm_multireg_updt('UpdateBiasField',dat,mu,sett);
        dat      = spm_multireg_updt('UpdateIntensity',dat, sett);
        E3       = sum(sum(cat(2,dat.E),2),1) + te;
            
        % Update velocities
        dat      = spm_multireg_energ('VelocityEnergy',dat,sett);
        dat      = spm_multireg_updt('UpdateVelocities',dat,mu,sett);
        dat      = spm_multireg_energ('VelocityEnergy',dat,sett);
        E4old    = E4;
        E4       = sum(sum(cat(2,dat.E),2),1) + te;       

        if (iter==nit_zm) && zm>1
            oMmu     = sett.var.Mmu;
            sett.var = spm_multireg_io('CopyFields',sz(zm-1), sett.var);
            [dat,mu] = spm_multireg_util('ZoomVolumes',dat,mu,sett,oMmu);
        end

        % Update deformations
        dat = spm_multireg_updt('UpdateWarps',dat,sett);  
        
        % Print stuff
        fprintf('zm=%i it=%i\t%g\t%g\t%g\t%g\t%g\n', zm, iter, E0, E1, E2, E3, E4);        
        Objective = [Objective; E4];
                
        % Save stuff
        save(fullfile(dir_res,'results_Groupwise.mat'),'dat','mu','sett')
        
        % Show stuff
        spm_multireg_show('ShowAll',dat,mu,Objective,N,sett);
    end
    
    fprintf('%g seconds\n\n', toc); tic;
end

% Final mean update
[mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);

% Save stuff
dat = spm_multireg_io('SaveImages',dat,mu,sett);

% Print total runtime
spm_multireg_show('Speak','Finished',toc(t0));

end
%==========================================================================

%==========================================================================
function dat = Register(F,mu,sett)
if nargin < 3, sett = struct; end

t0 = tic;

%------------------
% Get algorithm settings
%------------------

sett                 = spm_multireg_par('Settings',sett);
sett.model.groupwise = false;

%------------------
% Get template
%------------------

mu0 = spm_multireg_io('GetData',mu);
Mmu = spm_multireg_io('GetMat',mu);
K   = size(mu0,4);

%------------------
% Get zoom (multi-scale) settings
%------------------

d        = [size(mu0,1) size(mu0,2) size(mu0,3)];
nz       = max(ceil(log2(min(d(d ~= 1))) - log2(8)),1);
if ~sett.do.zoom, nz = 1; end
sz       = spm_multireg_par('ZoomSettings',d,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_multireg_io('CopyFields',sz(end), sett.var);

%------------------
% Init dat, deformation, bias field and GMM
%------------------

dat = spm_multireg_init('InitDat',F,sett); clear F
dat = spm_multireg_init('InitDef',dat,sett);
dat = spm_multireg_init('InitBiasField',dat,sett);
dat = spm_multireg_init('InitGMM',dat,K,sett);
N   = numel(dat); % Number of subjects

%------------------
% Shrink template
%------------------

mu = spm_multireg_util('ShrinkTemplate',mu0,Mmu,sett);

%------------------
% Start algorithm (Register)
%------------------

spm_multireg_show('Speak','Register',N,K);

Objective = [];
E         = Inf;
prevt     = Inf;
d         = spm_multireg_io('GetSize',dat(1).f); % dim of one input image

if sett.do.updt_aff
    spm_multireg_show('Speak','Init',sett.nit.init);
    for iter=1:sett.nit.init

        %------------------
        % Updates affine and GMM parameters (at largest template resolution)    
        %------------------

        Eold      = E; tic;
        dat       = spm_multireg_updt('UpdateSimpleAffines',dat,mu,sett); % GMM parameters are updated in here
        E         = sum(sum(cat(2,dat.E),2),1);                           % Cost function after affine update
        t         = toc;
        Objective = [Objective; E];
        
        % Print stuff
        fprintf('it=%i q \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);        
        prevt = t;

%         % Update bias field
%         Eold = E; tic;                
%         dat  = spm_multireg_updt('UpdateBiasField',dat,mu,sett);                        
%         E    = sum(sum(cat(2,dat.E),2),1);
%         t    = toc;
        
        % Print stuff
        fprintf('it=%i bf \t%g\t%g\t%g\n', iter, E, t, (Eold-E)/prevt);
        prevt = t;
        Objective = [Objective; E];
        
        % Finished?
        if (Eold-E)/prod(d) < 1e-4, break; end        

        % Show stuff
        spm_multireg_show('ShowAll',dat,mu,Objective,N,sett);
    end
end

%------------------
% Iteratively decrease the template resolution (Register)
%------------------

spm_multireg_show('Speak','Iter',numel(sz)); tic;
for zm=numel(sz):-1:1 % loop over zoom levels
    
    %------------------    
    % Updates affine, velocity and GMM parameters
    %------------------
        
    % Resize template
    mu = spm_multireg_util('ShrinkTemplate',mu0,Mmu,sett);
     
    % Update velocity energy
    dat = spm_multireg_energ('VelocityEnergy',dat,sett); 
    
    niter = sett.nit.zm;
    niter = niter + (zm-1);  
    for iter=1:niter

        if iter == 2 %iter<niter && zm~=numel(sz) && iter~=1
            % Update affine parameters            
            Eold   = E; tic;
            dat    = spm_multireg_updt('UpdateAffines',dat,mu,sett);
            E      = sum(sum(cat(2,dat.E),2),1); % Cost function after affine update
            t      = toc;
            
            % Print stuff
            fprintf('zm=%i it=%i q \t%g\t%g\t%g\n', zm, iter, E, t, (Eold-E)/prevt);
            prevt = t;
            
            Objective = [Objective; E];
        end
        
        Eold = E; tic;
        
        % Update velocities
        dat  = spm_multireg_energ('VelocityEnergy',dat,sett);
        dat  = spm_multireg_updt('UpdateVelocities',dat,mu,sett);
        E    = sum(sum(cat(2,dat.E),2),1); % Cost function after rigid update
        if (iter == niter || done) && zm > 1
            oMmu      = sett.var.Mmu;
            sett.var  = spm_multireg_io('CopyFields',sz(zm-1), sett.var);
            [dat,mu]  = spm_multireg_util('ZoomVolumes',dat,mu,sett,oMmu);
        end        
        dat = spm_multireg_energ('VelocityEnergy',dat,sett);
        
        % Update deformations
        dat = spm_multireg_updt('UpdateWarps',dat,sett);
        
        % Update bias field        
        dat  = spm_multireg_updt('UpdateBiasField',dat,mu,sett);                        
        E    = sum(sum(cat(2,dat.E),2),1);          
        t    = toc;
        
        % Print stuff
        fprintf('zm=%i it=%i v \t%g\t%g\t%g\n', zm, iter, E, t, (Eold-E)/prevt);
        prevt = t;
        
        Objective = [Objective; E];

        % Save stuff
        save(fullfile(sett.write.dir_res,'results_Register.mat'),'dat','mu','sett')        

        % Show stuff
        spm_multireg_show('ShowAll',dat,mu,Objective,N,sett);
    end
    
    fprintf('%g seconds\n\n', toc); tic;
end

% Print total runtime
spm_multireg_show('Speak','Finished',toc(t0));

end
%==========================================================================