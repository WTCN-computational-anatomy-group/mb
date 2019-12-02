function [dat,mu,sett] = spm_mb_fit(F,sett)
%__________________________________________________________________________
%
% Multi-Brain - Groupwise normalisation and segmentation of images
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin < 2, sett = struct; end

% Set boundary conditions and path
spm_mb_io('SetBoundCond');
spm_mb_io('SetPath');

t0 = tic;

%------------------
% Get algorithm settings
%------------------

sett = spm_mb_param('Settings',sett);

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

% Clear figures
spm_mb_show('Clear',sett);

%------------------
% Init dat
%------------------

dat = spm_mb_io('InitDat',F,sett); 
F   = [];
N   = numel(dat); % Number of subjects

% Get number of template classes (if not using GMM)
if ~do_gmm, [~,K] = spm_mb_io('GetSize',dat(1).f); end

%------------------
% Get template size and orientation
%------------------

[Mmu, d] = spm_mb_shape('SpecifyMean',dat,vx);

%------------------
% Get zoom (multi-scale) settings
%------------------

nz       = max(ceil(log2(min(d(d~=1))) - log2(8)),1);
if ~do_zoom, nz = 1; end
sz       = spm_mb_param('ZoomSettings',d,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_mb_io('CopyFields',sz(end), sett.var);

%------------------
% Init deformation, bias field and GMM
%------------------

dat = spm_mb_shape('Init',dat,sett);
dat = spm_mb_appearance('Init',dat,K,sett);

%------------------
% Start algorithm (Groupwise)
%------------------

spm_mb_show('Speak','Groupwise',N,K);

Objective = [];
E         = Inf;
prevt     = Inf;
mu        = zeros([sett.var.d K],'single'); % Initial template (uniform)

if do_updt_aff
    spm_mb_show('Speak','Init',sett.nit.init);
    for iter=1:nit_init

        %------------------
        % Updates template, affine and GMM parameters (at largest template resolution)    
        %------------------

        for subit=1:nit_init_mu
            % Update template, bias field and intensity model
            Eold     = E; tic;                        
            [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
            te       = spm_mb_shape('TemplateEnergy',mu,sett);
            dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
            E        = sum(sum(cat(2,dat.E),2),1) + te;
            t        = toc;
                               
            % Print stuff
            fprintf('it=%i mu \t%g\t%g\t%g\n', iter, E, t, (Eold - E)/prevt);
            prevt     = t;
            Objective = [Objective; E];
            
            % Show stuff
            spm_mb_show('All',dat,mu,Objective,N,sett);
        end

        % Update affine
        Eold = E; tic;
        dat  = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett);
        dat  = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E    = sum(sum(cat(2,dat.E),2),1) + te;
        t    = toc;
        
        % Print stuff
        fprintf('it=%i q  \t%g\t%g\t%g\n', iter, E, t, (Eold - E)/prevt);
        prevt = t;
        Objective = [Objective; E];

        % Save stuff
        save(fullfile(dir_res,'results_Groupwise.mat'),'dat','mu','sett')
        
        % Show stuff
        spm_mb_show('All',dat,mu,Objective,N,sett);
    end
end

%------------------
% Iteratively decrease the template resolution (Groupwise)
%------------------

spm_mb_show('Speak','Iter',numel(sz)); tic;
for zm=numel(sz):-1:1 % loop over zoom levels
    
    E0 = 0;
    if zm ~= numel(sz) || zm == 1
        % Runs only at finest resolution
        for i=1:nit_init_mu
            % Update template, bias field and intensity model                        
            [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
            dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
            
            % Show stuff
            spm_mb_show('All',dat,mu,Objective,N,sett);
        end
        te = spm_mb_shape('TemplateEnergy',mu,sett);
        E0 = sum(sum(cat(2,dat.E),2),1) + te;
    end    
        
    E4 = Inf;
    for iter=1:nit_zm

        % Update template, bias field and intensity model
        % Might be an idea to run this multiple times                
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
        te       = spm_mb_shape('TemplateEnergy',mu,sett);
        dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E1       = sum(sum(cat(2,dat.E),2),1) + te;        
                           
        % Update affine
        % (Might be an idea to run this less often - currently slow)
        dat      = spm_mb_shape('UpdateAffines',dat,mu,sett);
        dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E2       = sum(sum(cat(2,dat.E),2),1) + te;

        % Update template, bias field and intensity model
        % (Might be an idea to run this multiple times)                
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett); % An extra mean iteration
        te       = spm_mb_shape('TemplateEnergy',mu,sett);
        dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
                
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
        te       = spm_mb_shape('TemplateEnergy',mu,sett);
        dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E3       = sum(sum(cat(2,dat.E),2),1) + te;
            
        % Update velocities
        dat      = spm_mb_shape('VelocityEnergy',dat,sett);
        dat      = spm_mb_shape('UpdateVelocities',dat,mu,sett);
        dat      = spm_mb_shape('VelocityEnergy',dat,sett);
        dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E4old    = E4;
        E4       = sum(sum(cat(2,dat.E),2),1) + te;       

        if (iter==nit_zm) && zm>1
            oMmu     = sett.var.Mmu;
            sett.var = spm_mb_io('CopyFields',sz(zm-1), sett.var);
            [dat,mu] = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
        end

        % Update deformations
        dat = spm_mb_shape('UpdateWarps',dat,sett);  
        
        % Print stuff
        fprintf('zm=%i it=%i\t%g\t%g\t%g\t%g\t%g\n', zm, iter, E0, E1, E2, E3, E4);        
        Objective = [Objective; E4];
                
        % Save stuff
        save(fullfile(dir_res,'results_Groupwise.mat'),'dat','mu','sett')
        
        % Show stuff
        spm_mb_show('All',dat,mu,Objective,N,sett);
    end
    
    fprintf('%g seconds\n\n', toc); tic;
end

% Final mean update
[mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);

% Save template
dat = spm_mb_io('SaveTemplate',dat,mu,sett);

% Print total runtime
spm_mb_show('Speak','Finished',toc(t0));

if sett.show.level >= 1
    % Show stuff
    spm_mb_show('Model',mu,Objective,N,sett);
    spm_mb_show('Subjects',dat,mu,sett);
    spm_mb_show('Parameters',dat,mu,sett);
    spm_mb_show('BiasField',dat,sett);
    spm_mb_show('IntensityPrior',dat,sett);
end
end
%==========================================================================