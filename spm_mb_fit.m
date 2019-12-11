function [dat,model,sett] = spm_mb_fit(data,varargin)
%__________________________________________________________________________
%
% Multi-Brain - Groupwise normalisation and segmentation of images
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Parse input
p              = inputParser;
p.FunctionName = 'spm_mb_fit';
p.addParameter('model', struct(), @isstruct);
p.addParameter('sett',  struct(), @isstruct);
p.parse(varargin{:});
model = p.Results.model;
sett  = p.Results.sett;

% Set boundary conditions and path
spm_mb_io('SetBoundCond');
spm_mb_io('SetPath');

t0 = tic;

%------------------
% Get algorithm settings
%------------------

sett        = spm_mb_param('Settings',sett);
dir_res     = sett.write.dir_res;
do_gmm      = sett.do.gmm;
do_updt_aff = sett.do.updt_aff;
do_zoom     = sett.do.zoom;
K           = sett.model.K; 
nit_init    = sett.nit.init;
nit_init_mu = sett.nit.init_mu;
nit_zm0     = sett.nit.zm;
show_level  = sett.show.level;
vx          = sett.model.vx;

spm_mb_show('Clear',sett); % Clear figures

%------------------
% Decide what to learn
%------------------

N = numel(data); % Number of subjects

[sett,template_given] = spm_mb_param('SetFit',model,sett,N);

do_updt_int      = sett.do.updt_int;
do_updt_template = sett.do.updt_template;

%------------------
% Init dat
%------------------

dat  = spm_mb_io('InitDat',data,sett); 
data = [];

% Get number of template classes (if not using GMM)
if ~do_gmm, [~,K] = spm_mb_io('GetSize',dat(1).f); end
if template_given
    [~,K]        = spm_mb_io('GetSize',model.shape.template);
    sett.model.K = K;
end

%------------------
% Get template size and orientation
%------------------

if template_given    
    d         = spm_mb_io('GetSize',model.shape.template);
    [mu0,Mmu] = spm_mb_io('GetData',model.shape.template);   
else
    [Mmu,d] = spm_mb_shape('SpecifyMean',dat,vx);
end

%------------------
% Get zoom (multi-scale) settings
%------------------

nz       = max(ceil(log2(min(d(d~=1))) - log2(8)),1);
if ~do_zoom, nz = 1; end
sz       = spm_mb_param('ZoomSettings',d,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_mb_io('CopyFields',sz(end), sett.var);

%------------------
% Init shape and apperance model parameters
%------------------

dat = spm_mb_shape('Init',dat,sett);
dat = spm_mb_appearance('Init',dat,model,K,sett);

spm_mb_show('Speak','Groupwise',N,K);

%------------------
% Init template
%------------------

if template_given    
    % Shrink given template
    mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
else
    % Initial template
    [dat,mu,sett] = spm_mb_shape('InitMu',dat,K,sett);
end

spm_mb_show('All',dat,mu,[],N,sett);

%------------------
% Start algorithm
%------------------

Objective = [];
E         = Inf;
prevt     = Inf;
te        = spm_mb_shape('TemplateEnergy',mu,sett);

if do_updt_aff
    
    %------------------
    % Update shape (only affine) and appearance, on coarsest resolution
    %------------------
        
    spm_mb_show('Speak','InitAff',sett.nit.init);
    for it_init=1:nit_init
                               
        if do_updt_template
            for subit=1:nit_init_mu
                % Update template and intensity prior
                oE       = E; tic;                        
                [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
                te       = spm_mb_shape('TemplateEnergy',mu,sett);
                dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);
                E        = sum(sum(cat(2,dat.E),2),1) + te;
                t        = toc;

                % Print stuff
                fprintf('it=%i mu \t%g\t%g\t%g\t%g\n', it_init, E, t, (oE - E)/prevt, (E - oE)/abs(E));
                prevt     = t;
                Objective = [Objective; E];               
            end
        end                
            
        % Update affine
        oE = E; tic;
        dat  = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett);
        dat  = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E    = sum(sum(cat(2,dat.E),2),1) + te;
        t    = toc;        
        
        fprintf('it=%i q  \t%g\t%g\t%g\t%g\n', it_init, E, t, (oE - E)/prevt, (E - oE)/abs(E));
        prevt     = t;
        Objective = [Objective; E];
        
        if do_updt_template || do_updt_int
            % Save stuff
            save(fullfile(dir_res,'results_Groupwise.mat'),'dat','mu','sett')
        end
        
        % Show stuff
        spm_mb_show('All',dat,mu,Objective,N,sett);
        
        if it_init > 1 && (E - oErig)/abs(E) > -eps('single')*1e4
           % Finished rigid alignment
           break
        end
        oErig = E;
    end
end

%------------------
% Iteratively decrease the template resolution
%------------------

spm_mb_show('Speak','Iter',numel(sz)); tic;
for zm=numel(sz):-1:1 % loop over zoom levels
    
    if template_given && ~do_updt_template
        % Resize template
        mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
    end
    
    E0 = 0;
    if do_updt_template && (zm ~= numel(sz) || zm == 1)
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
        
    E4     = Inf;
    nit_zm = nit_zm0 + (zm - 1);
    %nit_zm = nit_zm0;
    for it_zm=1:nit_zm

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

        if (it_zm == nit_zm) && zm>1
            oMmu     = sett.var.Mmu;
            sett.var = spm_mb_io('CopyFields',sz(zm-1), sett.var);
            [dat,mu] = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
        end

        % Update deformations
        dat = spm_mb_shape('UpdateWarps',dat,sett);  
        
        % Print stuff
        fprintf('zm=%i it=%i\t%g\t%g\t%g\t%g\t%g\n', zm, it_zm, E0, E1, E2, E3, E4);        
        Objective = [Objective; E4];
                
        if do_updt_template || do_updt_int
            % Save stuff
            save(fullfile(dir_res,'results_Groupwise.mat'),'dat','mu','sett')
        end
        
        % Show stuff
        spm_mb_show('All',dat,mu,Objective,N,sett);
    end
    
    fprintf('%g seconds\n\n', toc); tic;
end

% Final mean update
[mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);

% Save template
dat = spm_mb_io('SaveTemplate',dat,mu,sett);

% Make model
model = spm_mb_io('MakeModel',dat,model,sett);

% Print total runtime
spm_mb_show('Speak','Finished',toc(t0));
end
%==========================================================================