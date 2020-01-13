function [dat,model,sett] = spm_mb_fit(data,varargin)
% Multi-Brain - Groupwise normalisation and segmentation of images
%
% FORMAT [dat,model,sett] = spm_mb_fit(data,varargin)
%
% INPUT
% data - input subjects data (see spm_mb_io('InitDat',data,sett))
%
% OUTPUT
% dat                 - struct of length N storing each subject's information (see spm_mb_io('InitDat',data,sett))
% model (inputParser) - struct storing shape and appearance model (see spm_mb_io('MakeModel',dat,model,sett))
% sett  (inputParser) - struct storing final algorithm settings (see spm_mb_param('Settings'))
%
%__________________________________________________________________________
%
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

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

% Repeatable random numbers
rng('default'); rng(1);

t0 = tic;

%------------------
% Get algorithm settings
%------------------

sett         = spm_mb_param('Settings',sett);
dir_res      = sett.write.dir_res;
do_gmm       = sett.do.gmm;
do_updt_aff  = sett.do.updt_aff;
do_zoom      = sett.do.zoom;
init_mu_dm   = sett.model.init_mu_dm;
K            = sett.model.K; 
nit_init     = sett.nit.init;
nit_init_mu  = sett.nit.init_mu;
nit_zm0      = sett.nit.zm;
print2screen = sett.show.print2screen;
vx           = sett.model.vx;
write_ws     = sett.write.workspace;

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
clear data

if isempty(dir_res) 
    pth     = fileparts(dat(1).f(1).dat.fname);
    dir_res = pth; 
end

% Get number of template classes (if not using GMM)
if ~do_gmm, [~,K] = spm_mb_io('GetSize',dat(1).f); end
if template_given
    [~,K]        = spm_mb_io('GetSize',model.shape.template);
    sett.model.K = K;
end
if isscalar(sett.model.mg_ix)
    sett.model.mg_ix = repelem(1:K + 1,sett.model.mg_ix);
end

%------------------
% Get template size and orientation
%------------------

if template_given    
    dmu       = spm_mb_io('GetSize',model.shape.template);
    [mu0,Mmu] = spm_mb_io('GetData',model.shape.template);
    sett      = spm_mb_shape('MuValOutsideFOV',mu0,sett); % For dealing with voxels outside of template's FOV (adds field sett.model.mu_bg)
else
    [Mmu,dmu] = spm_mb_shape('SpecifyMean',dat,vx,sett);
end
vxmu     = sqrt(sum(Mmu(1:3,1:3).^2));
sett.Mmu = Mmu;
sett.dmu = dmu;

%------------------
% Set affine bases
%------------------

if dmu(3) == 1 % 2D
    sett.registr.B      = spm_mb_shape('AffineBases','SE(2)'); % three parameter rigid transform
    denom_aff_tol       = N*100^3;                             % smaller convergence threshold
    sett.var.v_settings = 3/2*sett.var.v_settings;             % more v regularisation
else           % 3D
    sett.registr.B = spm_mb_shape('AffineBases','SE(3)');      % six parameter rigid transform
    denom_aff_tol  = N*100^4;
end

%------------------
% Get zoom (multi-scale) settings
%------------------

nz       = max(ceil(log2(min(dmu(dmu~=1))) - log2(init_mu_dm)),1);
if ~do_zoom, nz = 1; end
sz       = spm_mb_param('ZoomSettings',dmu,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_mb_io('CopyFields',sz(end), sett.var);

%------------------
% Init shape and apperance model parameters
%------------------

dat        = spm_mb_shape('InitDef',dat,sett);
[dat,sett] = spm_mb_appearance('Init',dat,model,K,sett);

spm_mb_show('Speak','Start',sett,N,K);

%------------------
% Init template
%------------------

if template_given    
    % Shrink given template
    mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);    
else
    % Uninformative template
    mu = zeros([sett.var.d K],'single');
    
    ix_mri = [];
    ix_ct  = [];
    for n=1:N
        if any(dat(n).is_ct == true), ix_ct  = [ix_ct n];
        else,                         ix_mri = [ix_mri n];
        end
    end

    if ~isempty(ix_ct)        
        [mu,dat(ix_mri)] = spm_mb_shape('UpdateSimpleMean',dat(ix_mri), mu, sett);
        dat              = spm_mb_appearance('UpdatePrior',dat, sett);
        [mu,dat]         = spm_mb_shape('UpdateSimpleMean',dat, mu, sett);
        dat              = spm_mb_appearance('UpdatePrior',dat, sett);
    end
end

% Save template
spm_mb_io('SaveTemplate',dat,mu,sett);
        
% Show stuff
spm_mb_show('All',dat,mu,[],N,sett);

%------------------
% Start algorithm
%------------------

Objective = [];
E         = Inf;
prevt     = Inf;
te        = 0;
if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end

if do_updt_aff
    
    %------------------
    % Update shape (only affine) and appearance, on coarsest resolution
    %------------------
        
    sett.gen.samp = min(max(vxmu(1),numel(sz)),5); % coarse-to-fine sampling of observed data
    
    spm_mb_show('Speak','InitAff',sett);
    for it_init=1:nit_init
                               
        if do_updt_template
            for subit=1:nit_init_mu
                % Update template and intensity prior
                oE       = E; tic;                        
                [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
                te       = spm_mb_shape('TemplateEnergy',mu,sett);
                dat      = spm_mb_appearance('UpdatePrior',dat, sett);
                E        = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after rigid update
                t        = toc;

                % Print stuff
                if print2screen > 0, fprintf('it=%i mu \t%g\t%g\t%g\n', it_init, E, t, (oE - E)/prevt); end
                prevt     = t;
                Objective = [Objective; E];               
            end
        end         
        
        if it_init > 1 && (oE - E)/denom_aff_tol < 1e-4
            % Finished rigid alignment
            break; 
        end        
        
        % Update affine
        oE  = E; tic;
        dat = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett);
        dat = spm_mb_appearance('UpdatePrior',dat, sett); 
        E   = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after mean update
        t   = toc;                        
        
        if print2screen > 0, fprintf('it=%i q  \t%g\t%g\t%g\n', it_init, E, t, (oE - E)/prevt); end
        prevt     = t;
        Objective = [Objective; E];        
        
        if write_ws && (do_updt_template || do_updt_int)
            % Save workspace (except template - saved as nifti separately) 
            save(fullfile(dir_res,'fit_spm_mb.mat'), '-regexp', '^(?!(mu)$).');
        end          
        
        % Save template
        spm_mb_io('SaveTemplate',dat,mu,sett);                
        
        if dmu(3) == 1, spm_mb_show('All',dat,mu,Objective,N,sett); end % If 2D, show stuff
    end            
    
    % Show stuff
    spm_mb_show('All',dat,mu,Objective,N,sett);
end

%------------------
% Iteratively decrease the template resolution
%------------------

spm_mb_show('Speak','Iter',sett,numel(sz)); 
if print2screen > 0, tic; end
for zm=numel(sz):-1:1 % loop over zoom levels
    
    % coarse-to-fine sampling of observed data    
    sett.gen.samp = min(max(vxmu(1),zm),5);     
    
    if template_given && ~do_updt_template
        mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
    end
    
    E0 = 0;
    if do_updt_template && (zm ~= numel(sz) || zm == 1)
        % If not largest zoom level
        for i=1:nit_init_mu                    
            [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
            dat      = spm_mb_appearance('UpdatePrior',dat, sett);
        end
        te = spm_mb_shape('TemplateEnergy',mu,sett);
        E0 = sum(sum(cat(2,dat.E),2),1) + te;
    end    
        
    nit_zm = nit_zm0 + (zm - 1);
    for it_zm=1:nit_zm

        % Update template                  
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end
        dat      = spm_mb_appearance('UpdatePrior',dat, sett);
        E1       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after diffeo update         
                           
        % Update affine
        dat      = spm_mb_shape('UpdateAffines',dat,mu,sett);
        dat      = spm_mb_appearance('UpdatePrior',dat, sett);
        E2       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean and int prior update

        % Update template and intensity prior (twice)
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);  
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end
        dat      = spm_mb_appearance('UpdatePrior',dat, sett);        
        E3       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after rigid update
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
        dat      = spm_mb_appearance('UpdatePrior',dat, sett);
        E4       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean and int prior update
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end
        
        % Update velocities
        dat      = spm_mb_shape('VelocityEnergy',dat,sett);
        dat      = spm_mb_shape('UpdateVelocities',dat,mu,sett);        
        dat      = spm_mb_appearance('UpdatePrior',dat, sett);        
        E5       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean and int prior update      
        dat      = spm_mb_shape('VelocityEnergy',dat,sett);

        Objective = [Objective; E5];
        
        if (it_zm == nit_zm) && zm>1
            oMmu     = sett.var.Mmu;
            sett.var = spm_mb_io('CopyFields',sz(zm-1), sett.var);
            if do_updt_template, [dat,mu] = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
            else,                dat      = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
            end
            dat      = spm_mb_shape('VelocityEnergy',dat,sett);
        end

        % Save template
        spm_mb_io('SaveTemplate',dat,mu,sett);   
        
        % Update deformations
        dat = spm_mb_shape('UpdateWarps',dat,sett);  
        
        % Print stuff
        if print2screen > 0, fprintf('zm=%i it=%i\t%g\t%g\t%g\t%g\t%g\t%g\n', zm, it_zm, E0, E1, E2, E3, E4, E5); end               
                
        if write_ws && (do_updt_template || do_updt_int)
            % Save workspace (except template - saved as nifti separately) 
            save(fullfile(dir_res,'fit_spm_mb.mat'), '-regexp', '^(?!(mu)$).');
        end                                           
        
        if dmu(3) == 1, spm_mb_show('All',dat,mu,Objective,N,sett); end % If 2D, show stuff
    end              
    
    % Show stuff
    spm_mb_show('All',dat,mu,Objective,N,sett);
        
    if print2screen > 0, fprintf('%g seconds\n\n', toc); tic; end               
end

% Final mean and intensity prior update
[mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
dat      = spm_mb_appearance('UpdatePrior',dat, sett);  

% Save template
spm_mb_io('SaveTemplate',dat,mu,sett);

% Make model
model = spm_mb_io('MakeModel',dat,model,sett);

% Print total runtime
spm_mb_show('Speak','Finished',sett,toc(t0));
end
%==========================================================================