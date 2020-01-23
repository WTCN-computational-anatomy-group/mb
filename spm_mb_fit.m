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
nit_aff      = sett.nit.init;
nit_zm0      = sett.nit.zm;
print2screen = sett.show.print2screen;
vx           = sett.model.vx;
write_ws     = sett.write.workspace;
nit_init_mu  = sett.nit.init_mu;

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
    sett.var.v_settings = 3/2*sett.var.v_settings;             % more v regularisation
else           % 3D
    sett.registr.B = spm_mb_shape('AffineBases','SE(3)');      % six parameter rigid transform
end

%------------------
% Get zoom (multi-scale) settings
%------------------

nz       = max(ceil(log2(min(dmu(dmu~=1))) - log2(init_mu_dm)),1);
if ~do_zoom, nz = 1; end
sz       = spm_mb_param('ZoomSettings',dmu,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_mb_io('CopyFields',sz(end), sett.var);

%------------------
% Init shape model parameters
%------------------

dat = spm_mb_shape('InitDef',dat,sett);

%------------------
% Initial alignment of template
%------------------

if template_given  
    for n=1:N
        % Image params
        Vf = spm_vol(dat(n).f(1).dat.fname);
        Mn = dat(n).f(1).mat;

        % Register atlas to image to get get R (so that Mmu\R*Mf)
        mu               = spm_load_priors8(spm_vol('/scratch/Results/diffeo-segment/20200120-K11-T1w/mu_softmax_spm_mb.nii'));    
        c                = (Vf(1).dim+1)/2;
        Vf(1).mat(1:3,4) = -Mn(1:3,1:3)*c(:);
        [Affine1,ll1]    = spm_maff8(Vf(1),8,(0+1)*16,mu,[],'mni'); % Closer to rigid
        Affine1          = Affine1*(Vf(1).mat/Mn);

        % Run using the origin from the header
        Vf(1).mat     = Mn;
        [Affine2,ll2] = spm_maff8(Vf(1),8,(0+1)*16,mu,[],'mni'); % Closer to rigid

        % Pick the result with the best fit
        if ll1>ll2, R  = Affine1; else R  = Affine2; end

        % Fit final
        R = spm_maff8(dat(n).f(1).dat.fname,8,32,mu,R,'mni');
        R = spm_maff8(dat(n).f(1).dat.fname,8,1,mu,R,'mni');
        
        e  = eig(R);
        if isreal(e) && any(e<=0), disp('Possible problem!'); disp(eig(R)); end
        B1 = reshape(sett.registr.B,[16 size(sett.registr.B,3)]);
        % B1 = reshape(B,[16 size(B,3)]);
        q  = B1\reshape(real(logm(R)),[16 1]);
        
        dat(n).q = q;
    end
end

%------------------
% Init apperance model parameters
%------------------

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
    
    % Show stuff
    spm_mb_show('All',dat,mu,[],N,sett);

    % Init template with one population then use that template to init GMM
    % parameters of other populations
    [dat,mu] = spm_mb_shape('PropagateTemplate',dat,mu,sz,sett);        
end

% Save template
spm_mb_io('SaveTemplate',dat,mu,sett);
        
% Show stuff
spm_mb_show('All',dat,mu,[],N,sett);

%------------------
% Start algorithm
%------------------

Objective = [];
E1        = Inf;
te        = 0;
if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end

if do_updt_aff
    
    %------------------
    % Update shape (only affine) and appearance, on coarsest resolution
    %------------------
    
    spm_mb_show('Speak','InitAff',sett);
    
    sett.gen.samp = numel(sz); % coarse-to-fine sampling of observed data
        
    for it_init=1:nit_aff
                             
        % UPDATE: mean        
        oE1      = E1;
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);        
        E1       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after rigid update         
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end  
                
        Objective = [Objective; E1]; % Append to objective function vector
        
        if write_ws && (do_updt_template || do_updt_int)
            % Save workspace (except template - saved as nifti separately) 
            save(fullfile(dir_res,'fit_spm_mb.mat'), '-regexp', '^(?!(mu)$).');
        end          
        
        % Save template
        spm_mb_io('SaveTemplate',dat,mu,sett);                
        
        if dmu(3) == 1, spm_mb_show('All',dat,mu,Objective,N,sett); end % If 2D, show stuff    
        
        if (oE1 - E1)/(N*prod(dmu(1:3))) < 1e-4 || it_init == nit_aff
            % Finished rigid alignment
            break; 
        end   
        
        % UPDATE: mean                  
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);        
        E2       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean update         
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end                  
        
        % UPDATE: intensity prior
        dat = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E3  = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after mean update
                             
        % UPDATE: rigid
        dat = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett);        
        E4  = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after intensity prior update                

        if print2screen > 0, fprintf('it=%i\tE1=%g\tE2=%g\tE3=%g\tE4=%g\n', it_init, E1, E2, E3, E4); end                                 
    end            
        
    spm_mb_show('All',dat,mu,Objective,N,sett); % Show stuff
end

%------------------
% Iteratively decrease the template resolution
%------------------

spm_mb_show('Speak','Iter',sett,numel(sz)); 
if print2screen > 0, tic; end
for zm=numel(sz):-1:1 % loop over zoom levels
    
    sett.gen.samp = zm; % coarse-to-fine sampling of observed data
    
    if template_given && ~do_updt_template
        mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
    end
    
    E0 = 0;
    if (zm ~= numel(sz) || zm == 1)        
        for i=1:nit_init_mu                   
            [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
        end        
        E0 = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean update
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end
    end    
        
    nit_zm = nit_zm0 + (zm - 1);
    for it_zm=1:nit_zm

        % UPDATE: mean                  
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);        
        E1       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean update         
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end        
                           
        % UPDATE: intensity prior
        dat = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E2  = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean update    
        
        % UPDATE: rigid
        dat = spm_mb_shape('UpdateAffines',dat,mu,sett);
        E3  = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after intensity prior update

        % UPDATE: mean
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);              
        E4       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after rigid update
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end
        
        % UPDATE: intensity prior
        dat = spm_mb_appearance('UpdatePrior',dat, mu, sett);
        E5  = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean update    
        
        % UPDATE: velocities        
        dat = spm_mb_shape('UpdateVelocities',dat,mu,sett);                      
        E6  = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after intensity prior update      
        dat = spm_mb_shape('VelocityEnergy',dat,sett);
        
        % UPDATE: mean
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);               
        E7       = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after diffeo update
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end
        
        Objective = [Objective; E7]; % Append to objective function vector
        
        % Print stuff
        if print2screen > 0, fprintf('zm=%i it=%i\tE0=%g\tE1=%g\tE2=%g\tE3=%g\tE4=%g\tE5=%g\tE6=%g\tE7=%g\n', zm, it_zm, E0, E1, E2, E3, E4, E5, E6, E7); end               
        
        if (it_zm == nit_zm) && zm>1
            oMmu     = sett.var.Mmu;
            sett.var = spm_mb_io('CopyFields',sz(zm-1), sett.var);
            if do_updt_template
                [dat,mu] = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);                                
                te       = spm_mb_shape('TemplateEnergy',mu,sett); % Compute template energy
            else
                dat = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
            end                        
            dat = spm_mb_shape('VelocityEnergy',dat,sett); % Compute velocity energy
        end

        % Save template
        spm_mb_io('SaveTemplate',dat,mu,sett);   
        
        % Shoot new deformations
        dat = spm_mb_shape('UpdateWarps',dat,sett);                  
                
        if write_ws && (do_updt_template || do_updt_int)
            % Save workspace (except template - saved as nifti separately) 
            save(fullfile(dir_res,'fit_spm_mb.mat'), '-regexp', '^(?!(mu)$).');
        end                                           
        
        if dmu(3) == 1, spm_mb_show('All',dat,mu,Objective,N,sett); end % If 2D, show stuff
    end              
        
    spm_mb_show('All',dat,mu,Objective,N,sett); % Show stuff
        
    if print2screen > 0, fprintf('%g seconds\n\n', toc); tic; end               
end

% Final mean and intensity prior update
[mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
dat      = spm_mb_appearance('UpdatePrior',dat, mu, sett);  

% Save template
spm_mb_io('SaveTemplate',dat,mu,sett);

% Make model
model = spm_mb_io('MakeModel',dat,model,sett);

% Print total runtime
spm_mb_show('Speak','Finished',sett,toc(t0));
end
%==========================================================================