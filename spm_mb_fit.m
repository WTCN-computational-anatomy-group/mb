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
vx           = sett.model.vx;
write_ws     = sett.write.workspace;
nit_mu       = sett.nit.init_mu;
print2screen = sett.show.print2screen;
tol          = sett.model.tol;

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
% Initial alignment of mean
%------------------

if template_given  
    % TODO: Refactor into its own function and remove hardcoding
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
te        = 0;
if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end

E  = inf(1,3);
oE = E;
if do_updt_aff    
%------------------
% Update affine only
%------------------

spm_mb_show('Speak','Affine',sett);

sett.gen.samp = numel(sz); % coarse-to-fine sampling of observed data

for it=1:nit_aff

    t = tic; % Start timer
    
    % UPDATE: mean    
    for i=1:nit_mu  
        [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett); oE(1) = E(1);
        E(1)     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after rigid update         
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end  
    end
    
    % Append to objective function vector
    Objective = [Objective; E(1)]; 
              
    % If 2D, show stuff    
    if dmu(3) == 1, spm_mb_show('All',dat,mu,Objective,N,sett); end
    
    % Check convergence
    done = (oE(1) - E(1))/E(1);
    if it > 1 && (done < tol || it == nit_aff) % Finished rigid alignment        
        spm_mb_show('PrintProgress',it,E,oE,toc(t),done,sett); % Print to command window 
        break;
    end                 

    % UPDATE: intensity prior
    dat  = spm_mb_appearance('UpdatePrior',dat, mu, sett); oE(2) = E(2);
    E(2) = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after mean update

    % UPDATE: rigid
    dat  = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett); oE(2) = E(2);       
    E(3) = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after intensity prior update                

    % Print to command window
    spm_mb_show('PrintProgress',it,E,oE,toc(t),done,sett);  
end            

if write_ws && (do_updt_template || do_updt_int)
    % Save workspace (except template - saved as nifti separately) 
    save(fullfile(dir_res,'fit_spm_mb.mat'), '-regexp', '^(?!(mu)$).');
end          

% Save template
spm_mb_io('SaveTemplate',dat,mu,sett);  
    
% Show stuff
spm_mb_show('All',dat,mu,Objective,N,sett);
end

%------------------
% Update affine and diffeo (iteratively decreases the template resolution)
%------------------

spm_mb_show('Speak','AffineDiffeo',sett,numel(sz)); 

E  = [E(1) inf(1,5)];
oE = E;
for zm=numel(sz):-1:1 % loop over zoom levels
    
    sett.gen.samp = zm; % coarse-to-fine sampling of observed data
    
    if template_given && ~do_updt_template
        mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
    end    
        
    nit_zm = nit_zm0 + (zm - 1);
    for it=1:nit_zm

        t = tic; % Start timer

        if ~(it == 1 && zm == numel(sz))
            % UPDATE: mean    
            for i=1:nit_mu  
                [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett); oE(1) = E(1);
                E(1)     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after diffeo update         
                if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end  
            end
        end        
    
        % Append to objective function vector                           
        Objective = [Objective; E(1)]; 
        
        % If 2D, show stuff
        if dmu(3) == 1, spm_mb_show('All',dat,mu,Objective,N,sett); end 
        
        % Check convergence
        done = (oE(1) - E(1))/E(1);
        if it > 1 && (done < tol || it == nit_zm) % Finished diffeo alignment            
            spm_mb_show('PrintProgress',[zm it],E,oE,toc(t),done,sett); % Print to command window
            break;
        end    
    
        % UPDATE: intensity prior
        dat  = spm_mb_appearance('UpdatePrior',dat, mu, sett); oE(2) = E(2);
        E(2) = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean update    
        
        % UPDATE: rigid
        dat  = spm_mb_shape('UpdateAffines',dat,mu,sett); oE(3) = E(3);
        E(3) = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after intensity prior update
        
        % UPDATE: mean    
        for i=1:nit_mu  
            [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett); oE(4) = E(4);
            E(4)     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after rigid update         
            if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end  
        end
        
        % UPDATE: intensity prior
        dat  = spm_mb_appearance('UpdatePrior',dat, mu, sett); oE(5) = E(5);
        E(5) = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after mean update    
        
        % UPDATE: diffeo        
        dat  = spm_mb_shape('UpdateVelocities',dat,mu,sett); oE(6) = E(6);                     
        E(6) = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after intensity prior update      
        dat  = spm_mb_shape('VelocityEnergy',dat,sett);                                
        
        % Shoot new deformations
        dat = spm_mb_shape('UpdateWarps',dat,sett);        
        
        % Print to command window
        spm_mb_show('PrintProgress',[zm it],E,oE,toc(t),done,sett);                                            
    end              
    if print2screen, fprintf('\n'); end
    
    % Show stuff
    spm_mb_show('All',dat,mu,Objective,N,sett);             
    
    if zm > 1
        oMmu     = sett.var.Mmu;
        sett.var = spm_mb_io('CopyFields',sz(zm-1), sett.var);
        [dat,mu] = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);                                
        if do_updt_template, te = spm_mb_shape('TemplateEnergy',mu,sett); end % Compute template energy
        dat      = spm_mb_shape('VelocityEnergy',dat,sett); % Compute velocity energy                
        dat      = spm_mb_shape('UpdateWarps',dat,sett);    % Shoot new deformations
    end
    
    if write_ws && (do_updt_template || do_updt_int)
        % Save workspace (except template - saved as nifti separately) 
        save(fullfile(dir_res,'fit_spm_mb.mat'), '-regexp', '^(?!(mu)$).');
    end                                           
        
    % Save template
    spm_mb_io('SaveTemplate',dat,mu,sett);  
end

% Make model
model = spm_mb_io('MakeModel',dat,model,sett);

% Print total runtime
spm_mb_show('Speak','Finished',sett,toc(t0));
end
%==========================================================================