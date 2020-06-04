function [dat,mu,sett,model] = spm_mb_fit(data,varargin)
% Multi-Brain - Groupwise normalisation and segmentation of images
%
% FORMAT [dat,mu,sett,model] = spm_mb_fit(data,varargin)
%
% INPUT
% data - input subjects data (see spm_mb_io('InitDat',data,sett))
%
% OUTPUT
% dat                 - struct of length N storing each subject's information (see spm_mb_io('InitDat',data,sett))
% mu                  - array with template data
% sett  (inputParser) - struct storing final algorithm settings (see spm_mb_param('Settings'))
% model (inputParser) - struct storing shape and appearance model (see spm_mb_io('MakeModel',dat,model,sett))
%
%__________________________________________________________________________
%
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% Parse input
p              = inputParser;
p.FunctionName = 'spm_mb_fit';
p.addParameter('PthModel', '',       @ischar);
p.addParameter('sett',     struct(), @isstruct);
p.parse(varargin{:});
PthModel       = p.Results.PthModel;
sett           = p.Results.sett;

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
init_mu_dm   = sett.model.init_mu_dm;
K            = sett.model.K;
nit_aff      = sett.nit.init;
nit_zm0      = sett.nit.zm;
vx           = sett.model.vx;
write_ws     = sett.write.workspace;
nit_mu       = sett.nit.init_mu;
print2screen = sett.show.print2screen;
tol_aff      = sett.model.tol_aff;
tol_diffeo   = sett.model.tol_diffeo;
updt_diff    = sett.do.updt_vel;
updt_aff     = sett.do.updt_aff;
samp_min     = sett.gen.samp_min;

spm_mb_show('Clear',sett); % Clear figures

N = numel(data); % Total number of subjects

%------------------
% Load model (if given)
%------------------

model = spm_mb_io('LoadModel',PthModel,sett);

%------------------
% Decide what to learn
%------------------

[template_given,~,sett] = spm_mb_param('SetFit',model,sett);
updt_mu                 = sett.do.updt_template;
updt_intp               = sett.do.updt_int;

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

%------------------
% Get template size and orientation
%------------------

if template_given
    dmu               = spm_mb_io('GetSize',model.shape.template);
    [mu0,Mmu]         = spm_mb_io('GetData',model.shape.template);
    sett.model.pth_mu = model.shape.template;  % Store path to template
else
    [Mmu,dmu]         = spm_mb_shape('SpecifyMean',dat,vx,sett);
    sett.model.pth_mu = fullfile(sett.write.dir_res,'spm_mb_mu.nii');  % Store path to template
end
if dmu(3) == 1 % for 2D data
    vx  = sqrt(sum(Mmu(1:3,1:3).^2));
    Mmu = [diag(vx) zeros(3,1); 0 0 0 1];
end
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
sz       = spm_mb_shape('ZoomSettings',dmu,Mmu,sett.var.v_settings,sett.var.mu_settings,nz);
sett.var = spm_mb_io('CopyFields',sz(end), sett.var);

%------------------
% Init shape model parameters
%------------------

dat = spm_mb_shape('InitDef',dat,model,sett);

%------------------
% Init apperance model parameters
%------------------

[dat,sett] = spm_mb_init('Init',dat,model,K,sett);

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

    % Init template with one population then use that template to init GMM
    % parameters of other populations
    [dat,mu] = spm_mb_init('PropagateTemplate',dat,mu,sz,sett);
end

%------------------
% Start algorithm
%------------------

Objective = [];
te        = 0;
if updt_mu, te = spm_mb_shape('TemplateEnergy',mu,sett); end
E         = inf(1,max(1,sum([updt_mu, updt_aff]))); % For tracking objfun
oE        = E;

if updt_aff
%------------------
% Update affine only
%------------------

spm_mb_show('Speak','Affine',sett);
spm_mb_show('PrintDateTime',sett);

sett.gen.samp = numel(sz); % coarse-to-fine sampling of observed data

for it0=1:nit_aff

    t = tic; % Start timer
    i = 1;   % For tracking objfun

    if ~(updt_mu && it0 == 1)
        % UPDATE: rigid
        dat   = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett);
        oE(i) = E(i);
        E(i)  = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
        if updt_intp, dat = spm_mb_appearance('UpdatePrior',dat, sett); end
    end
    i = i + 1;

    if updt_mu
        % UPDATE: mean
        for it1=1:max(nit_mu - 1,1)
            [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
            oE(i)    = E(i);
            E(i)     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
            if updt_intp, dat = spm_mb_appearance('UpdatePrior',dat, sett); end
            te       = spm_mb_shape('TemplateEnergy',mu,sett);
        end
    end

    % Check convergence
    ix_done   = min(numel(E),2);
    Objective = [Objective, E];
    done      = abs(oE(ix_done) - E(ix_done))/abs(E(ix_done));

    % Print to command window
    spm_mb_show('PrintProgress',it0,E,oE,toc(t),done,tol_aff,sett);

    % Show stuff
    spm_mb_show('All',dat,mu,Objective,N,sett);

    if do_gmm && it0 == 1 && (numel(sett.model.mg_ix_intro) ~= numel(sett.model.mg_ix))
        % Introduce multiple Gaussians per tissue
        [dat,sett] = spm_mb_appearance('IntroduceMG',dat,sett);
    end

    % Finished rigid alignment?
    if done < tol_aff, break; end
end

if write_ws
    % Save workspace (except template - saved as nifti separately)
    save(fullfile(dir_res,'spm_mb_workspace.mat'), '-regexp', '^(?!(mu)$).');
end

% Save template
spm_mb_io('SaveTemplate',dat,mu,sett);
end

%------------------
% Update affine and diffeo (iteratively decreases the template resolution)
%------------------

spm_mb_show('Speak','AffineDiffeo',sett,numel(sz));

E  = inf(1,sum([(updt_aff + updt_diff)*updt_mu, updt_aff, updt_diff])); % For tracking objfun
oE = E;

for zm=numel(sz):-1:1 % loop over zoom levels

    spm_mb_show('PrintDateTime',sett);
    
    sett.gen.samp = max(zm,samp_min); % coarse-to-fine sampling of observed data

    if template_given && ~updt_mu
        mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
    end

    nit_zm = nit_zm0 + (zm - 1); % use nit_zm0 only for zm = 1
    for it0=1:nit_zm

        t = tic; % Start timer
        i = 1;   % For tracking objfun

        if updt_diff
            % UPDATE: diffeo
            dat   = spm_mb_shape('UpdateVelocities',dat,mu,sett);
            oE(i) = E(i);
            E(i)  = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
            dat   = spm_mb_shape('VelocityEnergy',dat,sett);
            if updt_intp, dat = spm_mb_appearance('UpdatePrior',dat, sett); end

            % Shoot new deformations
            dat = spm_mb_shape('UpdateWarps',dat,sett);
            i   = i + 1;
        end

        if updt_diff && updt_mu
            % UPDATE: mean
            for it1=1:nit_mu
                [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
                oE(i)    = E(i);
                E(i)     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
                if updt_intp, dat = spm_mb_appearance('UpdatePrior',dat, sett); end
                te       = spm_mb_shape('TemplateEnergy',mu,sett); 
            end
            i = i + 1;
        end

        if updt_aff
            % UPDATE: rigid
            dat   = spm_mb_shape('UpdateAffines',dat,mu,sett);
            oE(i) = E(i);
            E(i)  = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
            if updt_intp, dat = spm_mb_appearance('UpdatePrior',dat, sett); end
            i     = i + 1;
        end

        if updt_aff && updt_mu
            % UPDATE: mean
            for it1=1:max(nit_mu - 1,1) % one less iteration for affine than diffeo
                [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
                oE(i)    = E(i);
                E(i)     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
                if updt_intp, dat = spm_mb_appearance('UpdatePrior',dat, sett); end
                te       = spm_mb_shape('TemplateEnergy',mu,sett);
            end
        end

        % Check convergence
        ix_done   = updt_aff + updt_diff;
        Objective = [Objective, E];
        done      = abs(oE(ix_done) - E(ix_done))/abs(E(ix_done));        

        % Print to command window
        spm_mb_show('PrintProgress',[zm it0],E,oE,toc(t),done,tol_diffeo,sett);

        % Show stuff
        spm_mb_show('All',dat,mu,Objective,N,sett);

        % Finished diffeo alignment?
        if done < tol_diffeo, break; end
    end
    if print2screen, fprintf('\n'); end

    if zm > 1
        oMmu           = sett.var.Mmu;
        sett.var       = spm_mb_io('CopyFields',sz(zm-1), sett.var);
        if      updt_mu &&  updt_diff
            [dat,mu]   = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
        elseif  updt_mu && ~updt_diff
            [~,mu]     = spm_mb_shape('ZoomVolumes',[],mu,sett,oMmu);
        elseif ~updt_mu &&  updt_diff
            dat        = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
        end
        if updt_mu, te = spm_mb_shape('TemplateEnergy',mu,sett); end % Compute template energy
        if updt_diff
            dat        = spm_mb_shape('VelocityEnergy',dat,sett);    % Compute velocity energy
            dat        = spm_mb_shape('UpdateWarps',dat,sett);       % Shoot new deformations
        end
    end

    if write_ws
        % Save workspace (except template - saved as nifti separately)
        save(fullfile(dir_res,'spm_mb_workspace.mat'), '-regexp', '^(?!(mu)$).');
    end

    % Save template
    spm_mb_io('SaveTemplate',dat,mu,sett);
end

% Make model
model = spm_mb_io('MakeModel',dat,model,sett);

% Clean-up temporary files
if exist(sett.show.dir_vis,'dir') == 7, rmdir(sett.show.dir_vis,'s'); end

% Print total runtime
spm_mb_show('Speak','Finished',sett,toc(t0));
end
%==========================================================================
