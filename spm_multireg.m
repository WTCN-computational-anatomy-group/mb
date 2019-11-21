function varargout = spm_multireg(varargin)
%__________________________________________________________________________
%
% Subject- or population-level image segmentation + normalisation.
%
% FORMAT [dat,mu,sett] = spm_multireg('Groupwise',F,sett)
% FORMAT [dat,mu,sett] = spm_multireg('Register',F,mu,sett)
% FORMAT                 spm_multireg('WriteNormalised',dat,mu,sett)
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
    case 'WriteNormalised'
        [varargout{1:nargout}] = WriteNormalised(varargin{:});        
    otherwise
        help spm_multireg
        error('Unknown function %s. Type ''help spm_multireg'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% Groupwise()
function [dat,mu,sett] = Groupwise(F,sett)
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

        % Save stuff
        save(fullfile(dir_res,'results_Groupwise.mat'),'dat','mu','sett')
        
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

% Save template
dat = spm_multireg_io('SaveTemplate',dat,mu,sett);

% Print total runtime
spm_multireg_show('Speak','Finished',toc(t0));

end
%==========================================================================

%==========================================================================
% Register()
function [dat,mu,sett] = Register(F,mu,sett)
if nargin < 3, sett = struct; end

t0 = tic;

%------------------
% Get algorithm settings
%------------------

sett                 = spm_multireg_par('Settings',sett);
sett.model.groupwise = false;
sett.gen.threads     = 0;

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

%==========================================================================
% WriteNormalised()
function WriteNormalised(dat,mu0,sett)

% Parse function settings
B        = sett.registr.B;
dmu      = sett.var.d;
dir_res  = sett.write.dir_res;
do_infer = sett.do.infer;
fwhm     = sett.bf.fwhm;
Mmu      = sett.var.Mmu;
reg      = sett.bf.reg;

for n=1:numel(dat)
    
    % Get parameters
    [df,C] = spm_multireg_io('GetSize',dat(n).f);
    q      = double(dat(n).q);
    Mr     = spm_dexpm(q,B);
    Mn     = dat(n).Mat;        
    if isa(dat(n).f(1),'nifti')
        [~,namn] = fileparts(dat(n).f(1).dat.fname);                
    else
        namn     = ['n' num2str(n)];
    end
           
    if isfield(dat(n),'mog')   
        % Get subject-space template (softmaxed K + 1)
        psi1 = spm_multireg_io('GetData',dat(n).psi);
        psi0 = spm_multireg_util('Affine',df,Mmu\Mr*Mn);
        psi  = spm_multireg_util('Compose',psi1,psi0);
        clear psi0 psi1
        
        mu = spm_multireg_util('Pull1',mu0,psi);
        clear psi
        
        mu = log(spm_multireg_util('softmaxmu',mu,4));
        mu = reshape(mu,[prod(df(1:3)) size(mu,4)]);

        % Get bias field
        chan = spm_multireg_io('GetBiasFieldStruct',C,df,Mn,reg,fwhm,[],dat(n).bf.T);
        bf   = spm_multireg_io('GetBiasField',chan,df);
        
        % Get image(s)
        fn   = spm_multireg_io('GetData',dat(n).f);
        fn   = reshape(fn,[prod(df(1:3)) C]);
        fn   = spm_multireg_util('MaskF',fn);
        code = spm_gmm_lib('obs2code', fn);
        
        % Get responsibilities
        zn = spm_multireg_io('ComputeResponsibilities',dat(n),bf.*fn,mu,code); clear mu
        K1 = size(zn,2);
                
        % Get bias field modulated image data
        fn = bf.*fn; clear bf
        if do_infer
            % Inger missing values
            sample_post = do_infer > 1;
            MU = dat(n).mog.po.m;    
            A  = bsxfun(@times, dat(n).mog.po.V, reshape(dat(n).mog.po.n, [1 1 K1]));            
            fn = spm_gmm_lib('InferMissing',fn,zn,{MU,A},{code,unique(code)},sample_post);        
        end
        clear code
        
         % Make 3D        
        zn = reshape(zn,[df(1:3) K1]);
        fn = reshape(fn,[df(1:3) C]);
        
        % Write segmentations
        descrip = 'Tissue (';
        for k=1:K1           
            nam  = ['c' num2str(k) namn '.nii'];
            fpth = fullfile(dir_res,nam);            
            spm_multireg_util('WriteNii',fpth,zn(:,:,:,k),Mn,[descrip 'k=' num2str(k) ')']);
        end  
    end
    
    % Get inverse deformation
    psi = spm_multireg_io('GetData',dat(n).psi);    
    psi = spm_diffeo('invdef',psi,dmu(1:3),eye(4),eye(4));    
%     psi = spm_extrapolate_def(psi,Mmu);
    M   = inv(Mmu\Mr*Mn);
    psi = reshape(reshape(psi,[prod(dmu) 3])*M(1:3,1:3)' + M(1:3,4)',[dmu 3]);        
    
    % Trilinear resampling in template space    
    if isfield(dat(n),'mog')    
        % Input data are intensity images        
        descrip = 'Normalised image (';
        for c=1:C
            nam = ['ic' num2str(c) namn '.nii'];
            fpth   = fullfile(dir_res,nam);            
            img = single(spm_diffeo('bsplins',fn(:,:,:,c),psi,[1 1 1 0 0 0]));
            spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(c) ')']);            
        end
        
        descrip = 'Normalised tissue (';
        for k=1:K1           
            nam = ['wc' num2str(k) namn '.nii'];
            fpth   = fullfile(dir_res,nam);            
            img = single(spm_diffeo('bsplins',zn(:,:,:,k),psi,[1 1 1 0 0 0]));
            spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(k) ')']);            
        end        
    else       
        % Input data are segmentations
        descrip = 'Normalised tissue (';
        fn      = spm_multireg_io('GetData',dat(n).f);
        for k=1:C            
            nam = ['wc' num2str(k) namn '.nii'];
            fpth   = fullfile(dir_res,nam);            
            img = single(spm_diffeo('bsplins',fn(:,:,:,k),psi,[1 1 1 0 0 0]));
            spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(k) ')']);            
        end        
        nam = ['wc' num2str(k + 1) namn '.nii'];
        fpth   = fullfile(dir_res,nam);
        img = 1 - sum(img,4);
        img = single(spm_diffeo('bsplins',img,psi,[1 1 1 0 0 0]));
        spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(k) ')']);
    end        
end
end
%==========================================================================