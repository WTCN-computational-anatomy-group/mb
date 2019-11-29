function varargout = spm_multireg(varargin)
%__________________________________________________________________________
%
% Subject- or population-level image segmentation + normalisation.
%
% FORMAT [dat,mu,sett] = spm_multireg('Groupwise',F,sett)
% FORMAT [dat,mu,sett] = spm_multireg('Register',F,mu,sett)
% FORMAT res             spm_multireg('WriteResults',dat,mu,sett)
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
    case 'WriteResults'
        [varargout{1:nargout}] = WriteResults(varargin{:});        
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

% Clear figures
spm_multireg_show('Clear',sett);

%------------------
% Init dat
%------------------

dat = spm_multireg_init('InitDat',F,sett); 
F   = [];
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
            dat      = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
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
        dat  = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
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
            dat      = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
            
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
        dat      = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
        E1       = sum(sum(cat(2,dat.E),2),1) + te;        
                           
        % Update affine
        % (Might be an idea to run this less often - currently slow)
        dat      = spm_multireg_updt('UpdateAffines',dat,mu,sett);
        dat      = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
        E2       = sum(sum(cat(2,dat.E),2),1) + te;

        % Update template, bias field and intensity model
        % (Might be an idea to run this multiple times)                
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett); % An extra mean iteration
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);
        dat      = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
                
        [mu,dat] = spm_multireg_updt('UpdateMean',dat, mu, sett);
        te       = spm_multireg_energ('TemplateEnergy',mu,sett);
        dat      = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
        E3       = sum(sum(cat(2,dat.E),2),1) + te;
            
        % Update velocities
        dat      = spm_multireg_energ('VelocityEnergy',dat,sett);
        dat      = spm_multireg_updt('UpdateVelocities',dat,mu,sett);
        dat      = spm_multireg_energ('VelocityEnergy',dat,sett);
        dat      = spm_multireg_updt('UpdateIntensityPrior',dat, mu, sett);
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

if sett.show.level >= 1
    % Show stuff
    spm_multireg_show('ShowModel',mu,Objective,N,sett);
    spm_multireg_show('ShowSubjects',dat,mu,sett);
    spm_multireg_show('ShowParameters',dat,mu,sett);
    spm_multireg_show('ShowBiasField',dat,sett);
    spm_multireg_show('ShowIntensityModel',dat,sett);
end
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

dat = spm_multireg_init('InitDat',F,sett); 
F   = [];
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
% WriteResults()
function res = WriteResults(dat,mu0,sett)

% Parse function settings
B        = sett.registr.B;
dmu      = sett.var.d;
dir_res  = sett.write.dir_res;
do_infer = sett.do.infer;
fwhm     = sett.bf.fwhm;
Mmu      = sett.var.Mmu;
reg      = sett.bf.reg;
write_bf = sett.write.bf; % field
write_df = sett.write.df; % forward, inverse
write_im = sett.write.im; % image, corrected, warped, warped corrected
write_tc = sett.write.tc; % native, warped, warped-mod

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('bf',cl,'im',cl,'imc',cl,'c',cl,'y',cl,'iy',cl,'wim',cl,'wimc',cl,'wc',cl,'mwc',cl);

for n=1:N % Loop over subjects
    
    % Get parameters
    [df,C] = spm_multireg_io('GetSize',dat(n).f);
    K      = size(mu0,4);
    K1     = K + 1;
    if isa(dat(n).f(1),'nifti'), [~,namn] = fileparts(dat(n).f(1).dat.fname);                
    else,                           namn  = ['n' num2str(n)];
    end            
    Mr = spm_dexpm(double(dat(n).q),B);
    Mn = dat(n).Mat;                
    
    % Integrate K1 and C into write settings
    if size(write_bf,1) == 1 && C  > 1, write_bf = repmat(write_bf,[C  1]); end    
    if size(write_im,1) == 1 && C  > 1, write_im = repmat(write_im,[C  1]); end   
    if size(write_tc,1) == 1 && K1 > 1, write_tc = repmat(write_tc,[K1 1]); end
    
    if any(write_bf(:) == true) || any(write_im(:) == true) || any(write_tc(:) == true)
        if isfield(dat(n),'mog')   
            % Input data were intensity images
            %------------------

            % Get subject-space template (softmaxed K + 1)
            psi1 = spm_multireg_io('GetData',dat(n).psi);
            psi0 = spm_multireg_util('Affine',df,Mmu\Mr*Mn);
            psi  = spm_multireg_util('Compose',psi1,psi0);
            psi0 = []; psi1 = [];        

            mu  = spm_multireg_util('Pull1',mu0,psi);
            psi = [];

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
            L    = unique(code);
            
            % GMM posterior
            m  = dat(n).mog.po.m;
            b  = dat(n).mog.po.b;
            V  = dat(n).mog.po.V;
            nu = dat(n).mog.po.n;

            % Get responsibilities
            zn = spm_multireg_io('ComputeResponsibilities',m,b,V,nu,bf.*fn,mu,L,code); 
            mu = [];     

            % Get bias field modulated image data
            fn = bf.*fn;
            if do_infer
                % Infer missing values
                sample_post = do_infer > 1;
                MU = dat(n).mog.po.m;    
                A  = bsxfun(@times, dat(n).mog.po.V, reshape(dat(n).mog.po.n, [1 1 K1]));            
                fn = spm_gmm_lib('InferMissing',fn,zn,{MU,A},{code,unique(code)},sample_post);        
            end

            % TODO: Possible post-processing (MRF + clean-up)


            % Make 3D        
            bf = reshape(bf,[df(1:3) C]);
            fn = reshape(fn,[df(1:3) C]);
            zn = reshape(zn,[df(1:3) K1]);

            if any(write_bf == true)
                % Write bias field
                descrip = 'Bias field (';
                pths    = {};
                for c=1:C
                    if ~write_bf(c,1), continue; end
                    nam  = ['bf' num2str(c) '_' namn '.nii'];
                    fpth = fullfile(dir_res,nam);            
                    spm_multireg_util('WriteNii',fpth,bf(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);                
                    pths{end + 1} = fpth;
                end
                res(n).bf = pths;
            end

            if any(write_im(:,1) == true)
                % Write image
                descrip = 'Image (';
                pths    = {};
                for c=1:C
                    if ~write_im(c,1), continue; end
                    nam  = ['im' num2str(c) '_' namn '.nii'];
                    fpth = fullfile(dir_res,nam);            
                    spm_multireg_util('WriteNii',fpth,fn(:,:,:,c)./bf(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);
                    pths{end + 1} = fpth;
                end
                res(n).im = pths;

                % Write image corrected
                descrip = 'Image corrected (';
                pths    = {};
                for c=1:C
                    if ~write_im(c,2), continue; end
                    nam  = ['imc' num2str(c) '_' namn '.nii'];
                    fpth = fullfile(dir_res,nam);            
                    spm_multireg_util('WriteNii',fpth,fn(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);
                    pths{end + 1} = fpth;
                end
                res(n).imc = pths;
            end

            if any(write_tc(:,1) == true)
                % Write segmentations
                descrip = 'Tissue (';
                pths    = {};
                for k=1:K1 
                    if ~write_tc(k,1), continue; end
                    nam  = ['c' num2str(k) '_' namn '.nii'];
                    fpth = fullfile(dir_res,nam);            
                    spm_multireg_util('WriteNii',fpth,zn(:,:,:,k),Mn,[descrip 'k=' num2str(k) ')']);
                    pths{end + 1} = fpth;
                end  
                res(n).c = pths;
            end
        else
            % Input data were segmentations
            %------------------
            
            zn = spm_multireg_io('GetData',dat(n).f);
            zn = cat(4,zn,1 - sum(zn,4));
        end    
    end

    if any(write_df == true) || any(reshape(write_tc(:,[2 3]),[],1) == true) ||  any(reshape(write_im(:,[3 4]),[],1) == true)
        % Write forward deformation and/or normalised images
        %------------------
        
        % For imporved push - subsampling density in each dimension
        sd = spm_multireg_util('SampDens',Mmu,Mn);
        
        % Get forward deformation
        psi1 = spm_multireg_io('GetData',dat(n).psi);
        psi  = spm_multireg_util('Compose',psi1,spm_multireg_util('Affine',df,Mmu\Mr*Mn));
        psi1 = [];

        if df(3) == 1, psi(:,:,:,3) = 1; end % 2D

        if write_df(1)
            % Write forward deformation
            descrip   = 'Forward deformation';
            nam       = ['y_' namn '.nii'];
            fpth      = fullfile(dir_res,nam);            
            spm_multireg_util('WriteNii',fpth,psi,Mn,descrip);
            res(n).y = fpth;
        end  

        if isfield(dat(n),'mog') && any(write_im(:,3) == true)
            % Write normalised image
            descrip = 'Normalised image (';
            pths    = {};
            for c=1:C
                if ~write_im(c,3), continue; end
                nam  = ['wim' num2str(c) '_' namn '.nii'];
                fpth = fullfile(dir_res,nam);            
                img   = spm_multireg_util('Push1',fn(:,:,:,c)./bf(:,:,:,c),psi,dmu,sd);
                spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'c=' num2str(c) ')']);            
                pths{end + 1} = fpth;
            end
            res(n).wim = pths;
        end

        if isfield(dat(n),'mog') && any(write_im(:,4) == true)
            % Write normalised image corrected
            descrip = 'Normalised image corrected (';
            pths    = {};
            for c=1:C
                if ~write_im(c,4), continue; end
                nam  = ['wimc' num2str(c) '_' namn '.nii'];
                fpth = fullfile(dir_res,nam);            
                img   = spm_multireg_util('Push1',fn(:,:,:,c),psi,dmu,sd);
                spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'c=' num2str(c) ')']);            
                pths{end + 1} = fpth;
            end
            res(n).wimc = pths;
        end

        if any(write_tc(:,2) == true)
            % Write normalised segmentations
            descrip = 'Normalised tissue (';
            pths    = {};
            for k=1:K1           
                if ~write_tc(k,2), continue; end
                nam  = ['wc' num2str(k) '_' namn '.nii'];
                fpth = fullfile(dir_res,nam);            
                img   = spm_multireg_util('Push1',zn(:,:,:,k),psi,dmu,sd);
                spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(k) ')']);            
                pths{end + 1} = fpth;
            end    
            res(n).wc = pths;
        end  
        
        if any(write_tc(:,3) == true)
            % Write normalised modulated segmentations (correct?)
            descrip = 'Normalised modulated tissue (';
            pths    = {};
            for k=1:K1           
                if ~write_tc(k,3), continue; end
                nam   = ['mwc' num2str(k) '_' namn '.nii'];
                fpth  = fullfile(dir_res,nam);
                img   = spm_multireg_util('Push1',zn(:,:,:,k),psi,dmu);
                img   = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
                spm_multireg_util('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(k) ')']);            
                pths{end + 1} = fpth;
            end    
            res(n).mwc = pths;
        end  
                
        if write_df(2)
            % Get inverse deformation (correct?)
            psi = spm_multireg_io('GetData',dat(n).psi);    
            psi = spm_diffeo('invdef',psi,dmu(1:3),eye(4),eye(4));    
            %psi = spm_extrapolate_def(psi,Mmu);
            M   = inv(Mmu\Mr*Mn);
            psi = reshape(reshape(psi,[prod(dmu) 3])*M(1:3,1:3)' + M(1:3,4)',[dmu 3]);        
        
            % Write inverse deformation
            descrip = 'Inverse deformation';
            nam     = ['iy_' namn '.nii'];
            fpth    = fullfile(dir_res,nam);            
            spm_multireg_util('WriteNii',fpth,psi,Mmu,descrip);
            res(n).iy = fpth;
        end       
    end
end
end
%==========================================================================
