function res = spm_mb_output(dat,model,sett)
%__________________________________________________________________________
%
% Write output from groupwise normalisation and segmentation of images.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('bf',cl,'im',cl,'imc',cl,'c',cl,'y',cl,'iy',cl,'wim',cl,'wimc',cl,'wc',cl,'mwc',cl);

% Get template
mu = spm_mb_io('GetData',model.shape.template);

for n=1:N % Loop over subjects
    res(n) = ProcessSubject(dat(n),res(n),mu,n,sett);
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% CleanGWC()
function zn = CleanGWC(zn,ixt,level)
if nargin < 2 || isempty(ixt)
    % Default SPM12 template ordering
    ixt = struct('gm',1,'wm',2,'csf',3);
end
if nargin < 3, level = 1; end

b = sum(zn(:,:,:,ixt.wm),4);

% Build a 3x3x3 seperable smoothing kernel
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Erosions and conditional dilations
th1 = 0.15;
if level==2, th1 = 0.2; end
niter  = 32;
niter2 = 32;
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = double(sum(zn(:,:,i,ixt.gm),4));
        wp       = double(sum(zn(:,:,i,ixt.wm),4));
        bp       = double(b(:,:,i));
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = bp;
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
end

% Also clean up the CSF.
if niter2 > 0
    c = b;
    for j=1:niter2
        for i=1:size(b,3)
            gp       = double(sum(zn(:,:,i,ixt.gm),4));
            wp       = double(sum(zn(:,:,i,ixt.wm),4));
            cp       = double(sum(zn(:,:,i,ixt.csf),4));
            bp       = double(c(:,:,i));
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = bp;
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
    end
end

th = 0.05;
for i=1:size(b,3)
    slices = cell(1,size(zn,4));
    for k1=1:size(zn,4)
        slices{k1} = double(zn(:,:,i,k1));
    end
    bp           = double(b(:,:,i));
    bp           = ((bp>th).*(sum(cat(3,slices{ixt.gm}),3)+sum(cat(3,slices{ixt.wm}),3)))>th;
    for i1=1:numel(ixt.gm)
        slices{ixt.gm(i1)} = slices{ixt.gm(i1)}.*bp;
    end
    for i1=1:numel(ixt.wm)
        slices{ixt.wm(i1)} = slices{ixt.wm(i1)}.*bp;
    end
    
    if niter2>0
        cp           = double(c(:,:,i));
        cp           = ((cp>th).*(sum(cat(3,slices{ixt.gm}),3)+sum(cat(3,slices{ixt.wm}),3)+sum(cat(3,slices{ixt.csf}),3)))>th;
        
        for i1=1:numel(ixt.csf)
            slices{ixt.csf(i1)} = slices{ixt.csf(i1)}.*cp;
        end        
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(zn,4)
        tot   = tot + slices{k1};
    end
    for k1=1:size(zn,4)
        zn(:,:,i,k1) = slices{k1}./tot;
    end 
end
end
%==========================================================================

%==========================================================================
% PostProcMRF()
function zn = PostProcMRF(zn,Mn,strength,nit)
P   = zeros(size(zn),'uint8');   
G   = ones([size(zn,4),1],'single')*strength;
vx  = sqrt(sum(Mn(1:3,1:3).^2));
vx2 = 1./single(vx);
for i=1:nit      
    spm_mrf(P,zn,G,vx2);        
end       
zn = single(P)/255; 
end
%==========================================================================

%==========================================================================
% ProcessSubject()
function resn = ProcessSubject(datn,resn,mun,ix,sett)

% Parse function settings
B          = sett.registr.B;
clean_def  = sett.write.clean_def;
clean_vel  = sett.write.clean_vel;
dmu        = sett.var.d;
dir_res    = sett.write.dir_res;
do_infer   = sett.do.infer;
fwhm       = sett.bf.fwhm;
gwc_level  = sett.clean_z.gwc_level;
gwc_tix    = sett.clean_z.gwc_tix;
mg_ix      = sett.model.mg_ix;
Mmu        = sett.var.Mmu;
mrf        = sett.clean_z.mrf;
nit_mrf    = sett.clean_z.nit_mrf;
reg        = sett.bf.reg;
write_bf   = sett.write.bf; % field
write_df   = sett.write.df; % forward, inverse
write_im   = sett.write.im; % image, corrected, warped, warped corrected
write_tc   = sett.write.tc; % native, warped, warped-mod

% Get parameters
[df,C] = spm_mb_io('GetSize',datn.f);
K      = size(mun,4);
K1     = K + 1;
Kmg    = numel(mg_ix);
if isa(datn.f(1),'nifti') 
    [pth,namn] = fileparts(datn.f(1).dat.fname);                
else
    namn  = ['n' num2str(ix)];
    pth   = '.';
end            
Mr    = spm_dexpm(double(datn.q),B);
Mn    = datn.Mat;                
do_bf = datn.do_bf;
is_ct = datn.is_ct;

% Set output path
if ~isempty(dir_res) && ~(exist(dir_res,'dir') == 7), mkdir(dir_res); end
if  isempty(dir_res), dir_res = pth; end
s       = what(dir_res); % Get absolute path
dir_res = s.path;

% Integrate K1 and C into write settings
if size(write_bf,1) == 1 && C  > 1, write_bf = repmat(write_bf,[C  1]); end    
if size(write_im,1) == 1 && C  > 1, write_im = repmat(write_im,[C  1]); end   
if size(write_tc,1) == 1 && K1 > 1, write_tc = repmat(write_tc,[K1 1]); end

if ~(all(write_bf(:) == false) && all(write_im(:) == false) && all(write_tc(:) == false))   
    psi0 = spm_mb_io('GetData',datn.psi);
end

if isfield(datn,'mog') && (any(write_bf(:) == true) || any(write_im(:) == true) || any(write_tc(:) == true))    
    % Input data were intensity images
    %------------------

    % Get subject-space template (softmaxed K + 1)
    psi = spm_mb_shape('Compose',psi0,spm_mb_shape('Affine',df,Mmu\Mr*Mn));    
    mun = spm_mb_shape('Pull1',mun,psi);
    clear psi
    
    % Make K + 1 template    
    mun = reshape(mun,[prod(df(1:3)) K]);
    mun = spm_mb_shape('TemplateK1',mun,2);

    if any(do_bf == true)
        % Get bias field
        chan = spm_mb_appearance('BiasFieldStruct',datn,C,df,reg,fwhm,[],datn.bf.T);
        bf   = spm_mb_appearance('BiasField',chan,df);
    else
        bf   = ones([1 C],'single');
    end
    
    % Get image(s)
    fn      = spm_mb_io('GetData',datn.f);
    fn      = reshape(fn,[prod(df(1:3)) C]);
    fn      = spm_mb_appearance('Mask',fn,is_ct);
    
    % Get labels
    labels = spm_mb_appearance('GetLabels',datn,sett);
    mun    = mun + labels;
    clear labels
    
    % Integrate use of multiple Gaussians per tissue
    mg_w = datn.mog.mg_w;
    mun  = mun(:,mg_ix);
    mun  = mun + log(mg_w);   
        
    % Format for spm_gmm
    [bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', bf.*fn);
    mun                       = spm_gmm_lib('obs2cell', mun, code_image, false);
    
    % GMM posterior
    m = datn.mog.po.m;
    b = datn.mog.po.b;
    W = datn.mog.po.W;
    n = datn.mog.po.n;

    % Get responsibilities
    zn  = spm_mb_appearance('Responsibility',m,b,W,n,bffn,mun,msk_chn); 
    zn  = spm_gmm_lib('cell2obs', zn, code_image, msk_chn);        
    clear mun msk_chn

    % Get bias field modulated image data
    fn = bf.*fn;
    if do_infer
        % Infer missing values
        code        = sum(bsxfun(@times, ~isnan(fn), 2.^(0:size(fn,2)-1)), 2);
        sample_post = do_infer > 1;
        MU          = datn.mog.po.m;    
        A           = bsxfun(@times, datn.mog.po.W, reshape(datn.mog.po.n, [1 1 Kmg]));            
        fn          = spm_gmm_lib('InferMissing',fn,zn,{MU,A},code,sample_post);  
        clear code
    end        
    
    % If using multiple Gaussians per tissue, collapse so that zn is of
    % size K1
    if Kmg > K1
        for k=1:K1, zn(:,k) = sum(zn(:,mg_ix==k),2); end
        zn(:,K1 + 1:end)    = [];
    end
    
    % Make 3D    
    if any(do_bf == true)
        bf = reshape(bf,[df(1:3) C]);
    else
        bf = reshape(bf,[1 1 1 C]);
    end
    fn = reshape(fn,[df(1:3) C]);
    zn = reshape(zn,[df(1:3) K1]);    
    
    if mrf > 0
        % Ad-hoc MRF clean-up of segmentation     
        zn = PostProcMRF(zn,Mn,mrf,nit_mrf);
    end
    
    if ~isempty(gwc_tix) && gwc_level > 0
        % Use an ad hoc brain cleanup procedure
        zn = CleanGWC(zn,gwc_tix,gwc_level);
    end        

    if any(write_bf == true) && any(do_bf == true)
        % Write bias field
        descrip = 'Bias field (';
        pths    = {};
        for c=1:C
            if ~write_bf(c,1), continue; end
            nam  = ['bf' num2str(c) '_' namn '.nii'];
            fpth = fullfile(dir_res,nam);            
            spm_mb_io('WriteNii',fpth,bf(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);                
            pths{end + 1} = fpth;
        end
        resn.bf = pths;
    end

    if any(write_im(:,1) == true)
        % Write image
        descrip = 'Image (';
        pths    = {};
        for c=1:C
            if ~write_im(c,1), continue; end
            nam  = ['im' num2str(c) '_' namn '.nii'];
            fpth = fullfile(dir_res,nam);            
            spm_mb_io('WriteNii',fpth,fn(:,:,:,c)./bf(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);
            pths{end + 1} = fpth;
        end
        resn.im = pths;

        % Write image corrected
        descrip = 'Image corrected (';
        pths    = {};
        for c=1:C
            if ~write_im(c,2), continue; end
            nam  = ['imc' num2str(c) '_' namn '.nii'];
            fpth = fullfile(dir_res,nam);            
            spm_mb_io('WriteNii',fpth,fn(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);
            pths{end + 1} = fpth;
        end
        resn.imc = pths;
    end

    if any(write_tc(:,1) == true)
        % Write segmentations
        descrip = 'Tissue (';
        pths    = {};
        for k=1:K1 
            if ~write_tc(k,1), continue; end
            nam  = ['c' num2str(k) '_' namn '.nii'];
            fpth = fullfile(dir_res,nam);            
            spm_mb_io('WriteNii',fpth,zn(:,:,:,k),Mn,[descrip 'k=' num2str(k) ')']);
            pths{end + 1} = fpth;
        end  
        resn.c = pths;
    end
else
    % Input data were segmentations
    %------------------

    zn = spm_mb_io('GetData',datn.f);
    zn = cat(4,zn,1 - sum(zn,4));
end

if any(write_df == true) || any(reshape(write_tc(:,[2 3]),[],1) == true) ||  any(reshape(write_im(:,[3 4]),[],1) == true)
    % Write forward deformation and/or normalised images
    %------------------

    % For imporved push - subsampling density in each dimension
    sd = spm_mb_shape('SampDens',Mmu,Mn);

    % Get forward deformation
    psi = spm_mb_shape('Compose',psi0,spm_mb_shape('Affine',df,Mmu\Mr*Mn));    

    if df(3) == 1, psi(:,:,:,3) = 1; end % 2D

    if write_df(1)
        % Write forward deformation
        descrip   = 'Forward deformation';
        nam       = ['y_' namn '.nii'];
        fpth      = fullfile(dir_res,nam);            
        spm_mb_io('WriteNii',fpth,psi,Mn,descrip);
        resn.y = fpth;
    end  

    if isfield(datn,'mog') && any(write_im(:,3) == true)
        % Write normalised image
        descrip = 'Normalised image (';
        pths    = {};
        for c=1:C
            if ~write_im(c,3), continue; end
            nam     = ['wim' num2str(c) '_' namn '.nii'];
            fpth    = fullfile(dir_res,nam);            
            [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c)./bf(:,:,:,c),psi,dmu,sd);
            spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),Mmu,[descrip 'c=' num2str(c) ')']);            
            pths{end + 1} = fpth;
        end
        resn.wim = pths;
    end

    if isfield(datn,'mog') && any(write_im(:,4) == true)
        % Write normalised image corrected
        descrip = 'Normalised image corrected (';
        pths    = {};
        for c=1:C
            if ~write_im(c,4), continue; end
            nam       = ['wimc' num2str(c) '_' namn '.nii'];
            fpth      = fullfile(dir_res,nam);            
            [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c),psi,dmu,sd);
            spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),Mmu,[descrip 'c=' num2str(c) ')']);            
            pths{end + 1} = fpth;
        end
        resn.wimc = pths;
    end

    if any(write_tc(:,2) == true)
        % Write normalised segmentations
        descrip = 'Normalised tissue (';
        pths    = {};
        for k=1:K1           
            if ~write_tc(k,2), continue; end
            nam       = ['wc' num2str(k) '_' namn '.nii'];
            fpth      = fullfile(dir_res,nam);            
            [img,cnt] = spm_mb_shape('Push1',zn(:,:,:,k),psi,dmu,sd);
            spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),Mmu,[descrip 'k=' num2str(k) ')']);            
            pths{end + 1} = fpth;
        end    
        resn.wc = pths;
    end  

    if any(write_tc(:,3) == true)
        % Write normalised modulated segmentations
        descrip = 'Normalised modulated tissue (';
        pths    = {};
        for k=1:K1           
            if ~write_tc(k,3), continue; end
            nam  = ['mwc' num2str(k) '_' namn '.nii'];
            fpth = fullfile(dir_res,nam);
            img  = spm_mb_shape('Push1',zn(:,:,:,k),psi,dmu);
            img  = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
            spm_mb_io('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(k) ')']);            
            pths{end + 1} = fpth;
        end    
        resn.mwc = pths;
    end  

    if write_df(2)
        % Get inverse deformation (correct?)
        psi = spm_diffeo('invdef',psi0,dmu(1:3),eye(4),eye(4));    
        %psi = spm_extrapolate_def(psi,Mmu);
        M   = inv(Mmu\Mr*Mn);
        psi = reshape(reshape(psi,[prod(dmu) 3])*M(1:3,1:3)' + M(1:3,4)',[dmu 3]);        

        % Write inverse deformation
        descrip = 'Inverse deformation';
        nam     = ['iy_' namn '.nii'];
        fpth    = fullfile(dir_res,nam);            
        spm_mb_io('WriteNii',fpth,psi,Mmu,descrip);
        resn.iy = fpth;
    end       
end
if clean_def && isa(datn.psi,'nifti') && isfile(datn.psi.dat.fname), delete(datn.psi.dat.fname); end
if clean_vel && isa(datn.v,'nifti') && isfile(datn.v.dat.fname),     delete(datn.v.dat.fname);   end
end
%==========================================================================