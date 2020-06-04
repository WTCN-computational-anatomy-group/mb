function res = spm_mb_output(dat,mu,sett)
%__________________________________________________________________________
%
% Write output from groupwise normalisation and segmentation of images.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('bf',cl,'im',cl,'imc',cl,'c',cl,'y',cl,'wim',cl, ...
             'wimc',cl,'wc',cl,'mwc',cl,'v',cl, 'sm', cl);

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
    if j>2
        th       = th1;
    else
        th       = 0.6;
    end  % Dilate after two its of erosion
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
function resn = ProcessSubject(datn,resn,mun0,ix,sett)

% Parse function settings
B          = sett.registr.B;
clean_def  = sett.write.clean_def;
clean_vel  = sett.write.clean_vel;
dmu        = sett.var.d;
dir_res    = sett.write.dir_res;
do_infer   = sett.do.infer;
gwc_level  = sett.clean_z.gwc_level;
gwc_tix    = sett.clean_z.gwc_tix;
mg_ix      = sett.model.mg_ix;
Mmu        = sett.var.Mmu;
mu_bg      = sett.model.mu_bg;
mrf        = sett.clean_z.mrf;
nit_mrf    = sett.clean_z.nit_mrf;
reg        = sett.bf.reg;
write_bf   = sett.write.bf; % field
write_df   = sett.write.df; % forward, inverse
write_im   = sett.write.im; % image, corrected, warped, warped corrected
write_tc   = sett.write.tc; % native, warped, warped-mod
write_vel  = sett.write.vel;
write_aff  = sett.write.affine;
has_spine  = sett.gen.has_spine;
write_sm   = sett.write.scal_mom;
pth_mu     = sett.model.pth_mu;
bb         = sett.write.bb;
vx_out     = sett.write.vx;

% Get parameters
[df,C] = spm_mb_io('GetSize',datn.f);
K      = size(mun0,4);
K1     = K + 1;
Kmg    = numel(mg_ix);
namn   = datn.nam;
Mr     = spm_dexpm(double(datn.q),B);
Mn     = datn.Mat;
do_bf  = datn.do_bf;

% Set output path
if ~isempty(dir_res) && ~(exist(dir_res,'dir') == 7), mkdir(dir_res); end
if  isempty(dir_res), dir_res = pth; end
s       = what(dir_res); % Get absolute path
dir_res = s.path;

% Integrate K1 and C into write settings
if size(write_bf,1) == 1 && C  > 1, write_bf = repmat(write_bf,[C  1]); end
if size(write_im,1) == 1 && C  > 1, write_im = repmat(write_im,[C  1]); end
if size(write_tc,1) == 1 && K1 > 1, write_tc = repmat(write_tc,[K1 1]); end

if ~(all(write_bf(:) == false) && all(write_im(:) == false) && all(write_tc(:) == false) && all(write_df(:) == false) && all(write_sm == true))
    psi0 = spm_mb_io('GetData',datn.psi);
end

if write_vel
    % Write initial velocity
    descrip = 'Velocity';
    nam     = ['v_' namn '.nii'];
    fpth    = fullfile(dir_res,nam);    
    if ~isfile(fpth)
        v       = spm_mb_io('GetData',datn.v);
        WriteNii(fpth,v,Mmu,descrip);        
        clear v
    end
    resn.v  = fpth;
end

if write_aff
    % Write affine matrix
    nam  = ['Affine_' namn '.mat'];
    fpth = fullfile(dir_res,nam);
    save(fpth,'Mr');
end

if any(write_df == true) || any(reshape(write_tc(:,[2 3]),[],1) == true) || any(reshape(write_im(:,[3 4]),[],1) == true) || any(write_sm == true) ...
   || (isfield(datn,'mog') && (any(write_bf(:) == true) || any(write_im(:) == true) || any(write_tc(:) == true) || any(write_sm == true)))
    % Get forward deformation (pulls template into subject space)
    psi = spm_mb_shape('Compose',psi0,spm_mb_shape('Affine',df,Mmu\Mr*Mn));    
    clear psi0
end

if isfield(datn,'mog') && (any(write_bf(:) == true) || any(write_im(:) == true) || any(write_tc(:) == true) || any(write_sm == true))
    % Input data were intensity images
    %------------------

    % Get subject-space template (softmaxed K + 1)    
    mun = spm_mb_shape('Pull1',mun0,psi,mu_bg);    

    % Make K + 1 template
    mun = reshape(mun,[prod(df(1:3)) K]);
    mun = spm_mb_shape('TemplateK1',mun,2);

    if any(do_bf == true)
        % Get bias field
        chan = spm_mb_appearance('BiasBasis',datn.T,df,datn.Mat,reg);
        bf   = spm_mb_appearance('BiasField',datn.T,chan);
    else
        bf   = ones([1 C],'single');
    end

    % Get image(s)
    fn = spm_mb_io('GetImage',datn,false);

    % Get labels
    labels = spm_mb_appearance('GetLabels',datn,sett);
    mun    = mun + labels;
    clear labels

    % Integrate use of multiple Gaussians per tissue
    mg_w = datn.mog.mg_w;
    mun  = mun(:,mg_ix);
    mun  = mun + log(mg_w);

    % Format for spm_gmm
    [bffn,code_image,msk_chn] = spm_gmm_lib('obs2cell', reshape(bf.*fn,[prod(df) C]));
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
        sample_post = do_infer > 1;
        MU          = datn.mog.po.m;
        A           = bsxfun(@times, datn.mog.po.W, reshape(datn.mog.po.n, [1 1 Kmg]));
        fn          = spm_gmm_lib('InferMissing',reshape(fn,[prod(df) C]),zn,{MU,A},code_image,sample_post);
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
        if has_spine
            zn = CleanGWCSpine(zn,gwc_tix);
        else
            zn = CleanGWC(zn,gwc_tix,gwc_level);
        end        
    end

    if any(write_bf == true) && any(do_bf == true)
        % Write bias field
        descrip = 'Bias field (';
        pths    = {};
        for c=1:C
            if ~write_bf(c,1), continue; end
            nam  = ['bf' num2str(c) '_' namn '.nii'];
            fpth = fullfile(dir_res,nam);
            WriteNii(fpth,bf(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);
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
            WriteNii(fpth,fn(:,:,:,c)./bf(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')'],'int16');
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
            WriteNii(fpth,fn(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')'],'int16');
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
            WriteNii(fpth,zn(:,:,:,k),Mn,[descrip 'k=' num2str(k) ')'],'uint8');
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
    
if write_df(1)
    % Write forward deformation (map to mm)             
    descrip = 'Forward deformation';
    nam     = ['y_' namn '.nii'];
    fpth    = fullfile(dir_res,nam);
    WriteNii(fpth,...
             reshape(reshape(psi,[prod(df) 3])*Mmu(1:3,1:3)' + Mmu(1:3,4)',[df 1 3]),...
             Mn,descrip);
    resn.y  = fpth;
end

if any(isfinite(bb(:))) || any(isfinite(vx_out))
    % If a bounding box is supplied, combine this with the closest
    % bounding box derived from the dimensions and orientations of
    % the tissue priors. (from spm_preproc8_write)
    [bb1,vx1]                    = spm_get_bbox(pth_mu, 'old');
    bb(~isfinite(bb))            = bb1(~isfinite(bb));
    if ~isfinite(vx_out), vx_out = abs(prod(vx1))^(1/3); end
    bb(1,:)                      = vx_out*round(bb(1,:)/vx_out);
    bb(2,:)                      = vx_out*round(bb(2,:)/vx_out);
    % Bounding box dimensions
    dbb = abs(round((bb(2,1:3)-bb(1,1:3))/vx_out)) + 1;
    % Bounding box orientation matrix
    mm  = [[bb(1,1) bb(1,2) bb(1,3)
            bb(2,1) bb(1,2) bb(1,3)
            bb(1,1) bb(2,2) bb(1,3)
            bb(2,1) bb(2,2) bb(1,3)
            bb(1,1) bb(1,2) bb(2,3)
            bb(2,1) bb(1,2) bb(2,3)
            bb(1,1) bb(2,2) bb(2,3)
            bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];
    vx3 = [[1      1      1
            dbb(1) 1      1
            1      dbb(2) 1
            dbb(1) dbb(2) 1
            1      1      dbb(3)
            dbb(1) 1      dbb(3)
            1      dbb(2) dbb(3)
            dbb(1) dbb(2) dbb(3)]'; ones(1,8)];
    Mbb = mm/vx3;
    % Apply bounding box to deformations
    M = Mbb\Mmu;
    psi = reshape(reshape(psi,[prod(df) 3])*M(1:3,1:3)' + M(1:3,4)',[df 3]);    
    % New normalised dimensions and orientation matrix
    dmu = dbb;
    Mmu = Mbb;
end

if any(write_df == true) || any(reshape(write_tc(:,[2 3]),[],1) == true) || any(reshape(write_im(:,[3 4]),[],1) == true) || any(write_sm == true)
    % Write forward deformation and/or normalised images
    %------------------

    % For imporved push - subsampling density in each dimension
    sd = spm_mb_shape('SampDens',Mmu,Mn);
    
    if isfield(datn,'mog') && any(write_im(:,3) == true)
        % Write normalised image
        descrip = 'Normalised image (';
        pths    = {};
        for c=1:C
            if ~write_im(c,3), continue; end
            nam     = ['wim' num2str(c) '_' namn '.nii'];
            fpth    = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c)./bf(:,:,:,c),psi,dmu,sd);
            WriteNii(fpth,img./(cnt + eps('single')),Mmu,[descrip 'c=' num2str(c) ')'],'int16');
            pths{end + 1} = fpth;
        end
        clear img cnt
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
            WriteNii(fpth,img./(cnt + eps('single')),Mmu,[descrip 'c=' num2str(c) ')'],'int16');
            pths{end + 1} = fpth;
        end
        clear img cnt
        resn.wimc = pths;
    end

    if any(write_tc(:,2) == true) || any(write_tc(:,3) == true)
        % Write normalised segmentations
        descrip_wc  = 'Normalised tissue (';
        pths_wc     = {};
        descrip_mwc = 'Normalised modulated tissue (';
        pths_mwc    = {};
        for k=1:K1
            if write_tc(k,2) || write_tc(k,3)
                % Push
                [img,cnt] = spm_mb_shape('Push1',zn(:,:,:,k),psi,dmu,sd);
            end            
            if write_tc(k,2)
                % Normalised
                nam              = ['wc' num2str(k) '_' namn '.nii'];
                fpth             = fullfile(dir_res,nam);                
                WriteNii(fpth,img./(cnt + eps('single')),Mmu,[descrip_wc 'k=' num2str(k) ')'],'uint8');
                pths_wc{end + 1} = fpth;
            end
            if write_tc(k,3)
                % Modulated normalised
                nam               = ['mwc' num2str(k) '_' namn '.nii'];
                fpth              = fullfile(dir_res,nam);                
                img               = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
                WriteNii(fpth,img,Mmu,[descrip_mwc 'k=' num2str(k) ')'],'int16');
                pths_mwc{end + 1} = fpth;
            end
        end
        clear img cnt
        resn.wc = pths_wc;
        resn.mwc = pths_mwc;
    end

    if any(write_sm == true)
        % Compute scalar momentum (SM)     
        % See appendix of:        
        % Monte-Rubio, Gemma C., et al. "A comparison of various MRI
        % feature types for characterizing whole brain anatomical
        % differences using linear pattern recognition methods." NeuroImage
        % 178 (2018): 753-768.
        
        Ka     = numel(write_sm);             % Number of classes to compute SM for
        fwhm_a = 10;                          % FWHM for smoothing of SM
        vx     = sqrt(sum(Mmu(1:3,1:3).^2));  % template voxel size
        if Ka == 1
            % Include all tissue classes
            write_sm = 1:K1;
            Ka       = numel(write_sm);             
        end
        
        A    = zeros([dmu(1:3) Ka],'single');      % Stores linear combination of momentum
        mun0 = spm_mb_shape('TemplateK1',mun0,4);  % All K + 1 classes of template are needed here
        for k=1:Ka  % Loop over classes
            [img,cnt]  = spm_mb_shape('Push1',zn(:,:,:,k),psi,dmu,sd);  % Warp z(k) to template space
            msk        = ~isfinite(img);                                % Identify missing data
            a          = cnt.*mun0(:,:,:,k) - img;                      % Compute SM
            a(msk)     = 0;                                             % Assume all values are zero outside FOV
            spm_smooth(a,a,fwhm_a./vx);                                 % Smooth SM
            A(:,:,:,k) = a;
            clear a
        end        
        clear img cnt
        
        % Write scalar momentum
        descrip = 'Scalar momentum (';
        pths    = {};
        for k=1:Ka
            if ~write_sm(k), continue; end
            nam           = ['sm' num2str(k) '_' namn '.nii'];
            fpth          = fullfile(dir_res,nam);
            WriteNii(fpth,A(:,:,:,k),Mmu,[descrip 'k=' num2str(k) ')']);
            pths{end + 1} = fpth;
        end
        resn.sm = pths;        
        clear A
    end    
end

% Clean-up
if clean_def  && isa(datn.psi,'nifti') && isfile(datn.psi.dat.fname),          delete(datn.psi.dat.fname); end
if ~write_vel && clean_vel && isa(datn.v,'nifti') && isfile(datn.v.dat.fname), delete(datn.v.dat.fname);   end
end
%==========================================================================

%==========================================================================
% CleanGWCSpine()
function zn = CleanGWCSpine(zn,ixt)
if nargin < 2 || isempty(ixt)
    % Default SPM12 template ordering
    ixt = struct('gm',1,'wm',2,'csf',3);
end

% Get GM and WM
gw = sum(zn(:,:,:,[ixt.gm ixt.wm]),4);

% Get brain mask (inc spine)
b = gw > 0.5;
% figure(111);montage(bmsk(:,:,1:8:176))

% Find largest connected component in brain mask
C    = bwlabeln(b);
unqC = unique(C);
vls  = zeros([1 numel(unqC)]);
for i=1:numel(unqC)
    vls(i) = sum(sum(sum(C == unqC(i))));
end
[~,ix] = sort(vls);
ix     = (ix(end - 1) - 1);
b   = C == ix;
% figure(111);montage(bmsk(:,:,1:8:176))

% Close holes
b = BinMorph('close',b,6); 
% bmsk1 = bmsk1 + bmsk;
% figure(111);montage(bmsk1(:,:,1:8:176));
% figure(111);clf;imagesc(bmsk1(:,:,8*11));

c = BinMorph('dilate',b,1);

% Modify resps
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
    
    cp = double(c(:,:,i));
    cp = ((cp>th).*(sum(cat(3,slices{ixt.gm}),3)+sum(cat(3,slices{ixt.wm}),3)+sum(cat(3,slices{ixt.csf}),3)))>th;

    for i1=1:numel(ixt.csf)
        slices{ixt.csf(i1)} = slices{ixt.csf(i1)}.*cp;
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
% BinMorph()
function vol = BinMorph(action,bim,n)

if nargin < 3, n = 1; end % Iterations
    
vol = uint8(bim);

kx = [1 1 1];
ky = [1 1 1];
kz = [1 1 1];

order = sum(kx(:) ~= 0)*sum(ky(:) ~= 0);

switch lower(action)
    
	case 'dilate'
	
        sz = size(vol);
        vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
        vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
        for i = 1:n
            spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
            vol2 = uint8(vol2~=0);
            
%             imagesc3d(vol2); axis off; drawnow
        end
        vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);

	case 'erode'

        for i = 1:n
            spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
            vol = uint8(vol>=order);
            
%             imagesc3d(vol); axis off; drawnow
        end
        
	case 'close'
        
        sz = size(vol);
        vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
        vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
        for i = 1:n
            spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
            vol2 = uint8(vol2~=0);
            
%             imagesc3d(vol2); axis off; drawnow
        end                
        
        for i = 1:n - 2
            spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
            vol2 = uint8(vol2>=order);
            
%             imagesc3d(vol2); axis off; drawnow
        end
        vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);
        clear vol2                        
end

vol = logical(vol);
end
%==========================================================================

%==========================================================================
function WriteNii(f,img,M,descrip,typ)
if nargin<5, typ = 'float32'; end
switch typ
case 'float32'
    fa = file_array(f,size(img),typ,0);
case 'uint8'
    mx = max(img(isfinite(img(:))));
    fa = file_array(f,size(img),typ,0,mx/255,0);
case 'int16'
    mx = max(img(isfinite(img(:))));
    mn = min(img(isfinite(img(:))));
    s  = max(mx/32767,-mn/32768);
    fa = file_array(f,size(img),typ,0,s,0);
otherwise
    error('Can''t do datatype "%s"', typ);
end
Nii         = nifti;
Nii.dat     = fa;
Nii.mat     = M;
Nii.mat0    = M;
Nii.descrip = descrip;
create(Nii);
Nii.dat(:,:,:,:,:,:) = img;
end
%==========================================================================