function res = spm_mb_output(cfg)
% Write output from groupwise normalisation and segmentation of images
% FORMAT res = spm_mb_output(cfg)
%
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Centre for Human Neuroimaging

% $Id: spm_mb_output.m 8164 2021-10-15 16:05:27Z john $

res  = cfg.result;
if iscell(res), res = res{1}; end
if ischar(res), res  = load(res); end
sett = res.sett;
dat  = res.dat;

if isfield(sett.mu,'exist')
    mu = sett.mu.exist.mu;
elseif isfield(sett.mu,'create')
    mu = sett.mu.create.mu;
end
mu = nifti(mu);
mu = single(mu.dat(:,:,:,:,:));

% If SPM has been compiled with OpenMP support then the number of threads
% are here set to speed up the algorithm
%--------------------------------------------------------------------------
nw   = spm_mb_shape('get_num_workers',sett,max(27,sett.K*5+17));
if sett.nworker > 1
    setenv('SPM_NUM_THREADS',sprintf('%d',0));
else
    setenv('SPM_NUM_THREADS',sprintf('%d',-1));
end

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('inu',cl,'i',cl,'mi',cl,'c',cl,'wi',cl, ...
             'wmi',cl,'wc',cl,'mwc',cl,'sm',cl);

write_tc = false(sett.K+1,4);
ind = cfg.c;   ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,1) = true;
ind = cfg.wc;  ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,2) = true;
ind = cfg.mwc; ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,3) = true;
ind = cfg.sm;  ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,4) = true;

opt = struct('write_inu',cfg.inu,...
             'write_im',[cfg.i cfg.mi cfg.wi cfg.wmi],...
             'write_tc',write_tc,...
             'mrf',cfg.mrf,...
             'vx',cfg.vox,...
             'bb',cfg.bb,...
             'odir',cfg.odir,...
             'fwhm',cfg.fwhm);
opt.clean_gwc = cfg.clean_gwc;

if nw > 1 && numel(dat) > 1 % PARFOR
    fprintf('Write output: ');
    parfor(n=1:N,nw)
        fprintf('.');
        res(n) = ProcessSubject(dat(n),res(n),mu,sett,opt);
    end
    fprintf(' done!\n');
else
    spm_progress_bar('Init',N,'Writing MB output','Subjects complete');
    for n=1:N % FOR
        res(n) = ProcessSubject(dat(n),res(n),mu,sett,opt);
        spm_progress_bar('Set',n);
    end
    spm_progress_bar('Clear');
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% mrf_apply()
function zn = mrf_apply(zn,Mn,strength,nit)
if nargin < 4, nit = 10; end
P   = zeros(size(zn),'uint8');
G   = ones([size(zn,4),1],'single')*strength;
vx  = sqrt(sum(Mn(1:3,1:3).^2));
vx2 = 1./single(vx);
for i=1:nit
    spm_mrf(P,zn,G,vx2);
end
zn = single(P)/255;
%==========================================================================

%==========================================================================
% ProcessSubject()
function resn = ProcessSubject(datn,resn,mu,sett,opt)

% Parse function settings
dmu         = sett.ms.d;
Mmu         = sett.ms.Mmu;
if isempty(opt.odir)
    dir_res = sett.odir;
else
    dir_res = opt.odir;
    if ~(exist(dir_res,'dir') == 7)
        mkdir(dir_res)
    end
end
do_infer    = true;
mrf         = opt.mrf;
write_inu   = opt.write_inu; % field
write_im    = opt.write_im;  % image, corrected, warped, warped corrected
write_tc    = opt.write_tc;  % native, warped, warped-mod, scalar momentum
fwhm        = opt.fwhm;      % FWHM for smoothing of warped tissues
vx          = opt.vx;        % Template space voxel size
bb          = opt.bb;        % Template space bounding box
clean_gwc   = opt.clean_gwc; % Settings for cleaning up tissue classes

cl   = cell(1,1);
resn = struct('inu',cl,'i',cl,'mi',cl,'c',cl,'wi',cl, ...
              'wmi',cl,'wc',cl,'mwc',cl,'sm',cl);

if ((~any(write_inu(:)) && ~any(write_im(:))) || ~isfield(datn.model,'gmm')) && ~any(write_tc(:))
    return;
end

% Get parameters
df     = datn.dm;
onam   = datn.onam;
Mn     = datn.Mat;
do_inu = true;

if isfield(datn.model,'gmm')

    % Input data were intensity images
    %------------------

    gmm    = datn.model.gmm;
    gmms   = sett.gmm(gmm.pop);
    C      = gmms.C;
    mg_ix  = gmms.mg_ix;
    K      = sett.K;
    K1     = K + 1;
    Kmg    = numel(mg_ix);

    % Integrate K1 and C into write settings
    if size(write_inu,1) == 1
        write_inu = repmat(write_inu,[C  1]);
    end
    if size(write_im, 1) == 1
        write_im  = repmat(write_im, [C  1]);
    end
    if size(write_tc, 1) == 1
        write_tc  = repmat(write_tc, [K1 1]);
    end

    if any(do_inu == true)
        % Get bias field
        chan = spm_mb_appearance('inu_basis',gmm.T,df,datn.Mat);
        inu  = spm_mb_appearance('inu_field',gmm.T,chan);
        clear chan

        if any(write_inu == true) && any(do_inu == true)
            % Write bias field
            inu      = reshape(inu,[df(1:3) C]);
            resn.inu = cell(1,sum(write_inu));
            c1       = 0;
            for c=1:C
                if ~write_inu(c,1), continue; end
                c1   = c1 + 1;
                nam  = sprintf('inu%d_%s.nii',c,onam);
                fpth = fullfile(dir_res,nam);
                write_nii(fpth,inu(:,:,:,c), Mn, sprintf('INU (%d)',c));
                c1           = c1 + 1;
                resn.inu{c1} = fpth;
            end
        end
    else
        inu   = ones([1 C],'single');
    end
end

if any(write_im(:)) || any(write_tc(:))
    Mat = Mmu\spm_dexpm(double(datn.q),sett.B)*datn.Mat;
    psi = spm_mb_io('get_data',datn.psi);
    psi = MatDefMul(psi,inv(Mmu));
    psi = spm_mb_shape('compose',psi,spm_mb_shape('affine',datn.dm,Mat));
end


if isfield(datn.model,'gmm') && (any(write_im(:)) || any(write_tc(:)))

    % Get image(s)
    [fn,msk] = spm_mb_io('get_image',gmm,false,do_infer,false);

    % Get warped tissue priors
    mun    = spm_mb_shape('pull1',mu,psi);
    mun    = spm_mb_classes('template_k1',mun,datn.delta,4);

    % Integrate use of multiple Gaussians per tissue
    gam  = gmm.gam;

    % Format for spm_gmm
    chan                     = spm_mb_appearance('inu_basis',gmm.T,df,datn.Mat,ones(1,C));
    [~,mf,vf]                = spm_mb_appearance('inu_recon',fn,[],chan,gmm.T,gmm.Sig);

    % Get label data
    if isa(datn.lab,'struct')
        label = uint8(spm_mb_io('get_data', datn.lab.f));
    else
        label = [];
    end
    lnP = dirichlet_logexpect(gmm.Alpha);

    % For improved push - subsampling density in each dimension
    sd  = spm_mb_shape('samp_dens',Mmu,Mn);

    zn  = spm_gmmlib('resp',gmm.m,gmm.b,gmm.W,gmm.nu,gmm.gam,...
                     uint64(mg_ix), mun,mf,vf, uint64(gmm.samp), label,lnP);
    zn  = cat(4,zn,1-sum(zn,4));

    if mrf > 0
        % Ad-hoc MRF clean-up of segmentation
        zn = mrf_apply(zn,Mn,mrf);
    end

    if clean_gwc.do == true
        % Ad-hoc clean-up of GM, WM and CSF
        zn = do_clean_gwc(zn, clean_gwc.gm, clean_gwc.wm, clean_gwc.csf, clean_gwc.level);
    end

    if any(write_tc(:,1) == true)
        % Write segmentations
        resn.c  = cell(1,sum(write_tc(:,1)));
        k1      = 0;
        for k=1:K1
            if ~write_tc(k,1), continue; end
            nam  = sprintf('c%.2d_%s.nii',k,onam);
            fpth = fullfile(dir_res,nam);
            write_nii(fpth,zn(:,:,:,k), Mn, sprintf('Tissue (%d)',k), 'uint8');
            k1         = k1 + 1;
            resn.c{k1} = fpth;
        end
    end


    if do_infer
        % Infer missing values
        mf(msk) = NaN;
        mf      = spm_gmmlib('infer', gmm.m,gmm.b,gmm.W,gmm.nu,gmm.gam, uint64(mg_ix), ...
                             mun,mf,vf, uint64([1 1 1]), label,lnP);
        clear msk
    end    
    
    if any(write_im(:,1))
        % Write image
        resn.i = cell(1,sum(write_im(:,1)));
        c1     = 0;
        for c=1:C
            if ~write_im(c,1), continue; end
            nam  = sprintf('i%d_%s.nii',c,onam);
            fpth = fullfile(dir_res,nam);
            write_nii(fpth,mf(:,:,:,c)./inu(:,:,:,c), Mn, sprintf('Image (%d)',c), 'int16');
            c1         = c1 + 1;
            resn.i{c1} = fpth;
        end
    end

    if any(write_im(:,2))
        % Write image corrected
        resn.mi = cell(1,sum(write_im(:,2)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,2), continue; end
            nam  = sprintf('mi%d_%s.nii',c,onam);
            fpth = fullfile(dir_res,nam);
            write_nii(fpth, mf(:,:,:,c), Mn, sprintf('INU corr. (%d)',c), 'int16');
            c1          = c1 + 1;
            resn.mi{c1} = fpth;
        end
    end

    if any(write_im(:,3))
        % Write normalised image
        resn.wi = cell(1,sum(write_im(:,3)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,3), continue; end
            nam       = sprintf('wi%d_%s.nii',c,onam);
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('push1',mf(:,:,:,c)./inu(:,:,:,c), psi,dmu,sd);
            write_nii(fpth,img./(cnt + eps('single')), Mmu, sprintf('Norm. (%d)',c), 'int16');
            clear img cnt
            c1           = c1 + 1;
            resn.wi{c1} = fpth;
        end
    end
    clear inu

    if any(write_im(:,4))
        % Write normalised image corrected
        resn.wmi = cell(1,sum(write_im(:,4)));
        c1       = 0;
        for c=1:C
            if ~write_im(c,4), continue; end
            nam       = sprintf('wmi%d_%s.nii',c,onam);
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('push1',mf(:,:,:,c),psi,dmu,sd);
            write_nii(fpth,img./(cnt + eps('single')), Mmu, sprintf('Norm. INU corr. (%d)',c),'int16');
            clear img cnt
            c1           = c1 + 1;
            resn.wmi{c1} = fpth;
        end
    end
    clear mf vf
end

if isfield(datn.model,'cat') && (any(write_tc(:,2)) || any(write_tc(:,3)))
    % Input data were segmentations
    %------------------
    zn = spm_mb_io('get_data',datn.model.cat.f);
    zn = cat(4,zn,1 - sum(zn,4));
    K1 = sett.K+1;
end

% For improved push - subsampling density in each dimension
sd    = spm_mb_shape('samp_dens',Mmu,Mn);
vx_mu = sqrt(sum(Mmu(1:3,1:3).^2));

if any(write_tc(:,2)) || any(write_tc(:,3)) || any(write_tc(:,4))
    if any(write_tc(:,2)), resn.wc  = cell(1,sum(write_tc(:,2))); end
    if any(write_tc(:,3)), resn.mwc = cell(1,sum(write_tc(:,3))); end
    if any(write_tc(:,4)), resn.sm  = cell(1,sum(write_tc(:,4))); end
    kwc  = 0;
    kmwc = 0;
    ksm  = 0;
    if any(write_tc(:,4))
        % The scalar momentum residuals are computed in native space and then
        % pushed. We therefore here compute the softmaxed K + 1 template in
        % native space.
        mun = spm_mb_shape('pull1',mu,psi);
        clear mu
        mun = spm_mb_shape('softmax0',mun,4);
        mun = cat(4,mun,max(1 - sum(mun,4),0));
    end

    % Possibly modify template space images' voxel size and FOV
    [Mmu,dmu,vx_mu,psi] = modify_fov(bb,vx,Mmu,dmu,vx_mu,psi,sett);

    % Write output
    for k=1:K1
        if write_tc(k,2) || write_tc(k,3) || write_tc(k,4)
            if write_tc(k,2) || write_tc(k,3)
                [img,cnt] = spm_mb_shape('push1',zn(:,:,:,k),psi,dmu,sd);
            end
            if write_tc(k,2)
                % Write normalised segmentation
                kwc          = kwc + 1;
                fpth         = fullfile(dir_res, sprintf('wc%.2d_%s.nii',k,onam));
                resn.wc{kwc} = fpth;
                wimg         = img./(cnt + eps('single'));
                spm_smooth(wimg,wimg,fwhm(1)./vx_mu);  % Smooth
                write_nii(fpth, wimg, Mmu, sprintf('Norm. tissue (%d)',k), 'uint8');
            end
            if write_tc(k,3)
                % Write normalised modulated segmentation
                kmwc           = kmwc + 1;
                fpth           = fullfile(dir_res,sprintf('mwc%.2d_%s.nii',k,onam));
                resn.mwc{kmwc} = fpth;
                wimg           = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
                spm_smooth(wimg,wimg,fwhm(2)./vx_mu);  % Smooth
                write_nii(fpth,wimg, Mmu, sprintf('Norm. mod. tissue (%d)',k), 'int16');
            end
            clear img cnt
            if write_tc(k,4)
                % Write scalar momentum, reference:
                % "A comparison of various MRI feature types for characterizing
                %  whole brain anatomical differences using linear pattern
                %  recognition methods." Monte-Rubio, et al. NeuroImage (2018)
                ksm          = ksm + 1;
                fpth         = fullfile(dir_res,sprintf('sm%.2d_%s.nii',k,onam));
                resn.sm{ksm} = fpth;
                % Compute scalar momentum
                wimg                   = spm_mb_shape('push1',zn(:,:,:,k) - mun(:,:,:,k),psi,dmu,sd);
                wimg(~isfinite(wimg))  = 0;           % Assume all values are zero outside FOV
                spm_smooth(wimg,wimg,fwhm(3)./vx_mu); % Smooth
                write_nii(fpth,wimg, Mmu, sprintf('Scalar momentum (%d)',k), 'int16');
            end
            clear wimg
        end
    end
end
%==========================================================================

%==========================================================================
function lnP = dirichlet_logexpect(Alpha)
% Expectation of parameter lobs from Dirichlet distribution
% Note: this is a separate function because there is a variable psi
lnP = bsxfun(@minus, psi(Alpha), psi(sum(Alpha,1)));
if ~isempty(Alpha), lnP(1,:) = 0; end
%==========================================================================

%==========================================================================
function [Mmu,dmu,vx_mu,psi] = modify_fov(bb_out,vx_out,Mmu,dmu,vx_mu,psi,sett)
if any(isfinite(bb_out(:))) || any(isfinite(vx_out))
    % Get bounding-box
    if isfield(sett.mu,'exist')
        [bb0,vox0] = spm_get_bbox(sett.mu.exist.mu,  'old');
    else
        [bb0,vox0] = spm_get_bbox(sett.mu.create.mu, 'old');
    end
    vx_out     = vx_out(1)*ones(1,3);
    msk    = ~isfinite(vx_out); vx_out(msk) = vox0(msk);
    msk    = ~isfinite(bb_out); bb_out(msk) =  bb0(msk);
    bb_out = sort(bb_out);
    vx_out = abs(vx_out);
    % Adjust bounding box slightly - so it rounds to closest voxel.
    bb_out(:,1) = round(bb_out(:,1)/vx_out(1))*vx_out(1);
    bb_out(:,2) = round(bb_out(:,2)/vx_out(2))*vx_out(2);
    bb_out(:,3) = round(bb_out(:,3)/vx_out(3))*vx_out(3);
    dim = round(diff(bb_out)./vx_out+1);
    of  = -vx_out.*(round(-bb_out(1,:)./vx_out)+1);
    mat = [vx_out(1) 0 0 of(1) ; 0 vx_out(2) 0 of(2) ; 0 0 vx_out(3) of(3) ; 0 0 0 1];
    if det(Mmu(1:3,1:3)) < 0
        mat = mat*[-1 0 0 dim(1)+1; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    end
    M0    = mat\Mmu;
    Mmu   = mat;
    vx_mu = sqrt(sum(Mmu(1:3,1:3).^2));
    dmu   = dim;
    % Modify deformation
    psi = MatDefMul(psi,M0);
end
%==========================================================================

%==========================================================================
function phi = MatDefMul(phi,M)
d   = size(phi);
phi = reshape(bsxfun(@plus,reshape(phi,[prod(d(1:3)),3])*M(1:3,1:3)',M(1:3,4)'),d);
%==========================================================================

%==========================================================================
function write_nii(f,img,M,descrip,typ)
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
%==========================================================================

%==========================================================================
function zn = do_clean_gwc(zn,gm,wm,csf,level)
if nargin < 2, gm = 1; end
if nargin < 3, wm = 2; end
if nargin < 4, csf = 3; end
if nargin < 5, level = 1; end

ixt = struct('gm',gm,'wm',wm,'csf',csf);
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
%==========================================================================