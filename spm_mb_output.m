function res = spm_mb_output(cfg)
% Write output from groupwise normalisation and segmentation of images
% FORMAT res = spm_mb_output(cfg)
%
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Centre for Human Neuroimaging

res  = load(char(cfg.result));
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
             'fwhm',cfg.fwhm);

spm_progress_bar('Init',N,'Writing MB output','Subjects complete');
for n=1:N % Loop over subjects
    res(n) = ProcessSubject(dat(n),res(n),mu,sett,opt);
    spm_progress_bar('Set',n);
end
spm_progress_bar('Clear');
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% PostProcMRF()
function zn = PostProcMRF(zn,Mn,strength,nit)
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
dmu        = sett.mu.d;
Mmu        = sett.mu.Mmu;
dir_res    = sett.odir;
do_infer   = true;
mrf        = opt.mrf;
write_inu  = opt.write_inu; % field
write_im   = opt.write_im;  % image, corrected, warped, warped corrected
write_tc   = opt.write_tc;  % native, warped, warped-mod, scalar momentum
fwhm       = opt.fwhm;      % FWHM for smoothing of warped tissues
vx         = opt.vx;        % Template space voxel size
bb         = opt.bb;        % Template space bounding box

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
        chan = spm_mb_appearance('inu_basis',gmm.T,df,datn.Mat,ones(1,C));
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
    fn     = spm_mb_io('get_image',gmm);

    % Get warped tissue priors
    mun    = spm_mb_shape('pull1',mu,psi);
    mun    = spm_mb_classes('template_k1',mun,4);
    mun    = bsxfun(@plus,mun,spm_mb_classes('get_labels',datn,size(mun,4)));
    mun    = reshape(mun,size(mun,1)*size(mun,2)*size(mun,3),size(mun,4));

    % Integrate use of multiple Gaussians per tissue
    mg_w = gmm.mg_w;
    mun  = mun(:,mg_ix);
    mun  = bsxfun(@plus, mun, log(mg_w));

    % Format for spm_gmm
    chan                   = spm_mb_appearance('inu_basis',gmm.T,df,datn.Mat,ones(1,C));
    [~,mf,vf]              = spm_mb_appearance('inu_recon',fn,chan,gmm.T,gmm.Sig);
    clear fn
    mf                     = reshape(mf,[prod(df) C]);
    vf                     = reshape(vf,[prod(df) C]);
    [~,code_image,msk_chn] = spm_gmm_lib('obs2cell', reshape(mf,[prod(df) C]));    

    % Get responsibilities, making sure that missing values are 'filled in'
    % by the template. For example, for CT, CSF can have intensity zero;
    % but we consider this value as missing as background values can also be
    % zero, which would bias the fitting of the GMM.
    const             = spm_gmm_lib('Normalisation', {gmm.m,gmm.b}, {gmm.V,gmm.n}, msk_chn);
    if ~isempty(vf)
        zn            = spm_gmm_lib('Marginal', mf, {gmm.m,gmm.V,gmm.n}, const, msk_chn, vf);
    else
        zn            = spm_gmm_lib('Marginal', mf, {gmm.m,gmm.V,gmm.n}, const, msk_chn);
    end
    zn(~isfinite(zn)) = min(zn(:));  % NaN assumed to have small (log) probability
    zn                = spm_gmm_lib('Responsibility', zn, mun);    
    clear mun msk_chn vf

    % Get bias field modulated image data
    if do_infer
        % Infer missing values
        sample_post = do_infer > 1;
        A           = bsxfun(@times, gmm.V, reshape(gmm.n, [1 1 Kmg]));
        mf          = spm_gmm_lib('InferMissing',reshape(mf,[prod(df) C]),...
                                  zn,{gmm.m,A},code_image,sample_post);
        clear code
    end
    clear code_image

    mf = reshape(mf,[df(1:3) C]);

    if any(write_im(:,1))
        % Write image
        resn.im = cell(1,sum(write_im(:,1)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,1), continue; end
            nam  = sprintf('i%d_%s.nii',c,onam);
            fpth = fullfile(dir_res,nam);
            write_nii(fpth,mf(:,:,:,c)./inu(:,:,:,c), Mn, sprintf('Image (%d)',c), 'int16');
            c1          = c1 + 1;
            resn.m{c1} = fpth;
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
            c1           = c1 + 1;
            resn.mi{c1} = fpth;
        end
    end

    % For improved push - subsampling density in each dimension
    sd = spm_mb_shape('samp_dens',Mmu,Mn);

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
    clear mf

    % If using multiple Gaussians per tissue, collapse so that zn is of size K1
    if Kmg > K1
        for k=1:K1
            zn(:,k) = sum(zn(:,mg_ix==k),2);
        end
        zn(:,K1 + 1:end)    = [];
    end
    zn = reshape(zn,[df(1:3) K1]);

    if mrf > 0
        % Ad-hoc MRF clean-up of segmentation
        zn = PostProcMRF(zn,Mn,mrf);
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
        mun = spm_mb_shape('softmax',mun,4);
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
function [Mmu,dmu,vx_mu,psi] = modify_fov(bb_out,vx_out,Mmu,dmu,vx_mu,psi,sett)
if any(isfinite(bb_out(:))) || any(isfinite(vx_out))
    % Get bounding-box   
    bb1 = spm_get_bbox(sett.mu.exist.mu, 'old');
    bb_out(~isfinite(bb_out)) = bb1(~isfinite(bb_out));
    if ~isfinite(vx_out), vx_out = abs(prod(vx_mu))^(1/3); end
    bb_out(1,:) = vx_out*round(bb_out(1,:)/vx_out);
    bb_out(2,:) = vx_out*round(bb_out(2,:)/vx_out);
    % New output dimensions
    dmu = abs(round((bb_out(2,1:3)-bb_out(1,1:3))/vx_out))+1;
    % New orientation matrix
    mm  = [[bb_out(1,1) bb_out(1,2) bb_out(1,3)
            bb_out(2,1) bb_out(1,2) bb_out(1,3)
            bb_out(1,1) bb_out(2,2) bb_out(1,3)
            bb_out(2,1) bb_out(2,2) bb_out(1,3)
            bb_out(1,1) bb_out(1,2) bb_out(2,3)
            bb_out(2,1) bb_out(1,2) bb_out(2,3)
            bb_out(1,1) bb_out(2,2) bb_out(2,3)
            bb_out(2,1) bb_out(2,2) bb_out(2,3)]'; ones(1,8)];
    vx3 = [[1       1       1
            dmu(1) 1       1
            1       dmu(2) 1
            dmu(1) dmu(2) 1
            1       1       dmu(3)
            dmu(1) 1       dmu(3)
            1       dmu(2) dmu(3)
            dmu(1) dmu(2) dmu(3)]'; ones(1,8)];
    M1    = mm/vx3;
    M0    = M1\Mmu;
    vx_mu = sqrt(sum(Mmu(1:3,1:3).^2));
    Mmu   = M1;
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