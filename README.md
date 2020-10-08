# Multi-Brain: Unified segmentation of population neuroimaging data

## Overview
This repository contains the Multi-Brain (MB) model, which has the general aim of integrating a number of disparate image analysis components within a single unified generative modelling framework (segmentation, nonlinear registration, image translation, etc.). The model is described in Brudfors et al [2020], and it builds on a number of previous works. Its objective is to achieve diffeomorphic alignment of a wide variaty of medical image modalities into a common anatomical space. This involves the ability to construct a "tissue probability template" from a population of scans through group-wise alignment [Ashburner & Friston, 2009; Blaiotta et al, 2018]. Diffeomorphic deformations are computed within a *geodesic shooting* framework [Ashburner & Friston, 2011], which is optimised with a Gauss-Newton strategy that uses a multi-grid approach to solve the system of linear equations [Ashburner, 2007]. Variability among image contrasts is modelled using a much more sophisticated version of the Gaussian mixture model with bias correction framework originally proposed by Ashburner & Friston [2005], and which has been extended to account for known variability of the intensity distributions of different tissues [Blaiotta et al, 2018]. This model has been shown to provide a good model of the intensity distributions of different imaging modalities [Brudfors et al, 2019]. Time permitting, additional registration accuracy through the use of shape variability priors [Balbastre et al, 2018] will also be incorporated.

## Dependencies
The algorithm is developed using MATLAB and relies on external functionality from the SPM12 software:
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.
* **Shoot toolbox:** Add Shoot folder from the toolbox directory of the SPM source code.
* **Longitudinal toolbox:** Add Longitudinal folder from the toolbox directory of the SPM source code.

## Example use cases
This section contains example code demonstrating how the MB toolbox can be used for nonlinear image registration, spatial normalisation, tissue segmentation and bias-field correction. **Example 1** learns (or fits) the MB model from a population of MRI scans. Learning the MB model results in a bunch of tissue segmentations, bias-field corrected scans and forward deformations, aligning a template to each subject's scan (a combined nonlinear and rigid registration). The deformations can be used to warp a subject image to template space, or aligning two subjects' images together by composing their deformations. These two operations are demonstrated in **Example 2**. Finally, **Example 3** fits an already learned MB model to two subject MRIs. This is equivalent to registering to a pre-existing common space. Warping can then be done as shown in Example 2.

### 1. Learning the MB model
```
% This script fits the MB model to a population of MRIs. It is a
% group-wise image registration where an optimal K class tissue template is 
% built, as well as optimal intensity parameters learnt. The resulting
% deformations (written to dir_out prefixed y_) can then be used to warp
% between subject scans. Also native and template space tissue
% segmentaitons are written to disk, as well as a bias-field corrected 
% versions of the input scans.

% Data directory, assumed to contain MRI nifti images of the same contrast.
dir_data = '/directory/with/some/MRIs';

% MB settings
K        = 9;           % Template classes (there are K + 1 Gaussians)
vx       = 1.5;         % Template voxel size (smaller = slower)
dir_out  = 'mb-output'; % Directory where to write results

% Get paths to data
pth_im  = spm_select('FPList',dir_data,'^.*\_3.nii$');

% Create MB cfg struct. There are a bunch of settings here, they should all
% work well by default.
% -----------------
% spm_mb_init
run = struct;
run.mu.create.K = K;
run.mu.create.vx = vx;
run.mu.create.mu_settings = [1e-05 0.5 0];
run.aff = 'SE(3)';
run.v_settings = [0.0001 0 0.4 0.1 0.4];
run.onam = 'mb';
run.odir = {dir_out};
run.cat = {{}};
run.accel = 0.8;
run.min_dim = 16;
run.tol = 0.001;
run.sampdens = 2;
run.save = true;
run.nworker = Inf;
run.gmm.chan.images = cellstr(pth_im);
run.gmm.chan.inu.inu_reg = 1e4;
run.gmm.chan.inu.inu_co = 40;
run.gmm.chan.modality = 1;
run.gmm.labels.false = [];
run.gmm.pr.file = {};
run.gmm.pr.hyperpriors = {};
run.gmm.mg_ix = 1:K + 1;
run.gmm.tol_gmm = 0.0005;
run.gmm.nit_gmm_miss = 32;
run.gmm.nit_gmm = 8;
run.gmm.nit_appear = 4;
% spm_mb_output
out = struct;
out.i = false;
out.mi = true; % writes bias-field corrected version
out.wi = false;
out.wmi = false;
out.inu = false;
out.mrf = 0;
out.c = true(1,K + 1); % writes tissue classes in native space
out.wc = true(1,K + 1); % writes tissue classes in template space
out.mwc = true(1,K + 1); % writes tissue classes in modulated template space

% Init MB
[dat,sett] = spm_mb_init(run);

if ~isempty(dat)
    % Fit MB
    fprintf('Fitting MB to %i subjects...\n', numel(dat))
    [dat,sett] = spm_mb_fit(dat,sett);
    
    % Save results
    pth_res    = fullfile(sett.odir,['mb_fit_' sett.onam '.mat']);
    save(pth_res,'dat','sett');
    out.result = pth_res;

    % Write output
    spm_mb_output(out);   
end
```

### 2. Warping with MB deformations
```
% This script demonstrates warping with the deformation generated from
% fitting the MB model. The template file generated by MB, plus two subject 
% scans with their corresponding forward deformations (prefixed 'y_') are
% needed. The following three warps are demonstrated:
% 1. Warp an image to template space.
% 2. Warp the template to image space.
% 3. Warp one image to another image.

% Path to a MB tissue template (this is dir_out in Example 1)
pth_mu = fullfile(dir_out,'mu_mb.nii');
% Paths to two subject MRIs to which MB has been fitted
pth_img1 = '/directory/with/some/MRIs/img1.nii';
pth_img2 = '/directory/with/some/MRIs/img2.nii';
% Paths to corresponding forward deformations of the above subjects
pth_y1 = fullfile(dir_out,'/y_img1.nii');
pth_y2 = fullfile(dir_out,'y_img2.nii');
% Where to write the warped images
dir_out = 'warped';

% 1. Use forward deformation (pth_y1) to warp image 1 (pth_img1) to 
%    template space (pth_mu)
% -----------------
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def     = {pth_y1};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space           = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_img1};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';
spm_jobman('run',matlabbatch);

% 2. Use forward deformation (pth_y1) to template (pth_mu) to  image 
%    space (pth_img1) 
% -----------------
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.comp{1}.def         = {pth_y1};
matlabbatch{1}.spm.util.defs.comp{1}.space               = {pth_img1};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';
spm_jobman('run',matlabbatch);

% 3. Compose both forward deformations (pth_y1, pth_y2) via template space
%    to register image 2 (pth_img2) to image 1 (pth_img1).
% -----------------
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def     = {pth_y2};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space           = {pth_mu};
matlabbatch{1}.spm.util.defs.comp{2}.def                 = {pth_y1};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_img2};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';
spm_jobman('run',matlabbatch);
```

### 3. Fitting a learned MB model
```
% This script demonstrates how to fit a learned MB model (template and 
% intensity prior) to new data. In practice, it only involves modifying
% two entires in Example 1.

% Path to results MAT-file generated by fitting MB
pth_fit = fullfile(dir_out,'mb_fit_mb.mat');

% Make intensity prior file by extracting information from pth_fit
load(pth_fit,'sett'); 
pr      = sett.gmm.pr; 
mg_ix   = sett.gmm.mg_ix; 
save('prior_mb.mat','pr','mg_ix');

% In Example 1, to use a learned template simply replace all run.mu entires with:
run.mu.exist    = {pth_mu};
% and to include the intensity prior, simply point to it at the following run entry:
run.gmm.pr.file = {pth_int_prior};
```

TODO

## References
* Ashburner J, Friston KJ. Unified segmentation. Neuroimage. 2005 Jul 1;26(3):839-51.
* Ashburner J. A fast diffeomorphic image registration algorithm. Neuroimage. 2007 Oct 15;38(1):95-113.
* Ashburner J, Friston KJ. Computing average shaped tissue probability templates. Neuroimage. 2009 Apr 1;45(2):333-41.
* Ashburner J, Friston KJ. Diffeomorphic registration using geodesic shooting and Gauss–Newton optimisation. NeuroImage. 2011 Apr 1;55(3):954-67.
* Blaiotta C, Cardoso MJ, Ashburner J. Variational inference for medical image segmentation. Computer Vision and Image Understanding. 2016 Oct 1;151:14-28.
* Blaiotta C, Freund P, Cardoso MJ, Ashburner J. Generative diffeomorphic modelling of large MRI data sets for probabilistic template construction. NeuroImage. 2018 Feb 1;166:117-34.
* Balbastre Y, Brudfors M, Bronik K, Ashburner J. Diffeomorphic brain shape modelling using Gauss-Newton optimisation. In International Conference on Medical Image Computing and Computer-Assisted Intervention 2018 Sep 16 (pp. 862-870). Springer, Cham.
* Brudfors M, Ashburner J, Nachev P, Balbastre Y. Empirical Bayesian Mixture Models for Medical Image Translation. In International Workshop on Simulation and Synthesis in Medical Imaging 2019 Oct 13 (pp. 1-12). Springer, Cham.
* Brudfors M, Balbastre Y, Flandin G, Nachev P, Ashburner J. Flexible Bayesian Modelling for Nonlinear Image Registration. To appear in: International Conference on Medical Image Computing and Computer-Assisted Intervention 2020

## Acknowledgements
This work was funded by the EU Human Brain Project’s Grant Agreement No 785907 (SGA2).
