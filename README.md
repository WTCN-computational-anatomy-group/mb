# Multi-Brain: Group-wise diffeomorphic registration and segmentation of medical images

## Overview
This repository contains the Multi-Brain (MB) model, which has the general aim of integrating a number of disparate image analysis components within a single unified generative modelling framework (segmentation, nonlinear registration, image translation, etc.). The model is described in Brudfors et al [2020], and it builds on a number of previous works. Its objective is to achieve diffeomorphic alignment of a wide variaty of medical image modalities into a common anatomical space. This involves the ability to construct a "tissue probability template" from a population of scans through group-wise alignment [Ashburner & Friston, 2009; Blaiotta et al, 2018]. Diffeomorphic deformations are computed within a *geodesic shooting* framework [Ashburner & Friston, 2011], which is optimised with a Gauss-Newton strategy that uses a multi-grid approach to solve the system of linear equations [Ashburner, 2007]. Variability among image contrasts is modelled using a much more sophisticated version of the Gaussian mixture model with bias correction framework originally proposed by Ashburner & Friston [2005], and which has been extended to account for known variability of the intensity distributions of different tissues [Blaiotta et al, 2018]. This model has been shown to provide a good model of the intensity distributions of different imaging modalities [Brudfors et al, 2019]. Time permitting, additional registration accuracy through the use of shape variability priors [Balbastre et al, 2018] will also be incorporated.

## Dependencies
The algorithm is developed using MATLAB and relies on external functionality from the SPM12 software. The following are therefore required downloads and need to be placed on the MATLAB search path (using `addpath`):
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/download/.
* **Shoot toolbox:** The Shoot folder from the toolbox directory of SPM12.
* **Longitudinal toolbox:** The Longitudinal folder from the toolbox directory of SPM12.

## Installation instructions
1. Download/clone this repository to the `toolbox` folder of SPM12.
2. Next, in a terminal, `cd` to this folder and execute the `make` command (this compiles required MEX-files). Note, if you are on MS Windows you should follow the instructions found here https://en.wikibooks.org/wiki/SPM/SPM12_MEX_Compilation_on_Windows (up to the compilation process), to enable the `make` command. On Windows, a compilation error has been reported where the `make` command can fail, then instead try the explicit command `mex -O -largeArrayDims -DSPM_WIN32 spm_gmmlib.c gmmlib.c`.

### Troubleshooting
* If you get an error related to the MATLAB mex path, try giving it explicitly as:
```bash
make MEXBIN=/my/matlab/path/bin/mex
```

## Example use cases
This section contains example code demonstrating how the MB toolbox can be used for nonlinear image registration, spatial normalisation, tissue segmentation and bias-field correction. **Example 1** fits the MB model to a population of images. Fitting the MB model results in: learned spatial (the template) and intensity priors, and a bunch of segmentations and forward deformations. The deformations can be used to warp subject images to template space, or aligning two subjects' images together by composing their deformations. These two operations are demonstrated in **Example 2**. Finally, **Example 3** fits an already learned MB model to two subject images. This is equivalent to registering to a pre-existing template space. Warping can then be done as shown in Example 2. For a full description of the model settings, see `demo_mb.m`.

### 1. Learning the MB model
``` matlab
% This script fits the MB model to a population of images. It is a
% group-wise image registration where an optimal K class tissue template is 
% built, as well as optimal intensity parameters learnt. The resulting
% deformations (written to dir_out prefixed y_) can then be used to warp
% between subject scans. Also native and template space tissue
% segmentaitons are written to disk, as well as a bias-field corrected 
% versions of the input scans.

% Data directory, assumed to contain a bunch of images.
dir_data = '/directory/with/some/images';  % OBS: should be .nii (NIfTIs)

% Get paths to data
pth_im  = spm_select('FPList',dir_data,'^.*\.nii$'); % Selects all NIfTIs files in dir_data

% RUN module (fits  the model)
run              = struct;
run.mu.create.K  = 5;                        % Number of classes in the TPMs. The actual TPMs will contain one more (implicit) class
run.mu.create.vx = 1.5;                      % Voxel size of the template -> smaller == faster, but less precise
run.onam         = 'mb_test';                % A name for the model
run.odir         = {'mb-output'};            % Output directory (to write model files)
run.nworker      = Inf;                      % Number of parallel workers
run.save         = true;                     % Save model at each 'epoch'
run.gmm(1).chan(1).images = cellstr(pth_im); % Image files

% OUT module (writes different outputs)
out        = struct;
out.result = {fullfile(run.odir{1},  ['mb_fit_' run.onam '.mat'])};
out.c      = 1:run.mu.create.K + 1;  % write classes in native space
out.wc     = 1:run.mu.create.K + 1;  % write classes in template space
out.mwc    = 1:run.mu.create.K + 1;  % write classes in modulated template space
out.sm     = 1:run.mu.create.K + 1;  % write scalar momentum

% Run jobs
jobs{1}.spm.tools.mb.run = run;
jobs{2}.spm.tools.mb.out = out;
spm_jobman('run', jobs);
```

### 2. Warping with MB deformations
``` matlab
% This script demonstrates warping with the deformation generated from
% fitting the MB model. The template file generated by MB, plus two images 
% with their corresponding forward deformations (prefixed 'y_') are
% needed. The following three warps are demonstrated:
% 1. Warp an image to template space (by pulling).
% 2. Warp an image to template space (by pushing).
% 3. Warp the template to image space.
% 4. Warp one image to another image.

% Path to a MB tissue template (this is dir_out in Example 1)
pth_mu = fullfile(dir_out,'mu_mb_test.nii');
% Paths to twoimages, to which MB has been fitted
pth_img1 = '/directory/with/some/img1.nii';
pth_img2 = '/directory/with/some/img2.nii';
% Paths to corresponding forward deformations of the above subjects
pth_y1 = fullfile(dir_out,'/y_img1.nii');
pth_y2 = fullfile(dir_out,'y_img2.nii');
% Where to write the warped images
dir_out = 'warped';

% image-to-template (pull)
% 1. Use forward deformation (pth_y1) to warp image 1 (pth_img1) to 
%    template space (pth_mu) by pulling
% -----------------
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def     = {pth_y1};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space           = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_img1};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'wpull';
spm_jobman('run',matlabbatch);

% image-to-template (push)
% 2. Use forward deformation (pth_y1) to warp image 1 (pth_img1) to 
%    template space (pth_mu) by pushing (with modulation and smoothing)
% -----------------
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.def                 = {pth_y1};
matlabbatch{1}.spm.util.defs.out{1}.push.fnames          = {pth_img1};
matlabbatch{1}.spm.util.defs.out{1}.push.weight          = {''};
matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.push.fov.file        = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.push.preserve        = 1;
matlabbatch{1}.spm.util.defs.out{1}.push.fwhm            = [10 10 10];
matlabbatch{1}.spm.util.defs.out{1}.push.prefix          = 'wpush';
spm_jobman('run',matlabbatch);

% 3. Use forward deformation (pth_y1) to template (pth_mu) to  image 
%    space (pth_img1) 
% -----------------
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.def                 = {pth_y1};
matlabbatch{1}.spm.util.defs.comp{2}.id.space            = {pth_img1};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';
spm_jobman('run',matlabbatch);

% 4. Compose both forward deformations (pth_y1, pth_y2) via template space
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
``` matlab
% This script demonstrates how to fit a learned MB model (template and 
% intensity prior) to new (unseen/test) images. In practice, it only involves modifying
% some settings in Example 1.

% Path to results MAT-file generated by fitting MB
pth_fit = fullfile(dir_out,'mb_fit_mb_test.mat');

% Make intensity prior file by extracting information from pth_fit
load(pth_fit,'sett'); 
pr      = sett.gmm.pr; 
mg_ix   = sett.gmm.mg_ix; 
save('prior_mb.mat','pr','mg_ix');

% Add the following fields to the RUN module of Example 1:
run.mu.exist           = {pth_mu};
run.gmm.pr.file        = {pth_int_prior};
run.gmm.pr.hyperpriors = [];  % Avoids re-learning the intensity prior
% and remove run.mu.create.K and run.mu.create.vx
```

## References
* Brudfors M, Balbastre Y, Flandin G, Nachev P, Ashburner J. Flexible Bayesian Modelling for Nonlinear Image Registration. In: International Conference on Medical Image Computing and Computer-Assisted Intervention 2020
* Brudfors M, Ashburner J, Nachev P, Balbastre Y. Empirical Bayesian Mixture Models for Medical Image Translation. In: International Workshop on Simulation and Synthesis in Medical Imaging 2019 Oct 13 (pp. 1-12). Springer, Cham.
* Balbastre Y, Brudfors M, Bronik K, Ashburner J. Diffeomorphic brain shape modelling using Gauss-Newton optimisation. In: International Conference on Medical Image Computing and Computer-Assisted Intervention 2018 Sep 16 (pp. 862-870). Springer, Cham.
* Blaiotta C, Freund P, Cardoso MJ, Ashburner J. Generative diffeomorphic modelling of large MRI data sets for probabilistic template construction. NeuroImage. 2018 Feb 1;166:117-34.
* Blaiotta C, Cardoso MJ, Ashburner J. Variational inference for medical image segmentation. Computer Vision and Image Understanding. 2016 Oct 1;151:14-28.
* Ashburner J, Friston KJ. Diffeomorphic registration using geodesic shooting and Gauss–Newton optimisation. NeuroImage. 2011 Apr 1;55(3):954-67.
* Ashburner J, Friston KJ. Computing average shaped tissue probability templates. Neuroimage. 2009 Apr 1;45(2):333-41.
* Ashburner J. A fast diffeomorphic image registration algorithm. Neuroimage. 2007 Oct 15;38(1):95-113.
* Ashburner J, Friston KJ. Unified segmentation. Neuroimage. 2005 Jul 1;26(3):839-51.

## Acknowledgements
This work was funded by the EU Human Brain Project’s Grant Agreement No 785907 (SGA2).
