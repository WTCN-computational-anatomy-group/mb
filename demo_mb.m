% =========================================================================
%
% A (sort of) documentation for scripting SPM Multi Brain
%
% This file defines all possible fields of a MB batch.
% In your own script, you do not need to define them all, the jobman will
% fill all missing fields with default values for you.
%
% The only mandatory fields are at least one of:
% * run.cat
% * run.gmm.chan.images
%
% =========================================================================

% Yael Balbastre, 2020

% -------------------------------------------------------------------------
%  (1) Training the model
% -------------------------------------------------------------------------

% --- Describe the model

run.mu.create.K  = 9;                                % Number of classes, without background (the final template will have K+1 classes)
run.mu.create.vx = 1;                                % Voxel size of the template
run.onam         = 'mb';                             % A name for the model
run.odir         = {'model_dir'};                    % Output directory for model files
run.save         = true;                             % Save the model at each epoch
run.nworker      = 0;                                % Number of parallel workers

% --- First possible type of input: categorical maps 
%     (e.g., from spm unified segmentation)

run.cat = {                                          % List of categorical images (e.g., c1*, c2*, c3* from spm unified segmentation)
    {'c1_sub1.nii' 'c2_sub1.nii'}                    %  >> we have two classes (g, wm) for each subject
    {'c1_sub2.nii' 'c2_sub2.nii'}}';                 %  >> leave empty if no categorical images available

% --- Second possible type of input: intensity images (MRI, CT) 
%     They are organised in populations (or datasets)
%     Images from a population should be similar. In particular, they
%     should have the same channels. Population-specific priors are 
%     learned during training.

for i=1:nb_populations
    
% - Images - 
for j=1:nb_channels(i)
    run.gmm(i).chan(j).images      = {               % List of images for channel j of population i
        'popi_chanj_sub1.nii' 
        'popi_chanj_sub2.nii'}';                
    run.gmm(i).chan(j).inu.inu_reg = 10000;          % Regularisation of the intensity non-uniformity (bias) field
    run.gmm(i).chan(j).inu.inu_co  = 40;             % [EXPERT] Number of basis functions of the intensity non-uniformity (bias) field
    run.gmm(i).chan(j).modality    = 1;              % Modality type: 1 = MRI, 2 = CT
end % channel loop

% - Parameters -
run.gmm(i).pr.file        = {'priors.mat'};          % File that contains pre-trained priors. Leave empty if none.
run.gmm(i).pr.hyperpriors = {               
    'b0_priors' {0.1 0.01}};                         % [EXPERT] Prior on the degrees of freedom of the means: Gamma(alpha, beta)
run.gmm(i).tol_gmm        = 0.0005;                  % [EXPERT] Tolerance for early stopping when fitting the GMM
run.gmm(i).nit_gmm_miss   = 32;                      % [EXPERT] Max number of iterations when inferring missing values
run.gmm(i).nit_gmm        = 8;                       % [EXPERT] Max number of iterations when fitting the GMM
run.gmm(i).nit_appear     = 4;                       % [EXPERT] Max number of iterations when fitting the bias field (?)

% --- Semi supervision
%     If manual labels are available in some populations, they can be 
%     used to supervise the automated clustering.
%     The user should map each manual label to one or more learned classes.

if we_have_manual_labels
    run.gmm(i).labels.true.images  = {               % Manual labels
        'popi_sub1_labels.nii'
        'popi_sub2_labels.nii'
        }';
    run.gmm(i).labels.true.cm_map  = {               % Confusion matrix
        [1, 2]         % label 1                     % = mapping from each manual label to possible classes
        [3, 4]         % label 2                     % There should be L+1 rows, where L is the maximum manual label
        [5, 6]         % label 3                     % >> possible manual labels are [1, 2, 3], so we have 4 rows
        [7, 8, 9]}';   % label 0                     % >> the last row corresponds to unlabelled voxels
    run.gmm(i).labels.true.w       = 0.99;           % [EXPERT] Confidence in the labels: 1 = full, 0 = none
else
    run.gmm(i).labels.false        = [];             % No manual labels 
end

end % population loop

% --- Expert parameters

run.mu.create.mu_settings = [1e-05 0.5 0];           % [EXPERT] Regularisation of the template (absolute values, membrane, bending)
run.aff                   = 'SE(3)';                 % [EXPERT] Affine model, from  {'SE(3)', 'SO(3)', 'T(3)'}
run.v_settings            = [0.0001 0 0.4 0.1 0.4];  % [EXPERT] Regularisation of the deformation (absolute displacement, membrane, bending, shears, zooms)
run.accel                 = 0.8;                     % [EXPERT] Stabilisation of the optimisation: 0 = slower/more stable, 1 = faster/less stable
run.min_dim               = 16;                      % [EXPERT] Template dimension at the coarsest pyramid level
run.tol                   = 0.001;                   % [EXPERT] Tolerance for early stopping
run.sampdens              = 2;                       % [EXPERT] Sampling density when fitting the GMM


% -------------------------------------------------------------------------
%  (2) Writing subject-level data
% -------------------------------------------------------------------------
% By default, the fitting procedure only write model files (tissue 
% probability maps, population priors, etc.). The `out` module allows
% the user to write a selection of maps at the subject-level
% (segmentation, normalised images, etc.)

model_file = fullfile(run.odir{1},  ['mb_fit_' run.onam '.mat']);

out.result = {model_file};                           % Path to the model file generated by run. 
                                                     %    In theory we could use spm's dependency system 
                                                     %    but it is a bit compilcated/ugly in a script.
out.i      = false;                                  % Write images with missing values filled by the model.
out.mi     = false;                                  % Write bias-corrected images in native space
out.wi     = false;                                  % Write images in template space
out.wmi    = false;                                  % Write bias-corrected images in template space
out.inu    = false;                                  % Write estimated bias field
out.c      = [];                                     % List of classes to write, in native space
out.wc     = [];                                     % List of classes to write, in template space
out.mwc    = [];                                     % List of classes to write, jacobian-modulated, in template space
out.mrf    = 1;                                      % Use Markov Random Field cleaning


% -------------------------------------------------------------------------
%  (3) Execute
% -------------------------------------------------------------------------

jobs{1}.spm.tools.mb.run = run;
jobs{2}.spm.tools.mb.out = out;
spm_jobman('run', jobs);