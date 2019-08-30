clear;

%------------------
% Input data parameters
%------------------

dir_data = fullfile(fileparts(mfilename('fullpath')),'data-ixi-t1-32'); % Directory with MRIs
num_subj = 16;                                                          % Max number of subjects
samp     = 2;                                                           % Sample distance (> 1 does subsampling)
do_2d    = 2;                                                           % 0 = 3D, 1 - 3 = different axes in 2D

%------------------
% Get test data if not present
%------------------

if ~(exist(dir_data,'dir') == 7)
    fprintf('Getting test data...')
    mkdir(dir_data);
    url         = 'https://www.dropbox.com/s/awo45fl3e5in0pb/data-ixi-t1-32.zip?dl=1';
    filename    = fullfile(dir_data,'data-ixi-t1-32.zip');
    outfilename = websave(filename,url);
    unzip(outfilename,dir_data);
    delete(outfilename);
    fprintf('done!\n')
end

%------------------
% Read data
%------------------

Nii = nifti(spm_select('FPList',dir_data,'^vx.*\.nii$'));

%------------------
% Make input data (F)
%------------------

S = min(num_subj,numel(Nii));
F = cell(1,S);
for s=1:S
    
    if samp > 1
        
        V       = spm_vol(Nii(s).dat.fname);    
        dm0     = V(1).dim;
        dm0     = [dm0 1];
        vx      = sqrt(sum(V(1).mat(1:3,1:3).^2));
        sk      = max([1 1 1],round(samp*[1 1 1]./vx));
        [x,y,z] = ndgrid(1:sk(1):dm0(1),1:sk(2):dm0(2),1:sk(3):dm0(3));
        dm      = size(x);
        
        nim = [];
        for k=1:dm0(4)
            imk = spm_sample_vol(V,x,y,z,0);
            
            nim = cat(4,nim,reshape(imk,dm));
        end        
        
        im = nim;
    else
        im = single(Nii(s).dat());            
    end
    
    if do_2d
        dm = size(im);
        dm = [dm 1];
        ix = ceil(dm(min(do_2d,3))*0.55);
        if do_2d == 1
            im = squeeze(im(ix,:,:,:));
        elseif do_2d == 2
            im = squeeze(im(:,ix,:,:));
        else
            im = squeeze(im(:,:,ix,:));
        end
    end
    
    F{s} = im;
end
clear im nim

%------------------
% Run algorithm
%------------------

[dat,mu] = multireg3D(F);