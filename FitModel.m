function FitModel(model,varargin)
%__________________________________________________________________________
%
% Some code for testing the diffeo_segment algorithm.
%__________________________________________________________________________

if strcmp(model,'groupwise')
    % Groupwise
    do_gw = true;
    P     = varargin{1};
    if numel(varargin) >= 2, sett = varargin{2};
    else,                    sett = struct; 
    end
elseif strcmp(model,'register')
    % Register
    do_gw = false;
    P     = varargin{1};
    model = varargin{2};
    if numel(varargin) >= 3, sett = varargin{3};
    else,                    sett = struct; 
    end
else
    error('Undefined model!')
end

% Get data and (default) settings (for the given test_case)
[dat,sett] = pop2dat(P,sett);
    
% Get input data
in  = dat2in(dat); 
clear dat

% If SPM has been compiled with OpenMP support then the number of threads
% can here be set to speed up the algorithm.
setenv('SPM_NUM_THREADS',sprintf('%d',-1));

if do_gw
    % Run Groupwise
    [dat,model,sett] = spm_mb_fit(in,'sett',sett);
else
    % Run Register
    [dat,model,sett] = spm_mb_fit(in,'model',model,'sett',sett);
end

% Write results in normalised space
spm_mb_output(dat,model,sett);

end
%==========================================================================

%==========================================================================
% dat2in()
function in = dat2in(dat)
Npop = numel(dat);
s    = struct('F',[],'do_bf',[],'ix_pop',[],'is_ct',[],'labels',[]);
for p=1:Npop % loop over populations
    pths_im = dat(p).pths_im;
    N       = size(pths_im{1},1);
    C       = numel(pths_im);
    for n=1:N % Loop over subjects        
        im = nifti;      
        for c=1:C % Loop over channels
            im(c) = nifti(deblank(pths_im{c}(n,:)));
        end
        if ~isempty(dat(p).pths_lab) 
            pths_lab = deblank(dat(p).pths_lab(n,:));
            lab      = {nifti(pths_lab),dat(p).cm_map};
        else
            lab = cell(1,2);
        end
        if p == 1 && n == 1
            in        = s;
            in.F      = im;
            in.do_bf  = dat(p).do_bf;
            in.ix_pop = dat(p).ix_pop;
            in.is_ct  = dat(p).is_ct;
            in.labels = lab;
        else
            ins        = s;
            ins.F      = im;
            ins.do_bf  = dat(p).do_bf;
            ins.ix_pop = dat(p).ix_pop;
            ins.is_ct  = dat(p).is_ct;
            ins.labels = lab;
            in         = [in ins];
        end
    end
end
end
%==========================================================================

%==========================================================================
% pop2dat()
function [dat,sett] = pop2dat(P,sett)

Npop = numel(P);
cl   = cell(Npop,1);
dat  = struct('pths_im',cl,'pths_lab',cl,'do_bf',cl,'ix_pop',cl,'is_ct',cl,'cm_map',cl);
for p=1:Npop % loop over populations
    
    % Defaults
    ix_pop = []; is_ct = false; do_bf = true; cm_map = []; Nsubj = Inf;
    
    dir_data                    = P{p}{1};
    modality                    = P{p}{2};    
    if numel(P{p}) >= 3, Nsubj  = P{p}{3}; end
    if numel(P{p}) >= 4, ix_pop = P{p}{4}; end
    if numel(P{p}) >= 5, cm_map = P{p}{5}; end
    if numel(P{p}) >= 6, is_ct  = P{p}{6}; end
    if numel(P{p}) >= 7, do_bf  = P{p}{7}; end
    
    if isempty(P{p}{4}), ix_pop = p; end
    
    dat(p).ix_pop = ix_pop;    
    dat(p).cm_map = cm_map;    
    dat(p).is_ct  = is_ct;    
    if is_ct, dat(p).do_bf = false;
    else,     dat(p).do_bf = do_bf;
    end
     
    % Get name of population (assumed to be folder name)
    name = strsplit(dir_data,filesep);
    name = name{end};
    
    if isempty(modality), C = 1;
    else,                 C = numel(modality);
    end
    
    pths_im  = cell(1,C); % Images
    pths_lab = {};        % Labels
    switch name % switch over defined populations
        
        case 'ATLAS'
            pths_im{1} = spm_select('FPList',dir_data,'^((?!labels).)*\.nii$');
            pths_lab   = spm_select('FPList',dir_data,'labels.*\.nii$');
        case 'BALGRIST'
            for c=1:C
                if strcmp(modality{c},'T1')
                    pths_im{c} = spm_select('FPList',dir_data,'^((?!spine_labels|_PDw).)*\.nii$');
                elseif strcmp(modality{c},'PD')
                    pths_im{c} = spm_select('FPList',dir_data,'_PDw.*\.nii$');
                end
            end            
            pths_lab = spm_select('FPList',dir_data,'spine_labels.*\.nii$');            
        case 'CROMISLABELS'
            pths_im{1} = spm_select('FPList',dir_data,'^((?!_smask).)*\.nii$');
            pths_lab   = spm_select('FPList',dir_data,'_smask.*\.nii$');
        case {'CROMIS','DELIRIUM'}
            pths_im{1} = spm_select('FPList',dir_data,'^.*\.nii$');
        case {'IXI','IXIMISALIGN'}
            for c=1:C
                pths_im{c} = spm_select('FPList',dir_data,['-' modality{c} '-.*\.nii$']);                                
            end
        case {'IXIC','IXIRC'}
            sett.do.gmm = false;    
            sett.do.gmm = false; 
            map         = containers.Map;
            map('GM')   = '1'; map('WM') = '2';  map('CSF') = '3'; 
            for i=1:C
                if strcmp(name,'IXIC')
                    pths_im{i} = spm_select('FPList',dir_data,['^c' map(modality{i}) '.*\.nii$']); 
                elseif strcmp(name,'IXIRC')
                    sett.model.vx    = 1.5;  
                    sett.do.updt_aff = false;
                    pths_im{i}       = spm_select('FPList',dir_data,['^rc' map(modality{i}) '.*\.nii$']);                     
                end
            end   
        case 'MICCAI2012'
            pths_im{1} = spm_select('FPList',dir_data,'^((?!_glm).)*\.nii$');
            pths_lab   = spm_select('FPList',dir_data,'_glm.*\.nii$');
        case 'MRBRAINS18'       
            for c=1:C
                pths_im{c} = spm_select('FPList',dir_data,['-' modality{c} '.*\.nii$']);
            end
            pths_lab = spm_select('FPList',dir_data,'segm.*\.nii$');  
        otherwise            
            error('Unknown population!')
    end
    
    if iscell(Nsubj)
        for c=1:C,        pths_im{c} = pths_im{c}(cell2mat(Nsubj),:); end
        if ~isempty(pths_lab), pths_lab = pths_lab(cell2mat(Nsubj),:); end
    else
        for c=1:C,        pths_im{c} = pths_im{c}(1:min(Nsubj,size(pths_im{c},1)),:); end
        if ~isempty(pths_lab), pths_lab = pths_lab(1:min(Nsubj,size(pths_lab,1)),:); end
    end
    
    % Set images and labels
    dat(p).pths_im  = pths_im;
    dat(p).pths_lab = pths_lab;
    
    if strcmp(name,'IXIC') || strcmp(name,'IXIRC')
        % Population requested has tissue segmentations, return only this
        % population
        dat = dat(p);
        break
    end
end
end
%==========================================================================