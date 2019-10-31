function varargout = spm_multireg_energ(varargin)
%__________________________________________________________________________
%
% Energy functions for spm_multireg.
%
% FORMAT E   = spm_multireg_energ('TemplateEnergy',mu,sett)
% FORMAT dat = spm_multireg_energ('VelocityEnergy',dat,sett)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_energ
    error('Not enough argument. Type ''help spm_multireg_energ'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id 
    case 'TemplateEnergy'
        [varargout{1:nargout}] = TemplateEnergy(varargin{:});         
    case 'VelocityEnergy'
        [varargout{1:nargout}] = VelocityEnergy(varargin{:});                 
    otherwise
        help spm_multireg_energ
        error('Unknown function %s. Type ''help spm_multireg_energ'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% TemplateEnergy()
function E = TemplateEnergy(mu,sett)
spm_multireg_util('SetBoundCond');
g    = spm_field('vel2mom', mu, sett.var.mu_settings);
E    = 0.5*mu(:)'*g(:);
end
%==========================================================================

%==========================================================================
% VelocityEnergy()
function dat = VelocityEnergy(dat,sett)
spm_multireg_util('SetBoundCond');
v_settings = sett.var.v_settings;
if sett.gen.threads>1 && numel(dat)>1
    % Memory  = 4*2*3*prod(GetSize(dat(1).v));
    % nthread = min(sett.gen.threads,floor(sett.maxmem/Memory));
    parfor(n=1:numel(dat),sett.gen.threads) % PARFOR
        spm_multireg_util('SetBoundCond');
        v           = spm_multireg_io('GetData',dat(n).v);
        u0          = spm_diffeo('vel2mom', v, v_settings); % Initial momentum
        dat(n).E(2) = 0.5*sum(u0(:).*v(:));                 % Prior term
    end
else
    for n=1:numel(dat)
        v           = spm_multireg_io('GetData',dat(n).v);
        u0          = spm_diffeo('vel2mom', v, v_settings); % Initial momentum
        dat(n).E(2) = 0.5*sum(u0(:).*v(:));                 % Prior term
    end
end
end
%==========================================================================