function varargout = spm_multireg_energ(varargin)
%__________________________________________________________________________
%
% Energy functions for spm_multireg.
%
% FORMAT lb  = spm_multireg_energ('LowerBound',type,varargin)
% FORMAT lb  = spm_multireg_energ('SumLowerBound',lb)
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
    case 'LowerBound'
        [varargout{1:nargout}] = LowerBound(varargin{:});       
    case 'SumLowerBound'
        [varargout{1:nargout}] = SumLowerBound(varargin{:});            
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
% LowerBound()
function lb = LowerBound(type,varargin)
if strcmpi(type,'XB')    
    bf = varargin{1};
    
    bf(isnan(bf)) = 1;        
    bf            = log(prod(bf,2)); 
    lb            = sum(double(bf));
elseif strcmpi(type,'X')    
    fn   = varargin{1};
    zn   = varargin{2};
    code = varargin{3};
    mean = varargin{4};
    prec = varargin{5};
    L    = unique(code);
    
    [lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', fn, zn, 1, {code,L});
    lb               = spm_gmm_lib('MarginalSum', lSS0, lSS1, lSS2, mean, prec, L);
else
    error('Undefined type!');
end
end
%==========================================================================

%==========================================================================
% SumLowerBound()
function lb = SumLowerBound(lb)
fields          = fieldnames(lb);
lb.sum(end + 1) = 0;
for i=1:numel(fields)
    field = fields{i};
    if ~any(strcmpi(field, {'sum' 'last'})) && ~isempty(lb.(field)) && ~isnan(lb.(field)(end))
        lb.sum(end) = lb.sum(end) + sum(lb.(field)(:,end));
    end
end
end
%==========================================================================

%==========================================================================
% TemplateEnergy()
function E = TemplateEnergy(mu,sett)

% Parse function settings
mu_settings = sett.var.mu_settings;

spm_multireg_util('SetBoundCond');
g    = spm_field('vel2mom', mu, mu_settings);
E    = 0.5*mu(:)'*g(:);
end
%==========================================================================

%==========================================================================
% VelocityEnergy()
function dat = VelocityEnergy(dat,sett)

% Parse function settings
threads    = sett.gen.threads;
v_settings = sett.var.v_settings;

spm_multireg_util('SetBoundCond');
if threads>1 && numel(dat)>1
    % Memory  = 4*2*3*prod(GetSize(dat(1).v));
    % nthread = min(sett.gen.threads,floor(sett.maxmem/Memory));
    parfor(n=1:numel(dat),threads) % PARFOR
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