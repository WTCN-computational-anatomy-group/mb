function varargout = spm_multireg_der(varargin)
%__________________________________________________________________________
%
% Derivative functions for spm_multireg.
%
% FORMAT [H,g] = spm_multireg_der('AffineHessian',mu,G,a,w,accel)
% FORMAT H     = spm_multireg_der('AppearanceHessian',mu,accel)
% FORMAT [H,g] = spm_multireg_der('SimpleAffineHessian',mu,G,H0,a,w)
% FORMAT H     = spm_multireg_der('VelocityHessian',mu,G,accel)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_multireg_der
    error('Not enough argument. Type ''help spm_multireg_der'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id 
    case 'AffineHessian'
        [varargout{1:nargout}] = AffineHessian(varargin{:});
    case 'AppearanceHessian'
        [varargout{1:nargout}] = AppearanceHessian(varargin{:});        
    case 'SimpleAffineHessian'
        [varargout{1:nargout}] = SimpleAffineHessian(varargin{:});        
    case 'VelocityHessian'
        [varargout{1:nargout}] = VelocityHessian(varargin{:});                
    otherwise
        help spm_multireg_der
        error('Unknown function %s. Type ''help spm_multireg_der'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% AffineHessian()
function [H,g] = AffineHessian(mu,G,a,w,accel)
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = Horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(1:d(1),1:d(2),1,1);
for i=1:d(3)
    x{3} = x{3}*0+i;
    gv   = reshape(sum(a(:,:,i,:).*G(:,:,i,:,:),4),[d(1:2) 1 3]);
    Hv   = w(:,:,i).*VelocityHessian(mu(:,:,i,:),G(:,:,i,:,:),accel);
    for i1=1:12
        k1g   = rem(i1-1,3)+1;
        k1x   = floor((i1-1)/3)+1;
        g(i1) = g(i1) + sum(sum(sum(x{k1x}.*gv(:,:,:,k1g))));
        for i2=1:12
            k2g      = rem(i2-1,3)+1;
            k2x      = floor((i2-1)/3)+1;
            H(i1,i2) = H(i1,i2) + sum(sum(sum(x{k1x}.*Hv(:,:,:,I(k1g,k2g)).*x{k2x})));
        end
    end
end
end
%========================================================================== 

%==========================================================================
% AppearanceHessian()
function H = AppearanceHessian(mu,accel)
M  = size(mu,4);
d  = [size(mu,1) size(mu,2) size(mu,3)];
if accel>0, s  = spm_multireg_util('softmax',mu); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
I  = Horder(M);
H  = zeros([d (M*(M+1))/2],'single');
for m1=1:M
    for m2=m1:M
        if accel==0
            tmp = Ab(m1,m2)*ones(d,'single');
        else
            if m2~=m1
                tmp = accel*(-s(:,:,:,m1).*s(:,:,:,m2))           + (1-accel)*Ab(m1,m2);
            else
                tmp = accel*(max(s(:,:,:,m1).*(1-s(:,:,:,m1)),0)) + (1-accel)*Ab(m1,m2);
            end
        end
        H(:,:,:,I(m1,m2)) = tmp;
    end
end
end
%==========================================================================

%==========================================================================
% SimpleAffineHessian()
function [H,g] = SimpleAffineHessian(mu,G,H0,a,w)
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = Horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(1:d(1),1:d(2),1,1);
for i=1:d(3)
    x{3} = x{3}*0+i;
    gv   = reshape(sum(a(:,:,i,:).*G(:,:,i,:,:),4),[d(1:2) 1 3]);
    Hv   = w(:,:,i).*H0(:,:,i,:);
    for i1=1:12
        k1g   = rem(i1-1,3)+1;
        k1x   = floor((i1-1)/3)+1;
        g(i1) = g(i1) + sum(sum(sum(x{k1x}.*gv(:,:,:,k1g))));
        for i2=1:12
            k2g      = rem(i2-1,3)+1;
            k2x      = floor((i2-1)/3)+1;
            H(i1,i2) = H(i1,i2) + sum(sum(sum(x{k1x}.*Hv(:,:,:,I(k1g,k2g)).*x{k2x})));
        end
    end
end
end
%==========================================================================

%==========================================================================
% VelocityHessian()
function H = VelocityHessian(mu,G,accel)
d  = [size(mu,1),size(mu,2),size(mu,3)];
M  = size(mu,4);
if accel>0, s  = spm_multireg_util('softmax',mu); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
H  = zeros([d 6],'single');
for m1=1:M
    for m2=1:M
        if accel==0
            tmp = Ab(m1,m2)*ones(d,'single');
        else
            if m2~=m1
                tmp = (-s(:,:,:,m1).*s(:,:,:,m2))*accel           + (1-accel)*Ab(m1,m2);
            else
                tmp = (max(s(:,:,:,m1).*(1-s(:,:,:,m1)),0))*accel + (1-accel)*Ab(m1,m2);
            end
        end
        H(:,:,:,1) = H(:,:,:,1) + tmp.*G(:,:,:,m1,1).*G(:,:,:,m2,1);
        H(:,:,:,2) = H(:,:,:,2) + tmp.*G(:,:,:,m1,2).*G(:,:,:,m2,2);
        H(:,:,:,3) = H(:,:,:,3) + tmp.*G(:,:,:,m1,3).*G(:,:,:,m2,3);
        H(:,:,:,4) = H(:,:,:,4) + tmp.*G(:,:,:,m1,1).*G(:,:,:,m2,2);
        H(:,:,:,5) = H(:,:,:,5) + tmp.*G(:,:,:,m1,1).*G(:,:,:,m2,3);
        H(:,:,:,6) = H(:,:,:,6) + tmp.*G(:,:,:,m1,2).*G(:,:,:,m2,3);
    end
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% Horder()
function I = Horder(d)
I = diag(1:d);
l = d;
for i1=1:d
    for i2=(i1+1):d
        l = l + 1;
        I(i1,i2) = l;
        I(i2,i1) = l;
    end
end
end
%==========================================================================