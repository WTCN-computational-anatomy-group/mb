function a1 = Pull1(a0,psi,r)
% Resample an image or set of images
% FORMAT a1 = Pull1(a0,psi,r)
%
% a0  - Input image(s)
% psi - Deformation
% r   - subsampling density in each dimension (default: [1 1 1])
%
% a1  - Output image(s)
%
% There are also a couple of Push1 and Pull1 functions, which might be of
% interest.  The idea is to reduce aliasing effects in the pushed images,
% which might also be useful for dealing with thick-sliced images.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$
if nargin<3, r=[1 1 1]; end

set_bound;

if isempty(a0)
    a1 = a0;
elseif isempty(psi)
    a1 = a0;
else
    if r==1
        a1 = spm_diffeo('pullc',a0,psi);
        return
    end
    d  = [size(a0) 1 1];
    if d(3)>1, zrange = range(r(3)); else, zrange = 0; end
    if d(2)>1, yrange = range(r(2)); else, yrange = 0; end
    if d(1)>1, xrange = range(r(1)); else, xrange = 0; end
    id = Identity(size(psi));
    a1 = zeros([size(psi,1),size(psi,2),size(psi,3),size(a0,4)],'single');
    for l=1:d(4)
        tmp = single(0);
        al  = single(a0(:,:,:,l));
        for dz=zrange
            for dy=yrange
                for dx=xrange
                    ids  = id  + cat(4,dx,dy,dz);
                    psi1 = spm_diffeo('pull',psi-id,    ids)+ids;
                    as   = spm_diffeo('pull',al,psi1);
                   %ids  = id  - cat(4,dx,dy,dz);
                    tmp  = tmp + spm_diffeo('push',as,  ids);
                end
            end
        end
        a1(:,:,:,l) = tmp/(numel(zrange)*numel(yrange)*numel(xrange));
    end
end

function r = range(n)
r = (-floor((n-1)/2):ceil((n-1)/2))/n;

