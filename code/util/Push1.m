function [f1,w1] = Push1(f,psi,d,r)
% Push an image (or set of images) accorging to a spatial transform
% FORMAT [f1,w1] = Push1(f,psi,d,r)
%
% f   - Image (3D or 4D)
% psi - Spatial transform
% d   - dimensions of output (default: size of f)
% r   - subsampling density in each dimension (default: [1 1 1])
%
% f1  - "Pushed" image
%
% There are also a couple of Push1 and Pull1 functions, which might be of
% interest.  The idea is to reduce aliasing effects in the pushed images,
% which might also be useful for dealing with thick-sliced images.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$
if nargin<4, r = [1 1 1]; end
if nargin<3, d = [size(f,1) size(f,2) size(f,3)]; end

set_bound;

%msk    = isfinite(f);
%f(msk) = 0;
if ~isempty(psi)
    if r==1
        if nargout==1
            f1      = spm_diffeo('pushc',single(f),psi,d);
        else
            [f1,w1] = spm_diffeo('pushc',single(f),psi,d);
        end
        return
    end

    if d(3)>1, zrange = range(r(3)); else, zrange = 0; end
    if d(2)>1, yrange = range(r(2)); else, yrange = 0; end
    if d(1)>1, xrange = range(r(1)); else, xrange = 0; end

    id    = Identity(size(psi));
    f1    = single(0);
    w1    = single(0);
    for dz=zrange
        for dy=yrange
            for dx=xrange
                ids       = id + cat(4,dx,dy,dz);
                psi1      = spm_diffeo('pull',psi-id,    ids)+ids;
                fs        = spm_diffeo('pull',single(f), ids);
               %fs=single(f);
                if nargout==1
                    fs        = spm_diffeo('push',fs,        psi1,d);
                    f1        = f1  + fs;
                else
                    [fs,ws]   = spm_diffeo('push',fs,        psi1,d);
                    f1        = f1  + fs;
                    w1        = w1  + ws;
                end
            end
        end
    end
    scale = 1/(numel(zrange)*numel(yrange)*numel(xrange));
    f1    = f1*scale;
    w1    = w1*scale;
else
    msk      = isfinite(f);
    f1       = f;
    f1(~msk) = 0;
    w1       = single(all(msk,4));
end

function r = range(n)
r = (-floor((n-1)/2):ceil((n-1)/2))/n;

