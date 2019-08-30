function [d,M] = GetSize(fin)
d = [GetDimensions(fin) 1 1 1];
M = d(4);
d = d(1:3);
end

function d = GetDimensions(fin)
if isnumeric(fin)
    d = size(fin);
    return
end
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    M    = numel(fin);
    d    = size(fin(1).dat,[1 2 3 4 5]);
    if M>1
        d(4) = M;
    else
        if numel(d)>4 && d(4)==1
            d = [d(1:3) d(5)];
        end
    end
    return
end
end
%==========================================================================