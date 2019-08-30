function fout = GetData(fin)
if isnumeric(fin)
    fout = single(fin);
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
        fout = zeros(d,'single');
        for m=1:M
            fout(:,:,:,m) = single(fin(m).dat(:,:,:,:,:));
        end
    else
        fout = single(fin.dat(:,:,:,:,:));
        if numel(d)>4 && d(4)==1
            fout = reshape(fout,[d(1:3) d(5)]);
        end
    end
    return
end
error('Unknown datatype.');
end
%==========================================================================