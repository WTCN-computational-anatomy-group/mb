function Mat = GetMat(fin)
if isnumeric(fin)
    Mat = eye(4);
    return;
end
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    Mat  = fin(1).mat;
    return
end
error('Unknown datatype.');
end
%==========================================================================