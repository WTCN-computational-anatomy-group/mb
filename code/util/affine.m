function psi0 = affine(d,Mat)
id    = identity(d);
psi0  = reshape(reshape(id,[prod(d) 3])*Mat(1:3,1:3)' + Mat(1:3,4)',[d 3]);
end
%==========================================================================