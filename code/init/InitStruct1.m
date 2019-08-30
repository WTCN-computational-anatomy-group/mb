function dat = InitStruct1(F,sett)

% Need a better way of initialising these parameters
K1    = sett.K+1;
mog   = struct('mu',(1:-(1/K1):(1/K1))*2000,'sig2',ones(1,K1)/3*sqrt(1000)); % Random
%mog  = struct('mu',ones(1,K1)*500,'sig2',ones(1,K1)*500^2); % Same (for existing TPM)

M0    = eye(4);
dat   = struct('f',F,'M',M0, 'q',zeros(6,1), 'v',[], 'psi',[], 'E',[0 0],'mog',mog);
for n=1:numel(dat)
    if isnumeric(dat(n).f)
        dat(n).Mat     = eye(4); % Should really do this better
    else
        if isa(dat(n).f,'char')
            dat(n).f   = nifti(dat(n).f);
        end
        if isa(dat(n).f,'nifti')
            dat(n).Mat = dat(n).f(1).mat;
        end
    end
end
end
%==========================================================================