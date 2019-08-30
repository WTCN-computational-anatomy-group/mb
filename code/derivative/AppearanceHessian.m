function H = AppearanceHessian(mu,accel)
M  = size(mu,4);
d  = [size(mu,1) size(mu,2) size(mu,3)];
if accel>0, s  = softmax(mu); end
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