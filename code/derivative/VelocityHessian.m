function H = VelocityHessian(mu,G,accel)
d  = [size(mu,1),size(mu,2),size(mu,3)];
M  = size(mu,4);
if accel>0, s  = softmax(mu); end
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