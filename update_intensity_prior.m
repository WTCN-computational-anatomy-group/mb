function pr = update_intensity_prior(pr,po)
m0  = pr{1};
b0  = pr{2};
W0  = pr{3};
n0  = pr{4};

N = numel(po);
D = size(m0,1);
K = size(m0,2);

for k=1:K
    
%______________________________________________________________________________
%
% Compute m_0

g = zeros(D,1);
H = zeros(D,D);
for i=1:N    
    [m,b,W,n] = get_po(po,i);
    
    g = g + b0(k)*n(k)*W(:,:,k)*(m(:,k)-m0(:,k));
    H = H + b0(k)*n(k)*W(:,:,k);
end
m0(:,k) = m0(:,k) + H\g;
%______________________________________________________________________________


%______________________________________________________________________________
%
% Compute \beta_0

g_const = 0;
for i=1:N
    [m,b,W,n] = get_po(po,i);
    
    g_const = g_const - 0.5*(D/b(k) + n(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)));
end
for subit=1:100,
    % Diff w.r.t. b0
    g  = 0.5*N*D/b0(k) + g_const; % Gradient
    H  = 0.5*N*D/b0(k)^2;         % Hessian
    b0(k) = max(b0(k) + H\g,1e-5);
    if norm(g)==0, break; end
end
%______________________________________________________________________________


%______________________________________________________________________________
%
% Set up some constants

nW_const = zeros(D);
for i=1:N
    [m,b,W,n] = get_po(po,i);
    
    nW_const = nW_const + n(k)*W(:,:,k);
end

ElogLam = N*D*log(2);
for i=1:N
    [m,b,W,n] = get_po(po,i);
    
    ElogLam = ElogLam + 2*sum(log(diag(chol(W(:,:,k)))));
    for j=1:D
        ElogLam = ElogLam + psi((n(k)+1-j)/2);
    end
end

% convergence = [];
E = -realmax;

for it=1:1000
    %______________________________________________________________________________
    %
    % Compute objective function (Equation 10.74 of Bishop)

    oE = E;
    logB = -n0(k)*sum(log(diag(chol(W0(:,:,k))))) - n0(k)*D/2*log(2) - D*(D-1)/4*log(pi);
    for j=1:D, 
        logB = logB - gammaln((n0(k)+1-j)/2); 
    end
    E = (0.5*D*log(b0(k)/(2*pi)) + logB)*N + 0.5*(n0(k)-D-1)*ElogLam;
    for i=1:N
        [m,b,W,n] = get_po(po,i);
        
        e = 0.5*(-D*b0(k)/b(k) - b0(k)*n(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)))...
          - 0.5*n(k)*trace(W0(:,:,k)\W(:,:,k));
        E = E + e;
    end
    %if E-oE<abs(E)*eps*D^2, break; end
    if E-oE==0, break; end

%     convergence = [convergence E];
%     plot(convergence,'.-'); drawnow;


    %______________________________________________________________________________
    %
    % Compute \nu_0

    % Objective function terms containing n0:
    % NlogB = -n0*N*(sum(log(diag(chol(W0)))) + D/2*log(2));
    % for j=1:D, NlogB = NlogB - N*gammaln((n0+1-j)/2); end
    % E = NlogB + n0*0.5*ElogLam

    g = (sum(log(diag(chol(W0(:,:,k))))) + D/2*log(2))*N - 0.5*ElogLam;
    H = 0;
    for j=1:D
        g = g + N*psi(  (n0(k)+1-j)/2)/2;
        H = H + N*psi(1,(n0(k)+1-j)/2)/4;
    end
    n0(k) = max(n0(k) - H\g,D-0.99999);
    %______________________________________________________________________________



    %______________________________________________________________________________
    %
    % Compute W_0

    % Objective function terms containing W0:
    % E = -n0*N*sum(log(diag(chol(W0))));
    % for i=1:N
    %    E = E - 0.5*n(i)*trace(W0\W(:,:,i));
    % end

    C = inv(chol(W0(:,:,k)));

    % Objective function terms containing W0, after
    % re-expressing using C = inv(chol(W0)):
    % E = n0*N*sum(log(diag(C)));
    % for i=1:N
    %    E = E - 0.5*n(i)*trace(C'*W(:,:,i)*C);
    % end

    G  = -n0(k)*N*diag(1./diag(C)) + nW_const*C;
    for d=1:D,
        c        = C(1:d,d);
        g        = G(1:d,d);
        H        = nW_const(1:d,1:d);
        H(d,d)   = H(d,d) + n0(k)*N/c(d)^2;
        C(1:d,d) = c - H\g;
    end
    C         = inv(C);
    W0(:,:,k) = C'*C;

end

end

pr{1} = m0;
pr{2} = b0;
pr{3} = W0;
pr{4} = n0;
%==========================================================================

