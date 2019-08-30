function P = softmax(mu)
d   = 4;
mx  = max(mu,[],d);
E   = exp(mu-mx);
den = sum(E,d)+exp(-mx);
P   = E./den;
end
%==========================================================================