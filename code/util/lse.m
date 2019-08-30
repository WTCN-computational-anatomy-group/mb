function l = lse(mu)
d   = 4;
mx  = max(mu,[],d);
l   = log(exp(-mx)+sum(exp(mu-mx),d))+mx;
end
%==========================================================================