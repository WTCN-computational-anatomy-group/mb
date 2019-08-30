function dat = UpdateWarps(dat,sett)
if sett.groupwise
    % Total initial velocity should be zero (Khan & Beg)
    avg_v = single(0);
    for n=1:numel(dat)
        avg_v = avg_v + GetData(dat(n).v); % For mean correcting initial velocities
    end
    avg_v = avg_v/numel(dat);
    d     = [size(avg_v,1) size(avg_v,2) size(avg_v,3)];
else
    avg_v = [];
    d     = GetSize(dat(1).v);
end
kernel = shoot(d,sett.v_settings);
if sett.threads>1 && numel(dat)>1
    % Memory = loads
    parfor(n=1:numel(dat),sett.threads) % PARFOR
        dat(n) = UpdateWarpsSub(dat(n),avg_v,sett,kernel);
    end
else
    for n=1:numel(dat)
        dat(n) = UpdateWarpsSub(dat(n),avg_v,sett,kernel);
    end
end
end
%==========================================================================


%==========================================================================
function datn = UpdateWarpsSub(datn,avg_v,sett,kernel)
v        = GetData(datn.v);
if ~isempty(avg_v)
    v    = v - avg_v;
end
datn.v   = SetData(datn.v,v);
psi1     = shoot(v, kernel, sett.args); % Geodesic shooting
datn.psi = SetData(datn.psi,psi1);
end
%==========================================================================