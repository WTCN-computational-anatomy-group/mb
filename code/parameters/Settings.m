function sett = Settings(is3d)
% Should really be able to pass settings to the functions but use these ones as defaults
sett.args        = 8;     % Timesteps for shooting integration
if is3d
    sett.B       = AffineBases('SE(3)');
else
    sett.B       = AffineBases('SE(2)');
end
sett.scal        = 1.0;   % Scaling of GN updates
sett.s_settings  = [3 2]; % Full multigrid settings
sett.accel       = 0;   % 0 <= accel <= 1: 0 -> slow & stable; 1 -> fast & unstable
sett.nits        = 3;     % Number of iterations at final resolution
sett.mu_settings = [1e-3 0.2 0];
sett.v_settings  = [0    0   0.2 0.05 0.2]*4;
sett.threads     = 8;
sett.is3d        = is3d;

% .model - model-specific
sett.model.K = 12;
if is3d
    sett.model.vx = 1.5;
else
    sett.model.vx = 1;
end

% .show - visualisation related
sett.show.mx_subjects      = 8;
sett.show.figname_subjects = '(SPM) Show subjects';
sett.show.figname_model    = '(SPM) Show model';
end
%==========================================================================