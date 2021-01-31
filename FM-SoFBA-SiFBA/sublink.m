function [FTarea,yo] = sublink(hour,daylength,sunrise,yi)

global v

% Photoperiod controls
ph = zeros(1, 4);
ph(1) = 24;         % Period
ph(2) = sunrise;    % Dawn
ph(3) = daylength;  % Photoperiod
ph(4) = 1;          % Amplitude

% Solvers options
accurate = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); 

% Solving ODE
[t2, y2] = ode15s(@flowering4, [hour-1 hour], yi, accurate, v, ph);

%output
FTarea = trapz(t2,y2(:,15));
yo = y2(end,:);
clear t2
clear y2

