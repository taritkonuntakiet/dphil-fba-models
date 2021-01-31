function [FTarea,yn] = link(daylength,sunrise)

global v

% Photoperiod controls
ph = zeros(1, 4);
ph(1) = 24;         % Period
ph(2) = sunrise;    % Dawn
ph(4) = 1;          % Amplitude

% Solvers options
ystart = ones(1, 15);
entrainment = 720; 
standard = [];
accurate = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); 
NMax = 1000;     % Maximum number of points in the limit cycle
bvopt = bvpset('BCJacobian', @flowering4bcjac, 'NMax', NMax);


% Solving ODE
if daylength <= 16
 
    ph(3) = daylength;

    % Entrainment
    [t0, y0] = ode15s(@flowering4, [0 entrainment], ystart, standard, v, ph);

    % Cycle used as seed
    [t1, y1] = ode15s(@flowering4, [0 24], y0(end, :), accurate, v, ph);

    % Searching for the limit cycle
    rinit.x = t1;
    rinit.y = y1';
    try
        soldl = bvp4c(@flowering4, @flowering4bc, rinit, bvopt, v, ph);
    catch
        soldl = rinit;
    end

    FTarea = trapz(t1,y1(:,15));
    yn = y1(end,:);
    
    
    
else
    ph(3) = 16;
    % Entrainment
    [t0, y0] = ode15s(@flowering4, [0 entrainment], ystart, standard, v, ph);
    % Cycle used as seed
    [t1, y1] = ode15s(@flowering4, [0 24], y0(end, :), accurate, v, ph);
    
    % Transfer to
    ph(3) = daylength;
    % Phase within the limit cycle
    LCphase = 0.8616*ph(3) + 3.2659; % for wt
      
    
    % Transferred for 3 days
    for i=1:3
        
        [t1, y1] = ode15s(@flowering4, [0 24], y1(end, :), accurate, v, ph);
        [Amp(i),Index(i)] = max(y1(:,5));
        Phase(i) = t1(Index(i));
        ErrorP(i) = (LCphase-Phase(i))^2;
        FTArea (i) = trapz(t1,y1(:,15));
        y(:,i) = y1(end,:);            
    end
    
    % Selecting the one closest to the limit cycle
    [MinP,Ind] = min(ErrorP);
    FTarea = FTArea(Ind);
    yn = y(:,Ind);
           
end 

%To run on an hourly basis:
start = 0;
for h = 1:24
    
    [t2, y2] = ode15s(@flowering4, [start h], yn, accurate, v, ph);
    FTarea(h) = trapz(t2,y2(:,15));
    yn = y2(end,:);
    clear t2
    clear y2
    start = start+1;
    
end    
        

