function dydt = flowering4(t, y, v, ph)
% 
%

% Named parameters
nu = v(54);
mu = v(55);
xi = v(56);
ylhy = v(57);
ytoc = v(58);

% Parameter values for CO protein

vCOm  = v(59);
vCOp  = v(60);
kCOp  = v(61);

%  Parameter values for FT
VCO  = v(62);
KCO  = v(63);
vFT  = v(64);
KFT  = v(65);
Bco = v(66);

% Photoperiod control
t = mod(t, ph(1));
light = (1 + tanh((t - ph(2))*6)) * (1 - tanh((t - (ph(2)+ph(3)))*6))/4;
dark = 1 - light;
force = ph(4) * light;


% RHS of the equations
dydt = [    % dydt(1) = Protein P levels which activates LHY
0.5*dark - light*y(1) - 1.2*y(1)/(1.2 + y(1));
	    
            % dydt(2) = LHY mRNA, activated by P and X
force*v(1)*y(1) + v(2)*y(10)^mu/(v(3) + y(10)^mu) - v(4)*y(2)/(v(5) + ...
						  y(2));
	    
            % dydt(3) = LHY Protein in Cytoplasm
v(6)*y(2) - v(7)*y(3) + v(8)*y(4) - v(9)*y(3)/(v(10) + y(3));
	    
            % dydt(4) = LHY Protein in Nucleus
v(7)*y(3) - v(8)*y(4) - v(11)*y(4)/(v(12) + y(4));
	    
            % dydt(5) = TOC1 mRNA, activated by Y and repressed by LHY
v(13)*y(13)^nu/(v(14) + y(13)^nu) * v(15)/(y(4)^nu + v(16)) - ...
v(17)*y(5)/(v(18) + y(5));
	    
            % dydt(6) = TOC1 Protein in Cytoplasm, extra-degraded by darkness
v(19)*y(5) - v(20)*y(6) + v(21)*y(7) - (v(22)*dark + v(23))*y(6)/(v(24) ...
						  + y(6));
	    
            % dydt(7) = TOC1 Protein in Nucleus, extra-degraded by darkness
v(20)*y(6) - v(21)*y(7) - (v(25)*dark + v(26))*y(7)/(v(27) + y(7));
	    
            % dydt(8) = Factor X mRNA, activated by TOC1
v(28)*y(7)^xi/(v(29) + y(7)^xi)- v(30)*y(8)/(v(31) + y(8));
	    
            % dydt(9) = X Protein in Cytoplasm
v(32)*y(8) - v(33)*y(9) + v(34)*y(10) - v(35)*y(9)/(v(36) + y(9));
	    
            % dydt(10) = X Protein in Nucleus
v(33)*y(9) - v(34)*y(10) - v(37)*y(10)/(v(38) + y(10));
	    
            % dydt(11) = Factor Y mRNA, activated by P and light,
            % represed by TOC1 and LHY
(force*v(53)*y(1) + (force*v(39) + v(40))/(v(41) + y(7)^ytoc)) * ...
v(51)/(y(4)^ylhy + v(52)) - v(42)*y(11)/(v(43) + y(11));
	    
            % dydt(12) = Y Protein in Cytoplasm
v(44)*y(11) - v(45)*y(12) + v(46)*y(13) - v(47)*y(12)/(v(48) + y(12));
	    
            % dydt(13) = Y Protein in Nucleus
v(45)*y(12) - v(46)*y(13) - v(49)*y(13)/(v(50) + y(13));
	    
	    % dydt(14) = CO protein, translated from nuclear TOC1 = mCO
	    % extra-degraded by darkness
vCOm*y(7) - dark*vCOp*y(14)/(kCOp + y(14)); 

        % dydt(15) = FT mRNA, which is activated by CO protein
(Bco + (VCO*y(14)/(KCO + y(14)))) - vFT*y(15)/(KFT + y(15));
];
