% MODEL EXERCISEFLAT
% Spatial redistribution of surface water.
% Stefan Dekker, Willem Bouten, Maarten Boerlijst en Max Rietkerk.
% Rietkerk et al. 2002. Self-organization of vegetation in arid ecosystems.
% The American Naturalist 160(4): 524-530.
% DeltaX, DeltaY, DifP, DifW and DifO have unrealistic values that differ from original
% publication to increase computation speed for educational purposes.

clear all

% System discretisation (see remark above)
DeltaX=101; % (m)
DeltaY=101; % (m)

% Diffusion constants for plants, soil water and surface water (see remark above)
DifP=100;   % (m2.d-1)
DifW=1;     % (m2.d-1)
DifO=1000;  % (m2.d-1)

% Initial fraction of grid cells with bare ground 
frac=0.90;  % (-)

% Parameter values
R       =   1.3;  % Rainfall (mm.d-1)
alpha   =   0.1;  % Proportion of surface water available for infiltration (d-1)
W0      =  0.15;  % Bare soil infiltration (-)
rw      =   0.1;  % Soil water loss rate due to seepage and evaporation (d-1)
c       =    10;  % Plant uptake constant (g.mm-1.m-2)
gmax    =  0.05;  % Plant growth constant (mm.g-1.m-2.d-1)   
d       =   0.4;  % Plant senescence rate and grazing rate (d-1)
k1      =     3;  % Half saturation constant for plant uptake and growth (mm)
k2      =     5;  % Half saturation constant for water infiltration (g.m-2)

% Number of grid cells
m=250;
NX=m;
NY=m;

% Timesteps
dT=1;     % timestep
Time=1;      % begin time       
EndTime=1000;    % end time
PlotStep=10; % (d)
PlotTime=PlotStep; % (d)

% Initialisation
popP = zeros(m,m);                
popW = zeros(m,m);                
popO = zeros(m,m);
dP=zeros(m,m);
dO=zeros(m,m);
dW=zeros(m,m);
NetP=zeros(m,m);
NetW=zeros(m,m);
NetO=zeros(m,m);

%Boundary conditions
FYP = zeros(NY+1,NX);   	% bound.con. no flow in/out to Y-direction	
FXP = zeros(NY,NX+1);		% bound.con. no flow in/out to X-direction
FYW = zeros(NY+1,NX);   	% bound.con. no flow in/out to Y-direction	
FXW = zeros(NY,NX+1);		% bound.con. no flow in/out to X-direction
FYO = zeros(NY+1,NX);   	% bound.con. no flow in/out to Y-direction	
FXO = zeros(NY,NX+1);		% bound.con. no flow in/out to X-direction

% Initial state
for i=1:m,
  for j=1:m,
    if (rand > frac)
      popO(i,j)=R/(alpha*W0); % Homogeneous equilibrium surface water in absence of plants
      popW(i,j)=R/rw; % Homogeneous equilibrium soil water in absence of plants
      popP(i,j)=90; % Initial plant biomass
    else
      popO(i,j)=R/(alpha*W0); % Homogeneous equilibrium surface water in absence of plants
      popW(i,j)=R/rw; % Homogeneous equilibrium soil water in absence of plants
      popP(i,j)=0; % Initial plant biomass
    end;
  end;
end;

% Timesteps
while Time<=EndTime,    
 
% Reaction
dO = R-alpha*popO*((popP+k2+W0)/(popP+k2);
dW = alpha*popO*((popP+k2+W0)/(popP+k2)-gmax*((W/(W+k1))*popP-rw*W;
dP = c*gmax*((W/(W+k1))*popP-d*popP;

% Diffusion
% calculate Flow in x-direction : Flow = -D * dpopP/dx;
FXP;	
FXW;	
FXO(1:i,2:j) = -DifO * (popO(1:i,2:j)- popO(1:i,1:j-1))/DeltaX;

% calculate Flow in y-direction: Flow = -D * dpopP/dy;
FYP;
FYW;
FYO;

% calculate netflow
NetP;	
NetW;	
NetO;	

% Update
   popW;
   popO;
   popP;

   Time=Time+dT;
   
   PlotTime=PlotTime-dT;
   if PlotTime<=0,
       imagesc (popP); title 'vegetation' 
       caxis([0 15]);
       colorbar
       drawnow;
       PlotTime=PlotStep;
   end
end
