% --------------------
% ID 109
% ASEN 2012-001
% Lab 2

% Purpose:
% erve as the defining set of differential equations to be solved during
% Phase 1 of rocket flight

% Inputs: t,z; global variables
% Outputs: dPhase1dt = [mDotRocket,accZ,accX,velZ,velX,volDotAir]';

% Last modified: 12/3 - ID 109 - Initial Release
% Last modified: 12/6 - ID 109 - Final Release. All requirments met
% --------------------

function dPhase1dt = PhaseOneODE(t,z)
% Sets up the system of ODEs for the first phase
% Rocket is thrusting with the water

% Pull global variables
global g gamma pAtm rhoWater rhoAir cDrag cDischarge nozzleArea volBottle volWaterInitial volAirInitial theta railLength p0 sectionBottle

% Determining the position of the rocket
posZ = z(4);
posX = z(5);

% Setting velocity in to be the inegral of acceleration
velZ = z(2);
velX = z(3);
vel = sqrt((velX^2)+(velZ^2));

% Setting rocket mass to be the integral of mDotRocket
mass = z(1);

% Checking if the heading should change or be constant
if railLength*cos(theta) > posZ
    % Rocket is still on the rail
%     headingZ = railLength;
%     headingX = railLength;
    headingZ = sin(theta);
    headingX = cos(theta);

else
    % Heading is determined by velocity
%     headingZ = velZ/(sqrt(velX^2 + velZ^2));
%     headingX = velX/(sqrt(velX^2 + velZ^2));
    headingZ = velZ/vel;
    headingX = velX/vel;
end

% Setting volume of air to the integral of volDot
volAir = z(6);

% Equation for determining the air pressure in the bottle
p = ((volAirInitial/volAir)^gamma)*p0;

% Equation for determining the exit velocity of water
vE = sqrt((2*(p-pAtm))/(rhoWater));

% Equation for thrust from water
forceThrust = (2*cDischarge*(p-pAtm)*nozzleArea);
forceZThrust = forceThrust*cos(atan(headingZ/headingX));
forceXThrust = forceThrust*sin(atan(headingZ/headingX));

% Equation for drag
forceDrag = -(rhoAir/2)*(vel^2)*cDrag*sectionBottle;
forceZDrag = -(rhoAir/2)*(velZ^2)/((cos(atan(headingZ/headingX)))^2)*cDrag*sectionBottle;
forceXDrag = -(rhoAir/2)*(velX^2)/((cos(atan(headingZ/headingX)))^2)*cDrag*sectionBottle;

% Defining the differential equation for mDot
mDotRocket = -rhoWater*cDischarge*nozzleArea*vE;

% Differential equation for rocket acceleration in Z
%accZ = forceZ/mass;
accZ = ((forceThrust+forceDrag)*headingZ/mass) - g;

% Differential equation for rocket acceleration in X
%accX = forceX/mass;
accX = (forceThrust+forceDrag)*headingX/mass;

% Differential equation for volDot
volDotAir = cDischarge*nozzleArea*sqrt((2/rhoWater)*(p0*((volAirInitial/volAir)^gamma)-pAtm));

% Combining in column vector for passing into ode45
dPhase1dt = [mDotRocket,accZ,accX,velZ,velX,volDotAir]';

end