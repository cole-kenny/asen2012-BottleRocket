% --------------------
% ID 109
% ASEN 2012-001
% Lab 2

% Purpose:
% erve as the defining set of differential equations to be solved during
% Phase 3 of rocket flight

% Inputs: t,z; global variables
% Outputs: dPhase1dt = [accZ,accX,velZ,velX]';

% Last modified: 12/4 - ID 109 - Initial Release
% Last modified: 12/6 - ID 109 - Final Release. All requirments met
% --------------------

function dPhase3dt = PhaseThreeODE(t,z)
% Sets up the system of ODEs for the third phase
% Rocket is ballistic

% Pull global variables
global g gamma pAtm rhoWater rhoAir cDrag cDischarge nozzleArea pAirInitial volBottle volWaterInitial volAirInitial theta railLength p0 sectionBottle mAirInitial mBottleEmpty R

% Setting velocity in to be the inegral of acceleration
velZ = z(1);
velX = z(2);
vel = sqrt((velX^2)+(velZ^2));

mass = mBottleEmpty;

% Heading is determined by velocity
headingZ = velZ/vel;
headingX = velX/vel;

% Equations for thrust
forceThrust = 0;

% Equation for drag
forceDrag = -(rhoAir/2)*(vel^2)*cDrag*sectionBottle;
forceZDrag = -(rhoAir/2)*(velZ^2)/((cos(atan(headingZ/headingX)))^2)*cDrag*sectionBottle;
forceXDrag = -(rhoAir/2)*(velX^2)/((cos(atan(headingZ/headingX)))^2)*cDrag*sectionBottle;

% Differential equation for rocket acceleration in Z
%accZ = forceZ/mass;
accZ = ((forceThrust+forceDrag)*headingZ/mass) - g;

% Differential equation for rocket acceleration in X
%accX = forceX/mass;
accX = (forceThrust+forceDrag)*headingX/mass;

% Combining in column vector for passing into ode45
dPhase3dt = [accZ,accX,velZ,velX]';

end