% --------------------
% ID 109
% ASEN 2012-001
% Lab 2

% Purpose:
% erve as the defining set of differential equations to be solved during
% Phase 2 of rocket flight

% Inputs: t,z; global variables
% Outputs: dPhase1dt = [mDotRocket,accZ,accX,velZ,velX]';

% Last modified: 12/4 - ID 109 - Initial Release
% Last modified: 12/6 - ID 109 - Final Release. All requirments met
% --------------------

function dPhase2dt = PhaseTwoODE(t,z)
% Sets up the system of ODEs for the second phase
% Rocket is thrusting with the air

% Pull global variables
global g gamma pAtm rhoWater rhoAir cDrag cDischarge nozzleArea pAirInitial volBottle volWaterInitial volAirInitial theta railLength p0 sectionBottle mAirInitial mBottleEmpty R

% Setting velocity in to be the inegral of acceleration
velZ = z(2);
velX = z(3);
vel = sqrt((velX^2)+(velZ^2));

% Setting rocket mass to be the integral of mDotRocket
mass = z(1);

% Equation for determining the air pressure in the bottle
pEnd = p0*((volAirInitial/volBottle)^gamma);
mAir = mass - mBottleEmpty;
p = ((mAir/mAirInitial)^gamma)*pEnd;

% Defining gas properties
rho = mAir/volBottle;
Temp = p/(rho*R);

% Heading is determined by velocity
headingZ = velZ/vel;
headingX = velX/vel;

% Calculate critical pressure
pCrit = p*((2/(gamma+1))^(gamma/(gamma-1)));

% Checking if the flow is choked
if pCrit > pAtm
    % Flow is choked
    % Computing the temperature at nozzle exit
    Te = (2/(gamma+1))*Temp;
    % Equation for determining exit velocity with a choked flow
    vE = sqrt(gamma*R*Te);
    % Defining pExit
    pExit = pCrit;
    % Defining rhoExit
    rhoExit = pExit/(R*Te);
else
    % The flow is not choked
    % Computing the exit Mach number
    machExit = sqrt((((p/pAtm)^((gamma-1)/gamma))-1)*(2/(gamma-1)));
    % Computing temperature at nozzle exit
    Te = Temp/(1+((gamma-1)/2)*(machExit^2));
    % Equation for determining exit velocity in a non-choked flow
    vE = machExit*sqrt(gamma*R*Te);
    % Defining pExit
    pExit = pAtm;
    % Defining rhoExit
    rhoExit = pAtm/(R*Te);
end

% Defining mDotAir
mDotAir = cDischarge*rhoExit*nozzleArea*vE;

% Equation for thrust from pressurized air
forceThrust = (mDotAir*vE)+((pExit-pAtm)*nozzleArea);
if forceThrust < 0
    forceThrust = 0;
end
forceZThrust = forceThrust*cos(atan(headingZ/headingX));
forceXThrust = forceThrust*sin(atan(headingZ/headingX));

% Equation for drag
forceDrag = -(rhoAir/2)*(vel^2)*cDrag*sectionBottle;
forceZDrag = -(rhoAir/2)*(velZ^2)/((cos(atan(headingZ/headingX)))^2)*cDrag*sectionBottle;
forceXDrag = -(rhoAir/2)*(velX^2)/((cos(atan(headingZ/headingX)))^2)*cDrag*sectionBottle;

% Defining the differential equation for mDot
mDotRocket = -mDotAir;

% Differential equation for rocket acceleration in Z
%accZ = forceZ/mass;
accZ = ((forceThrust+forceDrag)*headingZ/mass) - g;

% Differential equation for rocket acceleration in X
%accX = forceX/mass;
accX = (forceThrust-forceDrag)*headingX/mass;

% Combining in column vector for passing into ode45
dPhase2dt = [mDotRocket,accZ,accX,velZ,velX]';

end