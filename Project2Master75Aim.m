% --------------------
% ID 109
% ASEN 2012-001
% Lab 2

% Purpose:
% Model a water bottle rocket through 3 phases of flight. This lab will
% prepare students for a lab assigned in ASEN 2004 where we will expand on
% this code and actually build our model rockets.

% Inputs: None
% Outputs: No explicit outputs. Plots posX v. posZ and thrust vs. time.

% Last modified: 11/26 - ID 109 - Initial Release
% Last modified: 12/3 - ID 109 - Created PhaseOneODE
% Last modified: 12/4 - ID 109 - Created PhaseTwoODE and PhaseThreeODE
% Last modified: 12/5 - ID 109 - Sorted out returns from ODE function calls
% Last modified: 12/6 - ID 109 - Final Release. All requirments met
% --------------------

clear; close all; clc;

%% Defining variables and givens
global g gamma pAtm rhoWater rhoAir cDrag cDischarge nozzleArea pAirInitial volBottle volWaterInitial volAirInitial theta railLength p0 sectionBottle mAirInitial mBottleEmpty mWaterInitial R T

g = 9.81;                                       % m/s^2
gamma = 1.4;
pAtm = 83426.56;                                % Pa
rhoWater = 1000;                                % kg/m^3
rhoAir = 0.961;                                 % kg/m^3
cDrag = 0.38;
cDischarge = 0.8;
nozzleArea = ((2.1^2)/4)*pi/(100^2);            % m^2
pAirInitial = 400000;                           % Pa (gage)
volBottle = 0.002;                              % m^3
volWaterInitial = 0.001;                        % m^3
volAirInitial = volBottle-volWaterInitial;      % m^3
theta = pi/4;                                   % radians
railLength = 0.5;                               % m
p0 = pAtm + pAirInitial;                        % Pa (abs)
sectionBottle = ((10.5^2)/4)*pi/(100^2);        % m^2
mBottleEmpty = 0.15;                            % kg
mWaterInitial = volWaterInitial*rhoWater;       % kg
R = 287;                                        % J/kg*K
T = 300;                                        % Kelvin
mAirInitial = (p0*volAirInitial)/(R*T);         % kg   

mRocketInitial = mBottleEmpty+mWaterInitial+mAirInitial;    % kg
tSpanFull = [0 5];                              % s

%% Phase 1 (water thrust phase)
% Clearing integration variables and setting new initial conditions
clear t z;
mRocketI = mRocketInitial;
velZI = 0;
velXI = 0;
posZI = 0.01;
posXI = 0;
volAirI = volAirInitial;
tSpan = [tSpanFull(1) 1];

% Calling ode45
% output/IC order = [mRocket,velZ,velX,posZ,posX,volAir];
[t,z] = ode45('PhaseOneODE',tSpan,[mRocketI velZI velXI posZI posXI volAirI]);

% Sorting ode45 outputs
tP1 = t;
mRocketP1 = z(:,1);
velZP1 = z(:,2);
velXP1 = z(:,3);
posZP1 = z(:,4);
posXP1 = z(:,5);
volAirP1 = z(:,6);

% Finding where phase 1 is no longer applicible (water is exhausted)
%appTempIndex = volAirP1 < volBottle;
appTempIndex = mRocketP1 > mRocketInitial - mWaterInitial;
lastP1Index = sum(appTempIndex)+1;

% Removing data from phase 1 output that isn't valid
tP1(lastP1Index:end) = [];
mRocketP1(lastP1Index:end) = [];
velZP1(lastP1Index:end) = [];
velXP1(lastP1Index:end) = [];
posZP1(lastP1Index:end) = [];
posXP1(lastP1Index:end) = [];
volAirP1(lastP1Index:end) = [];

% Computing thrust from Phase 1
thrustP1 = zeros(length(tP1),1);
pressureP1 = zeros(length(tP1),1);
for i = 1:length(tP1)
    % Recomputing the pressure
    p = ((volAirInitial/volAirP1(i))^gamma)*p0;
    % Calculating the thrust
    forceThrust = (2*cDischarge*(p-pAtm)*nozzleArea);
    % Saving to a P1 thrust vector
    thrustP1(i) = forceThrust;
    pressureP1(i) = p;
end

%% Phase 2 (air thrust phase)
% Clearing integration variables and setting new initial conditions
clear t z tSpan appTempIndex;
mRocketI = mRocketInitial - mWaterInitial;
velZI = velZP1(end);
velXI = velXP1(end);
posZI = posZP1(end);
posXI = posXP1(end);
tSpan = [tP1(end) 0.6];

% Calling ode45
% output/IC order = [mRocket,velZ,velX,posZ,posX];
[t,z] = ode45('PhaseTwoODE',tSpan,[mRocketI velZI velXI posZI posXI]);

% Sorting ode45 outputs
tP2 = t;
mRocketP2 = z(:,1);
velZP2 = z(:,2);
velXP2 = z(:,3);
posZP2 = z(:,4);
posXP2 = z(:,5);

% Finding where phase 2 is no longer applicible (air pressure is exhausted)
appTempIndex = mRocketP2 > mBottleEmpty;
lastP2Index = sum(appTempIndex)+1;

% Removing data from phase 2 output that isn't valid
tP2(lastP2Index:end) = [];
mRocketP2(lastP2Index:end) = [];
velZP2(lastP2Index:end) = [];
velXP2(lastP2Index:end) = [];
posZP2(lastP2Index:end) = [];
posXP2(lastP2Index:end) = [];

% Computing thrust from Phase 2
thrustP2 = zeros(length(tP2),1);
pressureP2 = zeros(length(tP2),1);
for i = 1:length(tP2)
    % Recomputing the pressure
    pEnd = p0*((volAirInitial/volBottle)^gamma);
    mAir = mRocketP2(i) - mBottleEmpty;
    p = ((mAir/mAirInitial)^gamma)*pEnd;
    % Redefining gas parameters
    rho = mAir/volBottle;
    Temp = p/(rho*R);
    % Choked flow check
    pCrit = p*((2/(gamma+1))^(gamma/(gamma-1)));
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
    mDotAir = cDischarge*rhoExit*nozzleArea*vE;
    % Calculating the thrust
    forceThrust = (mDotAir*vE)+(pExit-pAtm)*nozzleArea;
    if forceThrust < 0
        forceThrust = 0;
    end
    % Saving to a P2 thrust vector
    thrustP2(i) = forceThrust;
    pressureP2(i) = p;
end

%% Phase 3 (ballistic phase)
% Clearing integration variables and setting new initial conditions
clear t z tSpan appTempIndex;
velZI = real(velZP2(end));
velXI = real(velXP2(end));
posZI = real(posZP2(end));
posXI = real(posXP2(end));
tSpan = [tP2(end) tSpanFull(2)];

% Calling ode45
% output/IC order = [velZ,velX,posZ,posX];
[t,z] = ode45('PhaseThreeODE',tSpan,[velZI velXI posZI posXI]);

% Sorting ode45 outputs
tP3 = t;
velZP3 = z(:,1);
velXP3 = z(:,2);
posZP3 = z(:,3);
posXP3 = z(:,4);

% Finding where phase 2 is no longer applicible (air pressure is exhausted)
appTempIndex = posZP3 > 0;
lastP3Index = sum(appTempIndex)+5;

% Removing data from phase 2 output that isn't valid
tP3(lastP3Index:end) = [];
velZP3(lastP3Index:end) = [];
velXP3(lastP3Index:end) = [];
posZP3(lastP3Index:end) = [];
posXP3(lastP3Index:end) = [];

% Thrust during ballistic phase is 0
thrustP3 = zeros(length(tP3(1:10)),1);

%% Plot: Position
posZFull = vertcat(posZP1,posZP2,posZP3);
posXFull = vertcat(posXP1,posXP2,posXP3);

figure
hold on
axis([0 80 0 30])
title('Rocket Altitude v. Horizontal Distance')
xlabel('Ground Distance [m]')
ylabel('Altitude [m]')
plot(posXP1,posZP1,'LineWidth',2)
plot(posXP2,posZP2,'LineWidth',2)
plot(posXP3,posZP3,'LineWidth',2)
legend('Phase 1','Phase 2','Phase 3')


%% Plot: Thrust

figure
hold on
axis([0 0.45 0 200])
title('Rocket Thrust vs. Time')
xlabel('Time [s]')
ylabel('Thrust [N]')
tFull = [tP1;tP2;tP3(1:10)];
thrustFull = [thrustP1; thrustP2; thrustP3];
plot(tFull,thrustFull,'LineWidth',2)