clear all ;close all; clc
% Author: Grant Lu
format("default")
%% INPUTS
% MTOW = Max Take-Off Weight, NLW = Normal Landing Weight
MTOW = 56185; % kg <- From weight estimation
NLW = 0.8 * MTOW; % Typically 80%

% Initial W/S Estimate Range
WS_range = [120, 160];

% Estimated Designed CLstar range
CLstar_range = [0.32 0.4];
Altitude = 35000; % ft
Mach = 0.78;

% Landing KTAS
LandSpeed = 140;
CL_max = [1.8 3];
stallfactor = 1.3; % minimum 1.3

%% CALCULATIONS
% DONT TOUCH BELOW FOR NOW
% Question about the Script, ask/message Grant
GW = 0:0.1:160000;

% Bound by WS range

WS_lower = 1/WS_range(1) * GW;
WS_upper = 1/WS_range(2) * GW;

figure(1)
grid on
hold on
plot(GW,WS_lower,GW,WS_upper);
xlabel('Gross Weight (lb)')
ylabel('Wing Area (sq ft)')
xline(MTOW*2.20462)

str1 = sprintf('W/S = %.f Limits',WS_range(1));
str2 = sprintf('W/S = %.f Limits',WS_range(2));

S2 = 1/WS_range(2) * MTOW*2.20462;
S1 = 1/WS_range(1) * MTOW*2.20462;

txtS2 = '  S2';
text(MTOW*2.20462,S2,txtS2)
txtS1 = '  S1';
text(MTOW*2.20462,S1,txtS1)


plot(MTOW*2.20462,S2,'o',MTOW*2.20462,S1,'o')

fprintf('\n')
fprintf('Wing Area Limit (S1): %10.4f ft2\n', S2)
fprintf('Wing Area Limit (S2): %10.4f ft2\n', S1)

T =@(h) 59 - 0.00356 * h; % F

gamma = 1.4;
P = @(h) 2116 * ((T(h) + 459.7)/518.6) ^ 5.256; % lbs/sq ft
q = @(M,h) gamma/2 .* P(h) .* M.^2; % lbs/sq ft (NEED CHECK)

% Bound by CL Star range

S3 = MTOW / ( CLstar_range(1) * q(Mach,Altitude));
S4 = MTOW / ( CLstar_range(2) * q(Mach,Altitude));

plot(MTOW*2.20462,S3,'o',MTOW*2.20462,S4,'o')

txtS3 = '  S3';
text(MTOW*2.20462,S3,txtS3)
txtS4 = '  S4';
text(MTOW*2.20462,S4,txtS4)

fprintf('Wing Area Limit (S3): %10.4f ft2\n', S3)
fprintf('Wing Area Limit (S4): %10.4f ft2\n', S4)

% Bound by Landing CL max
CLstall = CL_max;
CLapp = CLstall./stallfactor^2;
qapp = 53.2; % lbf/ft2 (from lecture, assumed SSL)

rhoSSL = 1.225; % kg/m3
q_SSL = 1/2 * rhoSSL * 0.00194 *(LandSpeed * 1.68781)^2; %lbs/sq ft (NEED CHECK)

S5 = NLW/(CLapp(1) * q_SSL);
S6 = NLW/(CLapp(2) * q_SSL);

xline(NLW*2.20462)
plot(NLW*2.20462,S5,'o',NLW*2.20462,S6,'o')

txtS5 = '  S5';
text(NLW*2.20462,S5,txtS5)
txtS6 = '  S6';
text(NLW*2.20462,S6,txtS6)

fprintf('Wing Area Limit (S5): %10.4f ft2\n', S5)
fprintf('Wing Area Limit (S6): %10.4f ft2\n', S6)

legend(str1,str2,'MTOW','S1','S2','S3','S4','NLW','S5','S6','Location','eastoutside')



