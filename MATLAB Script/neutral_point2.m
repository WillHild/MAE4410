%neutral point - Toby's notes

clear all; clc; close all; 

kf = 0.01; %empirical relation dependent on %% fuselage length 
wf = 3.60; %width of fuselage
lf = 30; %fuselage length 
CLa = 0.5383; %CL due to 2D airfoil of wing
lp = 5.07; %distance from CG to propulsion system
c = 3.168; %MAC 
ARw = 12.6277; %wing aspect ratio
ARh = 5.94776; %stabilitser aspect ratio
ARt = 8.76924; %truss aspect ratio
sweep_w = 27 ; %sweep angle of wing
sweep_h = 30; %sweep angle for horizontal stabiliser
sweep_t = 19; %sweep angle for truss
func = 1; %empirical value for transforming from elliptical to current shape of wing
func_h = 1;

lh = 18.27; %length from cg to HS 
Sh = 22.235; %stabiliser area
Sw = 102.6; %wing area
St = 45.8216; %truss area
CLah = 0.694; %CL due to horizontal stabiliter 2d aerofoil
CLa_t2 = 1.1649; %CL due to truss 2d aerofoil
lt = 0.0126; %distance from MAC of truss to aircraft MAC 
x_ac = 2.448; 
eta = 0.85; %efficiency of horizontal stabiliser 

Cmaf = (kf*wf^2*lf)/(Sw*c); %per degree: pitching moment of fuselage
Cmap = 0.02*CLa*lp/c; %moment coefficient due to propulsion system
CLa_w = (CLa*cosd(sweep_w))/(sqrt(1 + ((cosd(sweep_w))/(pi*ARw))^2) + (CLa*cosd(sweep_w))/(pi*ARw));

CLaw3 = func*CLa_w; %lifit coefficient for 3d wing 

Cla_h = (CLah*cosd(sweep_h))/(sqrt(1 + ((cosd(sweep_h))/(pi*ARh))^2) + (CLah*cosd(sweep_h))/(pi*ARh));

CLah3 = func_h*Cla_h;

CLa_t = (CLa_t2*cosd(sweep_t))/(sqrt(1 + ((cosd(sweep_t))/(pi*ARt))^2) + (CLa_t2*cosd(sweep_t))/(pi*ARt));

deda = (2*CLaw3)/(pi*ARw);
V_h = (lh*Sh)/(c*Sw);

xnp = c*(x_ac/c - Cmaf/CLaw3 - Cmap/CLaw3 - (eta*V_h*CLah3)*(1-deda)/CLaw3 - (CLa_t/CLaw3)*lt*St/Sw);

xnp_nose = 17.6 + 4.49197 + xnp;
xcg_nose = 19.643;

sm = xnp_nose/c - xcg_nose/c;