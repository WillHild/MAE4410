%tank weights and volume
%Author: Leon Phillips

clear all; close all; clc; 

%% calculating tank volume 
M_h = 1.05*5847.5281 + 200; %mass of hydrogen - this is a value calculated from initial estimates 
%M_h = 1660;
Vi = 0.072; %7.2% extra volume to maintain const pressure and allow for boil-off
rho_H = 72; %kg/m^3 is the density of LH2

Vt = (M_h*(1 + Vi))/rho_H; %volume required for tanks

L = 3; %tank height 

diff = 1; %initialise 
r = 1; %initial radius guess
k = 0.01; %increment to perform secant method

num_tanks = 4;

Vt = Vt/num_tanks;

%iterate to find radius required of cylindrical tank 
while diff > 0.000001 

%     Vt_guess1 = (4*pi*r^3)/3 + pi*r^2*L - Vt; %hemisphere caps
%     Vt_guess1 = pi*r^2*L + 2*pi*r^2 - Vt; % exact cylinder caps
%ellipsoid caps were chosen with a ratio, a/b = 1.6 from literature, where
% a = c = r, and V_e = 4/3*pi*a*b*c
    Vt_guess1 = (4/(3*1.6))*pi*r^3 + pi*r^2*L - Vt; %ellipsoid caps 


%     Vt_guess2 = (4*pi*(r+k)^3)/3 + pi*(r+k)^2*L - Vt;
%     Vt_guess2 =  pi*(r+k)^2*L + 2*pi*(r+k)^2 - Vt;
    Vt_guess2 = (4/(3*1.6))*pi*(r+k)^3 + pi*(r+k)^2*L - Vt; %ellipsoid caps

    r = r+k - ((Vt_guess2)*(r + k - (r)))/(Vt_guess2 - Vt_guess1); %secant method
    
%     diff = abs((4*pi*r^3)/3 + pi*r^2*L - Vt); %absolute difference
%     diff = abs(pi*r^2*L + 2*pi*r^2 - Vt); %absolute difference
    diff = abs((4/(3*1.6))*pi*r^3 + pi*r^2*L - Vt); %absolute difference


end

P_t = 1.4*10^5; %Pa - pressure of tank according to various sources
%P_a = 75.3*10^3; %Pa - external pressure outside of tank, equivalent to 8000 feet as per regulations 
P_a = 29.647*10^3; %Pa - external pressure at 30k feet
%P_a = 101325;

FoS = 5;

oy = 276*10^6; %Pa - yield strength of Al6061-T6

t = (FoS*(P_t-P_a)*r)/oy; %hoop stress calculation

rho_a = 2710; %kg/m^3 

%mass of the aluminium for one tank
% m_t = (rho_a*(4/3)*pi*(r + t)^3 + rho_a*pi*(r + t)^2*L - rho_a*Vt); %for
% hemsiphere tank
% m_t = (rho_a*pi*(r + t)^2*L + 2*rho_a*pi*(r + t)^2 - rho_a*Vt);
m_t = (rho_a*(4/(3*1.6))*pi*(r + t)^3 + rho_a*pi*(r + t)^2*L - rho_a*Vt); %for an ellipsoid 

%% Now, need to calculate adequate thickness of the tank and resultant mass
%calculating surface temperature of the insulant 
T_ext = 298.15; %room temp 
T_LH2 = 13.15; %storage temp of LH2 

T_s1 = 283.15; %initial guess 
T_inc = 0.01; %temperature increment


diff_T = 1;

a = -3.119*10^(-6) + 3.541*10^(-8)*T_ext + 1.679*10^(-10)*T_ext^2; %constants
v = -2.079*10^(-6) + 2.777*10^(-8)*T_ext + 1.077*10^(-10)*T_ext^2; %constants

k_air = 25.87; %thermal conductivity of air
k_i = [0.0096,0.0112,0.0046,0.0064,0.00016,0.000017,0.00017]; %thermal conductivity of rigid closed cell polymethacrylimide
rho_lin = [35.3,32.1,49.8,64.2,40,120,160]; %kg/m^3 density of liner

L_i = linspace(0.001,0.045,1000); %m - thickness of insulating liner - change this to determine minimum thickness

%initialising 
m_hloss_store = zeros(length(L_i),length(rho_lin));
m_all_tank_store = zeros(length(L_i),length(rho_lin));
m_liner_store = zeros(length(L_i),length(rho_lin));
m_hloss_rate_store = zeros(length(L_i),length(rho_lin));


g = 9.81;

for j = 1:length(k_i) %iterate through all insulating materials
    for i = 1:length(L_i)
        T_s1 = 283.15; %initial guess
        
        while diff_T > 0.01 %iterate to find surface temp
            
            T_s2 = T_s1 + T_inc; %second initial guess
        
            %Rayleigh coeficcients 
            R1 = ((g/(T_ext))*(T_ext - T_s1)*(2*r)^3)/(v*a);
            R2 = ((g/(T_ext))*(T_ext - T_s2)*(2*r)^3)/(v*a);
        
            %Nusselt Numbers 
            N1 = ((0.60 + 0.387*R1^(1/6))/(1 + (0.559/(P_a*R1))^(9/16))^(8/27))^2; %note: if getting weird numbers r might be R1
            N2 = ((0.60 + 0.387*R2^(1/6))/(1 + (0.559/(P_a*R2))^(9/16))^(8/27))^2;
            
            %convection coefficients 
            h1 = N1*k_air/(2*r);
            h2 = N2*k_air/(2*r);
        
            %heat in 
            Qi1 = h1*(T_ext - T_s1);
            Qi2 = h2*(T_ext - T_s2);
        
            %heat out 
            Qo1 = (k_i(j)*(T_s1 - T_LH2))/L_i(i); 
            Qo2 = (k_i(j)*(T_s2 - T_LH2))/L_i(i);
        
            %Perform Secant Method
            dE1 = Qo1 - Qi1;
            dE2 = Qo2 - Qi2;
        
            T_s1 = T_s2 - ((dE2)*(T_s2 - (T_s1)))/(dE2 - dE1); %secant method
            
            diff_T = abs(dE1); %absolute difference
        
        end
        

%         surf_a = 2*pi*r*L + 4*pi*r^2*L; % hemisphere
        surf_a = 2*pi*r*L + 4*pi*(((r^2)^1.6 + ((r^2)/1.6)^1.6 + ((r^2)/1.6))/3)^(1/1.6); % ellipsoid

        h_fg = 446592; %J/kg - latent heat of vaporization of LH2
        m_hloss = (k_i(j)*surf_a*(T_s1 - T_LH2))/(L_i(i)*h_fg); %resultant mass rate lost of LH2
        
        a = 292.4;
        V = 0.78*a;
        range = 4630*10^3; %maximum range in metres 
        time = range/V + 30*60; %flight time 
        m_hloss_tot = m_hloss*time;
        
        %calculating mass of liner
        v_liner = (4/3)*pi*(r + t + L_i(i))^3 + pi*(r + t + L_i(i))^2*L - ((4/3)*pi*(r + t)^3 + pi*(r + t)^2*L);
        m_liner = v_liner*rho_lin(j);
        
        m_all_tank = num_tanks*(m_liner + m_t);
    
        %storing important values
        m_hloss_rate_store(i,j) = m_hloss;
        m_hloss_store(i,j) = m_hloss_tot;
        m_all_tank_store(i,j) = m_all_tank; 
        m_liner_store(i,j) = m_liner;
    
    end
end

diam_tot = (2*(r + t + L_i))';
% height_tot = (L + 2*(r + t + L_i))';
height_tot = (L + 2*(r/1.6 + t + L_i))';

% percentage rate loss by weight per hour 
loss_perc_hr = ((m_hloss_rate_store.*60.*60)./M_h)*100;

figure(1)
plot(L_i,m_hloss_store)
title('LH2 loss throughout flight')
xlabel('insulation thickness (m)')
ylabel('LH2 mass lost (kg)')
legend('Rigid closed cell polymethacrylimide','Rigid open cell polyurethane','Rigid closed cell polyvinalchloride','Rigid closed cell polyurethane and chopped glass fiber','Evacuated aluminium foil weperated with fluffy glass mats','Evacuated aluminium foil and glass paper laminate,','Evacuated silica powder')

figure(2)
plot(L_i,m_all_tank_store)
title('Total Tank Mass')
xlabel('insulation thickness (m)')
ylabel('total tank mass (kg)')
legend('Rigid closed cell polymethacrylimide','Rigid open cell polyurethane','Rigid closed cell polyvinalchloride','Rigid closed cell polyurethane and chopped glass fiber','Evacuated aluminium foil weperated with fluffy glass mats','Evacuated aluminium foil and glass paper laminate,','Evacuated silica powder')


figure(3) 
plot(L_i,m_all_tank_store + m_hloss_store)
title('Total Tank Mass and boil off mass')
xlabel('insulation thickness (m)')
ylabel('total tank and boil off mass (kg)')
legend('Rigid closed cell polymethacrylimide','Rigid open cell polyurethane','Rigid closed cell polyvinalchloride','Rigid closed cell polyurethane and chopped glass fiber','Evacuated aluminium foil weperated with fluffy glass mats','Evacuated aluminium foil and glass paper laminate,','Evacuated silica powder')

% Now we need to find the geometry and mass for each rigid foam material
% which allows a <0.5% boil off mass per 
mat_names = ["Rigid closed cell polymethacrylimide","Rigid open cell polyurethane","Rigid closed cell polyvinalchloride","Rigid closed cell polyurethane and chopped glass fiber","Evacuated aluminium foil weperated with fluffy glass mats","Evacuated aluminium foil and glass paper laminate","Evacuated silica powder"];

[~,col] = size(loss_perc_hr);
id_ins_config = zeros(1,col);
for j = 1:col %iterate through all materials 
    id_ins_config(j) = find(loss_perc_hr(:,j) < 0.5,1); %this will return first element which is <0.5
    
    %printing configurations for each material 
    fprintf('%s',mat_names(j))
    fprintf('\n Geometry: \n total diameter is: %3.3f m \n total height is %3.3f m',diam_tot(id_ins_config(j),1),height_tot(id_ins_config(j),1))
    fprintf('\n Masses: \n total mass of each tank: %3.3f kg \n total mass of all tanks: %3.3f kg \n total LH2 boiled off: %3.3f kg \n\n',m_all_tank_store(id_ins_config(j),j)/num_tanks,m_all_tank_store(id_ins_config(j),j),m_hloss_store(id_ins_config(j),j))
    

end