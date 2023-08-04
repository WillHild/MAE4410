clear all; close all; clc
format("default")
W_TO = 63.1 * 1000; % kg (guess)
passenger = 110 + 8; % number of passenger
range_min = 2200; % minimum range requirement
range_max = 2500; % maximum range requirement

range = range_max;
%range = 3000;% * 1.852

pass_weight = 94; % NA Average Adults weight
bag_weight = 20; % International Economy Baggage weight


%% Calculations

W_payload = passenger * (pass_weight + bag_weight)* 1.1;
diff = 1;
while diff > 0.0001
    % Empty Weight Fraction
    % Raymer Parameters
    A = 0.97; C = -0.06;

    W_OE_TO = A * W_TO ^ C;
    % fprintf('\nEmpty Weight Fraction: %.4f',W_OE_TO)

    % Cruise
    % Fuel Efficiency (lower end) Roskam
    Cj = 0.5;
    LD = 14;

    TAS = 470.3; % to be changed

    W54 = exp(- range * Cj/(TAS * LD));

    % Loitering
    Cj_loiter = 0.4;
    LD_loiter = 16;
    E_loiter = 0.5; % hr

    W65 = exp(-E_loiter * Cj_loiter/LD_loiter);


    % Full Fuel Fractions
    % Use of standards values
    W8TO = 0.99 * 0.99 * 0.995 * 0.98 * W54 * W65 * 0.99 * 0.992;

    W_fuel_TO = 1 - W8TO;

    W_TO_calc = W_payload /  (1 - W_fuel_TO - W_OE_TO);

    diff = abs(W_TO_calc - W_TO);

    W_TO = W_TO_calc;
end
fprintf('\n      Numer of Passenger: %12.f', passenger)
fprintf('\nAverage Passenger Weight: %12.f kg', pass_weight)
fprintf('\n      Average Bag Weight: %12.f kg', bag_weight)
fprintf('\n             Loiter Time: %12.4f hr', E_loiter)
fprintf('\n Loiter Fuel Consumption: %12.4f kg', Cj_loiter)
fprintf('\n         Take-Off Weight: %12.4f kg',W_TO)
fprintf('\n')

