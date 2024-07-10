% Project - Sounding Rocket Simulator
% Developer info - Syed Affan, B.Tech ECE Vellore Institute of Technology
% Project under - Dr. Hari Hablani

% Given parameters
g = 9.81; % m/s^2, acceleration due to gravity
isp = 122.96; % s, specific impulse
Rmo = 12.137; % kg, launch mass
thrust = 2006; % N, thrust force
Mmo = 5.844; % kg, motor hardware + propellant mass
Pmo = 3.34956; % kg, propellant mass
cd = 0.64; % drag coefficient for truncated tangent ogive nose cone
cd2 = 1.0;
cdp = 1.75;
rho = 1.17123; % kg/m^3, air density for pressure 29.7 inHg and 26C

% Calculated parameters
c = isp * g; % m/s, effective exhaust velocity
A = 0.00723456; % m^2, cross-sectional area of airframe
A2 = 0.22752;
Ap = 1.76625;
mbo = Rmo - Pmo; % kg, burnout mass
mdot = thrust / (isp * g); % kg/s, mass burn rate
tb = Pmo / mdot; % s, time at burnout

k = cd * A * 0.5 * rho;
k2 = cd2 * A2 * 0.5 * rho;
kp = cdp * Ap * 0.5 * rho;



A1 = [];
A2 = [];
A3 = [];
% Burn phase calculations
for dt = 0:0.03:tb
    W = Rmo - (mdot * dt);
    v = (c * log(Rmo / W)) - (g * dt); 
    D = k * v^2;
   
    cosh_arg = ((dt / W) * sqrt(k * (thrust - (W * g))));
    sb = (W / k) * log(cosh(cosh_arg));
    
    A1(end+1) = sb;
    A2(end+1) = v;
    A3(end+1) = D;
end

% Coasting phase calculations
tc = sqrt(mbo / (g * k)) * atan(v * sqrt(k / (mbo * g)));
psi = k / mbo;

for dt1 = 0:0.03:tc
    vc = sqrt(g / psi) * tan(atan(sqrt(psi / g) * v) - sqrt(psi * g) * dt1);
    sc = (mbo / (2 * k)) * log(((mbo * g) + (k * (v - vc)^2)) / (mbo * g)); 
    Dc = k * vc^2;
    smax = sc + sb;
    
    A1(end+1) = smax; 
    A2(end+1) = vc;
    A3(end+1) = Dc;
end

% Descent phase calculations
tt = smax * sqrt(kp / (mbo * g));

for dt2 = 0:0.03:tt
    exponent = sqrt(2 * cdp * Ap * rho * g / mbo) * dt2;
    numerator = sqrt(2 * g * mbo) * (exp(exponent) - 1);
    denominator = sqrt(cdp * Ap * rho) * (exp(exponent) + 1);
    vf = numerator / denominator;
    vt = sqrt(mbo * g / kp);
    Df = kp * vf^2;
    sf = smax - vt * dt2 + ((vt * mbo) / kp) * (exp(-kp * dt2 / mbo) - 1);
    
    if sf > 0 
        A1(end+1) = sf;  
        A2(end+1) = vf;
        A3(end+1) = Df;
    end
end

% Time vectors for plotting
t = 0:0.03:(length(A1) - 1) * 0.03;
t2 = 0:0.03:(length(A2) - 1) * 0.03;
t3 = 0:0.03:(length(A3) - 1) * 0.03;

% Plotting results
figure;

subplot(3,1,1)
plot(t, A1);
xlabel('Time (s)');
ylabel('Altitude (m)');
title('Plot of Altitude vs Time');
grid on;

subplot(3,1,2)
plot(t2, A2);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Plot of Total Velocity vs Time');
grid on;

subplot(3,1,3)
plot(t3, A3);
xlabel('Time (s)');
ylabel('Drag Force (N)');
title('Plot of Drag Force vs Time');
grid on;
