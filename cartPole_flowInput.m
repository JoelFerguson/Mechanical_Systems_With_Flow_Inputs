%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Cart-pole system with flow input
% Description: Simulation files for the paper 'Port-Hamiltonian 
% representation of mechanical systems with flow inputs',
% submitted to CDC/LCSS
% Authours: Joel Ferguson
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialise
% Clear workspace
clear
close all
clc

% Set Figure default values
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize',11);
set(0,'DefaultLineLineWidth',2.0);
set(0,'DefaultAxesLineWidth',0.5);
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'defaultAxesNextPlot','add')

%% Simulation settings
% Simulation initial conditions
sim.tf = 40;             % Simulation end time
sim.q_init = [0.5; 0];  % Configuration initial conditions
sim.p_init = [0; 0];    % Canonical momentum initial conditions

%% Define cart-pole model
% Parameters
% kinematic parameters
sys.l = 1;
sys.G = [0; 1];

% Define degrees of freedom
sys.G_size = size(sys.G);
sys.N = sys.G_size(1);
sys.k = sys.G_size(2);
sys.Gp = [eye(sys.N-sys.k), zeros(sys.k)];

% Inertial parameters
sys.mp = 1;
sys.mc = 1;
sys.M = @(q) [sys.mp*sys.l^2, sys.mp*sys.l*cos(q(1));
              sys.mp*sys.l*cos(q(1)), sys.mc+sys.mp];
            
% Friction parameters. Note that the controller from Acosta 2005 assumes no
% friction.
sys.D = @(q) zeros(2);

% Potential parameters
sys.g = 9.8;
sys.V = @(q) sys.mp*sys.g*sys.l*cos(q(1));

% System energy
sys.H = @(q,p) 0.5*p.'*(sys.M(q)\p) + sys.V(q);

% Symbolic variables for gradient computations
syms q_sym [sys.N,1]
syms p_sym [sys.N,1]

% Energy gradients
sys.dHdq = matlabFunction(jacobian(sys.H(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);
sys.dHdp = @(q,p) sys.M(q)\p;

% System ODE
sys.dx = @(q,p,u) [zeros(sys.N) eye(sys.N); -eye(sys.N) -sys.D(q)]*[sys.dHdq(q,p); sys.dHdp(q,p)] + [zeros(size(sys.G)); sys.G]*u;

%% Constructe momentum transformation for velocity input model
% Partition mass matrix
sys.m11 = @(q) sys.Gp*sys.M(q)*sys.Gp.';
sys.m21 = @(q) sys.G.'*sys.M(q)*sys.Gp.';
sys.m22 = @(q) sys.G.'*sys.M(q)*sys.G;

% Construct momentum transformation matrix
tran.a22 = @(q) inv(sys.m22(q) - (sys.m21(q)/sys.m11(q))*sys.m21(q).');
tran.a21 = @(q) -tran.a22(q)*sys.m21(q)/sys.m11(q);
tran.A = @(q) [eye(sys.N-sys.k), zeros(sys.N-sys.k,sys.k);
                tran.a21(q), tran.a22(q)];

% Construct \bar C matrix
tran.dAipbdq = matlabFunction(jacobian(tran.A(q_sym)\p_sym,q_sym),'vars',[{q_sym}, {p_sym}]);
tran.Cbar = @(q,pb) tran.A(q)*(tran.dAipbdq(q,pb).' - tran.dAipbdq(q,pb))*tran.A(q).';
tran.Cb11 = @(q,p1,v) sys.Gp*tran.Cbar(q,[p1;v])*sys.Gp.';
tran.Cb21 = @(q,p1,v) sys.G.'*tran.Cbar(q,[p1;v])*sys.Gp.';
tran.Cb22 = @(q,p1,v) sys.G.'*tran.Cbar(q,[p1;v])*sys.G;

% Construct \bar D matrix
tran.Dbar = @(q) tran.A(q)*sys.D(q)*tran.A(q).';
tran.Db11 = @(q) sys.Gp*tran.Dbar(q)*sys.Gp.';
tran.Db21 = @(q) sys.G.'*tran.Dbar(q)*sys.Gp.';
tran.Db22 = @(q) sys.G.'*tran.Dbar(q)*sys.G;

% Define reduced energy function
tran.Tbar = @(q,p1,v) 0.5*p1.'*(sys.m11(q)\p1) + 0.5*v.'*(tran.a22(q)\v);
tran.Hbar = @(q,p1,v) tran.Tbar(q,p1,v) + sys.V(q);

% Symbolic variables for gradient computations
syms p1_sym [sys.N-sys.k,1]
syms v_sym [sys.k,1]

% Energy gradients of reduced Hamiltonian
tran.dHbdq = matlabFunction(jacobian(tran.Hbar(q_sym,p1_sym,v_sym),q_sym).','vars',[{q_sym}, {p1_sym}, {v_sym}]);
tran.dHbdq1 = matlabFunction(jacobian(tran.Hbar(q_sym,p1_sym,v_sym),q_sym1).','vars',[{q_sym}, {p1_sym}, {v_sym}]);
tran.dHbdq2 = matlabFunction(jacobian(tran.Hbar(q_sym,p1_sym,v_sym),q_sym2).','vars',[{q_sym}, {p1_sym}, {v_sym}]);
tran.dHbdp1 = matlabFunction(jacobian(tran.Hbar(q_sym,p1_sym,v_sym),p1_sym).','vars',[{q_sym}, {p1_sym}, {v_sym}]);
tran.dHdx = @(q,p1,v) [tran.dHbdq(q,p1,v); tran.dHbdp1(q,p1,v)];

% Define new input matrix
tran.Gb = @(q,p1,v) [-sys.m11(q)\sys.m21(q).';
                    eye(sys.k);
                    (-tran.Cb21(q,p1,v).' - tran.Db21(q).')/tran.a22(q)];

% Define flow input reduced order dynamics
tran.dx = @(q,p1,v) [zeros(sys.N), [eye(sys.N-sys.k); zeros(sys.k,sys.N-sys.k)];
                    [-eye(sys.N-sys.k); zeros(sys.k,sys.N-sys.k)].', tran.Cb11(q,p1,v)-tran.Db11(q)]*tran.dHdx(q,p1,v) ...
                    + tran.Gb(q,p1,v)*v;

% Define dynamic to convert effort input signal to flow input signal
tran.dv = @(q,p1,v,u) - tran.a21(q)*tran.dHbdq1(q,p1,v) - tran.a22(q)*tran.dHbdq2(q,p1,v) ...
                        + (tran.Cb21(q,p1,v) - tran.Db21(q))*(sys.m11(q)\p1) ...
                        + (tran.Cb22(q,p1,v) - tran.Db22(q))*(tran.a22(q)\v) + tran.a22(q)*u;

%% Define torque input control law from Acosta 2005
% Define the torque input controller tuning parameters
ctrl.q2d = 0;
ctrl.k = 1;
ctrl.m220 = 1;
ctrl.P = 1;
ctrl.kv = 5;

% The nonlinear control law requires a partial feedback linearisation
% control input. This is defined below
ctrl.u_fb = @(q,dq,w) (sys.mp^2*sys.l^2*sys.g*sin(q(1))*cos(q(1)) - sys.mp^2*sys.l^3*sin(q(1))*dq(1) ...
                        + det(sys.M(q))*w)/(sys.mp*sys.l^2);

% Model parameters for the partial feedback linearised system
ctrl.a = sys.g/sys.l;
ctrl.b = 1/sys.l;
ctrl.M_fb = eye(2);
ctrl.w_fb = @(q) ctrl.a*cos(q(1));
ctrl.G_fb = @(q) [-ctrl.b*cos(q(1)); 1];
ctrl.H_fb = @(q,p_fb) 0.5*p_fb.'*(ctrl.M_fb\p_fb) + ctrl.w_fb(q);
ctrl.dHdq_fb = matlabFunction(jacobian(ctrl.H_fb(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);

% Define the momentum vector used after feedback linearisation. Note that
% this is the velocity of the original system.
ctrl.p_fb = @(q,p) sys.M(q)\p;

% Define the closed-loop energy function and gradients
ctrl.Md = @(q) [(ctrl.k*ctrl.b^2/3)*cos(q(1))^3, -(ctrl.k*ctrl.b/2)*cos(q(1))^2;
                  -(ctrl.k*ctrl.b/2)*cos(q(1))^2, ctrl.k*cos(q(1))+ctrl.m220];
ctrl.wd = @(q) 3*ctrl.a/(ctrl.k*ctrl.b^2*cos(q(1))^2) + (ctrl.P/2)*(q(2) - ctrl.q2d + (3/ctrl.b)*log(sec(q(1)) + tan(q(1))) + (6*ctrl.m220/(ctrl.k*ctrl.b))*tan(q(1)))^2;
ctrl.Hd = @(q,p_fb) 0.5*p_fb.'*(ctrl.Md(q)\p_fb) + ctrl.wd(q);
ctrl.dHddp = @(q,p_fb) ctrl.Md(q)\p_fb;
ctrl.dHddq = matlabFunction(jacobian(ctrl.Hd(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);

% Define additional terms used in control input definition 
ctrl.gamma1 = @(q) -(ctrl.k*ctrl.b^2/6)*cos(q(1))^3;
ctrl.alpha = @(q) (ctrl.gamma1(q)*ctrl.k/2)*sin(q(1))*[-ctrl.b*cos(q(1)); 1];
ctrl.J2 = @(q,p_fb) p_fb.'*(ctrl.Md(q)\ctrl.alpha(q))*[0 1; -1 0];

% Define nonlinear torque input control law
ctrl.w = @(q,p_fb) (ctrl.G_fb(q).'*ctrl.G_fb(q))\ctrl.G_fb(q).'*(ctrl.dHdq_fb(q,p_fb) - (ctrl.Md(q)/ctrl.M_fb)*ctrl.dHddq(q,p_fb) + ctrl.J2(q,p_fb)*(ctrl.Md(q)\p_fb)) - ctrl.kv*ctrl.G_fb(q).'*(ctrl.Md(q)\p_fb);

% Define total control law including feedback linearisation
ctrl.u_cat = @(q,p_fb) ctrl.u_fb(q,p_fb,ctrl.w(q,p_fb));

%% Run simulation for torque input model
% Pack inital conditions
sim.ti = 0;
sim.x0 = [sim.q_init; sim.p_init];

% Concatenate plant and control functions
plant_cat = @(t,x) sys.dx(x(1:2),x(3:4),ctrl.u_cat(x(1:2),ctrl.p_fb(x(1:2),x(3:4))));

% Run ODE solver
[res_torque.tm,res_torque.xm] = ode45(plant_cat,[sim.ti sim.tf],sim.x0,odeset('RelTol',1e-6));

% Unpack and compute results
res_torque.q = res_torque.xm(:,1:2);
res_torque.p = res_torque.xm(:,3:4);
for i=1:length(res_torque.tm)
    res_torque.p_fb(i,:) = ctrl.p_fb(res_torque.q(i,:).',res_torque.p(i,:).').';
    res_torque.Hd(i) = ctrl.Hd(res_torque.q(i,:).',res_torque.p_fb(i,:).') - ctrl.Hd([0;0],[0;0]);
end

% Plot results
fig1 = figure("Name","Torque input simulation resutls");
plot(res_torque.tm,res_torque.q)
xlabel('time (s)')
ylabel('Configuration')
legend("$q_1$","$q_2$")

%% Run simulation for velocity input model
% Determine initial conditions from the canonical momentum
sim.ti = 0;
sim.pb0 = tran.A(sim.q_init)*sim.p_init;
sim.pb1_0 = sim.pb0(1);
sim.v_0 = sim.pb0(2);
sim.x0 = [sim.q_init; sim.pb1_0; sim.v_0];

% Define functino to reconstruct p_fb from p1, v for control implementation
tran.p_fb = @(q,p1,v) sys.M(q)\(tran.A(q)\[p1; v]);

% Concatenate plant and control functions
plant_cat = @(t,x) [tran.dx(x(1:2),x(3),x(4));
                    tran.dv(x(1:2),x(3),x(4),ctrl.u_cat(x(1:2),tran.p_fb(x(1:2),x(3),x(4))))];

% Run ODE solver
[res_velocity.tm,res_velocity.xm] = ode45(plant_cat,[sim.ti sim.tf],sim.x0,odeset('RelTol',1e-6));

% Unpack and compute results
res_velocity.q = res_velocity.xm(:,1:2);
res_velocity.p1 = res_velocity.xm(:,3);
res_velocity.v = res_velocity.xm(:,4);
for i=1:length(res_velocity.tm)
    res_velocity.u(i) = ctrl.u_cat(res_velocity.q(i,1:2).',tran.p_fb(res_velocity.q(i,1:2).',res_velocity.p1(i),res_velocity.v(i)));
end

% Plot results
fig2 = figure("Name","Velocity input simulation resutls");
plot(res_velocity.tm,res_velocity.q)
xlabel('time (s)')
ylabel('Configuration')
legend("$q_1$","$q_2$")

fig3 = figure("Name","Comparison of force and velocity input signals");
plot(res_velocity.tm,res_velocity.u,res_velocity.tm,res_velocity.v)
xlabel('time (s)')
ylabel('Input signals')
legend("Force input $u$","Velocity input $v$")
