%% Headland Term Project
% The following script and functions are my work for the JH Kinematics and
% Dynamics of Robots Course Project assignment. The two files that are not
% my own work are: odeabort.m and odeprog.m whose attributions are in the
% header of those files.

% The script generates torque profiles for each joint of a 6 degree of
% freedom UR5 robot and uses those to calculate joint angles as a function
% of time. TO do this, the equations of motion for the UR5 robot are
% generated using Matlab's symbolic toolbox and then are passed to ODE45
% along with the joint torque profiles.

% This script requires 3 essential files to run properly:
% -f_UR5_chriss_mat.m
% -f_UR5_i_matrix.m
% -f_UR5_dgrav.m

% These function files are generated via matlabFunction and are called from
% within the robot equations of motion ode function. They can be generated
% with the commented out code in the script below, but must be modifed to
% accept the correct input.

% There are 2 non-essential files used:
% -odeprog.m
% -odeabort.m

% These two functions are used to add progress and abort functionallity to
% ode45 by passing them in as options. To run this script without these
% files, simply comment out the options line in above the ode45 call.

clear all

%% Parameters
% In this section, robot parameters are defined. The intention with this
% organization was that the parameters of any robot, of arbitrary number
% joints could be specified.

% UR5 Robot Parameters
% Initial values for joint variables Ut1 through Ut8

% Joint position
Uti = [0,0,0,0,0,0];
% Joint velocity
Utdi = [0,0,0,0,0,0];

% This array discribes the joint types and is used with the Jacobian matrix
% generation.
UR5_joints = [1,1,1,1,1,1]; % 0 for prismatic, 1 for revolute

% Make joint variables symbolic (requires Maple toolbox)
syms Ut1 Ut2 Ut3 Ut4 Ut5 Ut6
q = [Ut1 Ut2 Ut3 Ut4 Ut5 Ut6];

% Dimensions are radians and meters.
UR5_DH_table = ...
    [Ut1, 0.08959 ,        0,  pi()/2;...
     Ut2, 0       ,   -0.425,       0;...
     Ut3, 0       , -0.39225,       0;...
     Ut4, 0.10915 ,        0,  pi()/2;...
     Ut5, 0.09465 ,        0, -pi()/2;...
     Ut6, 0.0823  ,        0,       0];
 
% UR5 arm masses taken from www.universal-robots.com
UR5_link_mass = [3.7,8.393,2.33,1.219,1.219,0.1879]; % in kg

% The UR5 link inertias are approximated as uniform cylinders by a function
% towards the bottom of this file. This program uses link inertia matrices
% as 3x3 matrices stored in a cell array.
UR5_link_inertias = UR5_MMOI_Cyl(UR5_link_mass,UR5_DH_table);


%% Script

% Create array of torq profiles data by calling a function defined towards
% the end of the file. The function generates the torq data based on a
% piece-wise expression.
ft = linspace(0,15,100); % Time data, here set as 0 to 15 seconds with 100 divisions
for i = 1:length(ft)
    torq_p(i,:) = Utp(ft(i));
end

% Set current robot DH table here. The idea being that everything below
% this point in the script would be general enough to be for any arbitrary
% robot. Because I could not get the output of matlabFunction to do what I
% wanted, this is not the case.
DHt = UR5_DH_table;

% Create cell array of symbolic A matrices
for i = 1:size(DHt,1)
    a_matrixes{i} = A(DHt(i,:));  
end

% Generate symbolic transformation matrices
h_matrixes = H_links(a_matrixes);

% Generate symbolic link midpoint jacobians
jacobian = Jac(h_matrixes,UR5_joints);


% Generate symbolic expression for inertia matrix
 i_matrix = D(UR5_link_mass,UR5_link_inertias,jacobian,h_matrixes);

% Generate symbolic Christoffel matrix and then convert to function file.
% Note that this line is commented out and the function file has been
% included with this main file. The reason for this is that the generated
% function file has been altered to work (in a way that I understand) and
% generating a new one will break the program.

%  fchris_mat = matlabFunction(C(UR5_link_mass,UR5_link_inertias,h_matrixes,jacobian), q, 'File','f_UR5_chris_mat','Optimize',false);


% Convert inertia matrix expression to function file. This is done so that
% the inertia matrix can be evaluated for a given vector or joint values
% inside the ode function. % Note that this line is commented out and the function file has been
% included with this main file. The reason for this is that the generated
% function file has been altered to work (in a way that I understand) and
% generating a new one will break the program.

% fi_matrix = matlabFunction(i_matrix, q,'FIle','f_UR5_i_matrix','Optimize',true);


% Convert change in potential expression to function file. This is done so that
% the function can be evaluated for a given vector or joint values
% inside the ode function. % Note that this line is commented out and the function file has been
% included with this main file. The reason for this is that the generated
% function file has been altered to work (in a way that I understand) and
% generating a new one will break the program.

% fpotent = matlabFunction(G(h_matrixes,UR5_link_mass), q,'File','f_UR5_dgrav','Optimize',false);


% Solve the set of odes that are the robot equations of motion
options=odeset('OutputFcn',@odeprog,'Events',@odeabort);
[t, q] = ode45(@(t,q) UR5_em(t,q, torq_p, ft), [0 15], [Uti ; Utdi],options);

%% Functions
% This sections includes all of the functions used to create and solve the
% equations of motion.

% This function creates an A matrix as a 4x4 symbolic matrix.
function a_matrix = A(DH_params)
% takes input in form of a single row of D-H parameters of the form:
% [theta, d, a, alpha]

p = DH_params;

a_matrix = [cos(p(1)), -sin(p(1))*cos(p(4)), sin(p(1))*sin(p(4)), p(3)*cos(p(1));...
            sin(p(1)),  cos(p(1))*cos(p(4)), cos(p(1))*sin(p(4)), p(3)*sin(p(1));...
                    0,            sin(p(4)),           cos(p(4)),           p(2);...
                    0,                    0,                   0,             1];


end

% This function creates all intermediary H matrices need for creating
% Jacobian. Good for any arbitrary robot.
function h_matrix = H_links(a_matrixes)
% takes a cell array of A matrixes as input

A = a_matrixes; % renamed input for convenience

% This is a lazy attempt at approximating the midpoint  of a link by
% halving the d and a parameters of the DH table. 
approx = [1 1 1 .5;...
          1 1 1 .5;...
          1 1 1 .5;
          1 1 1  1];

for i = 1:size(A,2)
    if i == 1
        h_matrix{i} = A{1}.*approx;
        h_intermediate{i} = A{1};
    else
        h_intermediate{i} = h_intermediate{i-1}*A{i};
        h_matrix{i} = h_intermediate{i-1}*(A{i}.*approx);
    end
end
end

% This function generates a symbolic expression for the end effector
% jacobian of a robot. Good for any arbitrary robot.
function h_matrix_end = H_end(a_matrixes)
% takes a cell array of A matrixes as input

A = a_matrixes; % renamed input for convenience
i = size(A,2);

for i = 1:size(A:2)
    if i == 1
        h_matrix_end = A{1};
    else
        h_matrix_end = h_matrix_end*A{i};
    end
end
end


% This function creates the Jacobian of all links of the robot. Good for
% any arbitrary robot.
function jac = Jac(h_matrix, joint_types)
% takes H matrix of all links as a cell array and array of joint types

% rename inputs for convenience
H = h_matrix; 
j = joint_types; 

z0 = [0;0;1];
for n = 1:size(H,2) % for each link
    for i = 1:size(H,2) % for each joint
        if i > n
            jac{n}(:,i) = [0;0;0;0;0;0];
        elseif i == 1
            jac{n}(:,i) = [cross(z0 , H{i}(1:3,4).^j(i)) ; z0.*j(i)];
        else
            jac{n}(:,i) = [cross(H{i-1}(1:3,3) , (H{n}(1:3,4)-H{i-1}(1:3,4)).^j(i)) ; H{i-1}(1:3,3).*j(i)];
        end
    end
end
end

% This function creates a symbolic inertia tensor for the robot, good for
% any arbitrary robot.
function inert_mat = D(link_mass,link_inertia,jacobian,h_matrixes)
% takes an array of link masses, a cell array of link inertias a cell array,
% of link jacobians, and a cell array of h matrixes

% rename inputs for convenience
M = link_mass;
J = jacobian;
H = h_matrixes;
I = link_inertia;

% for the first in the series
inert_mat = M(1).*J{1}(1:3,:).'*J{1}(1:3,:) + ...
            J{1}(4:6,:).'*H{1}(1:3,1:3)*I{1}*H{1}(1:3,1:3).'*J{1}(4:6,:);
        
for n = 2:size(M,2) % for every link except the first
    inert_mat = inert_mat + ...
            M(n).*J{n}(1:3,:).'*J{n}(1:3,:) + ...
            J{n}(4:6,:).'*H{n}(1:3,1:3)*I{n}*H{n}(1:3,1:3).'*J{n}(4:6,:);    
end
end



% This function generates symbolic Christoffel matrix, good for any
% arbitrary robot.
function c_mat = C(link_mass,link_inertias,h_matrixes,link_jacobians)
% The Coriolis matrix is generated from a factorization rather than by
% differentials. This is much more computationally efficient when using
% Matlab's symbolic toolbox. The proof is described in:
% Bjerkeng, Magnus, 'A new Coriolis matrix factorization', 2012, IEEE
% International Conference on Robotics and Automation.

% rename inputs for convenience
M = link_mass;
I = link_inertias;
H = h_matrixes;
J = link_jacobians;


% create first element of summation
c_mat = M(1).*J{1}(1:3,:).'*diff(J{1}(1:3,:)) +...
        J{1}(4:6,:).'*H{1}(1:3,1:3)*I{1}*H{1}(1:3,1:3).'*diff(J{1}(4:6,:)) +...
        J{1}(4:6,:).'*diff(H{1}(1:3,1:3))*I{1}*H{1}(1:3,1:3).'*J{1}(4:6,:);
    
for n = 2:size(M,2)
c_mat = c_mat + ...
        M(n).*J{n}(1:3,:).'*diff(J{n}(1:3,:)) +...
        J{n}(4:6,:).'*H{n}(1:3,1:3)*I{n}*H{n}(1:3,1:3).'*diff(J{n}(4:6,:)) +...
        J{n}(4:6,:).'*diff(H{n}(1:3,1:3))*I{n}*H{n}(1:3,1:3).'*J{n}(4:6,:);
end

end


% This function creates a symbolic expression for the change in potential
% energy of a robot. Good for any size robot.
function dgrav = G(h_matrixes,link_mass)

%rename inputs for convenience
H = h_matrixes;
M = link_mass;
g = 9.81; % m/s^2

    
for i = 1:size(M,2)
    dgrav(i) = M(i)*g*H{i}(3,4);
end

end


% This function describes the equations of motion of a UR5 robot as two
% first order ODEs. Only good for UR5 robot.
function dqdt = UR5_em(t,q, torq, torq_div)

% reassign q so I can use it in a more convenient way
q1 = q(1:6,1);
q2 = q(7:12,1);


% rename inputs for convenience
D = f_UR5_i_matrix(q1);
C = f_UR5_chris_mat(q1);
%G = -f_UR5_dgrav(q1);
G = [0 0 0 0 0 0];



% to get a value of torque for a given time, the torq profile data that was
% generated in the above script is interpolated.
T = interp1(torq_div,torq,t);

% Here in comments describes my method.

% D*q'' + C*q' + G = tau

% q'' = tau/D - C/D*q' - G/D

% q1 = q  ;   q1' = q'
% q2 = q1'
% q2' = D(q1)/tau - C(q1)/D(q1)*y2 - G(q1)/D(q1) 

% If I understand correctly, and I don't think that I do, the q input from
% the ode45 is the length of the number of initial conditions given to
% ode45 which is the number of joint variables times the number of first
% order ODEs ode45 is solving.


dqdt = [q2;...
        D\T.' - (D\C)*q2 - D\G.'];


end

% This function gives the torq of the UR5 robot joints for a given time
function tor = Utp(time)

ti = time;

% joint 1
if ti < 4
    tor(1) = 5*ti;
elseif ti >= 4 && ti < 8
    tor(1) = 20-5*(ti-4);
else
    tor(1) = 0;
end

% joint 2
if ti >= 5 && ti < 8.5
    tor(2) = 15/3.5*(ti-5);
elseif ti >= 8.5 && ti < 12
    tor(2) = 15 - 15/3.5*(ti-8.5);
else
    tor(2) = 0;
end

% joint 3
if ti >= 8 && ti < 11.5
    tor(3) = 10/3.5*(ti-8);
elseif ti >= 11.5 && ti < 15
    tor(3) = 10 - 10/3.5*(ti-11.5);
else
    tor(3) = 0;
end

% joint 4
if ti >= 10 && ti < 12.5
    tor(4) = 2*(ti-10);
elseif ti >= 12.5 && ti < 15
    tor(4) = 5 - 2*(ti-12.5);
else
    tor(4) = 0;
end

% joint 5
if ti >= 10 && ti < 12.5
    tor(5) = 2*(ti-10);
elseif ti >= 12.5 && ti < 15
    tor(5) = 5 - 2*(ti-12.5);
else
    tor(5) = 0;
end

% joint 6

tor(6) = 0;

end

function mmoi = UR5_MMOI_Cyl(link_mass,dh_table)
% This functoin approximates the moments of inertia of the links of the UR5
% robot as uniform solid cylinders. Finding acutal dimensions of the of the
% links was oddly difficult.

% rename input for convenience
m = link_mass;
dh = dh_table;
%dh(:,2:4) = double(dh(:,2:4)); % converts dh table to double from symbols

r = 0.075; % radius in meters

mmoi{1} = [1/12*m(1)*dh(1,2)^2 + 1/4*m(1)*r^2, 0, 0;
           0, 1/12*m(1)*dh(1,2)^2 + 1/4*m(1)*r^2, 0;
           0, 0, 1/2*m(1)*r^2];
       
mmoi{2} = [1/12*m(2)*dh(2,3)^2 + 1/4*m(2)*r^2, 0, 0;
           0, 1/12*m(2)*dh(2,3)^2 + 1/4*m(2)*r^2, 0;
           0, 0, 1/2*m(2)*r^2];
       
mmoi{3} = [1/12*m(3)*dh(3,3)^2 + 1/4*m(3)*r^2, 0, 0;
           0, 1/12*m(3)*dh(3,3)^2 + 1/4*m(3)*r^2, 0;
           0, 0, 1/2*m(3)*r^2];
       
mmoi{4} = [1/12*m(4)*dh(4,2)^2 + 1/4*m(4)*r^2, 0, 0;
           0, 1/12*m(4)*dh(4,2)^2 + 1/4*m(4)*r^2, 0;
           0, 0, 1/2*m(4)*r^2];
       
mmoi{5} = [1/12*m(5)*dh(5,2)^2 + 1/4*m(5)*r^2, 0, 0;
           0, 1/12*m(5)*dh(5,2)^2 + 1/4*m(5)*r^2, 0;
           0, 0, 1/2*m(5)*r^2];
       
mmoi{6} = [1/12*m(6)*dh(6,2)^2 + 1/4*m(6)*r^2, 0, 0;
           0, 1/12*m(6)*dh(6,2)^2 + 1/4*m(6)*r^2, 0;
           0, 0, 1/2*m(6)*r^2];
    

end