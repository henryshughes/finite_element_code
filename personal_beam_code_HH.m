% Henry Hughes code to compute vertical reactions and moments for beams.

clc; clear all; close all; 
syms x 

% Initialize variables
ndof = 2; % Number of degrees of freedom per node
nel = 4; % Number of elements
nnp = 5; % Number of nodes

neq = ndof*nnp; % Number of Equations

x_coordinates = [0, 3, 5.5, 8, 11]; % x-coordinates of nodes

   
conn = zeros(nel,2); % Initialize the connectivity matrix

conn(1,:) = [1,2]; % Insert node to node connectivity of element 1
conn(2,:) = [2,3]; % Insert node to node connectivity of element 2
conn(3,:) = [3,4]; % Insert node to node connectivity of element 3
conn(4,:) = [4,5]; % Insert node to node connectivity of element 4



I = []; % Element cross sectional area
E = []*30000; % Element modulus of elasticity

leng = zeros(1,nel); % Initialize element length vector


% For loop to load element lengths
for i = 1:nel
    x1 = x_coordinates(conn(i,1));
    x2 = x_coordinates(conn(i,2));
    leng(i) = (x2-x1);
 
end

% Indicate essential boundary conditions and settlements
fixed_dofs  = [1  ,2  ,3,  7];   
fixed_values = [0.0  ;0.0  ;0.0  ;0.0];

% Get free degrees of freedom
free_dofs = setxor(1:neq,fixed_dofs);

% Get vector of external forces

fext = zeros(neq,1);
fext(5) = 1;	 


% Initialize global stiffness matrix, force vector, displacement vector
K = zeros(neq,neq);
 

% Create global stiffness matrix
for e=1:nel
    % Get mapping sctr
    
    n1 = conn(e,1); % First node
    n2 = conn(e,2); % Second node
    
    sctr = zeros(4,1);
    sctr(1) = (n1-1)*2+1; % Local x1 dof
    sctr(2) = (n1-1)*2+2; % Local m1 dof
    sctr(3) = (n2-1)*2+1; % Local x2 dof
    sctr(4) = (n2-1)*2+2; % Local m2 dof 
    
    % Define variables to assemble local k matrix
    
    % Build the local k matrix
    % Consider adding E*I into the front of this expression. 
    Ke = 1/(leng(e)^3) *...
        [12, 6*leng(e), -12, 6*leng(e);
         6*leng(e), 4*(leng(e))^2, -6*leng(e), 2*(leng(e))^2;
        -12, - 6*leng(e), 12,  -6*leng(e);
         6*leng(e), 2*(leng(e))^2,  -6*leng(e), 4*(leng(e))^2];
    % Send to global k matrix
    K(sctr, sctr) = K(sctr, sctr) + Ke;
end



% Define submatrices
KE = K(fixed_dofs, fixed_dofs);
KF = K(free_dofs, free_dofs);
KEF = K(fixed_dofs, free_dofs);

% Define known force vectors
fE = fext(fixed_dofs);
fF = fext(free_dofs);

fW = zeros(neq,1);


fWE = fW(fixed_dofs);
fWF = fW(free_dofs);


% Define known displacement vector
dE = fixed_values;

% Solve unknown displacements
dF = inv(KF)*(fF + fWF - KEF'*dE);

% Solve unknown reactions
rE = KE*dE + KEF*dF - fE - fWE;

% Reassemble r and d to include all degrees of freedom

r = zeros(neq,1);
r(free_dofs) = 0;
r(fixed_dofs) = rE;

d = zeros(neq,1);
d(free_dofs) = dF;
d(fixed_dofs) = fixed_values;


 for e = 1:nel
     
     % Create the closed form N vector
     N = [1/(leng(e))^3*(2*x^3 - 3*x^2*leng(e) + (leng(e))^3),...
     1/(leng(e))^3*(leng(e)*x^3 - 2*x^2*(leng(e))^2 + x*(leng(e))^3),...
     1/(leng(e))^3*(-2*x^3 + 3*x^2*leng(e)),...
     1/(leng(e))^3*(leng(e)*(x)^3 - x^2*(leng(e))^2)];
 
     n1 = conn(e,1); % First node
     n2 = conn(e,2); % Second node
     sctr = zeros(4,1);  
     sctr(1) = (n1-1)*2+1; % Local x1 dof
     sctr(2) = (n1-1)*2+2; % Local m1 dof
     sctr(3) = (n2-1)*2+1; % Local x2 dof
     sctr(4) = (n2-1)*2+2; % Local m2 dof 
    
     delement = d(sctr);
    
     element_displacement = N*delement;
     pretty(element_displacement)
     
 
 end
 








