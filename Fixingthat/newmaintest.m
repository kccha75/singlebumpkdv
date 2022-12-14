clear;close
% 1D Fourier Newton Multigrid scheme on u_xx+au_x+bu+cu^2=f on domain 
% [-L/2,L/2] with periodic boundary conditions

% Uses vcycles to iteration (displays rms of residual after every v-cycle)

% Relaxation selections:
% SR - Stationary Richardson
% NSR - Non-Stationary Richardson (3 parameters)
% MRR - Minimum Residual Richardson (see Boyd)
% RSM - Residual smoothing method (see Canuto)
% cg_it - Conjugate Gradient

% Non-preconditioned relaxations ONLY!!!

% Important notes:
% SR and NSR assumes eigmin = 1, eigmax = N^2/4 (this is true for
% self-adjoint case only where a=1, b=0!!! with L=2*pi use with caution)

%--------------------------------------------------------------------------
global L

% Domain size
L = 30;

% Number of grids (total points being 2^(grids))
grids=9;

% Grid where solution is exactly solved at 2^(gridexact) points
gridexact=5; 

% Iterations on down cycle (on finest grid and after restriction)
Nd=3; 

% Iterations on up cycle (after prolongation)
Nu=3; 

% Number of v-cycles
num_vcycles=3; 

% Perform num_vcycles v-cycles
vcyclegrid=grids-gridexact+1; % number of vcycle grids 
% (since lower grids do not need storage as it is solved exactly)

% Choice of relaxation
relaxation=@MRR;

% Initial grid parameters
N = 2^grids; % initial grid size (Ideally 2^n sized)
k = 2*pi/L*[0:N/2-1 N/2 -N/2+1:-1]; % wave numbers
x = L*(-N/2:N/2-1)'/N; % Assumed to be periodic thus only going to N-1

% Define RHS 
RHS=0.0*sech(x).^2;

% a(x) function
a=zeros(N,1);

% b(x) function
b=1*ones(N,1);

% c(x) function
c=-3*ones(N,1);

% Preconditioner
preconditioner=@Hfd;

% Pseudo Arclength continuation step length
ds=0.01;

% Number of steps taken in pseudo arclength continuation
steps=5000;

%--------------------------------------------------------------------------
% SOLVE FOR INITIAL SOLUTION + LAMBDA USING NEWTON
%--------------------------------------------------------------------------

% Initial gamma (chosen)
gamma0=-4.6;

% Define RHS 
RHS=gamma0*sech(x).^2;

% Initial guess (near solution)
v0=-.75*sech(x).^2;

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k,a,b,c,RHS,v0);

% tic
cellv=Newton_vcycle(relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellc,cellRHS,cellv);
% toc

% Save solution
v0=cellv{1};

%--------------------------------------------------------------------------
% SOLVE FOR SECOND SOLUTION (NEAR INITIAL) USING NEWTON
%--------------------------------------------------------------------------

% Second gamma (chosen)
gamma1=gamma0+ds;

% Define RHS 
RHS=gamma1*sech(x).^2;

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k,a,b,c,RHS,v0);

% Initial guess uses previous solution (saved in cellv{1})
% tic
cellv=Newton_vcycle(relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellc,cellRHS,cellv);
% toc

% Save solution
v1=cellv{1};

%--------------------------------------------------------------------------
% FIND DV, DLAMBDA FOR PSEUDO ARCLENGTH INITIAL CONDITIONS
%--------------------------------------------------------------------------

gamma(2)=gamma1;
gamma(1)=gamma0;

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k, ...
    a,b,c,RHS,v0);

% Tangent approximation (found using secant approximation of two solutions)
dv=(v1-v0)/ds;
dgamma=(gamma(2)-gamma(1))/ds;

% Normalise
mag=sqrt(dot(dv,dv)+dgamma^2);
dv=dv/mag;
dgamma=dgamma/mag;

tic
% Pseudo arclength to trace solution family
[v,gamma]=pseudoarclength(v0,gamma,dv,dgamma,ds,steps,... 
    x,relaxation,vcyclegrid,cellN,cellk,Nd,Nu,cella,cellb,cellc,cellRHS);

toc