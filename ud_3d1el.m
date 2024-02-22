
function [DEFL,REACT,Ele_Forces,AFLAG] = ud_3d1el(...
		nnodes,coord,concen,fixity,nele,ends,A,Izz,Iyy,J,Cw,IsSym,Ysc,Zsc,Betay,Betaz,Betaw,Zzz,Zyy,Ayy,Azz,...
		E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,anatype)
%UD_3D1EL performs a user defined three-dimensional
% first-order elastic analysis of a structural system.
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Functions Called
%              < to be defined by the student >
%
%  Dictionary of Variables
%     Input Information:
%       nnodes         ==  total number of nodes
%       coord(i,1:3)   ==  node i's coordinates
%                            coord(i,1) = X coordinate
%                            coord(i,2) = Y coordinate
%                            coord(i,3) = Z coordinate
%       concen(i,1:6)  ==  concentrated loads for node i's 6 d.o.f.
%                            concen(i,1) = force in global X direction
%                            concen(i,2) = force in global Y direction
%                            concen(i,3) = force in global Z direction
%                            concen(i,4) = moment about global X axis
%                            concen(i,5) = moment about global Y axis
%                            concen(i,6) = moment about global Z axis
%       fixity(i,1:6)  ==  prescribed displacements for node i's 6 d.o.f.
%                          Note: A free d.o.f. will have a value of NaN
%                          and hence, you will find the Matlab function
%                          isnan very useful.
%                          Examples: If fixity(15,3) is set to NaN, then node 15's
%                                      Z-disp component is free;
%                                    If fixity(2,6) is set to 0.0, then node 2's
%                                      Z-rotation component is supported;
%                                    If fixity(5,2) is set to -2.1, then node 5's
%                                      Y-disp component is supported and defined
%                                      with a settlement of -2.1 units.
%                            fixity(i,1) = prescribed disp. in global X direction
%                            fixity(i,2) = prescribed disp. in global Y direction
%                            fixity(i,3) = prescribed disp. in global Z direction
%                            fixity(i,4) = prescribed rotation about global X axis
%                            fixity(i,5) = prescribed rotation about global Y axis
%                            fixity(i,6) = prescribed rotation about global Z axis
%       nele           ==  total number of elements
%       ends(i,1:14)   ==  element i's nodal information
%                            ends(i,1) = start node #
%                            ends(i,2) = finish node #
%                            ends(i,3) = flag to indicate whether or not flexural
%                            moments are released at start node.  ends(i,3)=0 both not
%                            released (rigid connection); ends(i,3)=1 both flexural
%                            moments are released (pinned connection); ends(i,3)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,4) = flag to indicate whether or not flexural
%                            moments are released at finish node.  ends(i,4)=0 both not
%                            released (rigid connection); ends(i,4)=1 both flexural
%                            moments are released (pinned connection); ends(i,4)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,5) = flag to indicate the degree of warping
%                            restraint at start node.  ends(i,5)=0 warping free;
%                            ends(i,5)=1 warping fixed; ends(i,5)=2 warping continuous
%                            ends(i,6) = flag to indicate the degree of warping
%                            restraint at finish node.  ends(i,6)=0 warping free;
%                            ends(i,6)=1 warping fixed; ends(i,6)=2 warping continuous
%                            ends(i,7) = rotational spring stiffness at the start
%                            node and about element i's local z-z axis.
%                            ends(i,8) = rotational spring stiffness at the start
%                            node and about element i's local y-y axis.
%                            ends(i,9) = rotational spring stiffness at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,10) = rotational spring stiffness at the finish
%                            node and about element i's local y-y axis.
%                            ends(i,11) = connection moment capacity Mpz at the start
%                            node and about element i's local z-z axis.
%                            ends(i,12) = connection moment capacity Mpy at the start
%                            node and about element i's local y-y axis.
%                            ends(i,13) = connection moment capacity Mpz at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,14) = connection moment capacity Mpy at the finish
%                            node and about element i's local y-y axis.
%       A(i)           ==  element i's cross sectional area
%       Izz(i)         ==  element i's moment of inertia about its local z-z axis
%       Iyy(i)         ==  element i's moment of inertia about its local y-y axis
%       J(i)           ==  element i's torsional constant
%       Cw(i)          ==  element i's warping constant
%       Zzz(i)         ==  element i's plastic section modulus about its local z-z axis
%       Zyy(i)         ==  element i's plastic section modulus about its local y-y axis
%       Ayy(i)         ==  element i's effective shear area along its local y-y axis
%       Azz(i)         ==  element i's effective shear area along its local z-z axis
%       E(i)           ==  element i's material elastic modulus, Young's Modulus
%       v(i)           ==  element i's material Poisson's ratio
%       Fy(i)          ==  element i's material yield strength
%       YldSurf(i)     ==  element i's yield surface maximum values
%                              YldSurf(i,1) = maximum P/Py value
%                              YldSurf(i,2) = maximum Mz/Mpz value
%                              YldSurf(i,3) = maximum My/Mpy value
%       Wt(i)          ==  element i's material weight density
%                          (Assume that gravity is directed in the negative global Y dir)
%       webdir(i,1:3)  ==  element i's unit web vector.  This is a unit vector
%                          that defines the element's local y-y axis with respect
%                          to the global coordinate system.  It is based on the
%                          structure's undeformed geometry.
%                              webdir(i,1) = x component of element's unit web vector
%                              webdir(i,2) = y component of element's unit web vector
%                              webdir(i,3) = z component of element's unit web vector
%                          NOTE: An element's 3x3 rotation matrix, [g], is constructed
%                          as follows: First, calculate a unit vector, x_vect, that
%                          describes the element's local x-axis. Second, take the
%                          cross product of x_vect and webdir(i,:) to obtain z_vect,
%                          i.e. z_vect = cross(x_vect,webdir(i,:)). Third, set z_vect 
%                          to a unit vector, i.e. z_vect = z_vect/norm(z_vect).
%                          Finally, the first row of [g] is x_vect, its second row is
%                          webdir(i,:), and its third row is z_vect.
%       beta_ang(i)    ==  element i's web rotation angle.  These values are
%                          provided for those students who are required to calculate
%                          their own unit web vectors (see above).  It is based
%                          on the structure's undeformed geometry.
%                          Note:  MASTAN2 uses the following convention for
%                                 defining a member's default web orientation:
%                                 A vector defing the element's local y-axis
%                                 with respect to the global coordinate system
%                                 will have a positive component in the global
%                                 Y direction.  If the element's local x-axis,
%                                 its length axis, is aligned with the global Y
%                                 axis, then element's local y-axis is aligned
%                                 with global negative X axis.  After this initial
%                                 orientation, element i may be rotated about
%                                 its local x-axis by the amount defined by
%                                 its web rotation angle, beta_ang(i).  The
%                                 angle is in radians and assumes a right-hand
%                                 convention about the local x-axis which runs from
%                                 the element's start node to its finish node.
%       w(i,1:3)         ==  element i's uniform load which references its
%                            local coordinate system
%                              w(i,1) = x component of uniform load
%                              w(i,2) = y component of uniform load
%                              w(i,3) = z component of uniform load
%       thermal(i,1:4)   ==  element i's thermal strain effects which reference its
%                            local coordinate system
%                              thermal(i,1) = coefficient of thermal expansion
%                              thermal(i,2) = change in temperature at centroid
%                              thermal(i,3) = linear temperature gradient in local y-dir
%                                           = (T_up_y - T_btm_y) / depth_y
%                              thermal(i,4) = linear temperature gradient in local z-dir
%                                           = (T_up_z - T_btm_z) / width_z
%       truss            ==  flag to indicate if structure is a truss or not
%                              truss = 0   System is not a truss
%                              truss = 1   System is a truss
%       anatype          ==  flag to indicate which type of analysis is requested
%                              anatype = 1  First-Order Elastic
%                              anatype = 2  Second-Order Elastic
%                              anatype = 3  First-Order Inelastic
%                              anatype = 4  Second-Order Inelastic
%                              anatype = 5  Elastic Buckling (Eigenvalue)
%                              anatype = 6  Inelastic Buckling (Eigenvalue)
%
%     Local Information:
%              < to be defined by the student >
%
%     Output Information:
%       DEFL(i,1:6)      ==  node i's calculated 6 d.o.f. deflections
%                              DEFL(i,1) = displacement in X direction
%                              DEFL(i,2) = displacement in Y direction
%                              DEFL(i,3) = displacement in Z direction
%                              DEFL(i,4) = rotation about X direction
%                              DEFL(i,5) = rotation about Y direction
%                              DEFL(i,6) = rotation about Z direction
%       REACT(i,1:6)     ==  reactions for supported node i's 6 d.o.f.
%                              REACT(i,1) = force in X direction
%                              REACT(i,2) = force in Y direction
%                              REACT(i,3) = force in Z direction
%                              REACT(i,4) = moment about X direction
%                              REACT(i,5) = moment about Y direction
%                              REACT(i,6) = moment about Z direction
%       ELE_FOR(i,1:1?)  ==  element i's internal forces and moments
%                            Note: All values reference the element's local
%                                  coordinate system.
%                              ELE_FOR(i,1)  = x-force at start node
%                              ELE_FOR(i,2)  = y-force at start node
%                              ELE_FOR(i,3)  = z-force at start node
%                              ELE_FOR(i,4)  = x-moment at start node
%                              ELE_FOR(i,5)  = y-moment at start node
%                              ELE_FOR(i,6)  = z-moment at start node
%                              ELE_FOR(i,7)  = x-force at end node
%                              ELE_FOR(i,8)  = y-force at end node
%                              ELE_FOR(i,9)  = z-force at end node
%                              ELE_FOR(i,10) = x-moment at end node
%                              ELE_FOR(i,11) = y-moment at end node
%                              ELE_FOR(i,12) = z-moment at end node
%                            If you are not programming warping torsion, the ELE_FOR
%                            array needs to contain only 12 columns, i.e. ELE_FOR(i,1:12)                           
%                            For those programming warping torsion, the bimoments and
%                            rates of twist should be stored as follows.
%                              ELE_FOR(i,13) = bimoment at start node
%                              ELE_FOR(i,14) = bimoment at end node
%                              ELE_FOR(i,15) = rate of twist at start node
%                              ELE_FOR(i,16) = rate of twist at end node
%       AFLAG            ==  logical flag to indicate if a successful
%                            analysis has been completed
%                              AFLAG = 1     Successful
%                              AFLAG = 0     Unstable Structure
%                              AFLAG = inf   No analysis code available
%
%
%       Version 1.0/Student's Initials/Date of Modification
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Start by defining all output arrays to be empty
%
	DEFL=[]; REACT=[]; Ele_Forces=[];
%  
nnodes;
coord;
concen;
fixity;
nele;
ends;
A;
Izz;
Iyy;
J;
Cw;
IsSym;
Ysc;
Zsc;
Betay;
Betaz;
Betaw;
Zzz;
Zyy;
Ayy;
Azz;
E;
v;
Fy;
YldSurf;
Wt;
webdir
beta_ang
w;
thermal;
truss;
anatype;
	AFLAG = 1;
%
%  STUDENT NOTE:
%     In order for this routine to become fully active AFLAG
%     must be changed.
%
%
%  Student's code starts here...
%
%
%
%  Good luck CE Student!!!
%Assign Nodal Degrees of Freedom to each Node
node_dofs=zeros(nnodes,6);   %node_dofs (creates a zeros matrix to store the dofs numberings at each node)
for i=1:nnodes
    for j=1:1:6
        node_dofs(i,j)=(i-1)*6+j;
    end
end
disp('node_dofs');disp(node_dofs);
% Generating array of element DOF numbers
memb_id=zeros(nele,12);  %memb_id (creates a zeros matrix for the dof for each element)
for i=1:nele
    startNodeNumber=ends(i,1);
    endNodeNumber=ends(i,2);
    for j=1:1:6
        memb_id(i,j)=(startNodeNumber-1)*6+j;
        memb_id(i,j+6)=(endNodeNumber-1)*6+j;
    end
end
disp('memb_id');disp(memb_id);

% Creating applied nodal load vector
appliedNodalLoads=zeros(nnodes*6,1);  %appliedNodalLoads (it's a column vector to store the applied loads at the nodes)
for i=1:nnodes
    appliedNodalLoads((i-1)*6+1:i*6)=concen(i,1:6);
end
disp('appliedNodalLoads');disp(appliedNodalLoads);

%length of the element/member
L = zeros(nele,1);  % L is a variable to store the length of each element
for i = 1:nele
    L(i,1) = norm(coord(ends(i,2),:)-coord(ends(i,1),:));
end
disp('L');disp(L);

%array for free, support and displaced 
fixity_tran = fixity';
D = fixity_tran(:); % D is a column vector (nx1) for the fixity matrix
freeDegreeOfFreedom = find(isnan(D));  %freeDegreeOfFreedom is a variable that finds the NaN values in the D matrix
SupportedDegreeOfFreedom = find(D==0);  %SupportedDegreeOfFreedom is a variable that finds the 0 values in the D matrix
DisplacedDegreeOfFreedom = find(D ~=0 & ~isnan(D)); %SupportedDegreeOfFreedom is a variable that finds not 0 and NaN values

disp('freeDegreeOfFreedom'),disp(freeDegreeOfFreedom)
disp('SupportedDegreeOfFreedom'),disp(SupportedDegreeOfFreedom)
disp('DisplacedDegreeOfFreedom'),disp(DisplacedDegreeOfFreedom)



%fixed end forces
fixed_end_forces = zeros(6*nnodes, 1);

for i = 1:nele
    %memberlocalFEF is a variable for the FEF in the local, and computedFEF
    %is a function that calculate the FEF for each element
    memberlocalFEF = computedFEF(w(i,:), L(i,1));
    memberlocalFEF = memberlocalFEF';

    % Calculate element transformation matrix
    etran_result = godfred_varun_etran(coord(ends(i,1),:), coord(ends(i,2),:), webdir(i,:));

    % Transform the member local FEF to global
    memberglobalFEF = etran_result' * memberlocalFEF;
    memberglobalFEF = memberglobalFEF';

    % Update fixed_end_forces using a loop for linear indexing
    for j = 1:length(memb_id(i,:))
        fixed_end_forces(memb_id(i,j), 1) = memberglobalFEF(j) + fixed_end_forces(memb_id(i,j), 1);
    end
end

disp('fixed_end_forces'),disp(fixed_end_forces);

%global stiffness matrix
k_global = zeros(6*nnodes, 6*nnodes); %k_global stores the global matrix for the structure
for i = 1:nele
    % Calculate element stiffness matrix
    estiff_result = godfred_varun_estiff(A(i,1), Izz(1,i), Iyy(1,i), J(i,1), Ayy(1,i), Azz(1,i), E(i,1), v(i,1), L(i,1));
    
    % Calculate element transformation matrix
    etran_result = godfred_varun_etran(coord(ends(i,1),:), coord(ends(i,2),:), webdir(i,:));

    %combined matrix
    k_combined = etran_result'*estiff_result*etran_result; %k_combined calculates the stiffness matrix 
    % for each element (T'KT - from notes)
    
    % Update the global matrices
    k_global(memb_id(i,:),memb_id(i,:)) = k_combined + k_global(memb_id(i,:),memb_id(i,:));
    
    
end

disp('k_global');disp(k_global);

%extracting the free, supported and displaced stiffness matrix
kff = k_global(freeDegreeOfFreedom,freeDegreeOfFreedom); %kff is the stiffness matrix for free dofs
kfn = k_global(freeDegreeOfFreedom,DisplacedDegreeOfFreedom); %kfn is the stiffness matrix for free and displaced dofs
ksf = k_global(SupportedDegreeOfFreedom,freeDegreeOfFreedom); %ksf is the stiffness matrix for the supported and free dofs
ksn = k_global(SupportedDegreeOfFreedom,DisplacedDegreeOfFreedom); %ksn is the stiffness matrix for the supported and displaced dofs

disp('kff');disp(kff);
disp('kfn');disp(kfn);
disp('ksf');disp(ksf);
disp('ksn');disp(ksn);

%fixed end forces for supported, free and displaced
FEFF = fixed_end_forces(freeDegreeOfFreedom); %FEFF is a variable that extracts the FEF for the free dofs
FEFN = fixed_end_forces(DisplacedDegreeOfFreedom); %FEFN is a variable that extracts the FEF for the displaced dofs
FEFS = fixed_end_forces(SupportedDegreeOfFreedom); %FEFS is a variable that extracts the FEF for the supported dofs

disp('FEFF');disp(FEFF);
disp('FEFN');disp(FEFN);
disp('FEFS');disp(FEFS);

%nodal forces for supported, free and displaced
NF = appliedNodalLoads(freeDegreeOfFreedom); %NF is a variable that extracts the appied nodal loads for free dofs
NN = appliedNodalLoads(DisplacedDegreeOfFreedom); %NN is a variable that extracts the appied nodal loads for displaced dofs
NS = appliedNodalLoads(SupportedDegreeOfFreedom); %NF is a variable that extracts the appied nodal loads for supported dofs

disp('NF');disp(NF);
disp('NN');disp(NN);
disp('NS');disp(NS);

%Displacement

delta_F = kff\(NF - FEFF - kfn*D(DisplacedDegreeOfFreedom,1)); %delta_F is the displacement at the free dofs

disp('delta_F');disp(delta_F);

assignin('base','Displacement_ff',delta_F)

delta_N = D(DisplacedDegreeOfFreedom); %delta_N is the displacement for the displaced dofs
disp('delta_N');disp(delta_N);

delta_S = D(SupportedDegreeOfFreedom); %delta_S is the displacement for the supported dofs
disp('delta_S');disp(delta_S);

delta_D = zeros(6*nnodes,1); %delta_D is column vector (nx1) that combines the displacement for free, supported and displaced
delta_D(freeDegreeOfFreedom,1)=delta_F;
delta_D(DisplacedDegreeOfFreedom,1)=delta_N;
disp('delta_D');disp(delta_D);

%reactions at the support
R = FEFS + ksf*delta_F + ksn*delta_N; %R is the variable for the reaction at the support
disp('R');disp(R);
assignin('base','reactions',R)

% Solving for member forces
Ele_Forces=zeros(nele,12);
for i=1:nele
    k_ele_local=godfred_varun_estiff(A(i),Izz(i),Iyy(i),J(i),Ayy(i),Azz(i),E(i),v(i),L(i)) ;
    gamma=godfred_varun_etran(coord(ends(i,1),:),coord(ends(i,2),:),webdir(i,:));

% Extract the vector of global displacements associated with element i
    Dele_global=delta_D(memb_id(i,:),1);
    Dele_local=gamma*Dele_global;
memberLocalFEF=computedFEF(w(i,:),L(i));
memberLocalFEF=memberLocalFEF';
localMemberForces=k_ele_local*Dele_local+memberLocalFEF;
Ele_Forces(i,:)=localMemberForces;
end

assignin('base','member_rxn',Ele_Forces)

re = k_global*delta_D; %re is a variable for the reaction matrix of the structure
assignin('base','react',re)
%Post-processing
DEFL = reshape(delta_D, 6, nnodes)'; %this code sends to MASTAN the displacements at each node
REACT = reshape(re, [6, nnodes])'; % this code sends to MASTAN the reactions at each node
