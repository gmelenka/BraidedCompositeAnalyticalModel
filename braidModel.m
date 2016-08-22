function [Ex, Ey, Ez, GxyCombined, GyzCombined, GzxCombined] = braidModel(S, Sm, angle, n, r0, a, b, beta, braidType)
%**************************************************************************
%% Tubular braided composite analytical model
%Author: Garrett Melenka
%Email: gmelenka@ualberta.ca
%date: August 12, 2016
%Description:
%This program calculates the elastic constants for tubular braided
%composites.  The model allows the user to input geometric parameters of
%the composite braid as well as the material properties of the matrix and
%fiber used.  Using these input values the elastic constants are computed
%while accounting for the braiding angle and undulations that occur during
%the braiding process.

%The following braid patterns can be analyzed using this method:
%1) diamond braid pattern (1/1)
%2) regular braid pattern (2/2)
%3) Hercules braid pattern (3/3)
%**************************************************************************

% %% Initial parameters for braid geometry
% %*************************************************************************
% a = 0.38; %yarn thickness
% b = 3.1; %yarn width
% D = 11.1; %mandrel diameter
% t = 2*a;
% braidType = 3; %type of braid see above
% 
% R = D/2; % mandrel size mm
% r0 = R+a; % nomial braid radius mm
% 
% Define the braid angle between 30 and 60 degrees
% angle = linspace(30,60,100);
% 
% %Define braiding machine parameters
% n = 18; %number of carriers
% Nc = 2*n; % total number of carriers
% beta = 2*pi / n; % braid shift angle (rad)
% 
% 
% **************************************************************************
% % Material Properties
% Matrix material properties
% Em = 3.5;
% Gm = 1.3;
% num = 0.3;
% 
% Vf = 0.6;
% Vv = 4.35 / 100;
% %Vv = 0;
% Vm = 1 - Vf - Vv
% *************************************************************************
% %Fiber, Matrix and Fiber+Matrix Material Properties from Carey/Ayranci
% %Thesis
% Ef1 = 130;
% Ef2 = 7.3;
% Ef3 = Ef2;
% Gf12 = 2.86;
% Gf13 = Gf12;
% nuf12 = 0.35;
% nuf13 = nuf12;
% nuf21 = nuf12*(Ef2/Ef1);
% nuf31 = nuf12*(Ef3/Ef1);
% nuf23 = 0.1;
% nuf32 = nuf23*(Ef3/Ef2);
% **************************************************************************
% %Source: Cagri model
% E1 = 79.7;
% E2 = 5.9;
% E3 = E2;
% G12 = 1.5;
% G13 = G12;
% eta23 = (3 - 4*num + (Gm / Gf12)) / (4*(1-num));
% G23 = (Gm*(Vf + eta23*(1-Vf))) / (eta23*(1-Vf) + Vf*(Gm/Gf12));
% nu12 = 0.3;
% nu13 = nu12;
% nu23 = (E2/(2*G23)) - 1;
% nu21 = nu12*(E2/E1);
% nu31 = nu13*(E3/E1);
% nu32 = nu23*(E3/E2);
% 
% initial transversely isotropic compliance matrix for yarns+epoxy
% S = [1/E1 -nu21/E2 -nu31/E3 0 0 0;...
%     -nu12/E1 1/E2 -nu32/E3 0 0 0;...
%     -nu13/E1 -nu23/E2 1/E3 0 0 0;...
%     0 0 0 1/G23 0 0;...
%     0 0 0 0 1/G13 0;...
%     0 0 0 0 0 1/G12];
% 
% 
% initial compliance matrix for epoxy
% Sm = [1/Em -num/Em -num/Em 0 0 0;...
%      -num/Em 1/Em -num/Em 0 0 0;...
%      -num/Em -num/Em 1/Em 0 0 0;...
%       0 0 0 1/Gm 0 0;...
%       0 0 0 0 1/Gm 0;...
%       0 0 0 0 0 1/Gm];
  
C = inv(S); Cm = inv(Sm);  

  

    

%% Calculate effective stiffness for undulating fibers
for int = 1:length(angle)
    
    %Diamond Braid Pattern = 1
    %Regular Braid Pattern = 2
    %Hercules Braid Pattern = 3
     
    %calls yarn stiffness function to account for yarn unduation of braid
    %patterns
    [Sout(:,:,int), Cout(:,:,int), braidName] = YarnStiffness(angle(int), r0, a, beta, braidType, S, C, int);
    

end

%% Transform stiffness matrix to global coordinate system XYZ
for i = 1:length(angle)
    
           
      %determine fiber volume within braid unit cell
      Vfuc = b * n ./ (8 * r0 * cosd(angle));
      
      %once unit cell has become jammed fiber volume is constant  
        for k = 1:length(Vfuc)
   
            if Vfuc(k) >=1
                Vfuc(k) =1;
            end
        end
    VucF = Vfuc(i);
 
    theta = (angle(i)) * (pi/180);
    
    %transforms braid geometry to account for braid angle
    [SxyzPos, CxyzPos] = transformXYZ(theta, Sout(:,:,i), Cout(:,:,i));
    [SxyzNeg, CxyzNeg] = transformXYZ(-theta, Sout(:,:,i), Cout(:,:,i));

    VucM = 1 - VucF;
    
    %Stiffness averaging as per Quek
    Ctotal = VucF*0.5*CxyzPos + VucF*0.5*CxyzNeg + VucM*Cm;
    
    %Compliance averaging
    CtotalCompliance = VucF*0.5*inv(SxyzPos) + VucF*0.5*inv(SxyzNeg) + VucM*Cm;

    Stotal(:,:,i) = inv(Ctotal);
    
    StotalByun(:,:,i) = inv(CtotalCompliance);
    
   

end

%% Determine Elastic Constants in the global coordinate system
Ex = squeeze(1./Stotal(1,1,:));
Ey = squeeze(1./Stotal(2,2,:));
Ez = squeeze(1./Stotal(3,3,:));
Gxy = squeeze(1./Stotal(6,6,:));
Gyz = squeeze(1./Stotal(5,5,:));
Gzx = squeeze(1./Stotal(4,4,:));
vxy = squeeze(-Stotal(1,2,:)./Stotal(1,1,:));

ExCompliance = squeeze(1./StotalByun(1,1,:));
EyCompliance = squeeze(1./StotalByun(2,2,:));
EzCompliance = squeeze(1./StotalByun(3,3,:));
GxyCompliance = squeeze(1./StotalByun(6,6,:));
GyzCompliance = squeeze(1./StotalByun(5,5,:));
GzxCompliance = squeeze(1./StotalByun(4,4,:));
vxyCompliance = squeeze(-StotalByun(1,2,:)./StotalByun(1,1,:));

GxyCombined = (Gxy + GxyCompliance)*0.5;
GyzCombined = (Gxy + GyzCompliance)*0.5;
GzxCombined = (Gxy + GzxCompliance)*0.5;

    

