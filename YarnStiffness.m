function [Syarn, Cyarn, braidName] = YarnStiffness(angle, r0, a, beta, braidType, S, C, int)

braidAngle = (angle)*(pi/180);
alpha = braidAngle;
Lund = r0 * beta / sin(braidAngle);

%% **************************************************************************
%Create a sinusoidal equation for each braid type using input from braid
%geometry.  The sinusoidal equation will be used to determine the
%undulations that occur to the braid.  Undulations will affect the in-plane
%and out of plane mechanical properties of the braid.
%**************************************************************************
%% Diamond Braid Pattern
    if braidType == 1
        
        
        Lund = r0 * beta / sin(braidAngle); % GM confirmed May 15, 2015
        braidName = 'Diamond';
        npts = 1000;
        x = linspace(0,Lund,npts);
        xSpacing = (Lund - 0 ) / (npts - 1);
        k = 2*pi / (Lund);
        z = a*0.5*sin(k*x);
        zprime = a*0.5*k*cos(k*x);       

    end

%% Regular Braid Pattern
%GM- changed undulation angle from cos to sin- braid angle measured from
%longitudinal axis- Feb 3, 2014
    if braidType == 2
        braidName = 'Regular';
        
        Lund = r0 * beta / sin(braidAngle); % GM confirmed May 15, 2015
        %1) Create a piecewise function for the regular braid
        npts = 1000;
        %regular braid total length = 4p == 2*Lund
        x = linspace(0, 2 * Lund, npts);
        xSpacing = (2 * Lund - 0 ) / (npts - 1);
        j=1;
        i = 1;
        k1 = 1.0;
        z = zeros(npts,1);
        while(x(i) < 2 * Lund)
        %while(x(i) < j*(k1+1)*beta*(r0/sin(alpha)))
            if abs(x(i)) >= (0 + ((j-1)*(k1+1)*beta)*(r0/sin(alpha))) && abs(x(i)) < (k1*beta*0.5 + ((j-1)*(k1+1)*beta))*(r0/sin(alpha))
                z(i) = 0.5*a;
                zprime(i) = 0;
            elseif abs(x(i)) >= (k1*beta*0.5 + ((j-1)*(k1+1)*beta))*(r0/sin(alpha)) && abs(x(i)) <= ((k1+1)*beta*0.5 + (j-1)*(k1+1)*beta)*(r0/sin(alpha))
                z(i) = -0.5*a*sin((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
                zprime(i) = -(0.5*a*2*pi*sin(alpha)/ (r0*beta)) * cos((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
            elseif abs(x(i)) >= ((k1+1)*beta*0.5 + ((j-1)*(k1+1)*beta))*(r0/sin(alpha)) && abs(x(i)) <= ((2*k1 + 1)*beta*0.5 + ((j-1)*(k1+1)*beta))*(r0/sin(alpha))
                z(i) = -0.5*a;
                zprime(i) = 0;
            %elseif abs(x(i)) >= ((2*k1 + 1)*beta*0.5 + ((j-1)*(k1+1)*beta))*(r0/sind(alpha))  && abs(x(i)) <= ((k1+1)*beta + ((j-1)*(k1+1)*beta))*(r0/sin(alpha))
            else %gm added May 19- ensures that x and z are the size length/same number of elements
                z(i) = 0.5*a*sin((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
                zprime(i) = (0.5*a*2*pi*sin(alpha)/ (r0*beta)) * cos((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
    
            end
            i = i+1;
        end 
    end
%% Hercules Braid Pattern
if braidType == 3
    braidName = 'Hercules';
    %1) Create a piecewise function for the regular braid
    npts = 1000;
    %hercules braid total length = 6p == 3*Lund
    x = linspace(0, 3 * Lund, npts);
    xSpacing = (3 * Lund - 0 ) / (npts - 1);
    i = 1;
    z = zeros(npts,1);
    while(x(i) < 3 * Lund)
        if abs(x(i)) >= 0 && abs(x(i)) < (beta*0.5)*(r0/sin(alpha))
            z(i) = 0.5*a;
            zprime(i) = 0;
        elseif abs(x(i)) >= (beta*0.5)*(r0/sin(alpha)) && abs(x(i)) < (2*beta*0.5)*(r0/sin(alpha))
            z(i) = 0.5*a;
            zprime(i) = 0;
        elseif abs(x(i)) >= (2*beta*0.5)*(r0/sin(alpha)) && abs(x(i)) < (3*0.5*beta)*(r0/sin(alpha))
            z(i) = 0.5*a*sin((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
            zprime(i) = (0.5*a*2*pi*sin(alpha)/ (r0*beta)) * cos((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
        elseif abs(x(i)) >= (3*0.5*beta)*(r0/sin(alpha)) && abs(x(i)) < (4*beta*0.5)*(r0/sin(alpha))
            z(i) = -0.5*a;
            zprime(i) = 0;
        elseif abs(x(i)) >= (4*beta*0.5)*(r0/sin(alpha)) && abs(x(i)) < (5*0.5*beta)*(r0/sin(alpha))
            z(i) = -0.5*a;
            zprime(i) = 0;
        %elseif abs(x(i)) >= (5*0.5*beta)*(r0/sin(alpha)) && abs(x(i)) < (6*beta*0.5)*(r0/sin(alpha))
        else %gm added May 19- ensures that x and z are the size length/same number of elements
             z(i) = 0.5*a*sin((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
             zprime(i) = (0.5*a*2*pi*sin(alpha)/ (r0*beta)) * cos((2*pi*sin(alpha)*x(i)) / (r0 * beta) + pi*0.5);
        end
        i = i+1;
        
    end
end
%% Diamond / Triaxial = 4
if braidType == 4
    braidName = 'Diamond / Triaxial';
    
    
    %a - braid yarn thickness
    %aA - axial yarn thickness
    aA = a; % assume axial yarns are same size as braid yarns for now...
    x = linspace(0,Lund,1000);
    k = 2*pi / (Lund);
    z = (a*0.5*aA*0.5)*sin(k*x);
    zprime = a*k*cos(k*x);
    
    
end
%% Regular / Triaxial = 5
if braidType == 5
    braidName = 'Regular / Triaxial';
    
    
    
end
%% No Undulation
if braidType == 6
    braidName = 'CLPT';
end

% figure; plot(x,z)
% xlabel('Position x (mm)');
% ylabel('Position z (mm)')
%% **************************************************************************
%Determine crimp angle along fiber strand path for each braiding type
%braid crimp angle is used to determine the effect of undulations on the
%braid effective stiffness
%**************************************************************************
%% Diamond Braid
if braidType == 1
tanB = a*k*0.5*cos(k*x);
cosB = 1 ./ sqrt(1 + tanB.^2);
sinB = tanB ./ sqrt(1 + tanB.^2);
end

%% Regular Braid
if braidType ==2
    tanB = zprime;
    cosB = 1 ./ sqrt(1 + tanB.^2);
    sinB = tanB ./ sqrt(1 + tanB.^2);
    
    x = x(2:end);
    
end
%% Hercules Braid
if braidType == 3
    tanB = zprime;
    cosB = 1 ./ sqrt(1 + tanB.^2);
    sinB = tanB ./ sqrt(1 + tanB.^2);
    x = x(2:end);
end

%% Diamond Triaxial Braid
if braidType == 4
    tanB = zprime;
    cosB = 1 ./ sqrt(1 + tanB.^2);
    sinB = tanB ./ sqrt(1 + tanB.^2);
    
end

%% Regular Triaxial Braid
if braidType == 5
end

%% No Undulations
if braidType == 6
end

%% Evaluate the effective stiffness due to undulations in the 123 coordinate system
%performs coordinate system transformation along yarn undulation and
%determines the average stiffness and compliance for one yarn undulation
Syarn = CS123Stiffness(cosB, sinB, tanB, x, xSpacing, S,int);
Cyarn = CS123Compliance(cosB, sinB, tanB, x, xSpacing, C,int);