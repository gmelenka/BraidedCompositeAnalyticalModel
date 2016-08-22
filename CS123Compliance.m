function Cout = CS123Compliance(cosB, sinB, tanB, x, xSpacing, C,int)

for i = 1:length(cosB)
m=cosB(i);
n=sinB(i);
  
T = [m^2 0 n^2 0 2*m*n 0;
     0 1 0 0 0 0;
     n^2 0 m^2 0 -2*m*n 0;
     0 0 0 m 0 -n;
     -m*n 0 m*n 0 m^2-n^2 0;
     0 0 0 n 0 m]; 
 
 R = [1 0 0 0 0 0; 
     0 1 0 0 0 0; 
     0 0 1 0 0 0; 
     0 0 0 2 0 0; 
     0 0 0 0 2 0; 
     0 0 0 0 0 2]; 

 
 

Cprime(:,:,i) = inv(T) * C * R * T * inv(R);


end



Cint11 = xSpacing * trapz(Cprime(1,1,:));
Cint12 = xSpacing * trapz(Cprime(1,2,:));
Cint13 = xSpacing * trapz(Cprime(1,3,:));
Cint14 = xSpacing * trapz(Cprime(1,4,:));
Cint15 = xSpacing * trapz(Cprime(1,5,:));
Cint16 = xSpacing * trapz(Cprime(1,6,:));

Cint21 = xSpacing * trapz(Cprime(2,1,:));
Cint22 = xSpacing * trapz(Cprime(2,2,:));
Cint23 = xSpacing * trapz(Cprime(2,3,:));
Cint24 = xSpacing * trapz(Cprime(2,4,:));
Cint25 = xSpacing * trapz(Cprime(2,5,:));
Cint26 = xSpacing * trapz(Cprime(2,6,:));

Cint31 = xSpacing * trapz(Cprime(3,1,:));
Cint32 = xSpacing * trapz(Cprime(3,2,:));
Cint33 = xSpacing * trapz(Cprime(3,3,:));
Cint34 = xSpacing * trapz(Cprime(3,4,:));
Cint35 = xSpacing * trapz(Cprime(3,5,:));
Cint36 = xSpacing * trapz(Cprime(3,6,:));

Cint41 = xSpacing * trapz(Cprime(4,1,:));
Cint42 = xSpacing * trapz(Cprime(4,2,:));
Cint43 = xSpacing * trapz(Cprime(4,3,:));
Cint44 = xSpacing * trapz(Cprime(4,4,:));
Cint45 = xSpacing * trapz(Cprime(4,5,:));
Cint46 = xSpacing * trapz(Cprime(4,6,:));

Cint51 = xSpacing * trapz(Cprime(5,1,:));
Cint52 = xSpacing * trapz(Cprime(5,2,:));
Cint53 = xSpacing * trapz(Cprime(5,3,:));
Cint54 = xSpacing * trapz(Cprime(5,4,:));
Cint55 = xSpacing * trapz(Cprime(5,5,:));
Cint56 = xSpacing * trapz(Cprime(5,6,:));

Cint61 = xSpacing * trapz(Cprime(6,1,:));
Cint62 = xSpacing * trapz(Cprime(6,2,:));
Cint63 = xSpacing * trapz(Cprime(6,3,:));
Cint64 = xSpacing * trapz(Cprime(6,4,:));
Cint65 = xSpacing * trapz(Cprime(6,5,:));
Cint66 = xSpacing * trapz(Cprime(6,6,:));


 Cout = (1/(x(end))) * [Cint11 Cint12 Cint13 Cint14 Cint15 Cint16;
                     Cint21 Cint22 Cint23 Cint24 Cint25 Cint26;
                     Cint31 Cint32 Cint33 Cint34 Cint35 Cint36;
                     Cint41 Cint42 Cint43 Cint44 Cint45 Cint46;
                     Cint51 Cint52 Cint53 Cint54 Cint55 Cint56;
                     Cint61 Cint62 Cint63 Cint64 Cint65 Cint66];



