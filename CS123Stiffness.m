function Sout = CS123Stiffness(cosB, sinB, tanB, x, xSpacing, S,int)


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

Sprime(:,:,i) = transpose(T) * S * T;

end
 

Sint11 = xSpacing * trapz(squeeze(Sprime(1,1,:)));
Sint12 = xSpacing * trapz(squeeze(Sprime(1,2,:)));
Sint13 = xSpacing * trapz(squeeze(Sprime(1,3,:)));
Sint14 = xSpacing * trapz(squeeze(Sprime(1,4,:)));
Sint15 = xSpacing * trapz(squeeze(Sprime(1,5,:)));
Sint16 = xSpacing * trapz(squeeze(Sprime(1,6,:)));

Sint21 = xSpacing * trapz(squeeze(Sprime(2,1,:)));
Sint22 = xSpacing * trapz(squeeze(Sprime(2,2,:)));
Sint23 = xSpacing * trapz(squeeze(Sprime(2,3,:)));
Sint24 = xSpacing * trapz(squeeze(Sprime(2,4,:)));
Sint25 = xSpacing * trapz(squeeze(Sprime(2,5,:)));
Sint26 = xSpacing * trapz(squeeze(Sprime(2,6,:)));

Sint31 = xSpacing * trapz(squeeze(Sprime(3,1,:)));
Sint32 = xSpacing * trapz(squeeze(Sprime(3,2,:)));
Sint33 = xSpacing * trapz(squeeze(Sprime(3,3,:)));
Sint34 = xSpacing * trapz(squeeze(Sprime(3,4,:)));
Sint35 = xSpacing * trapz(squeeze(Sprime(3,5,:)));
Sint36 = xSpacing * trapz(squeeze(Sprime(3,6,:)));

Sint41 = xSpacing * trapz(squeeze(Sprime(4,1,:)));
Sint42 = xSpacing * trapz(squeeze(Sprime(4,2,:)));
Sint43 = xSpacing * trapz(squeeze(Sprime(4,3,:)));
Sint44 = xSpacing * trapz(squeeze(Sprime(4,4,:)));
Sint45 = xSpacing * trapz(squeeze(Sprime(4,5,:)));
Sint46 = xSpacing * trapz(squeeze(Sprime(4,6,:)));

Sint51 = xSpacing * trapz(squeeze(Sprime(5,1,:)));
Sint52 = xSpacing * trapz(squeeze(Sprime(5,2,:)));
Sint53 = xSpacing * trapz(squeeze(Sprime(5,3,:)));
Sint54 = xSpacing * trapz(squeeze(Sprime(5,4,:)));
Sint55 = xSpacing * trapz(squeeze(Sprime(5,5,:)));
Sint56 = xSpacing * trapz(squeeze(Sprime(5,6,:)));

Sint61 = xSpacing * trapz(squeeze(Sprime(6,1,:)));
Sint62 = xSpacing * trapz(squeeze(Sprime(6,2,:)));
Sint63 = xSpacing * trapz(squeeze(Sprime(6,3,:)));
Sint64 = xSpacing * trapz(squeeze(Sprime(6,4,:)));
Sint65 = xSpacing * trapz(squeeze(Sprime(6,5,:)));
Sint66 = xSpacing * trapz(squeeze(Sprime(6,6,:)));


Sout = (1/(x(end))) * [Sint11 Sint12 Sint13 Sint14 Sint15 Sint16;
                     Sint21 Sint22 Sint23 Sint24 Sint25 Sint26;
                     Sint31 Sint32 Sint33 Sint34 Sint35 Sint36;
                     Sint41 Sint42 Sint43 Sint44 Sint45 Sint46;
                     Sint51 Sint52 Sint53 Sint54 Sint55 Sint56;
                     Sint61 Sint62 Sint63 Sint64 Sint65 Sint66];



