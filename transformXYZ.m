
function [Sout, Cout] = transformXYZ(alpha, S, C)

c = cos(alpha);
s = sin(alpha);
%Rotation matrix- transformation about the z axis
T = [c^2 s^2 0 0 0 2*c*s;...
     s^2 c^2 0 0 0 -2*c*s;...
     0   0   1 0 0  0;...
     0   0   0 c -s 0;...
     0   0   0 s  c 0;...
     -c*s c*s 0 0 0 (c^2-s^2)];
 
 
 R = [1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 2 0 0;
      0 0 0 0 2 0;
      0 0 0 0 0 2];
 
 %Byun transformation
 Sout = transpose(T) * S * T;
  
 %transform eqn and reuters matrix Jones pg 77
 Cout = inv(T) * C * R * T * inv(R);

end