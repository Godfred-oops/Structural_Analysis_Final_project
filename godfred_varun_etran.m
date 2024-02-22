%% ABABIO GODFRED OPOKU AND VARUN 
function [gamma] = godfred_varun_etran(coordi, coordj, webdir)
x_vect = coordj-coordi;  %difference betwen the i and j nodal coordinates for the element
x_vect = x_vect/norm(x_vect);  %unit vector of the element
z_vect = cross(x_vect,webdir);   %cross product of row one (unit vector) and row two (web direction vector)
z_vect = z_vect/norm(z_vect);  %unit vector of the cross product 
g = [x_vect;webdir;z_vect];  %forming the 3x3 matrix
d = zeros(12);
d(1:3,1:3) =g;
d(4:6,4:6) =g;
d(7:9,7:9) =g;
d(10:12,10:12) =g;

gamma = d;
end