function [Rot,minmax]=rotation(A,n)
%Rot is the rotation that minimizes the largest value of
%f_a(v_{s_xs_ys_z}) that is denoted as min. If this value is larger than 1
%try it again with a larger n (or use a more sophisticated optimization).

p=diag(sqrt(A*A'));

%X is the array of the eight vectors v=(s_x, s_y, s_z)
X= [1  1  1  1 -1 -1 -1 -1;
    1  1 -1 -1  1  1 -1 -1;
    1 -1  1 -1  1 -1  1 -1];

n=35;
minmax=5;

for i=0:n
    for j=0:n
        for k=0:n

%a, b, and c are the three euler angles
a=i/n*2*pi;
b=j/n*2*pi;
c=k/n*2*pi;

R1=[      1        0       0;
          0   cos(a)  sin(a);
          0  -sin(a)  cos(a)];

R2=[ cos(b)   sin(b)       0;
    -sin(b)   cos(b)       0;
          0        0       1];

R3=[      1        0       0;
          0   cos(c)  sin(c);
          0  -sin(c)  cos(c)];

R=R1*R2*R3;

%Theta is the array of values of p_i*\Theta(v_{s_xs_ys_z} \cdot a_i)
%the entries of f are the functional values of f_a(v_{s_xs_ys_z}) in this
%rotated frame
v=R*X;
Theta=(A*v+abs(A*v))/2;
f=sum(Theta,1);

if max(f)<minmax
    minmax=max(f);
    Rot=R;
end

        end
    end
end
end