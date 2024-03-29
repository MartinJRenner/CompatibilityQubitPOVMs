clear

%Write the normalized vectors \vec{a}_i in each row of the following
%array (the example is the SIC POVM from the paper)

a=[         0           0         1;
    sqrt(8)/3           0      -1/3;
   -sqrt(2)/3   sqrt(6)/3      -1/3;
   -sqrt(2)/3  -sqrt(6)/3      -1/3]

p=[0.5;
   0.5;
   0.5;
   0.5]

%The function randPOVM generates random points on the sphere and solves p
%according to \sum_i p_i \vec{a}_i=\vec{0} and \sum_i p_i=2. 4 is the number of outcomes
%here, you can change that number to 5 or more outcomes (less won't work).
[a,p]=randPOVM(4)

%You can do a quick check here that everything is well-define
%(it should be: sum(p)=2; p'*a=[0 0 0]; diag(a*a')=[1 1 1 1])
sum(p);
p'*a;
diag(a*a');

%"A" is defined for convenience:
A=diag(p)*a;

%The function "rotation" discretizes the three Euler angles with precision n
%and outputs the rotation in which the maximum of the eight values of f_a
%is the smallest. The smallest maximal value is denoted as minmax. If this
%value is smaller than 1, you found a suitable rotation, otherwise you can
%increase "n"
n=35;
[Rot,minmax]=rotation(A,n)

%The array X labels the eight octants of the rotated coordinate frame
X= [1  1  1  1 -1 -1 -1 -1;
    1  1 -1 -1  1  1 -1 -1;
    1 -1  1 -1  1 -1  1 -1];

%The array "v" are the eight vertices of the cube
%(the label corresponds to the column of X)
v=Rot*X

%Here, we calculate all the relevant functions of the paper.
%"Theta" is the array with values of p_i \Theta(v_{s_xs_ys_z} \cdot a_i)
Theta=(A*v+abs(A*v))/2;

%"f" is the function f_a for all eight vertices of the cube, the maximum
%of the eight values should be equal to "minmax"
f=sum(Theta,1);
alpha=p/2-sum(Theta,2)/8;

%"condprob" are the conditional probabilities for each octant of the rotated
%frame
condprob=Theta+alpha*(1-f)/sum(alpha)

%Now we can also compute the conditional
%probabilities p(i|a,lambda) for a given lambda (the outcome of the parent
%POVM). The algorithm below checks which of the eight cube vertices is the
%closest (this is equivalent to the octant of lambda in the rotated frame)
%and outputs the corresponding row of condprob.
lambda=[1 1 1];
lambda=lambda/sqrt(lambda*lambda');
[~,j]=max(lambda*v);

condproblam=condprob(:,j)

%Now we can also do a quick check that the conditional probabilities
%simulate the POVM

%The operators G_{s_xs_ys_z} from each octant of the rotatend frame are
%(the first row corresponds to the identity operator).
G(1,1:8)=1/8;
G(2:4,1:8)=v/16;
G

%noisy POVM we want to simulate
%(the first column corresponds to the identity operator)
Anoise(:,1)=p/2;
Anoise(:,2:4)=A/4;

Anoise
Simulation=condprob*G'
