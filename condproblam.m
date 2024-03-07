function condproblam=condproblam(A,Rot,lambda)
p=diag(sqrt(A*A'));

%The array X labels the eight octants of the rotated coordinate frame
X= [1  1  1  1 -1 -1 -1 -1;
    1  1 -1 -1  1  1 -1 -1;
    1 -1  1 -1  1 -1  1 -1];

%The array "v" are the eight vertices of the cube
%(the label corresponds to the column of X)
v=Rot*X;

%Theta is the array values of p_i \Theta(v_{s_xs_ys_z} \cdot a_i)
Theta=(A*v+abs(A*v))/2;
f=sum(Theta,1);
alpha=p/2-sum(Theta,2)/8;
condprob=Theta+alpha*(1-f)/sum(alpha);

%normalize lambda
lambda=lambda/sqrt(lambda*lambda');
lambdaprime=Rot'*lambda';

sgnlam=sign(lambdaprime);

help=sgnlam'*X;
min=0;

for i=1:8
    if help(i)>min
        min=help(i);
        condproblam=condprob(:,i);
    end
end

end