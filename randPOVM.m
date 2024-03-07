function [a,p]=randPOVM(n)
    a=randn(3,n);
    b=null(a);
    a=a*diag(b(:,1));
    t=real(trace(sqrt(a'*a)));
    a=a'*2/t;
    p=diag(sqrt(a*a'));
    for i=1:n
        pinv(i)=1/p(i);
    end
    a=diag(pinv)*a;
end