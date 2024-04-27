function PR=msefun(P,R)
[M,N]=size(P);
P=P(:);
R=R(:);
PR=sum((P-R).^2)/(M*N);
end