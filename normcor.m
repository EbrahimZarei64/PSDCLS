function nc=normcor(a,b)
    a=double(a(:));
    b=double(b(:));
    a1=sum((a.*a));
    b1=sum((b.*b));
    c1=sum(a.*b);
    c2=sqrt(a1*b1);
    nc=c1/c2;
end