function [EncImg,H1]=Encryption(PlainImg,ExternalKey,H1)
    N=size(PlainImg,1);
    H2=HashFunction(ExternalKey,'SHA-256');
    KeySHA =HashFunction([H1 H2],'SHA-256');
   
    n=size(KeySHA,2)/2;
    HexBin=Hex2Bin(KeySHA,4);
    HexBin=HexBin(:)';
        
    KeyDec=zeros(1,n);
    for i=1:1:n
        KeyDec(i)=bin2dec(HexBin((i-1)*8+1:i*8));
    end
    
    a=bitxor(KeyDec(1),KeyDec(2));
    for i=3:16
        a=bitxor(a,KeyDec(i));
    end
    
    b=bitxor(KeyDec(17),KeyDec(18));
    for i=19:32
        b=bitxor(b,KeyDec(i));
    end
    
    HexBin = HexBin(:)-'0';
    x0 = 0;
    B = HexBin(1:128);
    for i=1:128
        x0 = x0 + (B(i)*(2^(-i) ));   
    end
    
    y0 = 0;
    B = HexBin(129:256);
    for i=1:128
        y0 = y0 + (B(i)*(2^(-i) ));   
    end
    
    Len=N+a+b;
    X=zeros(1,Len);
    X(1)=x0;
    Y=zeros(1,Len);
    Y(1)=y0;
    
    for i=2:Len
        X(i)=mod(1-a*(sin(X(i-1))^2)+Y(i-1),1);
        Y(i)=mod( b*X(i-1),1);
    end
    
    X(1:a+b)=[];
    Y(1:a+b)=[];
    [~,Px]=sort(X);
    [~,Py]=sort(Y);
    
    
    
    Px=(Px-1);
    Px= repmat(Px,N,1);
    Py=(Py-1)';
    Py= repmat(Py,1,N);
    Latin=bitxor(Px,Py);


    LatinKey=mod(Latin,256);
    LatinPer=Latin+1;
    %%
    VDC=zeros(N);
    for i=1:N
        VDC(:,i)=bitxor(PlainImg(LatinPer(:,i),i),LatinKey(LatinPer(:,i),i));
    end
     
    TVHD=VDC';
    LatinKey=LatinKey';
    LatinPer=LatinPer';
    
    EncImg=zeros(N);
    for i=1:N
        EncImg(:,i)=bitxor(TVHD(LatinPer(:,i),i),LatinKey(LatinPer(:,i),i));
    end
end