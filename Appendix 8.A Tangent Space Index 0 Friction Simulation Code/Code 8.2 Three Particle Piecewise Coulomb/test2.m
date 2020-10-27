n=1;

sx1err(1)=0;
sx1derr(1)=0;
sx1dderr(1)=0;
sx2err(1)=0;
sx2derr(1)=0;
sx2dderr(1)=0;
syerr(1)=0;
syderr(1)=0;
sydderr(1)=0;

N=39900;

while n<N
    n=n+1;
x1err(n)=norm(x1c(n)-x1(n));
sx1err(n)=sx1err(n-1)+x1err(n);
x1derr(n)=norm(x1dc(n)-x1d(n));
sx1derr(n)=sx1derr(n-1)+x1derr(n);
x1dderr(n)=norm(x1ddc(n)-x1dd(n));
sx1dderr(n)=sx1dderr(n-1)+x1dderr(n);
x2err(n)=norm(x2c(n)-x2(n));
sx2err(n)=sx2err(n-1)+x2err(n);
x2derr(n)=norm(x2dc(n)-x2d(n));
sx2derr(n)=sx2derr(n-1)+x2derr(n);
x2dderr(n)=norm(x2ddc(n)-x2dd(n));
sx2dderr(n)=sx2dderr(n-1)+x2dderr(n);
yerr(n)=norm(yc(n)-y(n));
syerr(n)=syerr(n-1)+yerr(n);
yderr(n)=norm(ydc(n)-yd(n));
syderr(n)=syderr(n-1)+yderr(n);
ydderr(n)=norm(yddc(n)-ydd(n));
sydderr(n)=sydderr(n-1)+ydderr(n);
end

maxx1err=max(x1err);
maxx1derr=max(x1derr);
maxx1dderr=max(x1dderr);
maxx2err=max(x2err);
maxx2derr=max(x2derr);
maxx2dderr=max(x2dderr);
maxyerr=max(yerr);
maxyderr=max(yderr);
maxydderr=max(ydderr);

meanx1err=sx1err(N)/N;
meanx1derr=sx1derr(N)/N;
meanx1dderr=sx1dderr(N)/N;
meanx2err=sx2err(N)/N;
meanx2derr=sx2derr(N)/N;
meanx2dderr=sx2dderr(N)/N;
meanyerr=syerr(N)/N;
meanyderr=syderr(N)/N;
meanydderr=sydderr(N)/N;
