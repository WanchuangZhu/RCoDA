function y=ncintegnew(x,NCserial)
% ncinteg is a function which integates the spline result. x is the upper
% limit of the integration, the lower limit is set defaultly as the begining point of the
% spline.
% file is the datafile need to be import


 y=0;
fit=csaps(NCserial(:,1),NCserial(:,2));
coef=fit.coefs;
breaks=fit.breaks;breaks(length(breaks))=[];
dat=[coef,breaks'];
for i=1:size(dat,1)
    fun=@(z) dat(i,4)+dat(i,3)*(z-dat(i,5))+dat(i,2)*(z-dat(i,5)).^2+ ...
        dat(i,1)*(z-dat(i,5)).^3;
    func_cell{i}=fun;
end
breaks=fit.breaks;
sig=breaks>x;
[~,b]=find(sig==1);
stp=b(1);
if x<=breaks(2)
    y=integral(func_cell{1},breaks(1),x);
else
    for i=1:(stp-2)
        y=y+integral(func_cell{i},breaks(i),breaks(i+1));
    end
    y=y+integral(func_cell{stp-1},breaks(stp-1),x);
end