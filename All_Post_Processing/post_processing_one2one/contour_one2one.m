close all
clc
clear

out1name = 'LOO_contoured_resid_all';
title_opt = 'NEW Optimal Residuals'
% load up the residuals data, noting columns as follows:
% | meas | modeled
infile = 'LOO_MEAS_MOD.dat';
indata = log10(load(infile));


% find the data range
dmin = min(min(indata));
dmax = max(max(indata));
dmin = log10(0.05);
dmax = 0.2;

% discretize the results
% nvals indicates the number of increments to bin into
nvals = 31;
x = linspace(dmin,dmax,nvals);
y = x;
dx = (x(2)-x(1))/2;

[X,Y]  = meshgrid(x,y);
Zmeas  = zeros(size(X));
Zstart = Zmeas;
Zopt   = Zmeas;
% bin the values
for i = 1:numel(X)
    for j = 1:length(indata)
        if ((indata(j,2)>X(i)-dx) && (indata(j,2)<= X(i)+dx) && ...
            (indata(j,1)>Y(i)-dx) && (indata(j,1)<= Y(i)+dx));
            Zstart(i) = Zstart(i)+1;
        end
    end   
end
dlimits = [0:2:20];
limiting_val = 5;
out1 = figure(1);
hold on
contourf(X,Y,Zstart)%,dlimits);
axis([dmin dmax dmin dmax])
colormap jet
title('All Residuals')
plot([dmin,dmax],[dmin,dmax],'w')
colorbar
xlabel('Modeled Hg')
ylabel('Observed Hg')
print(out1,'-depsc',out1name)

