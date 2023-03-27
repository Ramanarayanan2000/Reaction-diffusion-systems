close all;

dt = 1;
C0 = 5;
field = 50;% size of the grid 
fluct1 = rand(field);
fluct2 = rand(field);

%initial conditions
a_init = fluct1; 
h_init = fluct2;

a = zeros(field,field);
h = zeros(field,field);


% initialisation
a = a_init;
h = h_init;

parameters = [0.005,0.01,0.01,0,0,0.2,0.02,0.02,0];%the model parameters

%solving the model
for i = 1:10000  
slopes = activator_inhibitor(a,h,parameters);
a = a + slopes{1}*dt;
h = h + slopes{2}*dt;  
disp(i)  
end


%plotting portion
figure
h = heatmap(a);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
grid off
title(['spots using Activator-inhibitor system']);
saveas(gcf,'Spot pattern using activator inhibitor system.png')



% A function that ensures periodic boundary conditions on a matrix 

function wrapd = WrapMatrix(x,i,j)
    %wrapd = x;
    m = size(x,1);
    n = size(x,2);
    
    %i = int8;
    %j = int8;

    if i > m;
        %i = rem(i,m);
        i = 1;
    elseif i<1;
        %i = m + i
       % i = abs(fix(i/m))*m + i + m + 1;
         i = m;
    end

    if j > n;
        %j = rem(j,n);
        j = 1;
    elseif j<1;
        %j = abs(fix(j/n))*n + j + n + 1; 
        j = n;
    end

    %wrapd(i,j) = x(i,j);
    wrapd = x(i,j);
end    

% the laplacian
function grad2 = Lap(x)

    grad2 = zeros(size(x));
    h = 1;
    for i = 1:size(x,1);
        for j = 1:size(x,2);

            grad2(i,j) = WrapMatrix(x,i-1,j) + WrapMatrix(x,i+1,j);
            grad2(i,j) = grad2(i,j) + WrapMatrix(x,i,j-1) + WrapMatrix(x,i,j+1);
            grad2(i,j) = grad2(i,j) - 4*(x(i,j));
            grad2(i,j) = grad2(i,j)/h^2;
        end
    end

end



% model for activator inhibitor
% parameters = [D_a,rho_a,mu_a,sig_a,k_a,D_h,rho_h,mu_h,sig_h]
function slopes = activator_inhibitor(a,h,para)

    da_dt = (para(1) * Lap(a)) + para(2) * ((a.^2)./(1. + (para(5) * a.^2)))./h - (para(3).*a) + para(4);
    dh_dt = (para(6) * Lap(h)) + (para(7) * (a.^2)) - para(8).*h  + para(9);
    slopes = {da_dt,dh_dt};
end

