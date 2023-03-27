close all;

% running the code


dt = 1;

%initial condition
points = 10;
index = randi([1,120*65],1,points)
a_init = zeros(120,65);
for pt = 1:points;
    a_init(index(pt)) = 5;

end

s_init = 3*ones(120,65);
y_init = zeros(120,65);


% initialisation
a = a_init;
s = s_init;
y = y_init;

parameters = [0.015,0.025,0,0,0.1;0.03,0.0025,0.00075,0.00225,20;0,0.03,0.003,0.00015,22];%giraffe coat parameters

for i = 1:10000  
slopes = animal_coats(a,s,y,parameters);

a = a + slopes{1}*dt;
s = s + slopes{2}*dt; 
y = y + slopes{3}*dt;

disp(i)  
end

%plotting portion
figure
h = heatmap(y)
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
grid off
title(['Giraffe Coat']);
saveas(gcf,'giraaffe.png')



% A function that ensures periodic boundary conditions on a matrix 

function wrapd = WrapMatrix(x,i,j)
  
    m = size(x,1);
    n = size(x,2);

    if i > m;
        i = 1;
    elseif i<1;        
         i = m;
    end

    if j > n;        
        j = 1;
    elseif j<1;        
        j = n;
    end

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






%params = [D_a,rho_a,mu_a,sig_a,k_a;D_s,rho_s,mu_s,sig_s,k_s;D_y,rho_y,mu_y,sig_y,k_y]

function slopes = animal_coats(a,s,y,para)

da_dt = (para(1,1) * Lap(a)) + para(1,2).*((((a.^2).*s)./(1 + para(1,5).*a.^2)) - a) ;
ds_dt = (para(2,1) * Lap(s)) + para(2,4)./(1 + (para(2,5).*y)) - (para(2,2).*(a.^2).*s)./(1 + (para(1,5).*(a.^2))) - para(2,3).*s;
dy_dt = para(3,2).*((y.^2)./(1 + (para(3,5).*(y.^2)))) - para(3,3).*y + para(3,4).*a;

slopes = {da_dt,ds_dt,dy_dt};
end
