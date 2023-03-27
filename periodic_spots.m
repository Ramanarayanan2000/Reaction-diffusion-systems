close all;

dt = 1;

field = 8;% size of the grid 
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

parameters = [0.006,0.01,0.01,0.001,0,0.2,0.02,0.02,0];%the model parameters

frame_repository = cell(1,44);
frame = zeros(52);
%solving the model
for i = 1:89000 

%updating the field after every 2000th iteration
if mod(i,2000) == 0;
    dimension = size(a);
    frame(1:dimension(1),1:dimension(2)) = a;
    frame_repository{i/2000} = frame;
    a = [a;rand(1,dimension(2))];
    a = [a,rand(dimension(1)+1,1)];
    
    h = [h;rand(1,dimension(2))];
    h = [h,rand(dimension(1)+1,1)];
end


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
title(['Periodic Spots']);
saveas(gcf,'Periodic spots.png')

% generating gif

for indx = 1:44;
    disp(indx)
    h = heatmap(frame_repository{indx});
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    colorbar off
    grid off
    exportgraphics(h,"Periodic spots with growing field.gif","Append",true)
end


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



% model for activator inhibitor
% parameters = [D_a,rho_a,mu_a,sig_a,k_a,D_h,rho_h,mu_h,sig_h]
function slopes = activator_inhibitor(a,h,para)

    da_dt = (para(1) * Lap(a)) + para(2) * ((a.^2)./(1. + (para(5) * a.^2)))./h - (para(3).*a) + para(4);
    dh_dt = (para(6) * Lap(h)) + (para(7) * (a.^2)) - para(8).*h  + para(9);
    slopes = {da_dt,dh_dt};
end
