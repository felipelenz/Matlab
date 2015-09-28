%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Felipe Lenz Carvalho
% This code models a 2-D dipole in with a 2-D Perfect Mathing Layers (PML) 
% in both ends (the goal is to attenuate the waves when they hit the PML)
% the top and bottom boundaries are PECs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we implement a 2-D version of the FDTD equations, as
% discretized by Tavlove and Hagness [2000] bounded by two PMLs along the 
% x axis and two PECs along the y axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We start by supplying our constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c = 1/sqrt(mu0*eps0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We then provide our S parameter, calculate our range for x, and use them
% to calculate delta_t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 1/sqrt(2);
f0 = 1e9;%Hz
lambda0 = c/f0;

delta_x = 1e-3;%meters
delta_y = delta_x;
N_PML = 100;
PML_sigma = 1;
PML_depth = N_PML*delta_x;

min_x = 0;%meters
max_x = PML_depth + 0.5 + PML_depth;%meters
x = (min_x:delta_x:max_x)';

min_y = 0;%meters
max_y = 0.1;%meters
y = (min_y:delta_y:max_y)';

[Y,X] = meshgrid(y,x);

source_x = PML_depth + 0.25;
source_y = 0.05;

delta_t = delta_x*S/c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We provide the permittivity, permeability, and conductivity as a function
% of space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = eps0.*ones(size(X));
mu = mu0.*ones(size(X));
sigma = zeros(size(X));
sigma_star = zeros(size(X));

%First PML
s = find((x >= 0) & (x < PML_depth));
eps(s) = eps0;
mu(s) = mu0;
sigma(s) = PML_sigma(1);
sigma_star(s) = sigma(s).*mu(s)./eps(s);

%Free Space
s = find((x >= PML_depth) & (x < (max_x-PML_depth)));
eps(s) = eps0;
mu(s) = mu0;
sigma(s) = 0;
sigma_star(s) = 0;

%Second PML
s = find((x >= (max_x-PML_depth)) & (x < max_x));
eps(s) = eps0;
mu(s) = mu0;
sigma(s) = PML_sigma(1);
sigma_star(s) = sigma(s).*mu(s)./eps(s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We calculate the "C" and "D" parameters, as given by Taflove and
% Hagness [2000]. We provide a half-spatial-step version as well
% for simplicity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ca = (1-sigma.*delta_t./2./eps)./(1+sigma.*delta_t./2./eps);
Cb = (delta_t./eps./delta_x)./(1+sigma.*delta_t./2./eps);
Da = (1-sigma_star.*delta_t./2./mu)./(1+sigma_star.*delta_t./2./mu);
Db = (delta_t./mu./delta_x)./(1+sigma_star.*delta_t./2./mu);
Ca_Ex = Ca(1:end-1,:);
Ca_Ey = Ca(:,1:end-1);
Ca_Ez = Ca;
Cb_Ex = Cb(1:end-1,:);
Cb_Ey = Cb(:,1:end-1);
Cb_Ez = Cb;
Da_Hx = Da(:,1:end-1);
Da_Hy = Da(1:end-1,:);
Da_Hz = Da(1:end-1,1:end-1);
Db_Hx = Db(:,1:end-1);
Db_Hy = Db(1:end-1,:);
Db_Hz = Db(1:end-1,1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We then set our initial conditions (E=H=J=D=0 everywhere). Based
% on the Yee grid, some vectors are longer than others. We keep 1
% past value of all vectors, except for J.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hx = zeros(length(x),length(y)-1,2);
Ey = zeros(length(x),length(y)-1,2);
Ez = zeros(length(x),length(y),2);
Jy = zeros(length(x),length(y)-1,1);
Jz = zeros(length(x),length(y),1);
Dy = zeros(length(x),length(y)-1,2);
Dz = zeros(length(x),length(y),2);

Hy = zeros(length(x)-1,length(y),2);
Hz = zeros(length(x)-1,length(y)-1,2);
Ex = zeros(length(x)-1,length(y),2);
Jx = zeros(length(x)-1,length(y),1);
Dx = zeros(length(x)-1,length(y),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we calculate the index of the source location.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, source_location_x] = min(abs(x - source_x));
[~, source_location_y] = min(abs(y - source_y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to calculate the reflection coefficient, transmission
% coefficient, attenuation rate, and phase constant in each medium,
% we save the values at certain locations for all time. Here we
% calculate the indices of those locations. These values change
% based on the parameters of medium 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, x_location1] = min(abs(x - (0.5)));
[~, x_location2] = min(abs(x - (0.59)));
[~, x_location3] = min(abs(x - (0.61)));
[~, x_location4] = min(abs(x - (0.69)));
[~, x_location5] = min(abs(x - (0.5)));

[~, y_location] = min(abs(y - (max_y/2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We identify the time steps at which we would like to make plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_t = 0;
max_t = 2*sqrt((max_x-PML_depth-source_x)^2+(max_y-PML_depth-source_y)^2)*sqrt(mu(1)*eps(1));
%Before reaching PML
t1 = (max_t-min_t)/6+min_t;
t2 = 2*(max_t-min_t)/6+min_t;
t3 = 3*(max_t-min_t)/6+min_t;

%After reaching PML
t4 = 4*(max_t-min_t)/6+min_t;
t5 = 5*(max_t-min_t)/6+min_t;
t6 = 6*(max_t-min_t)/6+min_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And finally we run the FDTD loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = min_t;
iteration = 0;
output_save_Ez = zeros(round(max_t/delta_t)+1,5);
eta1 = sqrt(sqrt(-1)*2*pi*f0*mu(1)/(sigma(1) +sqrt(-1)*2*pi*f0*eps(1)));
beta1 = 2*pi*f0*sqrt(mu(1)*eps(1))*sqrt(1/2*(sqrt(1+(sigma(1)/2/pi/f0/eps(1))^2)+1));

while (t <= max_t)
    iteration = iteration + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The crux of the FDTD routine are these next lines, followed
    % by perfect electric conductor boundary conditions.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hx(:,:,2) = Da_Hx.*Hx(:,:,1) + Db_Hx.*(Ez(:,1:end-1,1)-Ez(:,2:end,1));
    Hy(:,:,2) = Da_Hy.*Hy(:,:,1) + Db_Hy.*(Ez(2:end,:,1)-Ez(1:end-1,:,1));
    Hz(:,:,2) = Da_Hz.*Hz(:,:,1) + Db_Hz.*(Ex(:,2:end,1) - Ex(:,1:end-1,1) + Ey(1:end-1,:,1)-Ey(2:end,:,1));
    
    Dx(:,2:end-1,2) = Dx(:,2:end-1,1) + delta_t.*( (Hz(:,2:end,2)-Hz(:,1:end-1,2))./delta_y );
    Dx(:,[1 end],2) = 0;
    
    Ex(:,2:end-1,2) = Ex(:,2:end-1,1) + (1./eps(2:end,2:end-1)).*( Dx(:,2:end-1,2).*(1+sigma(2:end,2:end-1).*delta_t./2./eps(2:end,2:end-1)) - Dx(:,2:end-1,1).*(1-sigma(2:end,2:end-1).*delta_t./2./eps(2:end,2:end-1)) );
    Ex(:,[1 end],2) = 0;
    
    %     Dz(2:end-1,2:end-1,2) = Dz(2:end-1,2:end-1,1) + delta_t.*( (Hy(2:end,2:end-1,2)-Hy(1:end-1,2:end-1,2))./delta_x - (Hx(2:end-1,2:end,2)-Hx(2:end-1,1:end-1,2))./delta_y );
    %     Dz([1 end],[1 end],2) = 0;
    %
    %     Ez(2:end-1,2:end-1,2) = Ez(2:end-1,2:end-1,1) + (1./eps(2:end-1,2:end-1)).*( Dz(2:end-1,2:end-1,2).*(1+sigma(2:end-1,2:end-1).*delta_t./2./eps(2:end-1,2:end-1)) - Dz(2:end-1,2:end-1,1).*(1-sigma(2:end-1,2:end-1).*delta_t./2./eps(2:end-1,2:end-1)));
    %     Ez([1 end],[1 end],2) = 0;
    
    if (t <= 1e-9)
        Jz(source_location_x, source_location_y, 1) = 2.5.*sin(2*pi*f0.*t);
    else
        Jz(source_location_x, source_location_y, 1) = 0;
    end
    
    
    %     Ex(:,2:end-1,2) = Ca_Ex(:,2:end-1).*Ex(:,2:end-1,1) + Cb_Ex(:,2:end-1).*(Hz(:,2:end,2)-Hz(:,1:end-1,2) - Jx(:,2:end-1,1).*delta_x);
    %     Ex(:,[1 end],2) = 0;
    
    Ey(2:end-1,:,2) = Ca_Ey(2:end-1,:).*Ey(2:end-1,:,1) + Cb_Ey(2:end-1,:).*(Hz(1:end-1,:,2)-Hz(2:end,:,2) - Jy(2:end-1,:,1).*delta_x);
    Ey([1 end],:,2) = 0;
    
    Ez(2:end-1,2:end-1,2) = Ca_Ez(2:end-1,2:end-1).*Ez(2:end-1,2:end-1,1) + Cb_Ez(2:end-1,2:end-1).*(Hy(2:end,2:end-1,2)-Hy(1:end-1,2:end-1,2) + Hx(2:end-1,1:end-1,2)-Hx(2:end-1,2:end,2) - Jz(2:end-1,2:end-1,1).*delta_x);
    Ez([1 end],[1 end],2) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure-making code below. For E field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (iteration == round(t1/delta_t)),
        figure(2);
    elseif (iteration == round(t2/delta_t)),
        figure(3);
    elseif (iteration == round(t3/delta_t)),
        figure(4);
    elseif (iteration == round(t4/delta_t)),
        figure(5);
    elseif (iteration == round(t5/delta_t)),
        figure(6);
    elseif (iteration == round(t6/delta_t)),
        figure(7);
    else
        figure(1);
    end
    subplot(311);
    imagesc(x, y, Ex(:,:,2)');
    axis xy;
    axis equal;
    axis tight;
    caxis([-10 10]*1e-3);
    title('E_x');
    title(sprintf(['EEL6487 HW 4, Problem 1: Time: ' num2str(round(t*1e12)) ' picoseconds\nE_x']));
    
    subplot(312);
    imagesc(x, y, Ey(:,:,2)');
    axis xy;
    axis equal;
    axis tight;
    caxis([-10 10]*1e-3);
    title('E_y');
    ylabel('y (meters)');
    
    subplot(313);
    imagesc(x, y, Ez(:,:,2)');
    axis xy;
    axis equal;
    axis tight;
    caxis([-10 10]*1e-3);
    title('E_z');
    xlabel('x (meters)');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure-making code below. For H field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (iteration == round(t1/delta_t)),
        figure(9);
    elseif (iteration == round(t2/delta_t)),
        figure(10);
    elseif (iteration == round(t3/delta_t)),
        figure(11);
    elseif (iteration == round(t4/delta_t)),
        figure(12);
    elseif (iteration == round(t5/delta_t)),
        figure(13);
    elseif (iteration == round(t6/delta_t)),
        figure(14);
    else
        figure(8);
    end
    
    subplot(311);
    imagesc(x, y, real(eta1*Hx(:,:,2)'));
    axis xy;
    axis equal;
    axis tight;
    caxis([-20 20]*1e-4);
    title('H_x');
    title(sprintf(['EEL6487 HW 4, Problem 1: Time: ' num2str(round(t*1e12)) ' picoseconds\nH_x']));
    
    subplot(312);
    imagesc(x, y, real(eta1*Hy(:,:,2)'));
    axis xy;
    axis equal;
    axis tight;
    caxis([-20 20]*1e-4);
    title('H_y');
    ylabel('y (meters)');
    
    subplot(313);
    imagesc(x, y, real(eta1*Hz(:,:,2)'));
    axis xy;
    axis equal;
    axis tight;
    caxis([-20 20]*1e-4);
    title('H_z');
    xlabel('x (meters)');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the arrays/variables for the next iteration.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ex(:,:,1) = Ex(:,:,2);
    Ex(:,:,2) = 0;
    Ey(:,:,1) = Ey(:,:,2);
    Ey(:,:,2) = 0;
    Ez(:,:,1) = Ez(:,:,2);
    Ez(:,:,2) = 0;
    Hx(:,:,1) = Hx(:,:,2);
    Hx(:,:,2) = 0;
    Hy(:,:,1) = Hy(:,:,2);
    Hy(:,:,2) = 0;
    Hz(:,:,1) = Hz(:,:,2);
    Hz(:,:,2) = 0;
    t = t + delta_t;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the electric field values at the locations calculated
    % above.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_save_Ez(iteration,1) = Ez(x_location1,y_location,1);
    output_save_Ez(iteration,2) = Ez(x_location2,y_location,1);
    output_save_Ez(iteration,3) = Ez(x_location3,y_location,1);
    output_save_Ez(iteration,4) = Ez(x_location4,y_location,1);
    output_save_Ez(iteration,5) = Ez(x_location5,y_location,1);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate ideal reflection and transmission coefficients as well
% as ideal attenuation rates and phase constants in both media.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps1 = eps(1);
mu1 = mu(1);
sigma1 = sigma(2*N_PML);
sigmastar1 = sigma_star(2*N_PML);

eps2 = eps(end);
mu2 = mu(end);
sigma2 = sigma(end);
sigmastar2 = sigma_star(end);

eta1 = sqrt((sigmastar1 + sqrt(-1)*2*pi*f0*mu1)/(sigma1 +sqrt(-1)*2*pi*f0*eps1));
eta2 = sqrt((sigmastar2 + sqrt(-1)*2*pi*f0*mu2)/(sigma2 +sqrt(-1)*2*pi*f0*eps2));

ideal_alpha1 = 2*pi*f0*sqrt(mu1*eps1)*sqrt(1/2*(sqrt(1+(sigma1/2/pi/f0/eps1)^2)-1));
ideal_alpha2 = sigma2*eta1;

ideal_beta1 = 2*pi*f0*sqrt(mu1*eps1)*sqrt(1/2*(sqrt(1+(sigma1/2/pi/f0/eps1)^2)+1));
ideal_beta2 = ideal_beta1;

ideal_gamma = (eta2-eta1)/(eta2+eta1);
mag_ideal_gamma = abs(ideal_gamma)
angle_ideal_gamma = angle(ideal_gamma)*180/pi

ideal_T = 2*eta2/(eta2+eta1);
mag_ideal_T = abs(ideal_T);
angle_ideal_T = angle(ideal_T)*180/pi;

vp1 = 2*pi*f0/ideal_beta1;
vp2 = 2*pi*f0/ideal_beta2;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the chirp z-transform to calculate reflection and transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1_index = round((0.5-0.35)/vp1/delta_t)+1;
t2_index = t1_index + round(1e-9/delta_t);
t3_index = round(0.55/vp1/delta_t)+1;
t4_index = t3_index + round(1e-9/delta_t);

inc_field = output_save_Ez(t1_index:t2_index,5);
ref_field = output_save_Ez(t3_index:t4_index,5);

IFZ = czt(inc_field,1,1,exp(-sqrt(-1)*2*pi*f0*delta_t));
RFZ = czt(ref_field,1,1,exp(-sqrt(-1)*2*pi*f0*delta_t));


gamma = RFZ./IFZ;
T = gamma + 1;

mag_gamma = abs(gamma)
angle_gamma = angle(gamma)*180/pi

mag_T = abs(T)
angle_T = angle(T)*180/pi


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the chirp z-transform to calculate alpha and beta in medium 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1_index = round((0.61-0.35)/vp1/delta_t)+1;
t2_index = round((0.69-0.35)/vp1/delta_t);

t_field1 = output_save_Ez(t1_index:t2_index,1);
t_field2 = output_save_Ez(t1_index:t2_index,2);

TFZ1 = czt(t_field1,1,1,exp(-sqrt(-1)*2*pi*f0*delta_t));
TFZ2 = czt(t_field2,1,1,exp(-sqrt(-1)*2*pi*f0*delta_t));

alpha1=-log(abs(TFZ2./TFZ1))/0.08
beta1=mod(angle(TFZ2./TFZ1),2*pi)/0.08
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the chirp z-transform to calculate alpha and beta in medium 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1_index = round((0.61-0.35)/vp1/delta_t)+1;
t2_index = round((0.69-0.35)/vp1/delta_t);

t_field1 = output_save_Ez(t1_index:t2_index,3);
t_field2 = output_save_Ez(t1_index:t2_index,4);

TFZ1 = czt(t_field1,1,1,exp(-sqrt(-1)*2*pi*f0*delta_t));
TFZ2 = czt(t_field2,1,1,exp(-sqrt(-1)*2*pi*f0*delta_t));

alpha2=-log(abs(TFZ2./TFZ1))/0.08
beta2=mod(angle(TFZ2./TFZ1),2*pi)/0.08


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results for f=1 GHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mag_ideal_gamma = 0
% angle_ideal_gamma = 0
% mag_gamma = 0.2185
% angle_gamma = -8.3712
% mag_T = 3.1779
% angle_T = -5.7455
% alpha1 = 12.5167
% beta1 = 34.4385
% alpha2 = 64.9514
% beta2 = 41.8574

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results for f=2 GHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mag_ideal_gamma = 0
% angle_ideal_gamma = 0
% mag_gamma = 0.5770
% angle_gamma = -137.1767
% mag_T = 0.6975
% angle_T = -34.2171
% alpha1 = 12.9250
% beta1 = 39.8604
% alpha2 =  67.2456
% beta2 = 52.2199