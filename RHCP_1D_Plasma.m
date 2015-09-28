%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEL 6487
% Spring 2015
% Homework #4, Problem 4
% Original Code: Dr. Moore's solution for HW 2b, Problem 1
% Updated Code: Felipe Lenz Carvalho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c = 1/sqrt(mu0*eps0);
v=1e8;

qe=-1.602e-19;
me=9.11e-31;
Bo=3.12e-2;
Ne=0.5e16;

% wc=-1*qe*Bo/me;
% wp=sqrt( Ne*qe*qe/me/eps0 );
wc=2e9;
wp=6e9;

%B parallel to propagation
w_min=0;
delta_w=wp/100;
w_max =2*pi*4e9;

w=w_min:delta_w:w_max;

X=(wp^2)./(w.^2);
Y=wc./w;
Z=v./w;
U=1-sqrt(-1).*Z;

n_2_RH=1-X./(U-Y);
n_2_LH=1-X./(U+Y);
beta_2_RH=(w.^2).*(mu0*eps0).*n_2_RH;
beta_RH=sqrt(beta_2_RH);
x_RH=real(beta_RH);

beta_2_LH=(w.^2).*(mu0*eps0).*n_2_LH;
beta_LH=sqrt(beta_2_LH);
x_LH=real(beta_LH);
figure(4);
plot(x_RH,w,'b'); hold on;
plot(x_LH,w,'r');
xlabel('real(\beta)'); ylabel('\omega rad/s');
title('\omega k diagram (\nu non zero)');
ylim([0 5*2.9e9]);

figure(5);
plot(x_RH,w,'b'); hold on;
plot(x_LH,w,'r');
xlabel('real(\beta)'); ylabel('\omega rad/s');
title('Analytical \omega k diagram (\nu non zero)');
ylim([0 5*2.9e9]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We then provide our S parameter, calculate our range for x, and use them
% to calculate delta_t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = 1;
f_min=wp/5/2/pi;
f0 = [1 2 3 4 5 6 7 8 9 10].*f_min;%Hz

for a=1:length(f0),
    lambda0 = c/f0(a);
    
    min_z = 0;%meters
    max_z = 4;%meters
    
    delta_z = 1e-3;%meters
    z = (min_z:delta_z:max_z)';
    
    N_lambda = lambda0/delta_z;
    delta_t = delta_z*S/c;
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We provide the permittivity, permeability, and conductivity as a function
    % of space.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eps = zeros(size(z));
    mu = zeros(size(z));
    sigma = zeros(size(z));
    
    s = find((z >= 0) & (z < 4));
    eps(s) = eps0;
    mu(s) = mu0;
    sigma(s) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We calculate the "C" and "D" parameters, as given by Taflove and Hagness
    % [2000]. We provide a half-spatial-step version as well for simlicity.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ca = (1-sigma.*delta_t./2./eps)./(1+sigma.*delta_t./2./eps);
    Cb = (delta_t./eps./delta_z)./(1+sigma.*delta_t./2./eps);
    Da = ones(size(Ca));
    Db = delta_t./mu./delta_z;
    
    eps_half = eps(1:end-1);
    mu_half = mu(1:end-1);
    sigma_half = sigma(1:end-1);
    
    Ca_half = (1-sigma_half.*delta_t./2./eps_half)./(1+sigma_half.*delta_t./2./eps_half);
    Cb_half = (delta_t./eps_half./delta_z)./(1+sigma_half.*delta_t./2./eps_half);
    Da_half = ones(size(Ca_half));
    Db_half = delta_t./mu_half./delta_z;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We then set our initial conditions (E=H=J=0 everywhere). Based on the
    % Yee grid, some vectors are longer than others. We keep 1 past value of
    % all vectors, except for J.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Hx = zeros(length(z)-1,7);
    Ey = zeros(length(z),7);
    Dy = zeros(length(z),7);
    Jy = zeros(length(z),1);
    
    Hy = zeros(length(z)-1,7);
    Ex = zeros(length(z),7);
    Dx = zeros(length(z),7);
    Jx = zeros(length(z),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Here we calculate the index of the source location.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [val, source_location] = min(abs(z - 2.5));
    [val, boundary_location] = min(abs(z - 3));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In order to calculate the reflection coefficient, transmission
    % coefficient, attenuation rate, and phase constant in each medium, we save
    % the values at certain locations for all time. Here we calculate the
    % indices of those locations. These values change based on the parameters
    % of medium 2.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [val, save_location1] = min(abs(z - 2.60));
    [val, save_location2] = min(abs(z - 2.65));
    [val, save_location3] = min(abs(z - 3.01));
    [val, save_location4] = min(abs(z - 3.02));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We identify the time steps at which we would like to make plots.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min_t = 0;
    max_t = 2.0*sqrt(mu(1)*eps(1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the coefficients A,B,C,D1,D2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    v=1e8;
    wc=2e9;
    wp=6e9;
    
    A=[1, (4*v), (wp^2+2*wc^2+6*v^2), (3*(wp^2)*v+4*(wc^2)*v+4*v^3),...
        ((wp^2)*(wc^2)+3*(wp^2)*(v^2)+wc^4+2*(wc^2)*v^2+v^4),(wp^2)*(wc^2)*v+(wp^2)*v^3,0];
    B=[0,0,0,(wp*wc), 2*wp*wc*v, (wp*wc^3 + wp*wc*v^2),0];
    D1=[1, (4*v),(2*wp^2)+2*wc^2+6*v^2, 6*(wp^2)*v+4*(wc^2)*v+4*v^3,...
        2*(wp^2)*(wc^2)+6*(wp^2)*(v^2)+wp^4+wc^4+2*(wc^2)*(v^2)+v^4,...
        2*(wp^2)*(wc^2)*v+2*(wp^2)*(v^3)+2*(wp^4)*v,(wp^4)*v^2+(wp^4)*(wc^2)];
    C=[0,0,0,0,1,v,0];
    D2=[0,0,0,0,1,v,wp^2];
    
    del_6=(1/(delta_t^6)).*[1 -6 15 -20 15 -6 1];
    del_5=(1/2/(delta_t^5)).*[1 -4 5 0 -5 4 -1];
    del_4=(1/4/(delta_t^4)).*[1 -2 -1 4 -1 -2 1];
    del_3=(1/8/(delta_t^3)).*[1 0 -3 0 3 0 -1];
    del_2=(1/16/(delta_t^2)).*[1 2 -1 -4 -1 2 1];
    del_1=(1/32/(delta_t)).*[1 4 5 0 -5 -4 -1];
    del_0=(1/64).*[1 6 15 20 15 6 1];
    
    k_tensor=[del_6;del_5;del_4;del_3;del_2;del_1;del_0];
    A_1=A*k_tensor;
    B_1=B*k_tensor;
    D1_1=D1*k_tensor;
    C_1=C*k_tensor;
    D2_1=D2*k_tensor;
    %%
    t = min_t;
    iteration = 0;
    rangex = [-1 1];
    rangey = [-1.5 1.5];
    
    while (t <= max_t)
        iteration = iteration + 1;
        
        %In the free space region we implement a normal 2-D FDTD for a wave
        %propagating in the z direction
        Hx(1:boundary_location-1,7)=Da_half(1:boundary_location-1).*Hx(1:boundary_location-1,6)+...
            Db_half(1:boundary_location-1).*(Ey(2:boundary_location,6)-Ey(1:boundary_location-1,6));
        Hy(1:boundary_location-1,7)=Da_half(1:boundary_location-1).*Hy(1:boundary_location-1,6)+...
            Db_half(1:boundary_location-1).*(Ex(1:boundary_location-1,6)-Ex(2:boundary_location,6));
        
%         Jx(source_location) = 2.5*sin(2*pi*f0(a).*t);
%         Jy(source_location) = 2.5*sin(2*pi*f0(a).*t + pi/2); %LHCP 
        Jx(source_location) = 2.5*sin(2*pi*f0(a).*t);
        Jy(source_location) = 2.5*sin(2*pi*f0(a).*t - pi/2); %RHCP


        Ex(2:boundary_location-1,7) = Ca(2:boundary_location-1).*Ex(2:boundary_location-1,6)+...
            Cb(2:boundary_location-1).*( Hy(1:boundary_location-2,7)-Hy(2:boundary_location-1,7) - Jx(2:boundary_location-1).*delta_z );
        Ex(1,7)=0;
        
        Ey(2:boundary_location-1,7) = Ca(2:boundary_location-1).*Ey(2:boundary_location-1,6)+...
            Cb(2:boundary_location-1).*( Hx(2:boundary_location-1,7)-Hx(1:boundary_location-2,7) - Jy(2:boundary_location-1).*delta_z );
        Ey(1,7)=0;
        
        %Plasma Region
        Hx(boundary_location:end,7)=Hx(boundary_location:end,6)+(delta_t/mu0/delta_z).*(Ey(boundary_location+1:end,6)-Ey(boundary_location:end-1,6));
        Hy(boundary_location:end,7)=Hy(boundary_location:end,6)-(delta_t/mu0/delta_z).*(Ex(boundary_location+1:end,6)-Ex(boundary_location:end-1,6));
        
        Dx(boundary_location:end-1,7)=Dx(boundary_location:end-1,6)-(delta_t/delta_z).*(Hy(boundary_location:end,7)-Hy(boundary_location-1:end-1,7));
        Dx(end,7)=0;
        
        Dy(boundary_location:end-1,7)=Dy(boundary_location:end-1,6)+(delta_t/delta_z).*(Hx(boundary_location:end,7)-Hx(boundary_location-1:end-1,7));
        Dy(end,7)=0;
        
        Ex(boundary_location:end-1,7)=(1/D1_1(1,1)/eps0).*(A_1(1,1)*Dx(boundary_location:end-1,7) + A_1(1,2)*Dx(boundary_location:end-1,6) + A_1(1,3)*Dx(boundary_location:end-1,5)...
            + A_1(1,4)*Dx(boundary_location:end-1,4) + A_1(1,5)*Dx(boundary_location:end-1,3) + A_1(1,6)*Dx(boundary_location:end-1,2)+ A_1(1,7)*Dx(boundary_location:end-1,1)...
            - B_1(1,1)*Dy(boundary_location:end-1,7) - B_1(1,2)*Dy(boundary_location:end-1,6) - B_1(1,3)*Dy(boundary_location:end-1,5) - B_1(1,4)*Dy(boundary_location:end-1,4)...
            - B_1(1,5)*Dy(boundary_location:end-1,3) - B_1(1,6)*Dy(boundary_location:end-1,2) - B_1(1,7)*Dy(boundary_location:end-1,1))...
            - (1/D1_1(1,1)).*( D1_1(1,2)*Ex(boundary_location:end-1,6) + D1_1(1,3)*Ex(boundary_location:end-1,5) + D1_1(1,4)*Ex(boundary_location:end-1,4)...
            + D1_1(1,5)*Ex(boundary_location:end-1,3) + D1_1(1,6)*Ex(boundary_location:end-1,2) + D1_1(1,7)*Ex(boundary_location:end-1,1) );
        Ex(end,7)=0;
        
        Ey(boundary_location:end-1,7)=(1/D1_1(1,1)/eps0).*(B_1(1,1)*Dx(boundary_location:end-1,7) + B_1(1,2)*Dx(boundary_location:end-1,6) + B_1(1,3)*Dx(boundary_location:end-1,5)...
            +B_1(1,4)*Dx(boundary_location:end-1,4) + B_1(1,5)*Dx(boundary_location:end-1,3) + B_1(1,6)*Dx(boundary_location:end-1,2)+ B_1(1,7)*Dx(boundary_location:end-1,1)...
            + A_1(1,1)*Dy(boundary_location:end-1,7) + A_1(1,2)*Dy(boundary_location:end-1,6) + A_1(1,3)*Dy(boundary_location:end-1,5) + A_1(1,4)*Dy(boundary_location:end-1,4)...
            + A_1(1,5)*Dy(boundary_location:end-1,3) + A_1(1,6)*Dy(boundary_location:end-1,2) + A_1(1,7)*Dy(boundary_location:end-1,1))...
            - (1/D1_1(1,1)).*( D1_1(1,2)*Ey(boundary_location:end-1,6) + D1_1(1,3)*Ey(boundary_location:end-1,5) + D1_1(1,4)*Ey(boundary_location:end-1,4)...
            + D1_1(1,5)*Ey(boundary_location:end-1,3) + D1_1(1,6)*Ey(boundary_location:end-1,2) + D1_1(1,7)*Ey(boundary_location:end-1,1) );
        Ey(end,7)=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure-making code below.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1);
        subplot(211); plot(z, Ex(:,7));
        ylim(rangex);
        xlim([min_z max_z]);
        grid on;
        title(sprintf(['EEL6487 HW 4, Problem 4: Time: ' num2str(round(t*1e12))...
            ' picoseconds \nwave frequency =', num2str(f0(a)/1e9),' GHz, wc/2/pi=',...
            num2str(wc/2/pi/1e9),' GHz, wp/2/pi=',num2str(wp/2/pi/1e9), 'GHz']));
        h1 = line([3 3],rangex);
        set(h1,'color','r','LineWidth',2);
        xlabel('z (meters)'); ylabel('E_x (V/m)');
        
        subplot(212); plot(z, Ey(:,7));
        ylim(rangey);
        xlim([min_z max_z]);
        grid on;
        h1 = line([3 3],rangey);
        set(h1,'color','r','LineWidth',2);
        xlabel('z (meters)'); ylabel('E_y (V/m)');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save the electric field values at the locations calculated above.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        output_save_Ey(iteration,1) = Ey(save_location1,7);
        output_save_Ey(iteration,2) = Ey(save_location2,7);
        output_save_Ey(iteration,3) = Ey(save_location3,7);
        output_save_Ey(iteration,4) = Ey(save_location4,7);
        
        Hx(:,1:6)=Hx(:,2:7);
        Hx(:,7)=0;
        Ey(:,1:6) = Ey(:,2:7);
        Ey(:,7);
        Dy(:,1:6) = Dy(:,2:7);
        Dy(:,7);
        Hy(:,1:6) = Hy(:,2:7);
        Hy(:,7);
        Ex(:,1:6) = Ex(:,2:7);
        Ex(:,7);
        Dx(:,1:6) = Dx(:,2:7);
        Dx(:,7);
        t=t+delta_t;
               
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use the chirp z-transform to calculate the reflection and transmission
    % coefficients.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    t1_index = round(0.1/c/delta_t)+1;
    t2_index = t1_index + round(0.25/f0(a)/delta_t);
    t3_index = round(0.9/c/delta_t)+1;
    t4_index = t3_index + round(0.25/f0(a)/delta_t);
    
    incident_field = output_save_Ey(t1_index:t2_index,1);
    reflected_field = output_save_Ey(t3_index:t4_index,1);
    
    IFZ = czt(incident_field,1,1,exp(-sqrt(-1)*2*pi*f0(a)*delta_t));
    RFZ = czt(reflected_field,1,1,exp(-sqrt(-1)*2*pi*f0(a)*delta_t));
    
    gamma = RFZ./IFZ;
    T = gamma + 1;
    mag_gamma = abs(gamma)
    angle_gamma = angle(gamma)*180/pi
    mag_T = abs(T)
    angle_T = angle(T)*180/pi
    figure(2);
    subplot(211); plot(2*pi*f0(a),mag_gamma,'*'); hold all
    ylabel('|\Gamma|');
    title('\omega_{wave} vs. \Gamma diagram (Source is RHCP)');
    subplot(212); plot(2*pi*f0(a),angle_gamma,'*'); hold all
    xlabel('\omega_{wave} rad/s'); ylabel('\angle \Gamma degrees');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use the chirp z-transform to calculate the attenuation rate and phase
    % constant in medium 2.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t1_index = round(0.51/c/delta_t)+1;
    t2_index = t1_index + round(0.25/f0(a)/delta_t);
    t3_index = round(0.52/c/delta_t)+1;
    t4_index = t3_index + round(0.25/f0(a)/delta_t);
    
    transmitted_field1 = output_save_Ey(t1_index:t2_index,3);
    transmitted_field2 = output_save_Ey(t3_index:t4_index,4);
    TFZ1 = czt(transmitted_field1,1,1,exp(-sqrt(-1)*2*pi*f0(a)*delta_t));
    TFZ2 = czt(transmitted_field2,1,1,exp(-sqrt(-1)*2*pi*f0(a)*delta_t));
    
    alpha2 = -log(abs(TFZ2./TFZ1))/0.01
    transmitted_field1 = output_save_Ey(t3_index:t2_index,3);
    transmitted_field2 = output_save_Ey(t3_index:t2_index,4);
    TFZ1 = czt(transmitted_field1,1,1,exp(-sqrt(-1)*2*pi*f0(a)*delta_t));
    TFZ2 = czt(transmitted_field2,1,1,exp(-sqrt(-1)*2*pi*f0(a)*delta_t));
    beta2 = angle(TFZ2./TFZ1)/0.01
    
    figure(4);
    plot(real(beta2),2*pi*f0(a),'*'); hold all
    xlabel('real(\beta)'); ylabel('\omega_{wave} rad/s');
    title('\omega_{wave} k diagram (Source is RHCP)');
    
end