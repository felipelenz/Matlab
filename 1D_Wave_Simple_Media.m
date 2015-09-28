% Author: Felipe Lenz Carvalho
% This code models a 1-D triangular wave pulse propagating in simple 
% medium bounded by Perfect Electric Conductors (PEC) in both ends of 
% the slab

%%
miu=4*pi*1e-7;
epsilon0 = 8.854e-12; %F/m
c=1/sqrt(epsilon0*miu);
S=4/9;

x_0=1; %meters
deltax=0.0001; 
x=[0:deltax:x_0-deltax]; %x meter
N=1/deltax;
deltat=S*deltax/c;


Ez = zeros(N,3);
Ez(:,1)=triangularPulse(0.4,0.5,0.6,x);
Ez(:,2)=triangularPulse(0.4,0.5,0.6,x);
 
iteration = 0;

while Ez>-0.98,
     iteration = iteration + 1;  
        Ez(2:end-1,3)=((c*deltat)^2)*( (Ez(3:end,2)-2*Ez(2:end-1,2)+Ez(1:end-2,2))/( (deltax)^2 ) )+2*Ez(2:end-1,2)-Ez(2:end-1,1);
        Ez(1,3)=0;
        Ez(end,3)=0;
        
        figure(1); subplot(122); plot(x,Ez(:,3)); ylim([-1 1]); 
        title_string=['Elapsed Time=',num2str((iteration-1)*deltat),' seconds'];
        title(title_string); grid on; xlabel('Distance (meters)'); ylabel('E_z');
        
        if (iteration-1) == 0.20/(c*deltat),
            subplot(221); plot(x,Ez(:,3)); ylim([-1 1]);
            title_string=['Elapsed Time=',num2str((iteration-1)*deltat),' seconds'];
            title(title_string); grid on; xlabel('Distance (meters)'); ylabel('E_z');
        end
        if (iteration-1) == 0.80/(c*deltat),
            subplot(223); plot(x,Ez(:,3)); ylim([-1 1]);
            title_string=['Elapsed Time=',num2str((iteration-1)*deltat),' seconds'];
            title(title_string); grid on; xlabel('Distance (meters)'); ylabel('E_z');
        end

        Ez(:,1:2)=Ez(:,2:3);
        Ez(:,3)=0;
 end