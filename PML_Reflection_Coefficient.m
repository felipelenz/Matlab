% Author: Felipe Lenz Carvalho
% This code calculates analytically the magnitude of the reflection
% coefficient of an incident wave as a function of the PML (Perfect
% Matched Layer) conductivity

% Constants
mu0=4*pi*1e-7;
epsilon0=8.854e-12;
c=1/sqrt(mu0*epsilon0);
f=2e9;
w=2*pi*f;
deltax=0.001;

%Medium 1
mu1=mu0;
epsilon1=epsilon0;
eta1=sqrt(mu1/epsilon1);

%Medium 2
mu2=mu1;
epsilon2=epsilon1;
sigma=[0.001 0.01 0.1 1];
sigma_star=sigma.*mu1./epsilon1;

PML_depth=[10 20 50 100]*1e-2;
PML_samples=PML_depth./deltax;

for j=1:4,
    gamma_eff_out=zeros(1,4);
    for i=1:length(sigma),
        eta2=sqrt( (mu2*(1+sigma_star(i)/(sqrt(-1)*w*mu2)))/(epsilon2*(1+sigma(i)/(sqrt(-1)*w*epsilon2))) );
        %alpha2 = w*sqrt(mu2*epsilon2)*sqrt(1/2*(sqrt(1+(sigma(i)/w/epsilon2)^2)- 1));
        alpha2 = sigma(i)*eta1;
        
        gamma=(eta2-eta1)/(eta2+eta1);
        tau=2*eta2/(eta2+eta1);
        tau_2=2*eta1/(eta2+eta1);
        gamma_eff=gamma-tau*exp(-2*abs(alpha2)*PML_depth(j))*tau_2;
        gamma_eff_out(i)=gamma_eff;
    end
    gamma_eff_plot1=gamma_eff_out;
    loglog(sigma,abs(gamma_eff_plot1),'-v'); hold all;
    legend('PML depth of 10 cm','PML depth of 20 cm','PML depth of 50 cm','PML depth of 100 cm','location','southwest');
    xlabel('PML Conductivity (S/m)'); ylabel('|\Gamma_{eff}|'); grid on;
    title('HW 3a - Problem 1 - Felipe Lenz');
end
%legend('PML depth of 10cm','PML depth of 20cm','PML depth of 50cm','PML depth of 100cm');