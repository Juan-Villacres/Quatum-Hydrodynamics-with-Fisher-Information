clear
clc
%% definition of constants and rendering of initial value
%CONSTANTS & SYMPBOLS*****************************************************
fovx = 8; %field of view on the real line [-fovx, fovx]. 
numParticles = 51; %number of Bohmian particles. 
dt = 10^-2; %time step for numerical integration.
numSteps = 3000; %number of iterations of the numerical %integration. 
stencil = 3; %number of points left and right w.r.t. some %point for fitting 0.5*log(rho). 
degree = 6; %degree of polynomial for fitting 
            %please note: stencil=degree/2 means a 
            %precise polynomial fit of ’degree’ degree 
            %while stencil>degree/2 means a least square 
            %fit with a polynomial of degree ’degree’. 
s_stencil = 3; %log(sqrt(rho)) number of points left and 
               %right w.r.t. some point for fitting the 
               %phase. 
s_degree = 6; %degree of polynomial for fitting the phase.
              %the above note hols for the phase.
sym_x = sym('sym_x','real'); %dummy symbol for symbolic math.
tt= sym('tt','real');
%FUNCTIONS--------------------------------------------
% Analytic solution of an initial gaussian evolved according to the free 
% Schroedinger equation. 
% INPUT : t=time of evolution, x=position, sigma=initial width of the 
% gaussian 
% OUTPUT: complex value of the wavefunction 
gauss = @(t, x, sigma) (sigma/(pi*(sigma+1i*t)^2)).^(1/4).*exp(-x.^2/(2*(sigma^2+t^2))*(sigma-1i*t));
%----------------
% psi=analytic solution of a superposition of initial gaussians evolved 
% according to the free Schroedinger equation. 
% rho,S,vel=|psi|^2, the phase, the velocity field of psi 
% INPUT : t=time of evolution, x=position 
% OUTPUT: real resp. complex function value
psi=@(t, x) 1/sqrt(2)*(gauss(t,x-3,4)+gauss(t,x+3,4));
rho=@(t, x) psi(t,x).*conj(psi(t,x));
S = @(t, x) imag(log(psi(t,x)./abs(psi(t,x))));
vel = @(t, x) imag(subs(diff(psi(t,sym_x),sym_x)./psi(t,sym_x),sym_x,x));
C=@(t,x) 0.5.*log(rho(t,x));
dC=@(t,x) subs(diff(C(tt,x),tt),tt,x);
ddC=@(t,x) subs(diff(diff(C(tt,x),tt),tt),tt,t); 
%INITIAL DATA---
clear x; %initial distribution of the %Bohmian particles
x = (-fovx:(2*fovx/numParticles):fovx); %generate a uniform distribution
%plot the initial data: blue=initial rho and green=initial veleocity field 
plot(x, vel(0,x),'g',x,rho(0,x),'b.-'); 
axis([-fovx fovx -1 1]); 
grid on;
%% bohmian propagation
warning off;
%initialize the array that keeps track of the (x,0.5*log(rho),S,time) 
%data history. 
history ={x vel(0,x) 0.5*log(rho(0,x)) S(0,x) 0};
%loop over the number of time steps of numerical integration 
for step = 1:numSteps
    xlen = length(history{step,1});
    x_list = history{step,1}; 
    vel_list = history{step,2}; 
    logsqrtRho_list = history{step,3}; 
    s_list = history{step,4};
    %first loop over the all Bohmian particles updating the 
    %x_list and s_list 
    for c = 1:xlen
        %find stencil points for logsqrtRho 
        left_point=c - stencil; 
        right_point=c + stencil;
        deg = degree; 
        if (left_point < 1)
            left_point=1; 
            right_point= stencil*2+1 + round(numParticles/7); 
            %decrease degree at the boundary as dodgy fix to avoid 
            %bouadary problems 
            deg = 2;
        elseif (right_point > xlen)
            left_point=xlen -(2*stencil+1)-round(numParticles/7); 
            right_point = xlen; 
            %decrease degree at the boundary as dodgy fix to avoid 
            %bouadary problems 
            deg = 2;
        end
        %fit logsqrtRho_list
        fitdata_x=history{step,1}(left_point:right_point); 
        fitdata_y=history{step,3}(left_point:right_point);
        logsqrtRho_fit = polyfit(fitdata_x, fitdata_y, deg); 
        %compute derivatives 
        d_logsqrtRho_fit = polyder(logsqrtRho_fit); 
        dd_logsqrtRho_fit = polyder(d_logsqrtRho_fit); 
        %compute quantum potential
        Q=-1/2*(ddC(history{step,5},x_list(c)) +(dC(history{step,5},x_list(c))).^2); %Potencial cuántico que quiero investigar
        %Qc=-1/2*(polyval(dd_logsqrtRho_fit,x_list(c))+ polyval(d_logsqrtRho_fit,x_list(c)).^2);
        %La expresión de arriba que es el potencial cuántico usual para estas ecuaciones pienso habría que adicionarle 
        %al que yo plantéo y así probar mi hipótesis a velocidades cercanas a la de la luz.
        %update phase: S(t+dt)=S(t)+dt(1/2v(t)^2-Q(t)); 
        s_list(c) = s_list(c) + (1/2*vel_list(c).^2 - Q)*dt ;
        %update position: r(t+dt)=r(t)+v(t)*dt 
        x_list(c) = x_list(c) + vel_list(c)*dt;
    end
     %second loop over the all Bohmian particles updating 
     %the vel_list and logsqrtRho_list
     for c = 1:xlen 
         %find stecil points for S 
         s_left_point = c-s_stencil; 
         s_right_point= c+s_stencil; 
         deg=s_degree; 
         if (s_left_point < 1)
             s_left_point =1; 
             s_right_point=s_stencil*2+1 + round(numParticles/7); 
             %decrease degree at the boundary as dodgy fix to avoid 
             %bouadary problems 
             deg = 2;
         elseif (s_right_point > xlen)
             s_left_point=xlen -(2*s_stencil+1)-round(numParticles/7); 
             s_right_point=xlen; 
             %decrease degree at the boundary as dodgy fix to avoid 
             %bouadary problems 
             deg = 2;
         end
         %fit s_list 
         fitdata_x=x_list(s_left_point:s_right_point); 
         fitdata_y=s_list(s_left_point:s_right_point); 
         s_fit = polyfit(fitdata_x, fitdata_y, deg); 
         %compute derivatives 
         d_s_fit=polyder(s_fit); 
         dd_s_fit=polyder(d_s_fit); 
         %v(t+dt) = S’(t+dt) 
         vel_list(c)=polyval(d_s_fit,x_list(c)); 
         %R(t+dt) = R(t)*exp(-1/2 S’’(t+dt)*dt)Esta es la ecuación 6 de Lopreore  
         logsqrtRho_list(c)=logsqrtRho_list(c)-1/2*polyval(dd_s_fit,x_list(c))*dt;
     end
     %compile the tuple (x,vel,0.5*log(rho),time) 
     newSnapshot={x_list vel_list logsqrtRho_list s_list history{step,5}+dt}; 
     %and save the snapshot 
     history=[history; newSnapshot]; 
     %history=newSnapshot;
     %compute minimal distance of the Bohmian particles 
     dxmin = min(x_list(2:length(x_list))-x_list(1:length(x_list)-1));
     if (dxmin <= 0)
         sprintf('TRAJECTORY CROSSING OCCURED AT TIME: %f',history{step,5}+dt)
         break;
     end
     if (step == 1||mod(step/numSteps*100,10)== 0)
      %   %plot integrated velocity field 
       %  plot(x_list(1):0.1:x_list(xlen),vel(history{step,5}+dt,x_list(1):0.1:x_list(xlen)),'r',x_list,vel_list,'b');
       %  axis([-25 25 -5 5]); 
       %  grid on;
         sprintf('Progress: %.1f%%, timestep: %f fs, min deltax: %f',step/numSteps*100,history{step,5}+dt,dxmin) 
       %  pause(0.1);
     end 
end
