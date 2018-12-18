% Copyright (c) 2018, Wenxiao Wang All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
% IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%%
clear all
close all
clc
%% Define PM size, constraint image size
M1 = 40;    % PM size
M2 = M1*2;  % PSF image size

xx = 1:M1;
[xx,yy] = meshgrid(xx,xx);
cir = sqrt((xx-M1/2-0.5).^2 + (yy-M1/2-0.5).^2);
cir(cir<=M1/2) = 1;
cir(cir>M1/2) = 0;  %pupil function

%% Desired PSFs
load('desired_PSFs.mat');
u1 = rand(M1,M1);   
ran = 2;
load('non_linear_phase_delay.mat');
%% Gauss-Newton optimization
%load('H_high_resolution_40by40');
fprintf('Building Fourier transform matrix...\n');
H = zeros(M1*M1*4,M1*M1);
k1 = 1:M1;
k2 = 1:M1;
N1 = M2;
N2 = M2;
m = 1:M2;
n = 1:M2;
for i=1:M2*M2
    for j=1:M1*M1
        id1 = mod(j-1,M1)+1;
        id2 = ceil(j/M1);
        idn = ceil(i/M2);
        idm = mod(i-1,M2)+1;
        H(i,j) = 1/M2/M2*exp(2*pi*1i*(k1(id1)-1)*(m(idm)-1)/N1)*exp(2*pi*1i*(k2(id2)-1)*(n(idn)-1)/N2);
        
    end
end

x0=angle(u1)*0+rand(size(u1))*1;    %Initial guess
x0 = x0(:);
tol = 2e-19;
maxit = 15;     %Maximum iteration
iNT = 200;
debug_mode = 1;
if debug_mode
    vidObj = VideoWriter('f1.avi');
    open(vidObj);
end
fprintf('Calling Gauss-Newton optimizer...\n');
[x_recover,hfx,hcrit,iter] = myphase3(H,d,x0,tol,maxit,iNT,cir,debug_mode,ran,vidObj,phase);
if debug_mode
    close(vidObj);
end
%% Plotting PM recovery and convergence curve
figure;set(gcf,'position',[100 100 1200 400 ])
subplot(1,2,1);imagesc(cir.*reshape(x0,M1,M1));axis square;title('Initial Guess');colormap hot;
subplot(1,2,2);imagesc(cir.*reshape(angle(exp(1i*x_recover)),M1,M1));axis square;title('PM of WW');colormap hot;
s = 1;
figure;
 h = semilogy(1:numel(hcrit),hcrit(1:end),'-o');
set(h,'linewidth',2); 
            axis square; grid on
 xlabel('Iteration'); ylabel('scaled g-norm');
            grid on; drawnow; shg
set(gca,'Fontname','Aril','Fontsize',15)
xlabel('Iteration'); ylabel('scaled g-norm');
legend('Convergence Curve')


