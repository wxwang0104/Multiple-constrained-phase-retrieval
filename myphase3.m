    function [x,hfx,hcrit,iter,x_rec] = myphase3(H,d,x0,tol,maxit,iNT,cir,debug_mode,range,varargin)
    % Solve the phase retrieval least squares problem
    %        min 0.5*norm(|Hx|^2 - d)^2
    %
    % input:
    %       H, d = data
    %         x0 = initial guess
    %        tol = tolerance
    %      maxit = maximum iteration number
    %        iNT = switch to Newton at this iteration 
    % output:
    %          x = computed solution
    %        hfx = history vector storing fk = f(x_k)
    %      hcrit = history vector storing norm(gk)/(1+fk)
    %       iter = number of iteration taken
    %
    bs = 20;
    x = x0;
    iter = 0;
    M1 = sqrt(size(H,2));
    M1 = (M1);
    M2 = M1*2;
    
    m = size(H,1);
    crit = tol+1;
    while ((crit>tol)&&(iter<maxit))
        if iter<iNT
            
            iter = iter + 1;
            x_rec(:,iter) = x;
            [f,J,r,ind1,ind2]=func2(x,H,d,cir,range);
            delta = (J'*J)\(J'*r(ind2));
            x(ind1) = x(ind1) - delta*1;
           % x = x-floor(abs(x)).*sign(x);
            crit = norm(J'*r(ind2))/(1+abs(f));
            hcrit(iter) = crit;
            hfx(iter) = f;
            fprintf('iter: %3d f = %.4e, norm(g) = %.2e, crit = %.2e\n',iter,f,norm(J'*r(ind2)),crit);
       
        end
        if debug_mode
            pm = reshape(x,M1,M1);
            N=M1;
            fig = figure(101);
            set(gcf,'position',[100 100 1600 400 ])
            for i=-4:4
                phase=phase_ex_04022015(150/1.0,i*range*1e-7,N,N);
                u2 = H*reshape(cir.*exp(1i*pm+1i*phase),M1^2,1);
                iii = fftshift(reshape(d(:,i+5),M2,M2));
                subplot(2,9,i+5);imagesc(iii(M2/2-bs+1:M2/2+bs,M2/2-bs+1:M2/2+bs));axis square;caxis([0 3]*1e-3)
                set(gca,'Ytick',[]);
                set(gca,'Xtick',[]);
                iii = fftshift(reshape(abs(u2).^2,M2,M2));
                subplot(2,9,i+5+9);imagesc(iii(M2/2-bs+1:M2/2+bs,M2/2-bs+1:M2/2+bs));axis square;caxis([0 3]*1e-3)
                set(gca,'Ytick',[]);
                set(gca,'Xtick',[]);
            end
            text(-460/(40/bs)+2,-40*bs/10,['Iteration',num2str(iter)],'color','black','Fontsize',13);
            currFrame = getframe(fig);
            for rep = 1:1
                %writeVideo(varargin{1},currFrame);
            end
        end
    end
    %}
    
    end
    

    
function [f,J,r,ind1,ind2] = func2(x,H,d,cir,range)
    m = size(H,1);
    n = size(H,2);
    a = real (H);
    b = imag (H);
    r_r=[];
    f_r=[];
    J_r=[];
    N = sqrt(n);
    for i=-4:4
        phase=phase_ex_04022015(150/1.0,i*range*1e-7,N,N);
        r = abs(H*(cir(:).*exp(1i*x).*exp(1i*phase(:)))).^2-d(:,i+5);
               
        a_v1 = a*(cir(:).*cos(x+phase(:)))-b*(cir(:).*sin(x+phase(:)));
        a_t1 = repmat (a_v1,1,n);
        
        b_v1 = b*(cir(:).*cos(x+phase(:)))+a*(cir(:).*sin(x+phase(:)));
        b_t1 = repmat (b_v1,1,n);
        
        sin_v = (cir(:)').*sin(x'+phase(:)');
        sin_t = repmat (sin_v,m,1);
        
        cos_v = (cir(:)').*cos(x'+phase(:)');
        cos_t = repmat (cos_v,m,1);
        
        J = -a.*sin_t.*a_t1+a.*cos_t.*b_t1-b.*cos_t.*a_t1-b.*sin_t.*b_t1;
        J = J*2;
        J_r = [J_r;J];
        r_r = [r_r;r];
    end
    r = r_r;
    f = norm(r)^2/2;
    J = J_r;
    J_t = sum(J);
    ind1 = find(abs(J_t)>0);
    J = J(:,ind1);
    J_t = sum(J');
    ind2 = find(abs(J_t)>0);
    J = J(ind2,:);
   % figure;imagesc(J);
   % 1+1;
   %}
    
    
end


    