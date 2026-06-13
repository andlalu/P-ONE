%%%% 12.08.2014 V2.0  Beta for: any model;
function [ answ ] = quickcrit_true_states( X, dt, type, theta,lb,ub,theta1step)    
    

%%
  % Incercare de a reprograma functia de criteriu pentru a vedea daca optimizeaza cum trebuie: 
        %Criterion characteristics
        factor=100;
        [N,M]=size(X);
        [x,w]=nwspgr('KPN', M, 7);
        sig_r1=0.01;
        sig_r2=0.01;
        X2(:,1)=X(2:end,1)-X(1:end-1,1);
        X2(:,2)=X(2:end,2);
        N2=N-1;
        
        %Allocate parameters
        eta=theta(1);
        k=theta(2);
        vbar=theta(3);
        sigma=theta(4);
        rho=theta(5);
        
        %Define CCF system matrixes
        K0 = [ 0;      k*vbar;];
        K1 = [ 0, eta-0.5;
            0,  -k;];
        H00 = [ 0,      0;
            0,      0;];
        H1(:,:,1) = [ 0,0; 0,0;];
        H1(:,:,2) = [ 1,sigma*rho; sigma*rho,sigma^2;];
        L0=0;
        L1=[0;0;];
        JT=@(beta) (0);
        
        %Aici incepe partea cu calculul functiei cum vreau eu: 
         for i=1:length(w)
           %fourier part de la inceput
%            fourier=exp(i
           cf2=ccf(x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
           exp(1i*sum(repmat(x(i,:),[N2-1,1]).*X2))        
           exp(sum(repmat(cf2,[N2-1,1]).*X),2);
        end


%%
        %Calculata prin metoda veche:
        %%% Compute the FT part:
        X1=repmat(X(2:end,1)-X(1:end-1,1),[1,N-1]);
        DX1=X1-X1';
        X2=repmat(X(2:end,2),[1,N-1]);
        DX2=X2-X2';
        Ft=exp(-0.5*sig_r1^2*DX1.^2+-0.5*sig_r2^2*DX2.^2);
        
        H=zeros([N-1,length(w)]);
        for i=1:length(w)
            cf2=ccf(-x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
            H(:,i)=sqrt(w(i)).*factor.*(exp(-1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)-[X(1:end-1,1),zeros(N-1,1)]),2))-exp(sum(repmat(cf2,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),X(1:end-1,2:end)],2)));
        end
        p=H*H'/N;
        p=Ft.*p;
        p=sum(p,2);
%         if step2
%             Cm=H0*H0'/(N-5);
%             Cm=Ft.*Cm;
%             OptMatrix=inv(alpha*eye([N-1,N-1])+Cm^2);
%             out=real(p'*OptMatrix*p);
%         else
            answ=p'*eye*p;
%         end
%         if 2*k*vbar<sigma^2
%           out=out+1e4;
%         end
  
end
