
%Experimental implementation of the criterion function. 

function [ out ] = critstep1( theta, H0, X, dt, type, step2, factor, alpha)
%1st step estimator using C-GMM based on the characteristic function 


%Recognize the data-set
[N,M]=size(X);

%C-GMM Procedure controls
    [x,w]=nwspgr('KPN', M, 5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXPERIMENTAL REVIEW OF SV PART        
  %  case 'SV'
        sig_r1=1;
        sig_r2=1;
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
        
        %%%H0?
          %Compute the C-GMM criterion function 
        %%% Compute the H0 Matrix:
        H0=zeros([N-1,length(w)]);
        for i=1:length(w)
           cf=ccf(-x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
           H0(:,i)=sqrt(w(i)).*factor.*(exp(-1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)),2))-exp(sum(repmat(cf,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),X(1:end-1,2:end)],2)));
        end
          
        
        %%% Compute the FT part:
        X1=repmat(X(2:end,1)-X(1:end-1,1),[1,N-1]);
        DX1=X1-X1';
        X2=repmat(X(2:end,2),[1,N-1]);
        DX2=X2-X2';
        Ft=exp(-0.5*sig_r1^2*DX1.^2+-0.5*sig_r2^2*DX2.^2);
        H=zeros([N-1,length(w)]);
        for i=1:length(w)
            cf2=ccf(-x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
            H(:,i)=sqrt(w(i)).*factor.*(exp(-1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)),2))-exp(sum(repmat(cf2,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),X(1:end-1,2:end)],2)));
        end
        p=H0*H'/N;
        p=Ft.*p;
        p=sum(p,2);
        if step2
            Cm=H0*H0'/(N-5);
            Cm=Ft.*Cm;
            OptMatrix=inv(alpha*eye([N-1,N-1])+Cm^2);
            out=real(p'*OptMatrix*p);
        else
            out=real(p'*eye*p);
        end
        if 2*k*vbar<sigma^2 %artificially impose feller condition
          out=out+1e4;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXPERIMENTAL REVIEW OVER ABOVE ^^^^        

end

