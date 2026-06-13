function [ out ] = critstep2(theta, X, dt, type,H0,OptMatrix)
%1st step estimator using C-GMM based on the characteristic function 


%Recognize the data-set
[N,M]=size(X);

%Recognize which model to estimate
switch type 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXPERIMENTAL REVIEW OF SV PART        
    case 'SV'
        %C-GMM Procedure controls
%         [x,w]=nwspgr('KPN', 2, 5);
% [w,x]=fwtpts( 2, 5, 'Norm');
% x=x';w=w';
        [x,w]=nwspgr('GQU', M, 5);
        x=norminv(x,0,1);
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
        %%% Compute the FT part:
        sig_r1=10;
        sig_r2=100;
        X1=repmat(X(1:end-1,1),[1,N-1]);
        DX1=X1-X1';
        X2=repmat(X(1:end-1,2),[1,N-1]);
        DX2=X2-X2';
        Ft=exp(-0.5*sig_r1^2*DX1.^2+-0.5*sig_r2^2*DX2.^2);
        %%% Compute the criterion function part: 
        H=zeros([N-1,length(w)]);
            for i=1:length(w)
                cf2=ccf(x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
                H(:,i)=(exp(1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)),2))-exp(sum(repmat(cf2,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),X(1:end-1,2:end)],2)));
            end
        w=sum(Ft.*(H0.*(ones(N-1,1)*(w.'))*H'),2);   
        out=real(w).'*real(OptMatrix)*real(w)./(N-1)^2;
        
    case 'SVHJ'
   %C-GMM Procedure controls
% % % x=[1:5:100];
% % % x=[x, zeros(size(x)), flip(x);   zeros(size(x)), x, zeros(size(x));     zeros(size(x)),zeros(size(x)),x]';
% % % w=exp(-0.5*sum(x.*x,2)/1e4);
% % % % w=w./(sum(w));
[w,x]=fwtpts( 3, 5, 'Norm');
x=x';
w=w';

%         [x,w]=nwspgr('KPU', M,4);
%         x=norminv(x,0,5);
    %x=x/sqrt(10); %Narrowing the integration range to accommodate Hawkes Jumps
    %Allocate parameters    
    muj=theta(1);
    mujq=theta(2);
    sigmaj=theta(3);
    eta=theta(4);
    k=theta(5);
    vbar=theta(6);
    sigma=theta(7);
    rho=theta(8);
    kl=theta(9);
    lambda=theta(10);
    delta=theta(11);

    %Define CCF system matrixes
    c=100;
    K0 = [ 0;    k*vbar;   kl*lambda];
    K1 = [ 0,   eta-0.5, -c*(exp(mujq+0.5*sigmaj^2)-1);
           0,     -k,                           0;
           0,      0,                         -kl; ];
    H00=[ 0,   0,   0;
          0,   0,   0;
          0,   0,   0; ];
    H1(:,:,1) = [ 0,0,0; 0,0,0; 0,0,0;];
    H1(:,:,2) = [ 1,sigma*rho, 0; sigma*rho,sigma^2,0; 0,0,0;];
    H1(:,:,3) = [ 0,0,0; 0,0,0; 0,0,0;];
    L0 = 0;
    L1 = [0;0;c];
    JT = @(beta) (exp(beta(1)*muj+0.5*beta(1).^2*sigmaj^2+delta*beta(3))-1);
        
    %%% Compute the FT part:
        sig_r1=10;
        sig_r2=1;
        sig_r3=10;
        X1=repmat(X(2:end,1)-X(1:end-1,1),[1,N-1]);
        DX1=X1-X1';
        X2=repmat(X(1:end-1,2),[1,N-1]);
        DX2=X2-X2';
        X3=repmat(X(1:end-1,3),[1,N-1]);
        DX3=X3-X3';
        Ft=exp(-0.5*sig_r1^2*DX1.^2+-0.5*sig_r2^2*DX2.^2+-0.5*sig_r3^2*DX3.^2);
    %%% Compute the criterion function part: 
            H=zeros([N-1,length(w)]);
                for i=1:length(w)
                    cf2=ccf(x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
                    H(:,i)=(exp(1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)),2))-exp(sum(repmat(cf2,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),X(1:end-1,2:end)],2)));
                end
            w=sum(Ft.*(H0.*(ones(N-1,1)*(w.'))*H'),2);   
            out=real(w).'*real(OptMatrix)*real(w);%./(N-1)^2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% % % % % % % % %         case 'SVJ'
% % % % % % % % %         sig_r1=1;
% % % % % % % % %         sig_r2=1;
% % % % % % % % %         %Allocate parameters
% % % % % % % % %         muj=theta(1);
% % % % % % % % %         mujq=theta(2);
% % % % % % % % %         sigmaj=theta(3);
% % % % % % % % %         eta=theta(4);
% % % % % % % % %         k=theta(5);
% % % % % % % % %         vbar=theta(6);
% % % % % % % % %         sigma=theta(7);
% % % % % % % % %         rho=theta(8);
% % % % % % % % %         lambda=theta(9);
% % % % % % % % %         %Define CCF system matrixes
% % % % % % % % %         K0 = [ -(exp(mujq+0.5*sigmaj^2)-1)*lambda;      k*vbar;];
% % % % % % % % %         K1 = [ 0, eta-0.5;
% % % % % % % % %             0,  -k;];
% % % % % % % % %         H00 = [ 0,      0;
% % % % % % % % %             0,      0;];
% % % % % % % % %         H1(:,:,1) = [ 0,0; 0,0;];
% % % % % % % % %         H1(:,:,2) = [ 1,sigma*rho; sigma*rho,sigma^2;];
% % % % % % % % %         L0=lambda;
% % % % % % % % %         L1=[0;0;];
% % % % % % % % %         JT=@(beta) (exp(beta(1)*muj+0.5*beta(1)^2*sigmaj^2)-1);
% % % % % % % % %         %%% Compute the FT part:
% % % % % % % % %         X1=repmat(X(2:end,1)-X(1:end-1,1),[1,N-1]);
% % % % % % % % %         DX1=X1-X1';
% % % % % % % % %         X2=repmat(X(1:end-1,2),[1,N-1]);
% % % % % % % % %         DX2=X2-X2';
% % % % % % % % %         Ft=exp(-0.5*sig_r1^2*DX1.^2+-0.5*sig_r2^2*DX2.^2);
% % % % % % % % %         H=zeros([N-1,length(w)]);
% % % % % % % % %         for i=1:length(w)
% % % % % % % % %             cf2=ccf(-x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
% % % % % % % % %             H(:,i)=sqrt(w(i)).*factor.*(exp(-1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)-[X(1:end-1,1),zeros(N-1,2)]),2))-exp(sum(repmat(cf2,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),X(1:end-1,2:end)],2)));
% % % % % % % % %         end
% % % % % % % % %         p=H0*H'/N;
% % % % % % % % %         p=Ft.*p;
% % % % % % % % %         p=sum(p,2);
% % % % % % % % %         if step2
% % % % % % % % %             Cm=H0*H0'/(N-9);
% % % % % % % % %             Cm=Ft.*Cm;
% % % % % % % % %             OptMatrix=inv(alpha*eye([N-1,N-1])+Cm^2);
% % % % % % % % %             out=real(p'*OptMatrix*p);
% % % % % % % % %         else
% % % % % % % % %             out=real(p'*eye*p);
% % % % % % % % %         end
% % % % % % % % %         if 2*k*vbar<sigma^2
% % % % % % % % %             out=out+1e4;
% % % % % % % % %         end
% % % % % % % % %         
% % % % % % % % %         case 'PAN'
% % % % % % % % %         sig_r1=1;
% % % % % % % % %         sig_r2=1;
% % % % % % % % %         %Allocate parameters
% % % % % % % % %         muj=theta(1);
% % % % % % % % %         mujq=theta(2);
% % % % % % % % %         sigmaj=theta(3);
% % % % % % % % %         eta=theta(4);
% % % % % % % % %         k=theta(5);
% % % % % % % % %         vbar=theta(6);
% % % % % % % % %         sigma=theta(7);
% % % % % % % % %         rho=theta(8);
% % % % % % % % %         lambda1=theta(9);
% % % % % % % % %         lambda0=0;
% % % % % % % % %         %Define CCF system matrixes
% % % % % % % % %         K0 = [ -(exp(mujq+0.5*sigmaj^2)-1)*lambda0;      k*vbar;];
% % % % % % % % %         K1 = [ 0, eta-0.5-(exp(mujq+0.5*sigmaj^2)-1)*lambda1;
% % % % % % % % %             0,  -k;];
% % % % % % % % %         H00 = [ 0,      0;
% % % % % % % % %             0,      0;];
% % % % % % % % %         H1(:,:,1) = [ 0,0; 0,0;];
% % % % % % % % %         H1(:,:,2) = [ 1,sigma*rho; sigma*rho,sigma^2;];
% % % % % % % % %         L0=lambda0;
% % % % % % % % %         L1=[0;1];
% % % % % % % % %         JT=@(beta) (exp(beta(1)*muj+0.5*beta(1)^2*sigmaj^2)-1);
% % % % % % % % %         %%% Compute the FT part:
% % % % % % % % %         X1=repmat(X(2:end,1)-X(1:end-1,1),[1,N-1]);
% % % % % % % % %         DX1=X1-X1';
% % % % % % % % %         X2=repmat(X(1:end-1,2),[1,N-1]);
% % % % % % % % %         DX2=X2-X2';
% % % % % % % % %         Ft=exp(-0.5*sig_r1^2*DX1.^2+-0.5*sig_r2^2*DX2.^2);
% % % % % % % % %         H=zeros([N-1,length(w)]);
% % % % % % % % %         for i=1:length(w)
% % % % % % % % %             cf2=ccf(-x(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
% % % % % % % % %             H(:,i)=sqrt(w(i)).*factor.*(exp(-1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)-[X(1:end-1,1),zeros(N-1,2)]),2))-exp(sum(repmat(cf2,[N-1,1]).*[ones([N-1,1]),zeros([N-1,1]),X(1:end-1,2:end)],2)));
% % % % % % % % %         end
% % % % % % % % %         p=H0*H'/N;
% % % % % % % % %         p=Ft.*p;
% % % % % % % % %         p=sum(p,2);
% % % % % % % % %         if step2
% % % % % % % % %             Cm=H0*H0'/(N-9);
% % % % % % % % %             Cm=Ft.*Cm;
% % % % % % % % %             OptMatrix=inv(alpha*eye([N-1,N-1])+Cm^2);
% % % % % % % % %             out=real(p'*OptMatrix*p);
% % % % % % % % %         else
% % % % % % % % %             out=real(p'*eye*p);
% % % % % % % % %         end
% % % % % % % % %         if 2*k*vbar<sigma^2
% % % % % % % % %             out=out+1e4;
% % % % % % % % %         end
        
end
end

