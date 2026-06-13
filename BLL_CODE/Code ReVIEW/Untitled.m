c=100;
c2=1;
t=-100:.1:100;t=t';
t1=[t,zeros(size(t)),zeros(size(t))];
for i=1:length(t)
        cft1=ccf(t1(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
        test1(i)=exp(sum(cft1.*[1,0.01,0.15^2*c2,5/c],2));
%                 H0(:,i)=(exp(1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)),2))-
end

t2=-100:.1:100;t2=t2';
t2=[zeros(size(t2)),t2,zeros(size(t2))];
for i=1:length(t2)
        cft2=ccf(t2(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
        test2(i)=exp(sum(cft2.*[1,0.01,0.15^2*c2,5/c],2));
%                 H0(:,i)=(exp(1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)),2))-
end

t=-100:.1:100;t=t';
t3=[zeros(size(t)),zeros(size(t)),t];
for i=1:length(t)
        cft3=ccf(t3(i,:).',dt,K0,K1,H00,H1,L0,L1,JT);
        test3(i)=exp(sum(cft3.*[1,0.01,0.15^2*c2,5/c],2));
%                 H0(:,i)=(exp(1i*sum(repmat(x(i,:),[N-1,1]).*(X(2:end,:)),2))-
end



% subplot(1,3,1)
% plot(t,real(test1),t,imag(test1));
% 
% 
% subplot(1,3,2)
% plot(t,real(test2),t,imag(test2));
% 
% 
% subplot(1,3,3)
% plot(t,real(test3),t,imag(test3));
