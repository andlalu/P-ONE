function [ theta ] = freetobound(type,theta,lbound,ubound)
%boundtofree This function reverses a simple exponential transformation
switch type
    case 'SVHJ'
        c=1;
        theta=lbound+((ubound-lbound)/2).*(sin(theta)+1); %transform back from fminsearch output
        theta(9) = theta(11).*c + theta(9); %theta(9) stores the difference, theta(11) the level
        disp(theta)
%         theta(1) = theta(2) + theta(1); %theta(1) stores the difference, theta(2) the level
    case 'SV' 
%       theta=lbound+(ubound-lbound).*(1./(1+exp(-theta)));
        theta=lbound+((ubound-lbound)/2).*(sin(theta)+1);
    case 'SVJ' 
        theta=lbound+(ubound-lbound).*(1./(1+exp(-theta)));
        theta(1) = theta(2) + theta(1);
    case 'PAN'
        theta=lbound+(ubound-lbound).*(1./(1+exp(-theta)));
        theta(1) = theta(2) + theta(1);
    case 'COX'
        theta=lbound+(ubound-lbound).*(1./(1+exp(-theta)));
end
