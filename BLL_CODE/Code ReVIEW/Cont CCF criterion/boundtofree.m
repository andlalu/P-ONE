function [ free_theta ] = boundtofree(type,theta,lbound,ubound)
%boundtofree This function applies a simple exponential transformation
switch type
    case 'SVHJ'
        c=1;
        theta(9) = theta(9) - theta(11).*c; % delta (theta(11)) is level, kl = delta (theta(11)) + difference (stored for fminsearch in theta(9))
%         theta(1) = theta(1) - theta(2); %mujq (theta(2)) is level, mujp = mujq(theta(2)) + difference (stored for fminsearch in theta(1))
        free_theta = asin(2*((theta-lbound)./(ubound-lbound))-1);
    case 'SV' 
%         free_theta = -log((ubound-lbound)./(theta-lbound)-1);
          free_theta = asin(2*((theta-lbound)./(ubound-lbound))-1);
    case 'SVJ'
        theta(1) = theta(1) - theta(2); %mujq (theta(2)) is level, mujp = mujq(theta(2)) + difference (stored for fminsearch in theta(1))
        free_theta = -log((ubound-lbound)./(theta-lbound)-1);
    case 'PAN'
        theta(1) = theta(1) - theta(2); %mujq (theta(2)) is level, mujp = mujq(theta(2)) + difference (stored for fminsearch in theta(1))
        free_theta = -log((ubound-lbound)./(theta-lbound)-1);
    case 'COX'
        error('not programmed')
end

end