function Path_Plot(alpha,p2n,fig,view_az_el)
% PATH_PLOT Plot geodesic path
% -------------------------------------------------------------------------
% 
% Usage: Path_Plot(alpha,p2n,fig,col,view_az_el)
%
% Input:
% alpha: matrix (n,T,M) geodesic of curves
% p2n: registered final curve
% fig: figure number
% view_az_el: view of az el

% Output:

[n,T,k] = size(alpha);
stp = k-1;
if(n == 3)
    dt = 0.3;
else
    dt = 0.25;
end

figure(fig); clf;
axis ij;
f1 = zeros(3,T,stp+1);
for tau = 1:stp+1
    
    if (tau==stp+1)
        p = p2n;
    else
        p = q_to_curve(alpha(:,:,tau));
    end

    ft(1,:) = p(1,:) ;
    ft(1,:) = ft(1,:) - mean(ft(1,:)) +dt*tau;
    for i = 2:n
        ft(i,:) = p(i,:) - mean(p(i,:));
    end

    if(n == 2)
        z = plot(ft(1,:),ft(2,:),'Color',[(tau-1)/6 (tau-1)/12 0],'LineWidth',1); axis equal; hold on;axis equal; axis off; axis tight;  
        set(z,'LineWidth',3);
    else
        if(n == 3)
            f1(1,:,tau) = ft(1,:);
            f1(2,:,tau) = ft(2,:);
            f1(3,:,tau) = ft(3,:);
            plot3(ft(1,:),ft(2,:),ft(3,:),'Color',[(tau-1)/6 (tau-1)/12 0],'LineWidth',1.5); axis equal; hold on;axis equal; axis off; axis tight;             
            view(view_az_el(1),view_az_el(2));axis on;
            set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');     % Remove the tick labels
            set(gca,'TickLength',[0 0]);                  % Remove the ticks
        end
    end
end
axis ij;
