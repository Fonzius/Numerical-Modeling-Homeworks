function potential_surf(Dati,u)
    
    x_ax = linspace(0,Dati.domain(end),length(u(1,:)));
    y_ax = linspace(0,Dati.T,length(u(:,1)));
    
    figure(9);
    s = surf(x_ax,y_ax,u);
    s.EdgeColor = 'None';
    xlabel('Time');
    ylabel('Space');
    zlabel('Potential');
    view(0,90);
end