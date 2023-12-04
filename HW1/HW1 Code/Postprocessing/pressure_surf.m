function pressure_surf(Dati,u)

    for i = 1:length(u(:,1))
        for j = 2:length(u(1,:))
            u_t(i,j) = (u(i,j)-u(i,j-1))/Dati.dt;
        end
    end
    u_t = u_t/sqrt(Dati.c2);
    
    x_ax = linspace(0,Dati.domain(end),length(u_t(1,:)));
    y_ax = linspace(0,Dati.T,length(u_t(:,1)));

    figure(10);
    s = surf(x_ax,y_ax,u_t);
    s.EdgeColor = 'None';
    xlabel('Time');
    ylabel('Space');
    zlabel('Pressure');
    view(0,90);
end