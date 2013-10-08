function dt=delta_t(Nx,My,Xmin,Xmax,Ymin,Ymax,Re)
    dx=abs(Xmax-Xmin)/(Nx+1);
    dy=abs(Ymax-Ymin)/(My+1);
    if dx<dy
        dt=4/(Re);
    else
        dt=4*dy/(Re*dx);
    end
end