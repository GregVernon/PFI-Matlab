%%
close all
iter = 0;
for kk = 300000:length(DATA)
    iter = iter + 1;
    if iter == 1
        p = plot(DATA(1,:,kk),'k');
    else
        p.YData = DATA(1,:,kk);
    end
    if mod(kk,1000) == 0
        title(kk);
    end
    
    drawnow
end

%%
figure
load('SurfaceVelocity.mat');
DATA = squeeze(DATA);
p = surf(U);
p.LineStyle = 'none';
p.FaceColor = 'interp';
ax1 = gca();
ax1.CLim = [-max(max(abs(DATA))) max(max(abs(DATA)))];
%%
close all
U = load('U-VelocityProfile.mat'); U = squeeze(U.DATA);
V = load('V-VelocityProfile.mat'); V = squeeze(V.DATA);
for kk = 1:size(U,2);
    if kk == 1
        q = quiver(U(:,kk),V(:,kk),'ShowArrowHead','off','AutoScale','off');
    else
        q.UData = U(:,kk);
        q.VData = V(:,kk);
    end
    if mod(kk,100) == 0
        title(kk);
    end
    
    drawnow
end
