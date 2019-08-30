function ShowModel(mu,Objective,sett,N)
mu = softmax(mu);
if sett.is3d
    ShowCat(mu,1,2,3,1,sett.show.figname_model);
    ShowCat(mu,2,2,3,2,sett.show.figname_model);
    ShowCat(mu,3,2,3,3,sett.show.figname_model);
    subplot(2,1,2); plot(Objective,'.-');
    title(['K=' num2str(size(mu,4)) ', N=' num2str(N)]);
else
    ShowCat(mu,3,1,2,1,sett.show.figname_model);
    title(['K=' num2str(size(mu,4)) ', N=' num2str(N)]);
    subplot(1,2,2); plot(Objective,'.-');    
end
drawnow
end
%==========================================================================