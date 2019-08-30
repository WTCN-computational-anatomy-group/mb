function ShowSubjects(dat,mu,sett)

fg = findobj('Type', 'Figure', 'Name', sett.show.figname_subjects);
if isempty(fg)
    fg = figure('Name', sett.show.figname_subjects, 'NumberTitle', 'off');
else
    clf(fg);
end
set(0, 'CurrentFigure', fg);   

nd = min(numel(dat),sett.show.mx_subjects);
for n=1:nd
    d    = GetSize(dat(n).f);
    q    = double(dat(n).q);
    Mn   = dat(n).Mat;
    psi1 = GetData(dat(n).psi);
    psi  = compose(psi1,affine(d,sett.Mmu\spm_dexpm(q,sett.B)*Mn));
    mu1  = Pull1(mu,psi);
    
    % Get resp, image and template
    r   = GetClasses(dat(n),mu1,sett);
    f   = GetData(dat(n).f);
    mu1 = softmax(mu1);
    
    % Make plots    
    if sett.is3d
        nr = 9;
        for ax=1:3
            ShowIm(f,ax,nr,nd,   n + 3*(ax - 1)*nd,sett.show.figname_subjects);
            ShowCat(r,ax,nr,nd,  n + nd + 3*(ax - 1)*nd,sett.show.figname_subjects);
            ShowCat(mu1,ax,nr,nd,n + 2*nd + 3*(ax - 1)*nd,sett.show.figname_subjects);
        end
    else
        nr = 5;
        
        % Intensity histogram
        subplot(nr,nd,n)
        hist(f(isfinite(f)),40)                
        set(gca,'ytick',[])  
        
        % Image, segmentation, template
        ax = 3;
        ShowIm(f,ax,nr,nd,   n + nd,sett.show.figname_subjects);
        ShowCat(r,ax,nr,nd,  n + 2*nd,sett.show.figname_subjects);
        ShowCat(mu1,ax,nr,nd,n + 3*nd,sett.show.figname_subjects);
                        
        % Affine parameters
        M    = spm_dexpm(q,sett.B);
        q    = spm_imatrix(M);            
        q    = q([1 2 6]);
        q(3) = 180/pi*q(3);
        q    = abs(q);
        subplot(nr,nd,n + 4*nd) 
        hold on
        for k=1:numel(q)
            bar(k, q(k));
        end
        box on
        hold off            
        set(gca,'xtick',[])  
    end
end
drawnow
end
%==========================================================================