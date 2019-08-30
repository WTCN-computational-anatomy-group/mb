function ShowIm(in,ax,nr,nc,np,fn)

f  = findobj('Type', 'Figure', 'Name', fn);
if isempty(f)
    f = figure('Name', fn, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   

subplot(nr,nc,np);

dm  = size(in);
dm  = [dm 1];
mid = ceil(dm(1:3).*0.55);

if ax == 1
    in = permute(in(mid(1),:,:,:),[3 2 1 4]);
elseif ax == 2
    in = permute(in(:,mid(2),:,:),[1 3 2 4]);
else
    in = in(:,:,mid(3),:);
end

in = squeeze(in(:,:,:,:));
imagesc(in); axis off image xy;   
colormap(gray);

end 
%==========================================================================