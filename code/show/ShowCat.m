function ShowCat(in,ax,nr,nc,np,fn)

f  = findobj('Type', 'Figure', 'Name', fn);
if isempty(f)
    f = figure('Name', fn, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f); 

subplot(nr,nc,np);

dm  = size(in);
K   = dm(4);
pal = hsv(K);
mid = ceil(dm(1:3).*0.55);

if ax == 1
    in = permute(in(mid(1),:,:,:),[3 2 1 4]);
elseif ax == 2
    in = permute(in(:,mid(2),:,:),[1 3 2 4]);
else
    in = in(:,:,mid(3),:);
end

tri = false;
if numel(size(in)) == 4 && size(in, 3) == 1
    tri = true;
    dm  = [size(in) 1 1];
    in   = reshape(in, [dm(1:2) dm(4)]);
end
if isa(pal, 'function_handle')
    pal = pal(size(in,3));
end

dm = [size(in) 1 1];
c   = zeros([dm(1:2) 3]); % output RGB image
s   = zeros(dm(1:2));     % normalising term

for k=1:dm(3)
    s = s + in(:,:,k);
    color = reshape(pal(k,:), [1 1 3]);
    c = c + bsxfun(@times, in(:,:,k), color);
end
if dm(3) == 1
    c = c / max(1, max(s(:)));
else
    c = bsxfun(@rdivide, c, s);
end

if tri
    c = reshape(c, [size(c, 1) size(c, 2) 1 size(c, 3)]);
end

c = squeeze(c(:,:,:,:));
imagesc(c); axis off image xy;   

end 
%==========================================================================