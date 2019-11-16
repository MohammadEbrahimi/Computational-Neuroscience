
function axesx=pix2axes(dim, xdata, pixelx)

xfirst = xdata(1);
xlast = xdata(max(size(xdata)));

delta = xlast - xfirst;
if delta == 0
  xslope = 1;
else
  xslope = (dim - 1) / delta;
end

if ((xslope == 1) && (xfirst == 1))
  axesx = pixelx;
else
  axesx = (pixelx-1)/xslope + xfirst;
end
