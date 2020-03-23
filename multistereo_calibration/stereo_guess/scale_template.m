function pts_scaled = scale_template(pts2D, template)
%SCALE_TEMPLATE return 3D points in template frame with real scale

global TEMPLATE_SHORT_SIDE;

templ_im_side = min(size(template));

pts_scaled = (pts2D / templ_im_side) * TEMPLATE_SHORT_SIDE;
end

