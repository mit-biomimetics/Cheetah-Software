function fig = fn_open_figures(num_figure)

    original_x = 0;
    position_x = original_x;
    position_y = 500;
    width = 560;
    height = 420;
    increase_x = 560;
    increase_y = 420;
    
    x_limit = 2200;
    for i = 1:num_figure
      fig(i) = figure('position', [position_x, position_y, width, height]);
      position_x = position_x + width;
      if( position_x + width > x_limit)
        position_y = position_y -height;
        position_x = original_x;
      end
    end
end