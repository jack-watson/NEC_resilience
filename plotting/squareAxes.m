function squareAxes(ax)

% Make x/y ranges equal around their centers so the plot box can use full height
xl = xlim(ax);  yl = ylim(ax);
cx = mean(xl);  cy = mean(yl);
L  = max(diff(xl), diff(yl));
if L == 0 % guard against degenerate ranges
    L = 1;
end
xlim(ax, [cx - L/2, cx + L/2]);
ylim(ax, [cy - L/2, cy + L/2]);
pbaspect(ax, [1 1 1]); % proportional axes (X and Y commensurate)

end