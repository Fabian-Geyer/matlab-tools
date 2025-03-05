function animate(t, data, filepath, options)
    arguments
        t (1,:) double
        data (:,:) double {mustBeEqualSize(t, data)}
        filepath (1,1) string
        options.FigSize (1,2) double = [800, 600]
        options.LineWidth (1,1) double = 1.5
        options.animationDuration (1,1) double = 10
        options.XLabel (1,1) string = "Time"
        options.YLabel (1,1) string = "Value"
        options.Title (1,1) string = "Animated Plot"
        options.XLim (1,2) double = [min(t), max(t)]
        options.YLim (1,2) double = [
            min(data(:))-1.1/10*abs(min(data(:))-max(data(:))), ...
            max(data(:))+1.1/10*abs(min(data(:))-max(data(:)))]
        options.LatexStyle logical = true
        options.legend (1,:) string = ""
    end
    
    % Create figure
    fig = figure('Position', [100, 100, options.FigSize(1), options.FigSize(2)]);
    ax = axes(fig);
    set(gcf, 'Color', 'w'); % Sets the background to white
    hold on;
    
    % Initialize plot
    plotHandles = gobjects(1, size(data, 1));
    colors = lines(size(data, 1));
    for i = 1:size(data, 1)
        plotHandles(i) = plot(ax, t(1), data(i, 1), 'Color', colors(i, :), 'LineWidth', options.LineWidth);
    end
    
    % Apply LaTeX style if enabled
    if options.LatexStyle
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
    end
    
    xlabel(options.XLabel);
    ylabel(options.YLabel);
    title(options.Title);
    legend(options.legend, "Interpreter","latex")
    
    xlim(options.XLim);
    ylim(options.YLim);
    grid on;
    
    % Create video writer
    v = VideoWriter(filepath, 'MPEG-4');
    v.FrameRate = size(t,2)/options.animationDuration;
    open(v);
    
    % Animate
    for k = 1:length(t)
        for i = 1:size(data, 1)
            set(plotHandles(i), 'XData', t(1:k), 'YData', data(i, 1:k));
        end
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    % Close video
    close(v);
    close(fig);
end

function mustBeEqualSize(t, data)
    if size(t, 2) ~= size(data, 2)
        error("Time vector and data must have matching dimensions.");
    end
end
