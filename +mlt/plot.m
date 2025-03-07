function h = plot(x, y, varargin)
    % Create a new figure or use the current one
    if isempty(findobj('Type', 'figure'))
        h = hfigure('generic plot');
    else
        h = gcf;
    end

    % Plot the data
    plotHandle = plot(x, y, varargin{:});

    % Apply default settings
    set(plotHandle, 'LineWidth', 1); % Set default line width
    set(gca, 'TickLabelInterpreter', 'latex'); % Set interpreter for tick labels
    set(gca, 'FontSize', 12); % Set default font size
    title('', 'Interpreter', 'latex'); % Set title interpreter
    xlabel('', 'Interpreter', 'latex'); % Set x-label interpreter
    ylabel('', 'Interpreter', 'latex'); % Set y-label interpreter

    % Set legend interpreter to LaTeX
    if ~isempty(legend)
        legendHandle = legend;
        set(legendHandle, 'Interpreter', 'latex');
    end


    % Return the plot handle if needed
    if nargout > 0
        h = plotHandle;
    end
end