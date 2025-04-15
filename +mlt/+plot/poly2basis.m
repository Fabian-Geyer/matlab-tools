function figHandle = poly2basis(p,opts)
    arguments
        p {isa(p, 'casos.PD')}
        opts.figID char = 'Polynomial Coefficients'
        opts.title char = 'Polynomial Coefficients'
    end

    % Figure
    hfigure(opts.figID);

    % Coefficients
    c = full(poly2basis(p));

    % Monomial labels
    monoms = casos.PD(p.monomials.to_vector);
    labels = arrayfun(@(mon) mon.to_char, monoms, 'UniformOutput',false);

    % bar chart
    barh(c); % horizontal bar chart
    yticks(1:length(monoms)); % Set y-ticks to match number of monomials
    yticklabels(labels); % Assign monomial labels to y-axis
    xlabel('Coefficient Value');
    title(opts.title);
end