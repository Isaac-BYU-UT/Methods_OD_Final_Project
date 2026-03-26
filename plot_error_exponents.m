function [] = plot_error_exponents(relDiff_matrix)
    % --- SETTINGS ---
    fntSize = 16; 
    folderName = 'Figures'; 
    
    % *** SELECT MILK SCHEME ***
    % Option 1: 'WhiteToDarkBlue' (Light to Dark Sequential)
    % Option 2: 'BlueToWhiteToRed' (Diverging: White in middle)
    % Option 3: 'CustomCream' (Lighter, creamier version of jet)
    milkScheme = 'CustomCream'; 
    % ------------------------

    varName = inputname(1);
    if isempty(varName), varName = 'Matrix'; end

    % 1. Histogram (remains same)
    data = relDiff_matrix(:);
    logData = log10(abs(data));
    logData = logData(isfinite(logData)); 

    if ~isempty(logData)
        fig1 = figure('Visible', 'off'); % Set off to avoid pops
        edges = floor(min(logData)) : 1 : ceil(max(logData));
        histogram(logData, 'BinEdges', edges, 'FaceColor', [0.3 0.3 0.8]);
        xlabel('Order of Magnitude (10^x)'); ylabel('Number of Elements');
        title("Error Precision Distribution of " + strrep(varName,"_"," "));
        grid on; xticks(edges); set(gca, 'FontSize', fntSize);
        exportgraphics(fig1, fullfile(folderName, varName + "_Distribution.png"), 'Resolution', 300);
        close(fig1); % Auto-close
    end

    % 2. Heatmap Plot
    logMatrix = log10(abs(relDiff_matrix));
    maxVal = max(logData);
    minVal = min(logData);
    [rows, cols] = size(logMatrix);
    
    fig2 = figure('Visible', 'on');
    
    % --- Apply Milk Scheme and Contrast Colors ---
    switch milkScheme
        case 'WhiteToDarkBlue'
            % Option 1: White (lightest error) to Dark Blue (highest error)
            cm = parula; 
            cm = flipud(cm); % Reverse so dark is high
            colormap(gca, cm);
            bgRGB = [0.8 0.8 0.8]; % Light Gray for NaN
            
            infLabelColor = 'black';
            nanLabelColor = 'black';
            
            % Text Contrast Logic
            contrastLimit = minVal + 3; % Flip near the dark side
            flipToWhite = true; 

        case 'BlueToWhiteToRed'
            % Option 2: True milk scheme (white in the middle, blues and reds at edges)
            cm = flipud(colormap(gca, 'cool')); % Custom diverging
            cm = [ones(128,3); cm(129:end,:)]; % Inject white
            colormap(gca, cm);
            bgRGB = [0.2 0.2 0.2]; % Dark Gray for NaN (high contrast)
            
            infLabelColor = [0.8 0.8 0.8]; % Light gray
            nanLabelColor = [0.8 0.8 0.8];
            
            % Text Contrast Logic
            contrastLimit = (minVal + maxVal)/2; % Midpoint
            flipToWhite = false; % High values are red (needs black)

        case 'CustomCream'
            % Option 3: Retains jet, but applies a "milk wash" (lightened colors)
            cm = jet;
            cm = cm + (ones(size(cm)) - cm)*0.6; % Lighten all colors by 60%
            colormap(gca, cm);
            bgRGB = [0.95 0.95 0.95]; % Almost White for NaN
            
            infLabelColor = 'black';
            nanLabelColor = 'black';
            
            % Text Contrast Logic
            contrastLimit = minVal + 3; % Near the bottom
            flipToWhite = true;
            
        otherwise
            error('Invalid milkScheme selected');
    end
    
    set(gca, 'Color', bgRGB); % Apply NaN background
    hold on;

    % --- The -Inf Underlay (White/Cream for minimal error) ---
    % Use exact white to imply "clearest/best" precision
    hInf = image(ones(rows, cols, 3)); % RGB White block
    set(hInf, 'AlphaData', isinf(logMatrix)); 

    % --- The Main Data Layer ---
    hMain = imagesc(logMatrix); 
    set(hMain, 'AlphaData', isfinite(logMatrix)); % Only show finite numbers

    cb = colorbar;
    ylabel(cb, 'log10(Error)', 'FontSize', fntSize); 
    title("Error Precision Heatmap of " + strrep(varName,"_"," "), 'FontSize', fntSize+2);
    xlabel('Column Index'); ylabel('Row Index');

    % --- FLIP Y-AXIS TO MATCH MATRIX INDEXING ---
    axis ij; 
    
    set(gca, 'FontSize', fntSize);
    
    % Add text labels with dynamic contrast
    for i = 1:rows
        for j = 1:cols
            val = logMatrix(i,j);
            if isnan(val)
                text(j, i, 'NaN', 'HorizontalAlignment', 'center', 'Color', nanLabelColor, 'FontWeight', 'bold', 'FontSize', fntSize);
            elseif isinf(val)
                text(j, i, '0.0', 'HorizontalAlignment', 'center', 'Color', infLabelColor, 'FontWeight', 'bold', 'FontSize', fntSize);
            else
                % Use contrast logic defined in the switch statement
                txtColor = 'black';
                if flipToWhite
                    if val > contrastLimit, txtColor = 'white'; end
                else
                    if val < contrastLimit, txtColor = 'white'; end
                end
                
                text(j, i, sprintf('%.1f', val), 'HorizontalAlignment', 'center', 'Color', txtColor, 'FontWeight', 'bold', 'FontSize', fntSize);
            end
        end
    end
    
    exportgraphics(fig2, fullfile(folderName, varName + "_Heatmap_" + milkScheme + ".png"), 'Resolution', 900);
    fprintf('Figures saved to /%s folder using %s scheme.\n', folderName, milkScheme);
    close(fig2); % Auto-close heatmap
end