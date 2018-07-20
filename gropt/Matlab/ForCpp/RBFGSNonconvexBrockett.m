function PlotSoftICA(N)
    num = 1000;
    lr = 4;
    lb = 4;

    RBFGSARMIJO = load(['RBFGSARMIJO']);
    RBFGSWOLFE = load(['RBFGSWOLFE']);
    RBFGSARMIJO = reshape(RBFGSARMIJO, num, lr, 8);
    RBFGSWOLFE = reshape(RBFGSWOLFE, num, lb, 8);

    tab = zeros(5, 8);

    idx = [1, 2, 3, 5, 8];

    for i = 1 : 5
        if(i ~= 4)
            ave = sum(RBFGSARMIJO(:, :, idx(i))) / num;
        else
            ave = (sum(RBFGSARMIJO(:, :, idx(i))) + sum(RBFGSARMIJO(:, :, 6))) / num;
        end
        tab(i, 1) = ave(1);
        tab(i, 2) = ave(2);
        tab(i, 3) = ave(3);
        tab(i, 4) = ave(4);


        if(i ~= 4)
            ave = sum(RBFGSWOLFE(:, :, idx(i))) / num;
        else
            ave = (sum(RBFGSWOLFE(:, :, idx(i))) + sum(RBFGSWOLFE(:, :, 6))) / num;
        end
        tab(i, 5) = ave(1);
        tab(i, 6) = ave(2);
        tab(i, 7) = ave(3);
        tab(i, 8) = ave(4);
    end

    fout = fopen(['RBFGSNonconvex_ARMIJO_WOLFE_tab.txt'],'w');
    for i = 1 : 5
        for j = 1 : 8

            tabstr = outputfloat(tab(i, j));
            if(j < 8)
                fprintf(fout,['$' tabstr '$ & ']);
            else
                fprintf(fout,['$' tabstr '$']);
            end
        end
        fprintf(fout,' \\\\\n');
    end
    fclose(fout);
end


function str = outputfloat(x)
    if(x <= 0)
        str = '';
        return;
    end
    p = log(x)/log(10);
    p = - ceil(-p);
    x = round(x * 10^(-p) * 100);
    x = x / 100;
    strx = sprintf('%3.2f', x);
    if(p ~= 0)
        str = [strx '_{' num2str(p) '}'];
    else
        str = strx;
    end
end
