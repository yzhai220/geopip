function makepca(file)
    d = importdata (file);
    data = d.data;
    names = d.textdata;
    [u,s,v,f] = pca (data, 2);
    scatter (u(:,1), u(:,2), 40, 'filled');
    for i=1:size(u,1)
        text (u(i,1), u(i,2), names{i});
    end

    fprintf (2, 'Variance explained = %f\n', f);
return