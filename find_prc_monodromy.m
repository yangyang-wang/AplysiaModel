function M = find_prc_monodromy(model)
init_matrix=eye(6);
M = [];

for i = 1:6
    model.find_prc(init_matrix(i,:));
    M = [M; model.prc(end,:)];
end
end