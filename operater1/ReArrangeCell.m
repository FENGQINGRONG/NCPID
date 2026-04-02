function gamma = ReArrangeCell(alpha,D)

nD = length(D);
Level = length(alpha);

for ki = 1 : Level
    for ji = 1 : nD-1
        for jj = 1 : nD-1
            gamma{ki}{ji,jj} = alpha{ki};
        end
    end
end