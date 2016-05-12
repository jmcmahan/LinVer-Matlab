

N = 4000;

phi = 0.5;

Nlist = [10 20 40 100 200 400 1000 2000 4000]';
tstandard = zeros(size(Nlist));
tfast = zeros(size(Nlist));

if 1 == 0 
C1 = zeros(N);
for j = 1:N
    C1(j,:) = phi.^abs( (1:(N)) - j );
end




for j = 1:length(Nlist)
    Nj = Nlist(j);
    t0 = cputime;
    L1 = chol(C1(1:Nj,1:Nj));
    tstandard(j) = cputime - t0;

    t0 = cputime;
    L2 = fastchol_ar1(phi, Nj);
    tfast(j) = cputime - t0; 

    disp(100 * j / length(Nlist))
end

end

C2 = phi*ones(N) + (1 - phi)* eye(N);


for j = 1:length(Nlist)
    Nj = Nlist(j);
    t0 = cputime;
    L3 = chol(C2(1:Nj,1:Nj));
    tstandard(j) = cputime - t0;

    t0 = cputime;
    L4 = fastchol_eq(phi, Nj);
    tfast(j) = cputime - t0; 

    disp(100 * j / length(Nlist))
end
