function [Hermit, mono, Lagur] = LSMC(N, T, K, S0, r, sigma, num, k)
    rng(12345);
    paths = zeros(num, N+1);
    rnorm = randn(num, N);
    CF = zeros(num, N);
    CF2 = zeros(num, N);
    CF3 = zeros(num, N);
    dt = T/N;
    paths(:, 1) = S0;
    for i = 1:N
        paths(:,i+1) = paths(:, i).*exp( (r-sigma^2/2)*dt +  sigma* rnorm(:, i) * sqrt(dt));
    end
    CF(:, N) = max([K - paths(:, N+1) zeros(num,1)],[],2);
    CF2(:, N) = max([K - paths(:, N+1) zeros(num,1)],[],2);
    CF3(:, N) = max([K - paths(:, N+1) zeros(num,1)],[],2);
    for i = (N-1):-1:1
        temp = (paths(:,i+1) < K);
        x = paths(:, i+1).* temp;
        y = zeros(num, 1);
        y2 = zeros(num, 1);
        y3 = zeros(num, 1);
        for j = 1: num
            [M,I] = max(CF(j, i+1:N));
            [M2,I2] = max(CF2(j, i+1:N));
            [M3,I3] = max(CF3(j, i+1:N));
            if M > 0
                y(j) = CF(j, i+I) * exp(-r*dt*I);
            end
            if M2 > 0
                y2(j) = CF2(j, i+I2) * exp(-r*dt*I2);
            end
            if M3 > 0
                y3(j) = CF3(j, i+I3) * exp(-r*dt*I3);
            end
        end
        y = y.* temp;
        y2 = y2.* temp;
        y3 = y3.* temp;
        if k == 2
            X = [ones(length(x),1).* temp x];
            X2 = [ones(length(x),1).* temp 2*x];
            X3 = [exp(-x/2).* temp exp(-x/2).*(1-x).* temp];
        elseif k == 3
            X = [ones(length(x),1).* temp x x.^2];
            X2 = [ones(length(x),1).* temp 2*x 4*x.^2 - 2];
            X3 = [exp(-x/2).* temp exp(-x/2).*(1-x).* temp exp(-x/2).*(1-2*x+x.^2/2).* temp];
        elseif k == 4
            X = [ones(length(x),1).* temp x x.^2 x.^3];
            X2 = [ones(length(x),1).* temp 2*x 4*x.^2 - 2 8*x.^3 - 12*x];
            X3 = [exp(-x/2).* temp exp(-x/2).*(1-x).* temp exp(-x/2).*(1-2*x+x.^2/2).* temp exp(-x/2).*(1-3*x+3*x.^2/2 - x.^3/6).* temp];
        end
        b = X\y;
        b2 = X2\y2;
        b3 = X3\y3;
        continuoation = X*b;
        continuoation2 = X2*b2;
        continuoation3 = X3*b3;
        newCF = max([K - paths(:,i+1) zeros(num, 1)], [], 2);
        ngc = (continuoation < newCF);
        ngc2 = (continuoation2 < newCF);
        ngc3 = (continuoation3 < newCF);
        for j = 1: length(ngc)
            if ngc(j) == 1
                CF(j, : ) = zeros(1, N);
                CF(j, i) = newCF(j);
            end    
            if ngc2(j) == 1
                CF2(j, : ) = zeros(1, N);
                CF2(j, i) = newCF(j);
            end  
            if ngc3(j) == 1
                CF3(j, : ) = zeros(1, N);
                CF3(j, i) = newCF(j);
            end  
        end
    end
    
    B = transpose(exp(-r*dt:-r*dt:-N*r*dt));
    mono = mean(CF * B);
    Hermit = mean(CF2 * B);
    Lagur = mean(CF3 * B);

end