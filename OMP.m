function [x] = OMP(y, A, maxIteration)

Eps = 1e-10;

[numMeasurement, lengthX] = size(A);  % A: sensing matrix

indexCandidate = [];

ATest = A;

residual = y;

iTest = 1;

while ((norm(residual) > Eps) && (iTest < maxIteration))
    correlation = ATest'*residual;

    [~, currentCandidate] = max(abs(correlation));

    indexCandidate = [indexCandidate currentCandidate];

    ACandidate = A(:,indexCandidate);

    xCandidate = ACandidate\y;

    residual = y - ACandidate*xCandidate;

    ATest(:,currentCandidate) = zeros(numMeasurement,1);

    iTest = iTest + 1;
end

x = zeros(lengthX, 1);
x(indexCandidate) = xCandidate;

end