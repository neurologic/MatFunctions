function Ax = RidgeRegression(A,x,epsilon)

Ax = ((A' * A) + epsilon^2 * eye(size(A))) \ (A' * x);