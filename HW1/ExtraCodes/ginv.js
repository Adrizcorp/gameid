// Solve by Least Square method with Generalized Inverse Matrix

// numeric.js for general matrix computation (npm install numeric)
// see: http://numericjs.com/ 
var numeric = require("numeric");

// mechanism of moore-penrose psuedoinverse
// (Eigenvalues of tM*M should be positive real numbers)
// Some properties of ginv
// ginv(ginv(A)) = A
// A * ginv(A) * A = A
var ginv0 = function (M) {
    // compute ginv from t(M) * M square matrix
    // A: m rows * n cols matrix, tA: n*m matrix, tAA: n*n matrix
    // ginv(A): n*m matrix
    var A = (M.constructor == numeric.T) ? M : numeric.t(M);
    var tA = A.transpose();
    var tAA = tA.dot(A);
    var etAA = numeric.eig(tAA.x);
    if (etAA.lambda.y || etAA.lambda.x.some(function (v) {return v <= 0;}))
        return null; // would not be occurred

    // Singlar values of A: n values
    var S = etAA.lambda.x.map(function (v) {
        return Math.sqrt(v);
    });
    // V(Eigen vectors of tAA): n items of n-vector
    var V = etAA.E.transpose().x;
    // Unit vectors: n items of m-vector
    var U = S.map(function (s, i) {
        return A.dot(V[i]).mul(1/s).x;
    });
    var zero = tA.sub(tA);
    
    var ginv = S.reduce(function (r, s, i) {
        var tViUi = numeric.dot(numeric.transpose([V[i]]), [U[i]]);
        return r.add(numeric.mul(1/s, tViUi));
    }, zero).x;
    return ginv;
};

var ginv = function (M) {
    // ginv with using SVD (singlar value decomposition) function
    var svd = numeric.svd(M);
    // A = U * diag(S) * t(V)
    //console.log(numeric.dot(svd.U, numeric.dot(
    //    numeric.diag(svd.S), svd.V))); // ==> A
    // ginv(A) = V * diag(1/S) * t(U)
    var diS = numeric.diag(svd.S.map(function (s) {return 1/s;}));
    var ginv = numeric.dot(svd.V, numeric.dot(diS, numeric.transpose(svd.U)));
    return ginv;
};


// [main] Solving X as least square method with ginv
{
    console.log("<<", process.title, process.version, ">>");
    console.log("<< numeric", numeric.version, ">>");
    // example: find most fit value of x and y from 3 relationships:
    // 2x - y = 4
    // -x + 2y = 9
    // -x - y = 16
    var A = [
        [2, -1],
        [-1, 2],
        [-1, -1],
    ];
    var b = [4, 9, 16];
    console.log("\n[matrix A of A*x = b]");
    console.log(A);
    console.log("\n[vector b of A*x = b]");
    console.log(b);
    
    var ginvA0 = ginv0(A);
    console.log("\n[ginv(A) made from t(A)*A eigenvalues]");
    console.log(ginvA0);
    var ginvA = ginv(A);
    console.log("\n[ginv(A) made from SVD function]");
    console.log(ginvA);
    var x = numeric.dot(ginvA, b); // (almost same as ginvA0)
    console.log("\n[result x from ginv(A)*b]");
    console.log(x);
    console.log("\n[lengthes as |A[i]*x-b|]")
    console.log( A.map(function (a, i) {
        return Math.abs(numeric.dot(a, x) - b[i]);
    })); // least average value of |A[i]*x - b[i]|
    
    // [Summary of the Solution of x = ginv(A)*b]
    // [System of Equations]
    //   A*x = b
    // [Singlar Value Decomposition by svd(A)]
    //   A = U*diag(S)*t(V)
    //       (A: m*n matrix, U: m*n mat, diag(S): n*n mat, V: n*n mat)
    //   t(U)*U = t(V)*V = V*t(V) = I
    //   diag(1/S)*diag(S) = diag(S)*diag(1/S) = I
    // [Generalized Inverse Matrix]
    //   ginv(A) = V*diag(1/S)*t(U)
    //   ginv(A)*A  = {V*diag(1/S)*t(U)}*{U*diag(S)*t(V)} = I
    // [Least Square Method]
    //   ginv(A)*A*x = ginv(A)*b
    //   x = ginv(A)*b
    
    // [Decomposing ginv matrix to vector list]
    //   x = V*diag(1/S)*t(U)*b
    //   x = sum(V[i] * 1/S[i] * t(U[i]) * b)
    //     = sum(1/S[i] * V[i] * (U[i].b))
    //     = sum((U[i].b)/S[i] * V[i])
    var svd = numeric.svd(A);
    var tV = numeric.transpose(svd.V); // as vector list
    var tU = numeric.transpose(svd.U); // as vector list
    console.log("\n[A from U * diag(S) * t(V)]");
    console.log(numeric.dot(numeric.dot(svd.U, numeric.diag(svd.S)), tV));
    
    var x2 = tV.reduce(function (r, vi, i) {
        var e = numeric.mul(numeric.dot(tU[i], b) / svd.S[i], vi);
        return numeric.add(r, e);
    }, numeric.sub(svd.S, svd.S));
    console.log("\n[result x2 by non matrix computation]")
    console.log(x2); // = ginv(A)*b
    console.log("\n[lengthes as |A[i]*x2-b|]")
    console.log(A.map(function (a, i) {
        return Math.abs(numeric.dot(a, x2) - b[i]);
    }));
}
