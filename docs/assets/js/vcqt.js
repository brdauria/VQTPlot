// Utility functions
const kδ = bool => bool ? 1 : 0
const range = n => [...Array(n).keys()]
const make_array = (n, f) => math.map(range(n), f)
const make_diag = (n, f) => math.diag(math.map(range(n), f))
const make_matrix = (n, m, f) => math.map(range(n), i => math.map(range(m), j => f(i,j)))
const getRow = (M, i) => math.flatten(M.subset(math.index(i, math.range(0, M._size[0])))).toArray();
const getColumn = (M, i) => math.flatten(M.subset(math.index(math.range(0, M._size[0]),i))).toArray();
// apporxeq & chop functions
let approxeq = (v1, v2, epsilon = 0.001) => Math.abs(v1 - v2) <= epsilon;
let chop = (v, epsilon = 0.000001) => Math.abs(v) < epsilon ? 0 : v;
let capitalized = name => name.charAt(0).toUpperCase() + name.slice(1)
const null_space = M => math.multiply(math.ones(M._size[0]),math.inv(math.add(math.ones(M._size[0],M._size[0]), M)));
const normalize = v => math.multiply(1/math.norm(v),v);
const eigenvector = M => normalize(null_space(M.map(chop)));

function ErlangB(n,x) {
  var sum = num = den = 1
  for (i = 1; i <= n; i++) {
    num *= x
    den *= i
    sum += num/den
  }
  return  num / den / sum
}

function ErlangC(n,x) {
  return n*ErlangB(n,x)/(n-x*(1-ErlangB(n,x)))
}

function workloadMMc(c=2, λ=1, μ2=0.55) {
  const rho = λ/(c*μ2)
  const pWait = ErlangC(c,λ/μ2)
  const VQTcdf = x => 1 - pWait * math.exp(-(c * μ2 - λ)*x); 
  const VQTdf = x => pWait * math.exp(-(c * μ2 - λ)*x) * (c * μ2 - λ); 
  VQTmean = pWait / (c * μ2 - λ);

  return {"CDF": VQTcdf, "PDF": VQTdf, "Mean": VQTmean}
}

function workload(c=2, λ=1, μ1=0.45, μ2=0.55, k=0.1) {

  if (μ1 == μ2) return workloadMMc(c, λ, μ2)

  // Definitions
  const μ = (i, j) => i * μ1 + j * μ2;
  const Δ0 = make_diag(c, i => μ(i,c-1-i));
  const Δ1 = make_diag(c, i => μ(i+1,c-1-i));
  const Δ2 = make_diag(c,  i => μ(i,c-i));

  const μvec = n => make_array(n+1, i => μ(i,n-i));
  const hΔ = n => math.diag(μvec(n));
  const hB = n => make_matrix(n+2, n+1, (i,j) => kδ(i==j)*(n+1-i)*μ2+kδ(i==j+1)*i*μ1);
  const hI = n => make_matrix(n, n+1, (i,j) => kδ(i+1==j));

  const Q1 = x => math.expm(math.multiply(-x, Δ1))
  const Q2 = x => math.expm(math.multiply(-x, Δ2))
  const B1 = make_matrix(c, c, (i,j) => kδ(i==j)*(i+1)*μ1+kδ(i+1==j)*(c-1-i)*μ2);
  const B2 = make_matrix(c, c, (i,j) => kδ(i==j)*(c-i)*μ2+kδ(i-1==j)*i*μ1);

  const tΔ1 = math.multiply(math.inv(B1), Δ1, B1)
  const tΔ2 = math.multiply(math.inv(B2), Δ2, B2)
  const tQ1 = x => math.multiply(math.inv(B1), Q1(x), B1)
  const tQ2 = x => math.multiply(math.inv(B2), Q2(x), B2)

  // Solving homogeneous equations
  const D1 = x => math.add(
    math.multiply(x, x, math.identity(c)),
    math.multiply(-x, λ, math.identity(c)),
    math.multiply(x, tΔ1),
    math.multiply(λ, B1),
    math.multiply(-λ, tΔ1))
  const D2 = x => math.add(
    math.multiply(x, x, math.identity(c)),
    math.multiply(-x, λ, math.identity(c)),
    math.multiply(x, tΔ2),
    math.multiply(λ, B2),
    math.multiply(-λ, tΔ2))

  const θ1P = make_array(c, i => (λ-μ(i+1,c-1-i)+math.sqrt(math.square(λ-μ(i+1,c-1-i))+4*λ*μ(0,c-1-i)))/2)
  const θ1N = make_array(c, i => (λ-μ(i+1,c-1-i)-math.sqrt(math.square(λ-μ(i+1,c-1-i))+4*λ*μ(0,c-1-i)))/2)
  const θ2P = make_array(c, i => (λ-μ(i,c-i)+math.sqrt(math.square(λ-μ(i,c-i))+4*λ*μ(i,0)))/2)
  const θ2N = make_array(c, i => (λ-μ(i,c-i)-math.sqrt(math.square(λ-μ(i,c-i))+4*λ*μ(i,0)))/2)

  const ɸ1P = math.matrix(make_array(c, i => eigenvector(D1(θ1P[i]))));
  const ɸ1N = math.matrix(make_array(c, i => eigenvector(D1(θ1N[i]))));
  const ɸ2P = math.matrix(make_array(c, i => eigenvector(D2(θ2P[i]))));
  const ɸ2N = math.matrix(make_array(c, i => eigenvector(D2(θ2N[i]))));

  const U1P = math.multiply(math.inv(ɸ1P), math.diag(θ1P), ɸ1P);
  const U1N = math.multiply(math.inv(ɸ1N), math.diag(θ1N), ɸ1N);
  const U2P = math.multiply(math.inv(ɸ2P), math.diag(θ2P), ɸ2P);
  const U2N = math.multiply(math.inv(ɸ2N), math.diag(θ2N), ɸ2N);

  // Solving homogeneous equations
  const vη1 = eigenvector(D1(0));
  const Mη1 = math.add(
    math.multiply(λ, B1),
    math.multiply(-λ, tΔ1),
    math.multiply(-1, make_array(c, i => vη1)));
  const M0 = math.inv(Mη1);
  const vη2 = eigenvector(D2(0));
  const Mη2 = math.add(
    math.multiply(λ, B2),
    math.multiply(-λ, tΔ2),
    math.multiply(-1, make_array(c, i => vη2)));
  const M1 = math.inv(Mη2);
  const M2 = math.inv(math.add(
      math.multiply((c * μ1 + λ) * (c * μ1), math.identity(c)),
      math.multiply(-(c * μ1 + λ), tΔ2),
      math.multiply(λ, B2)));

  const invtΔ1 = math.inv(tΔ1);
  const invU2N = math.inv(U2N);
  const invU1PN = math.inv(math.subtract(U1P, U1N));
  const kU1P = math.expm(math.multiply(k, U1P));
  const kU1N = math.expm(math.multiply(k, U1N));
  const kU1PN = math.subtract(kU1P, kU1N);
  const kUU1PN = math.subtract(math.multiply(U1P, kU1P), math.multiply(U1N, kU1N));
  const tΔ12 = math.subtract(tΔ1, tΔ2);
  const BΔ12 = math.subtract(math.multiply(B1, math.inv(tΔ1), tΔ2), B2);

  ikexp = v => v.map(x => x == 0 ? k : (math.exp(k*x)-1) / x)
  const IkU1P = math.multiply(math.inv(ɸ1P), math.diag(ikexp(θ1P)), ɸ1P);
  const IkU1N = math.multiply(math.inv(ɸ1N), math.diag(ikexp(θ1N)), ɸ1N);

  const H1 = math.multiply(invU1PN, kU1PN);
  const H2 = math.multiply(M0, math.subtract(math.add(math.identity(c), math.multiply(U1N, H1)), kU1N));
  const H3 = math.add(H1, math.multiply(tΔ1, H2));
  const H4 = math.multiply(-λ, B1, H2);

  const H5 = math.multiply(invU1PN, kUU1PN);
  const H6 = math.subtract(math.multiply(M0, U1N, kU1N), math.multiply(M0, U1N, H5));
  const H7 = math.subtract(H5, math.multiply(tΔ1, H6));
  const H8 = math.multiply(λ, B1, H6);

  const H9 = math.add(math.multiply(tΔ12, M2, U2N), math.multiply(tΔ1, tΔ12, M2));
  const H10 = math.add(U2N, math.multiply(-λ, math.subtract(math.identity(c), math.multiply(B2, math.inv(tΔ2))), H9)); 
  const H11 = math.add(math.multiply(M1, U2N), math.multiply(math.inv(tΔ2), H9));
  const H12 = math.add(H10, math.multiply(λ, math.subtract(math.multiply(B1, math.inv(tΔ1), tΔ2), B2), H11));
  const H13 = math.multiply(tΔ2, H11);
  const H14 = math.multiply(λ, B1, math.inv(tΔ1), tΔ2, H11);

  const H15 = math.inv(math.add(H7, math.multiply(-1, H7, H9), math.multiply(-1, H3, H12), H13));
  const H16 = math.multiply(math.add(H14, math.multiply(H4, H12), math.multiply(-1, H8), math.multiply(H8, H9)), H15);

  const H17 = math.multiply(math.subtract(tΔ2, math.multiply(λ, H3, BΔ12)), M1);
  const H18 = math.add(math.multiply(λ, B1, math.inv(tΔ1), tΔ2, M1), math.multiply(λ, H4, BΔ12, M1));
  const H19 = math.subtract(math.identity(c), math.multiply(U2N, H15, H17));
  const H20 = math.subtract(math.multiply(H16, H17), H18);

  // Solving the final system of equations
  function hC(n) {
    if (n==0) return math.multiply(1/λ, hB(0));
    if (n<c-1) return math.multiply(hB(n), math.inv(math.add(math.multiply(λ, math.subtract(math.identity(n+1), math.multiply(hC(n-1), hI(n)))),hΔ(n))));
    if (n==c-1) return math.multiply(-1, U2N, H15, math.inv(math.add(math.multiply(λ, math.identity(c)), math.multiply(-λ, hC(c-2), hI(c-1)), math.subtract(hΔ(c-1),  H16))));
    return -1
  }

  function hH(n) {
    if (n==c-1) return hC(c-1);
    if (n<c-1) return math.multiply(hH(n+1), hC(n));
    return -1
  }

  const vη2sum = math.sum(make_array(c, n => math.sum(math.multiply(vη2, hH(n)))));
  const bcsol = 1/(math.sum(math.multiply(vη2, H19)) +math.sum(math.multiply(vη2, hH(c-1),H20)) + vη2sum)
  const δsol = n => math.multiply(bcsol, vη2, hH(n));

  // Constructing the workload function
  const vf0sol = math.subtract(math.multiply(δsol(c-1), H16), math.multiply(bcsol, vη2, U2N, H15));
  const vfksol = math.add(math.multiply(vf0sol, H7), math.multiply(δsol(c-1), H8));
  const vFksol = math.add(math.multiply(vf0sol, H3), math.multiply(δsol(c-1), H4));

  const 𝛼0sol = math.subtract(math.multiply(vf0sol, tΔ1), math.multiply(λ, δsol(c-1), B1));
  const 𝛼1sol = math.subtract(math.multiply(𝛼0sol, math.inv(tΔ1), tΔ2), math.multiply(λ, vFksol, BΔ12));
  const 𝛼2sol = math.add(math.multiply(𝛼1sol, math.inv(tΔ2)), math.multiply(-1,vfksol), math.multiply(λ, vFksol), math.multiply(-λ, vFksol, B2, math.inv(tΔ2)));

  const F2solInf = math.add(math.multiply(bcsol, vη2, H19), math.multiply(δsol(c-1), H20));
  const VQTpos = 1 - math.sum(F2solInf);

  const F1sol = x => math.add(math.multiply(math.add(vf0sol, math.multiply(𝛼0sol, M0, U1N)), invU1PN, math.subtract(math.expm(math.multiply(x,U1P)), math.expm(math.multiply(x,U1N)))), math.multiply(𝛼0sol, M0, math.subtract(math.identity(c), math.expm(math.multiply(x,U1N)))));
  const F2sol = x => math.add(
    math.multiply(vFksol, math.expm(math.multiply(x-k, U2N))),
    math.multiply(math.add(math.multiply(bcsol, vη2), math.multiply(𝛼1sol, M1)), math.subtract(math.identity(c), math.expm(math.multiply(x-k,U2N)))),
    math.multiply(-1, 𝛼2sol, tΔ12, M2, math.expm(math.multiply(x-k, U2N))),
    math.multiply(𝛼2sol, tQ1(x-k), tΔ12, M2));
  const Fsol = x => ((x < k) ? F1sol(x) : F2sol(x));
  const VQTcdf = x => VQTpos + math.sum(Fsol(x));

  const f1sol = x => math.subtract(math.multiply(math.add(vf0sol, math.multiply(𝛼0sol, M0, U1N)), invU1PN, math.subtract(math.multiply(U1P,math.expm(math.multiply(x,U1P))), math.multiply(U1N,math.expm(math.multiply(x,U1N))))), math.multiply(𝛼0sol, M0, U1N, math.expm(math.multiply(x,U1N))));
  const f2sol = x => math.subtract(
    math.multiply(math.subtract(vFksol, math.add(math.multiply(bcsol, vη2), math.multiply(𝛼1sol, M1), math.multiply(𝛼2sol, tΔ12, M2))), U2N , math.expm(math.multiply(x-k,U2N))), math.multiply(𝛼2sol, tQ1(x-k), tΔ1, tΔ12, M2));
  const fsol = x => ((x < k) ? f1sol(x) : f2sol(x));
  const VQTdf = x => math.sum(fsol(x));

  VQTmean = math.sum(math.multiply(math.add(vf0sol, math.multiply(𝛼0sol, M0, U1N)), invU1PN, math.subtract(math.multiply(k,math.expm(math.multiply(k,U1P))), IkU1P)));
  VQTmean -= math.sum(math.multiply(math.add(math.multiply(math.add(vf0sol, math.multiply(𝛼0sol, M0, U1N)), invU1PN),math.multiply(𝛼0sol, M0)),math.subtract(math.multiply(k,math.expm(math.multiply(k,U1N))), IkU1N))
    );
  VQTmean += math.sum(math.multiply(math.subtract(vFksol, math.add(math.multiply(bcsol, vη2), math.multiply(𝛼1sol, M1), math.multiply(𝛼2sol, tΔ12, M2))), math.subtract(invU2N, math.multiply(k, math.identity(c)))));
  VQTmean -= math.sum(math.multiply(𝛼2sol, math.add(invtΔ1, math.multiply(k, math.identity(c))), tΔ12, M2));

  return {"CDF": VQTcdf, "PDF": VQTdf, "Mean": VQTmean}
}