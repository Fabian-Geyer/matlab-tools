function basis = solid_harmonic_basis_S2(s, N)
%SOLID_HARMONIC_BASIS  Real solid spherical harmonic basis in Cartesian form.
%
%   basis = SOLID_HARMONIC_BASIS(s, N) returns the real solid spherical
%   harmonics of degree 0..N, expressed as polynomials in the components of
%   the direction vector s = (s1, s2, s3), and L2-normalized on the unit
%   sphere. On ||s|| = 1 these are the harmonic homogeneous polynomials that
%   form an orthonormal basis of the spherical harmonic space Y_n^3; they
%   remove the rank deficiency of the raw monomial basis on the sphere and
%   are therefore suitable as a polynomial basis for SOS programming.
%
%   INPUTS
%       s : 3x1 casos.PD vector, the (unit) flow-direction indeterminates
%           s = [s1; s2; s3]. Interpreted via the standard convention
%               s1 = sin(theta) cos(phi),
%               s2 = sin(theta) sin(phi),
%               s3 = cos(theta),
%           so that t := s3 = cos(theta).
%       N : nonnegative integer scalar, the maximum polynomial degree.
%
%   OUTPUT
%       basis : (N+1)^2 x 1 casos.PD vector of basis polynomials, ordered by
%               degree n = 0..N and within each degree by order
%               m = 0 (cosine), then (m, cos), (m, sin) for m = 1..n.
%
%   METHOD (all equation numbers refer to Atkinson & Han, 2012,
%   "Spherical Harmonics and Approximations on the Unit Sphere")
%       1. Legendre polynomials P_n(t) via the three-term recurrence,
%          Eq. (2.86) with d = 3.
%       2. Their m-th derivatives P_n^(m)(t) = d^m/dt^m P_n(t), obtained by
%          exact differentiation (the associated Legendre functions of
%          Prop. 2.43 without the (1-t^2)^{m/2} factor).
%       3. Assembly into the orthonormal 3-D basis of Example 2.48, using
%              (sin theta)^m cos(m phi) = Re[(s1 + i s2)^m],
%              (sin theta)^m sin(m phi) = Im[(s1 + i s2)^m],
%          and the normalization constant c_{n,m} of Example 2.48.
%
%   NOTE ON NORMALIZATION
%       The factor c_{n,m} makes the functions orthonormal in L2(S^2). It is
%       a pure change of basis: a subsequent least-squares / SOS fit returns
%       coefficients scaled by 1/c_{n,m}, and the reconstructed function is
%       unchanged. Its purpose is numerical conditioning, not correctness.
%
%   Reference:
%       K. Atkinson and W. Han, Spherical Harmonics and Approximations on
%       the Unit Sphere: An Introduction, Lecture Notes in Mathematics,
%       Springer, 2012. (Eq. (2.86); Prop. 2.43; Example 2.48.)

    arguments
        s (3,1) casos.PD          % must be a 3x1 vector of the toolbox type
        N (1,1) double {mustBeInteger, mustBeNonnegative}
    end

    t = s(3);                     % t = cos(theta) = s3

    % ---------------------------------------------------------------------
    % Step 1: Legendre polynomials P_n(t)              [A&H Eq. (2.86), d=3]
    %
    %   P_n(t) = (2n-1)/n * t * P_{n-1}(t) - (n-1)/n * P_{n-2}(t),
    %   P_0(t) = 1,  P_1(t) = t.
    %
    % Storage convention: P(k) holds degree k-1, i.e. P(1)=P_0, P(2)=P_1, ...
    % so the array index is (degree + 1).
    % ---------------------------------------------------------------------
    P = casos.PD.zeros(N+1,1);
    P(1) = 1;                     % P_0
    if N >= 1
        P(2) = t;                 % P_1
    end
    for d = 2:N                   % d = degree, exactly the n in Eq. (2.86)
        P(d+1) = (2*d-1)/d * t * P(d) - (d-1)/d * P(d-1);
    end

    % ---------------------------------------------------------------------
    % Step 2: derivatives P_n^(m)(t) = d^m/dt^m P_n(t)   [A&H Prop. 2.43,
    %                                                     differentiation]
    %
    % Column 1 holds the 0th derivative (= P itself); each further column is
    % the exact derivative of the previous one. CasADi/casos differentiates
    % symbolically, so nabla(.,t) is exact.
    %
    % Indexing after this block: P_der(k, m+1) = P_{k-1}^(m)(t), i.e.
    %   P_der(n+1, m+1) = P_n^(m)(t),   with derivative order m = 0..N.
    % ---------------------------------------------------------------------
    P_der = P;                    % order-0 column
    for m = 1:N
        P_der = [P_der, nabla(P_der(:,m), t)]; %#ok<AGROW>
    end

    % ---------------------------------------------------------------------
    % Step 3: horizontal (azimuthal) parts             [A&H Example 2.48]
    %
    %   A_m = Re[(s1 + i s2)^m] = (sin theta)^m cos(m phi)
    %   B_m = Im[(s1 + i s2)^m] = (sin theta)^m sin(m phi)
    %
    % Built by the complex-multiply recurrence (one factor of (s1 + i s2)
    % per step), starting from A_0 = 1, B_0 = 0. Same index convention:
    % A(m+1) = A_m, B(m+1) = B_m.
    % ---------------------------------------------------------------------
    A = casos.PD.zeros(N+1,1);  B = casos.PD.zeros(N+1,1);
    A(1) = 1;  B(1) = 0;          % m = 0
    for m = 1:N
        A(m+1) = s(1)*A(m) - s(2)*B(m);
        B(m+1) = s(1)*B(m) + s(2)*A(m);
    end

    % ---------------------------------------------------------------------
    % Step 4: assemble and normalize                    [A&H Example 2.48]
    %
    % Each basis function is
    %       c_{n,m} * P_n^(m)(t) * { A_m  (cosine) | B_m  (sine) },
    % with the Example 2.48 normalization constant
    %       c_{n,m} = sqrt( (2n+1) (n-m)! / (2*pi*(n+m)!) ).
    % Order: for each degree n, the m=0 term (cosine only), then for
    % m = 1..n the cosine term followed by the sine term.
    % ---------------------------------------------------------------------
    basis = {};
    for n = 0:N
        cn0 = sqrt( (2*n+1) / (2*pi) );                 % c_{n,0}
        basis{end+1} = cn0 * P_der(n+1,1) * A(1);       %#ok<AGROW> m = 0
        for m = 1:n
            cnm = sqrt( (2*n+1)*factorial(n-m) ...
                        / (2*pi*factorial(n+m)) );      % c_{n,m}
            basis{end+1} = cnm * P_der(n+1,m+1) * A(m+1); %#ok<AGROW> cos
            basis{end+1} = cnm * P_der(n+1,m+1) * B(m+1); %#ok<AGROW> sin
        end
    end
    basis = vertcat(basis{:});    % (N+1)^2 x 1 vector of polynomials in s
end