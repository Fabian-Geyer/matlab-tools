function basis = solid_harmonic_basis_S1(s, N)
%SOLID_HARMONIC_BASIS_S1  Real solid harmonic basis on the circle S^1 (R^2).
%
%   basis = SOLID_HARMONIC_BASIS_S1(s, N) returns the real solid harmonics
%   of degree 0..N on the unit circle S^1 = { s in R^2 : ||s|| = 1 },
%   expressed as polynomials in the components of the direction vector
%   s = (s1, s2), and L2-normalized on S^1. On ||s|| = 1 these are the
%   harmonic homogeneous polynomials that form an orthonormal basis of the
%   harmonic space Y_n^2; they remove the rank deficiency of the raw
%   monomial basis on the circle and are therefore suitable as a polynomial
%   basis for SOS programming.
%
%   This is the d = 2 member of the solid_harmonic_basis_* family. Compared
%   with the S^2 version, the "vertical" Legendre factor P_n^(m)(cos theta)
%   is absent: for d = 2 it collapses (P_{n,2}(cos theta) = cos(n theta)),
%   so the basis is exactly the horizontal / Fourier part. See A&H Sect. 2.2.
%
%   INPUTS
%       s : 2x1 casos.PD vector, the (unit) direction indeterminates
%           s = [s1; s2]. Interpreted via the standard convention
%               s1 = cos(alpha),
%               s2 = sin(alpha),
%           i.e. alpha is the polar angle on the circle.
%       N : nonnegative integer scalar, the maximum polynomial degree.
%
%   OUTPUT
%       basis : (2N+1) x 1 casos.PD vector of basis polynomials, ordered by
%               degree n = 0..N: the constant (n = 0), then for each n >= 1
%               the cosine term followed by the sine term.
%
%   METHOD (equation numbers refer to Atkinson & Han, 2012, "Spherical
%   Harmonics and Approximations on the Unit Sphere")
%       1. Horizontal (Fourier) parts, the harmonic homogeneous polynomials
%          in R^2 restricted to S^1:
%              A_m = Re[(s1 + i s2)^m] = cos(m alpha),
%              B_m = Im[(s1 + i s2)^m] = sin(m alpha).
%          These span Y_m^2 (A&H Sect. 2.2); (s1 + i s2)^m is homogeneous
%          and harmonic (cf. the (zeta . x)^n argument in Sect. 2.11).
%       2. Assembly with the Sect. 2.2 orthonormalization constants
%              c_0 = 1/sqrt(2*pi)   (constant term),
%              c_n = 1/sqrt(pi)     (cosine and sine terms, n >= 1).
%
%   NOTE ON NORMALIZATION
%       The factors c_n make the functions orthonormal in L2(S^1). This is a
%       pure change of basis: a subsequent least-squares / SOS fit returns
%       coefficients scaled by 1/c_n, and the reconstructed function is
%       unchanged. Its purpose is numerical conditioning, not correctness.
%
%   Reference:
%       K. Atkinson and W. Han, Spherical Harmonics and Approximations on
%       the Unit Sphere: An Introduction, Lecture Notes in Mathematics,
%       Springer, 2012. (Sect. 2.2; Example 2.48 for the d = 2 basis.)

    arguments
        s (2,1) casos.PD          % must be a 2x1 vector of the toolbox type
        N (1,1) double {mustBeInteger, mustBeNonnegative}
    end

    % ---------------------------------------------------------------------
    % Step 1: horizontal (Fourier) parts                 [A&H Sect. 2.2]
    %
    %   A_m = Re[(s1 + i s2)^m] = cos(m alpha)
    %   B_m = Im[(s1 + i s2)^m] = sin(m alpha)
    %
    % Built by the complex-multiply recurrence (one factor of (s1 + i s2)
    % per step), starting from A_0 = 1, B_0 = 0. Index convention:
    % A(m+1) = A_m, B(m+1) = B_m.
    %
    % For d = 2 there is no Legendre factor to multiply in: P_{n,2} = cos(n
    % theta) is already carried by A_n itself, so these ARE the basis
    % functions (up to normalization).
    % ---------------------------------------------------------------------
    A = casos.PD.zeros(N+1,1);  B = casos.PD.zeros(N+1,1);
    A(1) = 1;  B(1) = 0;          % m = 0
    for m = 1:N
        A(m+1) = s(1)*A(m) - s(2)*B(m);
        B(m+1) = s(1)*B(m) + s(2)*A(m);
    end

    % ---------------------------------------------------------------------
    % Step 2: assemble and normalize                     [A&H Sect. 2.2]
    %
    %   n = 0 : c_0 * A_0            with c_0 = 1/sqrt(2*pi)
    %   n >= 1: c_n * A_n (cosine),  c_n * B_n (sine),  c_n = 1/sqrt(pi)
    %
    % Degree n contributes 1 function for n = 0 and 2 functions for n >= 1,
    % giving 2N+1 basis functions in total (= dim of harmonics up to N).
    % ---------------------------------------------------------------------
    basis = {};
    basis{end+1} = (1/sqrt(2*pi)) * A(1);               % n = 0 (constant)
    for n = 1:N
        cn = 1/sqrt(pi);                                % c_n, n >= 1
        basis{end+1} = cn * A(n+1);                     %#ok<AGROW> cosine
        basis{end+1} = cn * B(n+1);                     %#ok<AGROW> sine
    end
    basis = vertcat(basis{:});    % (2N+1) x 1 vector of polynomials in s
end