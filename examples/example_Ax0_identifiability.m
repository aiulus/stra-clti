function example_Ax0_identifiability()
%PLOT_EXAMPLE3_1_LAYOUT
% Five-plot (1,2,2) visualization for A and Atilde using a real Jordan/modal basis.
% Row 1: invariant subspaces only (L1 line, L2 plane)
% Row 2: projections of x0^(a) [left] and x0^(b) [right] onto L1 and L2
% Row 3: trajectories from x0^(a) [left] and x0^(b) [right] under the same system
%
% Saves PNGs in ./figs and displays figures.
%
% Dependencies on path: JordanReal.m, rsolve.m

% ------------------------------ Problem data ------------------------------
A = [0 1 -1;
     2 0  0;
     3 1  0];

Atilde = [1 -1  0;
          0  4 -2;
          2  3 -1];

if ~exist('figs','dir'), mkdir figs; end

% ------------------------------ Real modal bases via JordanReal ----------
jrA   = realJordan(A);
Q_A   = jrA.Qmat;
assert(jrA.K1>=1 && jrA.K2>=1, 'Expected 1 real eigenvalue + 1 complex pair for A.');
q1_A  = Q_A(:,1);                         % L1 direction for A
U2_A  = Q_A(:, jrA.K1+1 : jrA.K1+2);      % L2 plane for A (spanned by Re/Im columns)

jrAt  = realJordan(Atilde);
Q_At  = jrAt.Qmat;
assert(jrAt.K1>=1 && jrAt.K2>=1, 'Expected 1 real eigenvalue + 1 complex pair for Atilde.');
q1_At = Q_At(:,1);                        % L1 direction for Atilde
U2_At = Q_At(:, jrAt.K1+1 : jrAt.K1+2);   % L2 plane for Atilde

% ------------------------------ Initial states from modal coords of A ----
% (As in the paper/snippet: z0 defined in A's modal coordinates; reused for both systems.)
z0_a = [ 2; -1; 0];
z0_b = [ 0; -2; 3];
x0_a = Q_A * z0_a;
x0_b = Q_A * z0_b;

% ------------------------------ Build both figures -----------------------
makeFigure_Ax0(A, q1_A,  U2_A,  x0_a, x0_b, 'A', fullfile('figs','A_layout.png'));
makeFigure_Ax0(Atilde, q1_At, U2_At, x0_a, x0_b, 'A~', fullfile('figs','Atilde_layout.png'));

end











