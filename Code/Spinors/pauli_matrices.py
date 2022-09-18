import itertools
from sympy import I, Matrix, Symbol, init_printing, latex, conjugate, trace
from itertools import permutations


# Pauli matrices
SIGMA_1 = Matrix([[0, 1], [1, 0]])
SIGMA_2 = Matrix([[0, -I], [I, 0]])
SIGMA_3 = Matrix([[1, 0], [0, -1]])
PAULI_MATRICES = [SIGMA_1, SIGMA_2, SIGMA_3]

# General SU(2) matrix
a = Symbol('a', real=False)
b = Symbol('b', real=False)
U = Matrix([[a, b], [-conjugate(b), conjugate(a)]])
U_dagger = Matrix([[conjugate(a), -b], [conjugate(b), a]])


# Adjoint action of an arbitrary SU(2) matrix on the input matrix 
def adjoint_action(matrix: Matrix) -> Matrix:
    return U * matrix * U_dagger


# Use to get all permutation of U \sigma_i U^\dagger \sigma^j
def get_all_permutations():

    adjoints = [adjoint_action(pauli_matrix) for pauli_matrix in PAULI_MATRICES]

    index_permutations = itertools.permutations(range(3), 2)

    return {f"{i+1}{j+1}": adjoints[i] * PAULI_MATRICES[j] for i, j in index_permutations}


# Construct the SO(3) matrix defined by \frac12 \mathrm{tr}(U \sigma_i U^\dagger \sigma^j)
def get_so3_matrix():

    rows = []
    for i in range(3):
        row = []
        for j in range(3):
            row.append(trace(adjoint_action(PAULI_MATRICES[i]) * PAULI_MATRICES[j])/2)
        rows.append(row)

    return Matrix(rows)


def main():

    perms = get_all_permutations()

    traced = {perm: trace(matrix)/2 for perm, matrix in perms.items()}

    so3_matrix = get_so3_matrix()

    so3_latex = latex(so3_matrix)
    
    determinant = so3_matrix.det()

    determinant_latex = latex(determinant)

    pass


if __name__ == "__main__":
    main()



