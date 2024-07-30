from pysat.formula import CNF
from pysat.solvers import Glucose4

# Global variable counter
var_counter = 0

def new_var():
    global var_counter
    var_counter += 1
    return var_counter

def create_hamiltonian_cycle_formula(n, edges):
    global var_counter
    formula = CNF()
    variables = {}

    def var(i, j):
        global var_counter
        if (i, j) not in variables:
            variables[(i, j)] = len(variables) + 1
        var_counter = max(var_counter, variables[(i, j)])
        return variables[(i, j)]

    # Add constraints for exactly one outgoing arc per vertex
    for i in range(1, n+1):
        formula.append([var(i, j) for j in range(1, n+1) if i != j and (i, j) in edges])
        for j in range(1, n+1):
            if i != j and (i, j) in edges:
                for k in range(j+1, n+1):
                    if i != k and (i, k) in edges:
                        formula.append([-var(i, j), -var(i, k)])

    # Add constraints for exactly one incoming arc per vertex
    for j in range(1, n+1):
        formula.append([var(i, j) for i in range(1, n+1) if i != j and (i, j) in edges])
        for i in range(1, n+1):
            if i != j and (i, j) in edges:
                for k in range(i+1, n+1):
                    if k != j and (k, j) in edges:
                        formula.append([-var(i, j), -var(k, j)])

    # Binary adder encoding for successor function
    # Using binary representation of positions
    m = (n-1).bit_length()  # Number of bits to represent positions

    position_vars = {}
    def pos_var(i, bit):
        if (i, bit) not in position_vars:
            position_vars[(i, bit)] = len(variables) + 1 + len(position_vars)
        return position_vars[(i, bit)]

    # Constraints for binary adder encoding
    for i in range(2, n+1):
        # P_i = 2 if H_1i
        formula.append([-var(1, i), pos_var(i, 1)])  # The 1st bit of P_i is 1 if H_1i
        formula.append([-var(1, i), -pos_var(i, 0)])  # The 0th bit of P_i is 0 if H_1i
        for bit in range(2, m):
            formula.append([-var(1, i), -pos_var(i, bit)])  # The other bits of P_i are 0 if H_1i

        # P_i = n if H_i1
        for bit in range(m):
            if (n >> bit) & 1:
                formula.append([-var(i, 1), pos_var(i, bit)])  # The bit-th bit of P_i is 1 if H_i1 and the bit-th bit of n is 1
            else:
                formula.append([-var(i, 1), -pos_var(i, bit)])  # The bit-th bit of P_i is 0 if H_i1 and the bit-th bit of n is 0

    for i in range(2, n+1):
        for j in range(2, n+1):
            if i != j:
                # P_j = P_i + 1 if H_ij
                carry = [var(i, j)]  # Initialize carry with H_ij
                for bit in range(m):
                    # Create new variables for the sum and carry
                    sum_var = new_var()
                    carry_var = new_var()

                    # Constraints for the binary adder
                    formula.append([-carry[bit], -pos_var(i, bit), sum_var])
                    formula.append([-carry[bit], pos_var(i, bit), -sum_var])
                    formula.append([carry[bit], -pos_var(i, bit), -sum_var])
                    formula.append([carry[bit], pos_var(i, bit), sum_var])

                    if bit < m - 1:
                        formula.append([-carry[bit], -pos_var(i, bit), -pos_var(i, bit+1), carry_var])
                        formula.append([-carry[bit], -pos_var(i, bit), pos_var(i, bit+1), -carry_var])
                        formula.append([-carry[bit], pos_var(i, bit), -pos_var(i, bit+1), -carry_var])
                        formula.append([-carry[bit], pos_var(i, bit), pos_var(i, bit+1), carry_var])
                        formula.append([carry[bit], -pos_var(i, bit), -pos_var(i, bit+1), -carry_var])
                        formula.append([carry[bit], -pos_var(i, bit), pos_var(i, bit+1), carry_var])
                        formula.append([carry[bit], pos_var(i, bit), -pos_var(i, bit+1), carry_var])
                        formula.append([carry[bit], pos_var(i, bit), pos_var(i, bit+1), -carry_var])

                    # The sum must be equal to P_j
                    formula.append([-sum_var, pos_var(j, bit)])
                    formula.append([sum_var, -pos_var(j, bit)])

                    # Update carry
                    if bit < m - 1:
                        carry.append(carry_var)

    return formula, variables

def read_graph_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Parse the first line to get the number of vertices and edges
    _, _, n, _ = lines[0].split()
    n = int(n)

    # Parse the remaining lines to get the edges
    edges = []
    for line in lines[1:]:
        _, i, j = line.split()
        edges.append((int(i), int(j)))

    return n, edges

def print_hamiltonian_cycle(model, variables):
    # Create a dictionary mapping variables to their values
    var_values = {var: (value > 0) for var, value in enumerate(model, start=1)}

    # Create a dictionary mapping vertices to their successors in the cycle
    successors = {i: j for (i, j), var in variables.items() if var_values[var]}

    # Print the cycle starting from vertex 1
    i = 1
    while True:
        print(f"({i}, {successors[i]})", end=' ')
        i = successors[i]
        if i == 1:
            break
    print()

# Example usage
filename = './input/hc-4.col'
n, edges = read_graph_from_file(filename)

formula, variables = create_hamiltonian_cycle_formula(n, edges)

solver = Glucose4()
solver.append_formula(formula)

if solver.solve():
    model = solver.get_model()
    print("Hamiltonian Cycle exists!")
    print_hamiltonian_cycle(model, variables)
else:
    print("No Hamiltonian Cycle exists.")
