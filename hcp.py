from pysat.formula import CNF
from pysat.solvers import Glucose4

def create_hamiltonian_cycle_formula(n, edges):
    formula = CNF()
    variables = {}

    def var(i, j):
        if (i, j) not in variables:
            variables[(i, j)] = len(variables) + 1
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
        formula.append([var(1, i), -pos_var(i, 0)])  # P_i = 2 if H_1i
        for bit in range(1, m):
            formula.append([var(1, i), -pos_var(i, bit)])

        formula.append([var(i, 1), -pos_var(i, 0)])  # P_i = n if H_i1
        for bit in range(1, m):
            formula.append([var(i, 1), -pos_var(i, bit)])

    for i in range(2, n+1):
        for j in range(2, n+1):
            if i != j:
                for bit in range(m):
                    formula.append([-var(i, j), -pos_var(i, bit), pos_var(j, bit)])

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
