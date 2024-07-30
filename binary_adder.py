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
    position_vars = {}

    # Hij
    def var(i, j):
        global var_counter
        if (i, j) not in variables:
            var_counter += 1
            variables[(i, j)] = var_counter
        return variables[(i, j)]

    # Pi (successor function)
    def pos_var(i, bit):
        global var_counter
        if (i, bit) not in position_vars:
            var_counter += 1
            position_vars[(i, bit)] = var_counter
        return position_vars[(i, bit)]

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
    m = n.bit_length()  # Number of bits to represent positions

    # Initial position constraint P1 = 1
    for bit in range(m):
        if bit == 0:
            formula.append([pos_var(1, bit)])
        else:
            formula.append([-pos_var(1, bit)])

    # H1i ⇒ Pi = 2 (00..010)
    for i in range(2, n+1):
        for bit in range(m):
            if bit == 1:
                formula.append([-var(1, i), pos_var(i, bit)])
            else:
                formula.append([-var(1, i), -pos_var(i, bit)])

    # Hi1 ⇒ Pi = n (10..00)
    for i in range(2, n+1):
        for bit in range(m):
            if bit == m-1:  # Set the most significant bit to 1
                formula.append([-var(i, 1), pos_var(i, bit)])
            else:  # Set all other bits to 0
                formula.append([-var(i, 1), -pos_var(i, bit)])
        
    for i in range(2, n+1):
        for j in range(2, n+1):
            if i != j and (i, j) in edges:
                for bit in range(m):
                    if bit == 0:
                        # For bit position 0: Y0 = ¬X0
                        formula.append([-var(i, j), pos_var(j, bit), pos_var(i, bit)])
                        formula.append([-var(i, j), -pos_var(j, bit), -pos_var(i, bit)])
                    elif bit == 1:
                        # For bit position 1: X0 ⇒ (Y1 = ¬X1) and ¬X0 ⇒ (Y1 = X1)
                        
                        # X0 ⇒ (Y1 = ¬X1)
                        # ¬X0 v ((Y1 v X1) ^ (¬Y1 v ¬X1))
                        # (¬X0 v (Y1 v X1)) ^ (¬X0 v (¬Y1 v ¬X1))
                        formula.append([-var(i, j), -pos_var(i, bit-1), pos_var(j, bit), pos_var(i, bit)])
                        formula.append([-var(i, j), -pos_var(i, bit-1), -pos_var(j, bit), -pos_var(i, bit)])

                        # ¬X0 ⇒ (Y1 = X1)
                        # X0 v ((Y1 v ¬X1) ^ (¬Y1 v X1))
                        # (X0 v (Y1 v ¬X1)) ^ (X0 v (¬Y1 v X1))
                        formula.append([-var(i, j), -pos_var(i, bit-1), pos_var(j, bit), -pos_var(i, bit)])
                        formula.append([-var(i, j), -pos_var(i, bit-1), -pos_var(j, bit), pos_var(i, bit)])
                    else:
                        # For bit positions i > 1, consider two bits at a time
                        # ¬Yi−1 ∧ Xi−1 ⇒ Yi = ¬Xi
                        # (Yi−1 v ¬Xi−1) v ((Yi v Xi) ^ (¬Yi v ¬Xi))
                        # (Yi−1 v ¬Xi−1 v Yi v Xi) ^ (Yi−1 v ¬Xi−1 v ¬Yi v ¬Xi)
                        
                        formula.append([-var(i, j), pos_var(j, bit-1), -pos_var(i, bit), pos_var(j, bit), pos_var(i, bit)])
                        formula.append([-var(i, j), pos_var(j, bit-1), -pos_var(i, bit), -pos_var(j, bit), -pos_var(i, bit)])

                        # ¬Yi−1 ∧ Xi−1 ∧ Xi ⇒ Yi+1 = ¬Xi+1
                        # (Yi−1 v ¬Xi−1 v ¬Xi) v ((Yi+1 v Xi+1) ^ (¬Yi+1 v ¬Xi+1))
                        # (Yi−1 v ¬Xi−1 v ¬Xi v Yi+1 v Xi+1) ^ (Yi−1 v ¬Xi−1 v ¬Xi v ¬Yi+1 v ¬Xi+1)
                        
                        formula.append([-var(i, j), pos_var(j, bit-1), -pos_var(i, bit-1), -pos_var(i, bit), pos_var(j, bit+1), pos_var(i, bit+1)])
                        formula.append([-var(i, j), pos_var(j, bit-1), -pos_var(i, bit-1), -pos_var(i, bit), -pos_var(j, bit+1), -pos_var(i, bit+1)])

                        # (Yi-1 v ¬Xi-1) ⇒ Yi = Xi ∧ Yi+1 = Xi+1
                        # (¬Yi−1 ∧ Xi−1) v ((Yi v Xi) ^ (Yi+1 v Xi+1) ^ (¬Yi v ¬Xi) ^ (¬Yi+1 v ¬Xi+1))
                        # (¬Yi−1 v Yi) ^ (¬Yi−1 v Xi) ^ (¬Yi−1 v Yi+1) ^ (¬Yi−1 v Xi+1) ^ (Xi−1 v Yi) ^ (Xi−1 v Xi) ^ (Xi−1 v Yi+1) ^ (Xi−1 v Xi+1)
                        
                        formula.append([-var(i, j), -pos_var(j, bit-1), pos_var(j, bit)])
                        formula.append([-var(i, j), -pos_var(j, bit-1), pos_var(i, bit)])
                        formula.append([-var(i, j), -pos_var(j, bit-1), pos_var(j, bit+1)])
                        formula.append([-var(i, j), -pos_var(j, bit-1), pos_var(i, bit+1)])
                        formula.append([-var(i, j), pos_var(j, bit-1), pos_var(j, bit)])
                        formula.append([-var(i, j), pos_var(j, bit-1), pos_var(i, bit)])
                        formula.append([-var(i, j), pos_var(j, bit-1), pos_var(j, bit+1)])
                        formula.append([-var(i, j), pos_var(j, bit-1), pos_var(i, bit+1)])
                        
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
filename = './input/hc-5.col'
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