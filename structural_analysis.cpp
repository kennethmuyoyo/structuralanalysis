#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
int main() {
    // Prompt the user to input the beam parameters
    double E, L, A, I;
    cout << "Enter Young's modulus (Pa): ";
    cin >> E;
    cout << "Enter length of the beam (m): ";
    cin >> L;
    cout << "Enter cross-sectional area of the beam (m^2): ";
    cin >> A;
    cout << "Enter moment of inertia of the beam (m^4): ";
    cin >> I;

    // Prompt the user to input the load parameters
    double P, q;
    cout << "Enter point load at the middle of the beam (N): ";
    cin >> P;
    cout << "Enter uniformly distributed load along the beam (N/m): ";
    cin >> q;

    // Prompt the user to input the finite element parameters
    int n;
    cout << "Enter the number of elements: ";
    cin >> n;
    int m = n + 1;
    double h = L / n;
}
// Define the global stiffness matrix
vector<vector<double> > K(m, vector<double>(m, 0.0));

// Define the global force vector
vector<double> F(m, 0.0);

// Define the global displacement vector
vector<double> U(m, 0.0);

// Function to assemble the global stiffness matrix from the local stiffness matrices
void assembleK() {
    // Loop over each element
    for (int e = 0; e < n; e++) {
        // Calculate the local stiffness matrix for each element
        vector<vector<double> > k(2, vector<double>(2));
        k[0][0] = k[1][1] = E * A / h + 12 * E * I / pow(h, 3);
        k[0][1] = k[1][0] = -E * A / h + 6 * E * I / pow(h, 2);
        k[1][1] += -12 * E * I / pow(h, 3);

        // Add the local stiffness matrix to the global stiffness matrix
        int i = e;
        int j = e + 1;
        K[i][i] += k[0][0];
        K[i][j] += k[0][1];
        K[j][i] += k[1][0];
        K[j][j] += k[1][1];
    }
}

// Function to assemble the global force vector from the nodal forces
void assembleF() {
    // Loop over each node
    for (int i = 0; i < m; i++) {
        // Calculate the nodal force due to the uniformly distributed load
        F[i] += q * h / 2;

        // Add the point load at the middle node
        if (i == m / 2) {
            F[i] += P;
        }
    }
}

// Function to apply the boundary conditions by modifying the global stiffness matrix and force vector
void applyBC() {
    // Fix the displacement at the first and last nodes to zero
    U[0] = U[m - 1] = 0.0;

    // Modify the global stiffness matrix and force vector accordingly
    for (int i = 1; i < m - 1; i++) {
        F[i] -= K[i][0] * U[0];
        F[i] -= K[i][m - 1] * U[m - 1];
        K[i][0] = K[0][i] = K[i][m - 1] = K[m - 1][i] = 0.0;
    }
    K[0][0] = K[m - 1][m - 1] = 1.0;
    F[0] = F[m - 1] = 0.0;
}

// Function to solve the system of equations using Gaussian elimination
void solve() {
    // Forward elimination
    for (int k = 0; k < m - 1; k++) {
        for (int i = k + 1; i < m; i++) {
            double factor = K[i][k] / K[k][k];
            for (int j = k + 1; j < m; j++) {
                K[i][j] -= factor * K[k][j];
            }
            F[i] -= factor * F[k];
        }
    }

    // Backward substitution
    for (int i = m - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < m; j++) {
            sum += K[i][j] * U[j];
        }
        U[i] = (F[i] - sum) / K[i][i];
    }
}

// Function to calculate and print the stress and strain at each node
void printResults() {
    cout << "Node\tDisplacement\tStrain\t\tStress" << endl;
    cout << fixed;
    cout.precision(6);
    for (int i = 0; i < m; i++) {
        double strain = U[i] / h;
        double stress = E * strain;
        cout << i + 1 << "\t" << U[i] << "\t\t" << strain << "\t\t" << stress << endl;
    }
}

// Main function
int main() {
    // Assemble the global stiffness matrix and force vector
    assembleK();
    assembleF();

    // Apply the boundary conditions
    applyBC();

    // Solve the system of equations
    solve();

    // Print the results
    printResults();

    return 0;
}
