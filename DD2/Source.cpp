#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <string>
#include <unordered_map>
#include <iomanip>

using namespace std;

struct Net {
    vector<int> components;
};

struct Netlist {
    int numCells;
    int numConnections;
    int numRows;
    int numColumns;
    vector<Net> nets;
};

Netlist parseNetlistFile(const string& filename) {
    ifstream file(filename);
    Netlist netlist;

    if (file.is_open()) {
        file >> netlist.numCells >> netlist.numConnections >> netlist.numRows >>
            netlist.numColumns;

        for (int i = 0; i < netlist.numConnections; ++i) {
            int numComponents;
            file >> numComponents;

            Net net;
            for (int j = 0; j < numComponents; ++j) {
                int component;
                file >> component;
                net.components.push_back(component);
            }

            netlist.nets.push_back(net);
        }

        file.close();
    }
    else {
        cout << "Unable to open file: " << filename << endl;
    }

    return netlist;
}

vector<vector<string>> createInitialGrid(const Netlist& netlist) {
    vector<vector<string>> grid(netlist.numRows,
        vector<string>(netlist.numColumns, "--"));
    return grid;
}

void placeCellsRandomly(vector<vector<string>>& grid, const Netlist& netlist) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> rowDist(0, netlist.numRows - 1);
    uniform_int_distribution<int> colDist(0, netlist.numColumns - 1);

    for (int cellId = 1; cellId <= netlist.numCells; ++cellId) {
        int row = rowDist(gen);
        int col = colDist(gen);

        while (grid[row][col] != "--") {
            row = rowDist(gen);
            col = colDist(gen);
        }

        grid[row][col] = to_string(cellId);
    }

}


double estimateWireLength(const vector<vector<string>>& grid, const Netlist& netlist) {
    double wireLength = 0.0; // HPWL Calculation
    unordered_map<string, pair<int, int>> componentPositions;

    // Build a map from component to its position
    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[0].size(); ++j) {
            componentPositions[grid[i][j]] = { i, j };
        }
    }

    for (const Net& net : netlist.nets) {
        int minRow = grid.size();    // Initialize with maximum possible value
        int maxRow = -1;             // Initialize with minimum possible value
        int minCol = grid[0].size(); // Initialize with maximum possible value
        int maxCol = -1;             // Initialize with minimum possible value

        for (int component : net.components) {
            auto it = componentPositions.find(to_string(component));
            if (it != componentPositions.end()) {
                minRow = min(minRow, it->second.first);
                maxRow = max(maxRow, it->second.first);
                minCol = min(minCol, it->second.second);
                maxCol = max(maxCol, it->second.second);
            }
        }

        wireLength += (maxRow - minRow) + (maxCol - minCol);
    }

    return wireLength;
}


double schedule_temp(double currentTemp, double coolingRate) {
    return currentTemp * coolingRate;
}

void simulatedAnnealing(const Netlist& netlist,
    const  double coolingRates) {
    // Create the initial placement
    vector<vector<string>> grid = createInitialGrid(netlist);
    placeCellsRandomly(grid, netlist);

    // Calculate the initial cost (wire length)
    double initialCost = estimateWireLength(grid, netlist);

    // Set initial temperature and final temperature based on initial cost
    double initialTemperature = 500 * initialCost;
    double finalTemperature = 5e-6 * initialCost / netlist.numConnections;
    double currentTemperature = initialTemperature;

    // Set number of moves per temperature based on number of cells
    int movesPerTemperature = 10 * netlist.numCells;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> probabilityDist(0.0, 1.0);

    while (currentTemperature > finalTemperature) {
        for (int i = 0; i < movesPerTemperature; ++i) {
            // Pick two random cells to swap
            int cell1Row = gen() % netlist.numRows;
            int cell1Col = gen() % netlist.numColumns;
            int cell2Row = gen() % netlist.numRows;
            int cell2Col = gen() % netlist.numColumns;

            // Check if both cells are empty before swapping
            if ((grid[cell1Row][cell1Col] != "--" && grid[cell2Row][cell2Col] != "--") ||
                (grid[cell1Row][cell1Col] != "--" && grid[cell2Row][cell2Col] == "--") ||
                (grid[cell1Row][cell1Col] == "--" && grid[cell2Row][cell2Col] != "--")) {
                // Swap the cells in the grid
                swap(grid[cell1Row][cell1Col], grid[cell2Row][cell2Col]);

                // Calculate the new cost (wire length)
                double newCost = estimateWireLength(grid, netlist);

                // Calculate the change in cost (delta L)
                double deltaL = newCost - initialCost;

                // Accept or reject the swap based on the Metropolis criterion
                if (deltaL < 0) {
                    // Accept the swap if it improves the cost
                    initialCost = newCost;
                }
                else {
                    // Reject the swap with a probability based on the temperature
                    double acceptanceProb = exp(-deltaL / currentTemperature);
                    if (probabilityDist(gen) > acceptanceProb) {
                        // Reject the swap and revert the cells back to their original positions
                        swap(grid[cell1Row][cell1Col], grid[cell2Row][cell2Col]);
                    }
                    else {
                        // Accept the swap
                        initialCost = newCost;
                    }
                }
            }
        }
        currentTemperature = schedule_temp(currentTemperature, coolingRates);
    }

    // Print the final placement and cost
    cout << "Final Placement:" << endl;
    for (const auto& row : grid) {
        for (const string& cell : row) {
            cout << setw(3) << cell << " ";
        }
        cout << endl;
    }

    cout << "Binary Placement:" << endl;
    for (const auto& row : grid) {
        for (const string& cell : row) {
            cout << (cell == "--" ? "0" : "1") << " ";
        }
        cout << endl;
    }

    cout << "Final Cost (Wire Length): " << initialCost << endl;
}


int main() {
    Netlist netlist = parseNetlistFile("netlist.txt");

    cout << "Number of cells: " << netlist.numCells << endl;
    cout << "Number of connections: " << netlist.numConnections << endl;
    cout << "Number of rows: " << netlist.numRows << endl;
    cout << "Number of columns: " << netlist.numColumns << endl;

    vector<vector<string>> grid = createInitialGrid(netlist);
    placeCellsRandomly(grid, netlist);

    cout << "Initial Grid:" << endl;
    for (int i = 0; i < netlist.numRows; ++i) {
        for (int j = 0; j < netlist.numColumns; ++j) {
            cout << setw(3) << grid[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Binary Grid:" << endl;
    for (const auto& row : grid) {
        for (const string& cell : row) {
            cout << (cell == "--" ? "0" : "1") << " ";
        }
        cout << endl;
    }

    cout << "Wirelength of initial placement: "
        << estimateWireLength(grid, netlist) << endl;


    double coolingRates = 0.75;

    auto startTime = chrono::high_resolution_clock::now();
    simulatedAnnealing(netlist, coolingRates);
    auto endTime = chrono::high_resolution_clock::now();

    double elapsedTime = chrono::duration<double>(endTime - startTime).count();

    cout << "Execution time: " << elapsedTime << "seconds";


    return 0;
}

