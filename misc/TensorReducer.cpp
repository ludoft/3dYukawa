#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <functional>

using namespace std;

int countLoops(const vector<pair<int, int>>& graph) {
    unordered_map<int, vector<int>> adjacencyList;
    unordered_map<int, bool> visited;
    int loops = 0;

    // Build the adjacency list
    for (const auto& edge : graph) {
        adjacencyList[edge.first].push_back(edge.second);
        adjacencyList[edge.second].push_back(edge.first);
    }

    // Initialize visited map
    for (const auto& edge : graph) {
        visited[edge.first] = false;
        visited[edge.second] = false;
    }

    // Recursive function to check loops
    function<void(int)> checkLoop = [&](int vertex) {
        visited[vertex] = true;
        for (int neighbor : adjacencyList[vertex]) {
            if (visited[neighbor] && neighbor != vertex) {
                loops++;
                return;
            }
            if (!visited[neighbor]) {
                checkLoop(neighbor);
            }
        }
    };

    // Iterate over the vertices and check loops
    for (const auto& edge : graph) {
        if (!visited[edge.first]) {
            checkLoop(edge.first);
        }
    }

    return loops;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: ./program input_file output_file" << endl;
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    ifstream input(inputFile);
    ofstream output(outputFile);

    if (!input) {
        cout << "Failed to open input file: " << inputFile << endl;
        return 1;
    }

    if (!output) {
        cout << "Failed to open output file: " << outputFile << endl;
        return 1;
    }

    string line;
    while (getline(input, line)) {
        vector<pair<int, int>> graph;
        istringstream iss(line);
        int vertex1, vertex2;
        while (iss >> vertex1 >> vertex2) {
            graph.push_back(make_pair(vertex1, vertex2));
        }

        int numLoops = countLoops(graph);
        output << "Number of loops: " << numLoops << endl;
    }

    input.close();
    output.close();

    return 0;
}

