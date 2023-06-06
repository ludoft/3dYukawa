#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

using namespace std;

//unordered_map<int, vector<int>> contraction_adjacency_dict;
vector<vector<int>> contraction_adjacency_dict;
const int num_colours = 3;

// Create a stringstream to store all results for this execution. This is output at the end.
stringstream resultStream;

pair<vector<vector<int>>, vector<vector<int>>> get_cycles_and_sequences_for_fixed_contraction(const vector<pair<int, int>>& graphWithoutContractions, const int & num_indices) {
    
    // Assumes the vertex IDs are consecutive integers starting from 1.
    vector<vector<int>> cycles;
    vector<vector<int>> sequences;

    vector<vector<int>> adjacency_dict(contraction_adjacency_dict); //Initialise it using the known data from the contraction_adjacency_dict.

    // graphWithoutContractions itself is never used again after this
    for (const auto& pair : graphWithoutContractions) { 
        int vertex1 = pair.first;
        int vertex2 = pair.second;
        adjacency_dict[vertex1].emplace_back(vertex2);
        adjacency_dict[vertex2].emplace_back(vertex1);
    }

    unordered_set<int> visited;

    /*
    unordered_map<int, int> lengths;
    for (const auto& entry : adjacency_dict) {
        lengths[entry.first] = entry.second.size();
    }
    */
    vector<pair<int, int>> lengths(num_indices + 1); //Reserve num_indices +1 
    // Iterate from 1, as that is the first vertex. We leave the zero entry empty.
    for (int i = 1; i < num_indices+1; ++i) {
        lengths[i] = make_pair(i,adjacency_dict[i].size());
    }

    vector<pair<int, int>> verticesAndValencies(lengths.begin() + 1, lengths.end()); //Once again, skip the first element, which we have left empty.
    sort(verticesAndValencies.begin(), verticesAndValencies.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    }); // Sort these in place


    // Iterating from smallest valency this way ensures that we always start in the beginning of sequences. This means we never start accidentally in the middle, then having to patch up the two ends afterwards, which would take time.
    for (const auto& vertexAndValency : verticesAndValencies) {
        int vertex = vertexAndValency.first;
        int valency = vertexAndValency.second;

        if (visited.find(vertex) == visited.end()) { //This does indeed check if an object is present.
            int previous = -1;
            int current = vertex;
            vector<int> path;
            path.push_back(current);

            const vector<int>& neighbors = adjacency_dict[current];
            if (neighbors[0] != previous) {
                previous = current;
                current = neighbors[0];
            } else {
                previous = current;
                current = neighbors[1];
            }

            while (true) {
                path.push_back(current);
                visited.insert(current);

                const vector<int>& neighbors = adjacency_dict[current];

                if (valency == 1 && neighbors.size() == 1) {
                    sequences.push_back(path);
                    //cout << "Seq" << previous << " " << current << endl;
                    break;
                }

                if (neighbors[0] != previous) { // Don't want to retrace our footsteps, so make sure we go to the neighbour not yet visited
                    previous = current;
                    current = neighbors[0];
                } else {
                    previous = current;
                    current = neighbors[1];
                }

                if (current == path[0]) {
                    //cout << "Cycl" << previous << " " << current << endl;
                    cycles.push_back(path);
                    break;
                }
            }
        }
    }

    return make_pair(cycles, sequences);
}

/*
pair<vector<vector<int>>, vector<vector<int>>> get_cycles_and_sequences(const vector<pair<int, int>>& graph, const int& num_vertices) {
    // Assumes the vertex IDs are consecutive integers starting from 1.
    vector<vector<int>> cycles;
    vector<vector<int>> sequences;

    unordered_map<int, vector<int>> adjacency_dict;
    
    //vector<vector<int>> adjacency_dict;
    //adjacency_list.reserve(graph.size()+1); // How many indices are there? Who knows.

    // graph itself is never used again after this
    for (const auto& pair : graph) { 
        int vertex1 = pair.first;
        int vertex2 = pair.second;
        adjacency_dict[vertex1].push_back(vertex2);
        adjacency_dict[vertex2].push_back(vertex1);
    }

    unordered_set<int> visited;

    unordered_map<int, int> lengths;
    for (const auto& entry : adjacency_dict) {
        lengths[entry.first] = entry.second.size();
    }
    //vector<pair<int,int>> lengths(adjacency_list.size());  // Using vector instead of unordered_map
    //for (int i = 1; i < adjacency_list.size(); ++i) { lengths[i] = make_pair(i,adjacency_list[i].size()); }

    vector<pair<int, int>> verticesAndValencies(lengths.begin(), lengths.end());
    sort(verticesAndValencies.begin(), verticesAndValencies.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    for (const auto& vertexAndValency : verticesAndValencies) {
        int vertex = vertexAndValency.first;
        int valency = vertexAndValency.second;

        if (visited.find(vertex) == visited.end()) { //This does indeed check if an object is present.
            int previous = -1;
            int current = vertex;
            vector<int> path;
            path.push_back(current);
            //cout << current << endl;

            const vector<int>& neighbors = adjacency_dict[current];
            if (neighbors[0] != previous) {
                previous = current;
                current = neighbors[0];
            } else {
                previous = current;
                current = neighbors[1];
            }

            while (true) {
                //cout << current << endl;
                path.push_back(current);
                visited.insert(current);

                const vector<int>& neighbors = adjacency_dict[current];

                if (valency == 1 && neighbors.size() == 1) {
                    sequences.push_back(path);
                    //cout << "Seq" << previous << " " << current << endl;
                    break;
                }

                if (neighbors[0] != previous) { // Don't want to retrace our footsteps, so make sure we go to the neighbour not yet visited
                    previous = current;
                    current = neighbors[0];
                } else {
                    previous = current;
                    current = neighbors[1];
                }

                if (current == path[0]) {
                    //cout << "Cycl" << previous << " " << current << endl;
                    cycles.push_back(path);
                    break;
                }
            }
        }
    }

    return make_pair(cycles, sequences);
}
*/


std::vector<std::pair<int, int>> parseMathematicaList(const std::string& input) {
    std::vector<std::pair<int, int>> result;
    std::stringstream ss(input);
    char c;
    int num1, num2;

     if (ss >> c && c == '{') {
        while (ss >> c && c != '}') {
            if (c == '{') {
                ss >> num1;
                if (ss >> c && c == ',') {
                    ss >> num2;
                    result.push_back(std::make_pair(num1, num2));
                    ss >> c; //Drop the } from the pair.
                }
            }
        }
    }

    return result;
}

std::vector<std::pair<int, int>> parseSpaceSeparatedList(const std::string& input) {
    std::vector<std::pair<int, int>> pairs;
    std::stringstream ss(input);
    int num;

    while (ss >> num) {
        int secondNum;
        ss >> secondNum;
        pairs.emplace_back(num, secondNum);
    }
    /*If the input is malformed, this will not notice.*/

    return pairs;
}

std::vector<int> flattenVector(const std::vector<std::vector<int>>& inputVector) {
    std::vector<int> flattened;
    for (const auto& subVector : inputVector) {
        flattened.insert(flattened.end(), subVector.begin(), subVector.end());
    }
    return flattened;
}

/*
 * This is just the following code:
out = {15, 3, 36, 33, 1, 13, 34, 31, 14, 2, 32, 35};
mapToNew = Thread[Sort[out] -> Range@Length[out]];
partitionedNew = Partition[out /. mapToNew, 2]

num_colours = 3;
thusResulting =
 Sort /@ SortBy[partitionedNew, Mod[# - 1, num_colours] &]
 */

std::vector<int> getTransposeForm(const std::vector<int>& out) {
    // Create a map to assign new values to sorted elements
    //for (auto const& c :fout)std::cout << c << ' ';
    //std::vector<int> out = {15, 3, 36, 33, 1, 13, 34, 31, 14, 2, 32, 35};
    std::unordered_map<int, int> mapToNew;
    std::vector<int> sortedOut = out;
    std::sort(sortedOut.begin(), sortedOut.end());
    for (size_t i = 0; i < sortedOut.size(); i++) {
        mapToNew[sortedOut[i]] = i + 1;
    }

    // Apply transformations and partition the resulting vector
    std::vector<int> partitionedNew;
    for (const auto& num : out) {
        partitionedNew.push_back(mapToNew[num]);
    }
    std::vector<std::vector<int>> partitioned;
    for (size_t i = 0; i < partitionedNew.size(); i += 2) {
        partitioned.push_back({partitionedNew[i], partitionedNew[i + 1]});
    }

    // Sort each pair within the resulting vector
    for (auto& pair : partitioned) {
        std::sort(pair.begin(), pair.end());
    }

    // Sort the transformed pairs by modulo operation, i.e. by colour. Break ties by whichever number is smaller.
    std::sort(partitioned.begin(), partitioned.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
            if (a[0] % num_colours == b[0] % num_colours){
                return a[0]<b[0];
            } else{
                return (a[0] - 1) % num_colours < (b[0] - 1) % num_colours;
            }

        } // comp: comparison function object (i.e. an object that satisfies the requirements of Compare) which returns true if the first argument is less than (i.e. is ordered before) the second. 
    );

    return flattenVector(partitioned);
} 

void analyseIndividualLine(const std::string & line, const int & num_indices) {
        std::vector<std::pair<int, int>> graph =parseSpaceSeparatedList(line);
        // for (const auto& pair : graph) { std::cout << "{" << pair.first << ", " << pair.second << "}" << std::endl; }
        
        //cout<< line<< endl;
        vector<vector<int>> cycles, sequences;
        tie(cycles, sequences) = get_cycles_and_sequences_for_fixed_contraction(graph, num_indices);

        /*
        cout << "Cycles:" << cycles.size() << endl;
        for (const auto& cycle : cycles) {
            for (int vertex : cycle) { cout << vertex << " "; resultStream << vertex << " ";}
            cout << " So colour is " << ((cycle.front()-1) % num_colours) + 1 << endl; }

        cout << "Sequences:" << endl;
        for (const auto& sequence : sequences) {
            for (int vertex : sequence) { cout << vertex << " "; resultStream << vertex << " ";}
            cout << endl; }*/

        std::vector<int> vectorOfColourLoopCounts(num_colours, 0); // Create a vector of integers with length num_colours, initialised to zero.

        for (const auto& cycle : cycles) {
            vectorOfColourLoopCounts[((cycle.front()-1) % num_colours)+ 1 - 1]++; // Subtract one to get 0-indexed object. The -1 +1 otherwise is related to the mod, to ensure that we get colour=3 rather than colour =0 when colour=num_colours;
            /*int whichColourInt =((cycle.front()-1) % num_colours) + 1;
            std::cout << whichColourInt <<  " " << vectorOfColourLoopCounts[whichColourInt] << std::endl;*/
        }
        
        resultStream << "{"<<  cycles.size() << ", n[";

         for (size_t i = 0; i < num_colours; ++i) {
            //std::cout << "Colour: " << (i+1) << ", Value: " << vectorOfColourLoopCounts[i] << std::endl;
            resultStream << vectorOfColourLoopCounts[i] << ",";
        }

        resultStream.seekp(-1,resultStream.cur); resultStream<< "], TensorTranspose[obj, {";

        std::vector<int> transposeFormToBe;
        transposeFormToBe.reserve(sequences.size()* 2);

        for (const auto& sequence : sequences) {
                //resultStream << sequence.front() << "," << sequence.back() << ",";
                transposeFormToBe.emplace_back(sequence.front());
                transposeFormToBe.emplace_back(sequence.back());
        }
        
        transposeFormToBe = getTransposeForm(transposeFormToBe);
        for (const auto& index: transposeFormToBe) {
                resultStream << index << ",";
        }
        resultStream.seekp(-1,resultStream.cur);
        resultStream<< "}]}" << endl;

}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cout << "Usage: ./get_cycles_and_sequences.out inputTrules inputContraction output_file" << endl;
        return 1;
    }
    const int num_colours=3;

    string inputFilename = argv[1];
    string contractionFilename = argv[2];
    string outputFilename = argv[3];

    ifstream inputFile(inputFilename);
    if (!inputFile) { cout << "Error opening input file: " << inputFilename << endl; return 1; }
    vector<string> inputLines;

    std::string line;
    while (std::getline(inputFile, line)) {
        inputLines.push_back(line);
    }
    inputFile.close();

    std::unordered_set<int> uniqueIntegers;
    vector<pair<int, int>> bla = parseSpaceSeparatedList(inputLines[0]);
    for (const auto& pair : bla) {
        uniqueIntegers.insert(pair.first);
        uniqueIntegers.insert(pair.second);
    }
    int num_indices = uniqueIntegers.size();
    // cout << num_indices << " is the graph size." << endl;

    ifstream contractionFile(contractionFilename);
    if (!contractionFile) { cout << "Error opening contraction file: " << contractionFilename << endl; return 1; }
    std::getline(contractionFile, line);
    contractionFile.close();

    std::vector<std::pair<int,int>> contractionRules= parseMathematicaList(line);
         // graph itself is never used again after this
    
     contraction_adjacency_dict.resize(num_indices+1);
     for (const auto& pair : contractionRules) {
         int vertex1 = pair.first;
         int vertex2 = pair.second;
         //cout << num_indices << " " <<  contraction_adjacency_dict.size() << endl;
         //cout << vertex1 << " " <<  vertex2 << endl;
         contraction_adjacency_dict[vertex1].push_back(vertex2);
         contraction_adjacency_dict[vertex2].push_back(vertex1);
     }

    //while (std::getline(inputFile, line)) {
    for (const auto& line : inputLines){
        analyseIndividualLine(line, num_indices);
    }

    //inputFile.close();

    ofstream outputFile(outputFilename);
    if (!outputFile) { cout << "Error opening output file: " << outputFilename << endl; return 1; }

    outputFile << resultStream.str();
    outputFile.close();


    //cout << "Results saved to: " << outputFilename << endl;

    return 0;
}
