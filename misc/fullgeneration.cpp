#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <functional>
#include <stack>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <omp.h>
//#include <boost/interprocess/file_mapping.hpp>
//#include <boost/interprocess/mapped_region.hpp>
//#include "mathlink.h"

using namespace std;

typedef int vert;
typedef vector<array<int,2>> adjList;
//typedef string analysisResultType;
typedef pair<uint, vector<vert>> analysisResultType;
typedef vector<size_t> IndexTuple;
typedef pair<IndexTuple, analysisResultType> indexedResult;

const int num_colours = 3;

vector<array<vert,2>>  contractionRules;
adjList totalAdjList;

pair<vector<vector<vert>>, vector<vector<vert>>> get_cycles_and_sequences_for_fixed_contraction(const adjList& adjacency_dict) {
    
    const int numVerticesTotalPlusOne =adjacency_dict.size();
    // Assumes the vertex IDs are consecutive integers starting from 1.
    vector<vector<vert>> cycles;
    vector<vector<vert>> sequences;

    /*vector<vector<int>> adjacency_dict(contraction_adjacency_dict); //Initialise it using the known data from the contraction_adjacency_dict.

    // graphWithoutContractions itself is never used again after this
    for (const auto& pair : graphWithoutContractions) { 
        int vertex1 = pair.first;
        int vertex2 = pair.second;
        adjacency_dict[vertex1].emplace_back(vertex2);
        adjacency_dict[vertex2].emplace_back(vertex1);
    }*/

    unordered_set<vert> visited;

    /*
    unordered_map<int, int> lengths;
    for (const auto& entry : adjacency_dict) {
        lengths[entry.first] = entry.second.size();
    }
    */
    vector<pair<vert, vert>> lengths(numVerticesTotalPlusOne); //Reserve num_indices +1 
    // Iterate from 1, as that is the first vertex. We leave the zero entry empty.
    
    for (int i = 1; i < numVerticesTotalPlusOne; ++i) {
        lengths[i] = make_pair(i, adjacency_dict[i][1]); //This is zero if we are on a valence one index, and nonzero otherwise. Perfect.
    }

    vector<pair<vert,vert>> verticesAndValencies(lengths.begin() + 1, lengths.end()); //Once again, skip the first element, which we have left empty.
    sort(verticesAndValencies.begin(), verticesAndValencies.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    }); // Sort these in place


    // Iterating from smallest valency this way ensures that we always start in the beginning of sequences. This means we never start accidentally in the middle, then having to patch up the two ends afterwards, which would take time.
    for (const auto& vertexAndValency : verticesAndValencies) {
        vert vertex = vertexAndValency.first;
        vert isValencyTwo = vertexAndValency.second; //Valency is 0 if started on vertex of valency one, otherwise nonzero.

        if (visited.find(vertex) == visited.end()) { //This does indeed check if an object is present.
            vert previous = 0;
            vert current = vertex;
            vector<vert> path;
            path.push_back(current);

            const array<vert,2>& neighbors = adjacency_dict[current];
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

                const array<vert,2>& neighbors = adjacency_dict[current];

                if (neighbors[1] == 0 && isValencyTwo == 0) { //isValencyTwo ==0 means that we started on a vertex of valency 1.
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

std::vector<std::pair<int, int>> parseMathematicaList(const std::string& input) {
    std::vector<std::pair<int, int>> result;
    std::istringstream ss(input);
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
    std::istringstream ss(input);
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




std::vector<adjList> processFile(const std::string& filename, int& totalCount) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return {};
    }

    std::string line;
    std::getline(inputFile, line);

    std::istringstream iss(line);
    std::vector<int> firstLine;
    int num;
    while (iss >> num) {
        firstLine.push_back(num);
    }



    std::vector<adjList> result;
    adjList firstAdjList;

    const int size=firstLine.size();
    firstAdjList.resize(size);
    for (size_t i = 0; i < size; ++i) {
        firstAdjList[i][0] = firstLine[i] + totalCount;
    }
    result.push_back(firstAdjList);
    // Output some data:
    //cout<<size<<endl;
    //cout<<firstAdjList.size()<<endl;

    std::vector<std::array<int, 2>> row(size);

    while (std::getline(inputFile, line)) {
        std::istringstream issLine(line);
        int value;
        std::array<int, 2> subvector;
        subvector[1] = 0;
        for (int i = 0; i < size && issLine >> value; ++i) {
            subvector[0] = value + totalCount;
            row[i] = subvector;
        }
        //adjList list; list.vertices = row;
        result.push_back(row);
    }
    inputFile.close();

    totalCount += size;
    return result;
}



//using IndexTuple = std::vector<size_t>;

template <typename Function, typename List, typename Result>
vector<pair<IndexTuple,Result>> mapTuples(Function&& function, const std::vector<List>& lists) {
    const size_t numLists = lists.size();
    std::vector<size_t> sizes(numLists);
    IndexTuple indices(numLists, 0);


    size_t totalTuples = 1;
    for (size_t i = 0; i < numLists; ++i) {
        sizes[i] = lists[i].size();
        //if (sizes[i] == 0) { totalTuples = 0; break; }
        totalTuples *= sizes[i];
    }
    //vector<pair<IndexTuple,Result>> results(totalTuples);
    //vector<Result> results(totalTuples);
    vector<pair<IndexTuple,Result>> results(totalTuples);
    std::vector<typename List::value_type> elements(numLists);

    
    for (size_t i = 0; i < totalTuples; ++i) {
        for (size_t j = 0; j < numLists; ++j) {
            elements[j] = lists[j][indices[j]];
        }

        //results[i] = make_pair(indices,function(elements));
        results[i] = make_pair(indices,function(elements));

        size_t k = numLists - 1;
        while (k != static_cast<size_t>(-1)) {
            if (++indices[k] < sizes[k])
                break;

            indices[k] = 0;
            --k;
        }
    }
    /*
    vector<IndexTuple> inputIndices(totalTuples);

    for (size_t i = 0; i < totalTuples; ++i) {
        inputIndices[i]=indices;

        size_t k = numLists - 1;
        while (k != static_cast<size_t>(-1)) {
            if (++indices[k] < sizes[k])
                break;

            indices[k] = 0;
            --k;
        }
    };
#pragma omp parallel shared(function, elements, indices) private(contractionRules, totalAdjList)
{
// Private copy of the global variable for each thread

    #pragma omp for
    for (size_t i = 0; i < totalTuples; ++i) {
        for (size_t j = 0; j < numLists; ++j) {
            elements[j] = lists[j][inputIndices[i][j]];
        }
        results[i] = function(elements);
    };
}
*/
    return results;
}

// Helper function to convert a vector to a string
string vectorToString(const vector<int>& vec) {
    stringstream ss;
    for (const auto& value : vec) {
        ss << value << ",";
    }
    return ss.str();
}

// Helper function to convert a vector to a string
string makeMMAStringOutput(const IndexTuple& indices, const int numLoops) {
    ostringstream ss;
    ss << "c[" << indices[0]+1;//Add one, as Mathematica indexes from one.
    for (size_t i=1; i<indices.size(); ++i) {
        ss <<","<< indices[i]+1; //Add one, as Mathematica indexes from one.
    }
    if(numLoops==0){
        ss << "]";
    } else {
        ss << "]"<<"N^" << numLoops;
    }
    /*
    ss << "c[";
    for (const auto& str : indices){ ss<<str<<",";}
    ss << "]"<<"N^" << numLoops;
    */
    
    return ss.str();
}

template <typename Function, typename List, typename Result>
map<string,vector<string>> mapTuplesAndGroup(Function&& function, const std::vector<List>& lists) {
    const size_t numLists = lists.size();
    std::vector<size_t> sizes(numLists);

    IndexTuple indices(numLists, 0);

    map<string,vector<string>> groups;

    size_t totalTuples = 1;
    for (size_t i = 0; i < numLists; ++i) {
        sizes[i] = lists[i].size();
        totalTuples *= sizes[i];
    }
    //vector<pair<IndexTuple,Result>> results(totalTuples);
    std::vector<typename List::value_type> elements(numLists);

    Result thisResult;
    string stringKey;
    
    for (size_t i = 0; i < totalTuples; ++i) {
        for (size_t j = 0; j < numLists; ++j) {
            elements[j] = lists[j][indices[j]];
        }

        //results[i] = make_pair(indices,function(elements));
        thisResult = function(elements);
        groups[vectorToString(thisResult.second)].push_back(makeMMAStringOutput(indices,thisResult.first));

        size_t k = numLists - 1;
        while (k != static_cast<size_t>(-1)) {
            if (++indices[k] < sizes[k])
                break;

            indices[k] = 0;
            --k;
        }
    }
    return groups;
}

//string analyseListOfAdjLists(const vector<adjList>& adjLists, adjList& totalAdjList, const int& numVerticesTotal) {
analysisResultType analyseListOfAdjLists(const vector<adjList>& adjLists) {
    //std::cout << "Elements: "; for (const auto& element :totalAdjList) { std::cout << element[0] << " "; } std::cout << std::endl;

    // Copy the vectors consecutively into the merged vector. Recall that we need to mess with the ordering to get the contractions to be correct.
    //auto currentPos = totalAdjList.begin()+ 1; //Not an integer but something else
    int currentPos =  1; //Not an integer but something else

    for (const auto& vector : adjLists) {
        //std::cout << "Elements: "; for (const auto& element :totalAdjList) { std::cout << element[0] << " "; } std::cout << std::endl;
        size_t vecSize = vector.size();
        for (size_t j = 0; j < vecSize; ++j) {
            totalAdjList[currentPos + j] =vector[j];
        }
        currentPos += vecSize;
    }
    //std::cout << "Elements: "; for (const auto& element :totalAdjList) { std::cout << "(" << element[0] << " " << element[1] << ") "; } std::cout << std::endl;


    //vector<array<int,2>>  contractionRules = {{7, 25}, {8, 26}, {9, 27}, {10, 28}, {11, 29}, {12, 30}};// {{7, 13}, {8, 14}, {9, 15}, {10, 16}, {11, 17}, {12, 18}, {19, 22}, {20, 23}, {21, 24}} ;

    // graphWithoutContractions itself is never used again after this
    for (const auto& pair : contractionRules) { 
        int vertex1 = pair[0];
        int vertex2 = pair[1];
        totalAdjList[vertex1][1] = vertex2;
        totalAdjList[vertex2][1] = vertex1; //First entry is always populated.
    }
    //std::cout << "Elements: "; for (const auto& element :totalAdjList) { std::cout << "(" << element[0] << " " << element[1] << ") "; } std::cout << std::endl;

    vector<vector<int>> cycles, sequences;
    tie(cycles, sequences) = get_cycles_and_sequences_for_fixed_contraction(totalAdjList);

    /*
    cout << "Cycles:" << cycles.size() << endl;
    for (const auto& cycle : cycles) {
        for (int vertex : cycle) { cout << vertex << " "; }
        cout << " So colour is " << ((cycle.front()-1) % num_colours) + 1 << endl; }

    cout << "Sequences:" << endl;
    for (const auto& sequence : sequences) {
        for (int vertex : sequence) { cout << vertex << " ";};
        cout << endl; }
        */

    std::vector<int> vectorOfColourLoopCounts(num_colours, 0); // Create a vector of integers with length num_colours, initialised to zero.

    for (const auto& cycle : cycles) {
        vectorOfColourLoopCounts[((cycle.front()-1) % num_colours)+ 1 - 1]++; // Subtract one to get 0-indexed object. The -1 +1 otherwise is related to the mod, to ensure that we get colour=3 rather than colour =0 when colour=num_colours;
        /*int whichColourInt =((cycle.front()-1) % num_colours) + 1;
        std::cout << whichColourInt <<  " " << vectorOfColourLoopCounts[whichColourInt] << std::endl;*/
    }
    
    /*
    // ostringstream thisResultStream;
    thisResultStream << "{"<<  cycles.size() << ", "<<"n[";

     for (size_t i = 0; i < num_colours; ++i) {
        //std::cout << "Colour: " << (i+1) << ", Value: " << vectorOfColourLoopCounts[i] << std::endl;
        thisResultStream << vectorOfColourLoopCounts[i] << ",";
    }

    thisResultStream.seekp(-1,thisResultStream.cur);
    thisResultStream<< "], "<<"TensorTranspose[obj, {";
    */

    //thisResultStream << "{"<<  cycles.size() << ", "<<"{";

    std::vector<int> transposeFormToBe;
    transposeFormToBe.reserve(sequences.size()* 2);
for (const auto& sequence : sequences) {
            //thisResultStream << sequence.front() << "," << sequence.back() << ",";
            transposeFormToBe.emplace_back(sequence.front());
            transposeFormToBe.emplace_back(sequence.back());
    }
    
    transposeFormToBe = getTransposeForm(transposeFormToBe);
    /*
    for (const auto& index: transposeFormToBe) {
            thisResultStream << index << ",";
    }
    thisResultStream.seekp(-1,thisResultStream.cur);
    //thisResultStream<< "}]}" << endl;
    thisResultStream<< "}}" << endl;
    */

    //Force segfault: transposeFormToBe[10000*currentPos] =10 ;
    //return thisResultStream.str();
    

    return make_pair(cycles.size(), transposeFormToBe);
}

const map<char, string> filePathMap = {
        {'g', "/tmp/gAdjList.al"},
        {'h', "/tmp/hAdjList.al"},
        {'l', "/tmp/lAdjList.al"}
};


vector<string> mapStringToFilePaths(const string& input) {
    vector<string> result;
    for (const auto& ch : input) {
        auto it = filePathMap.find(ch);
        if (it != filePathMap.end()) {
            result.push_back(it->second);
        }
    }

    return result;
}


vector<array<int,2>> parseStringContractionRules(const string& input) {
    vector<array<int,2>> result;
    istringstream ss(input);
    string token;

    while (getline(ss, token, '-')) {
        int first = stoi(token);
        getline(ss, token, '-');
        int second = stoi(token);
        result.push_back({first, second});
    }

    return result;
}


map<string,vector<string>> GroupBy22(const vector<indexedResult>&indexedResults) {
    //map<string, vector<pair<IndexTuple, int>>> groups;
    map<string,vector<string>> groups;

    // Group pairs based on the second element (string)
    for (const auto& pair : indexedResults) {
        string key = vectorToString(pair.second.second);
        groups[key].push_back(makeMMAStringOutput(pair.first, pair.second.first));
    }

    /*
    // Convert the map to a vector of vectors
    vector<vector<string>> result;
    for (const auto& group : groups) {
        result.push_back(group.second);
    }*/

    return groups;
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        cout << "Usage: ./get_cycles_and_sequences.out <contraction 1-4-5-8...> <vertices e.g. gghl> <outputfile>" << endl;
        return 1;
    }
    contractionRules = parseStringContractionRules(argv[1]);

    vector<string> fileNames =  mapStringToFilePaths(argv[2]);//{"/tmp/gAdjList.al", "/tmp/hAdjList.al"};//, "/tmp/gAdjList.al"};
    int totalCount = 0;  // Running total of totalCount
    vector<vector<adjList>> allNumbers;
     
    
     for (const auto& filename : fileNames) {
        std::vector<adjList> output = processFile(filename, totalCount);
         allNumbers.push_back(output);
         
         /*cout << "Total number of integers in " << filename << ": " << totalCount << endl;
        for (const auto& list : output) { for (const auto& subvector : list) { std::cout << "(" << subvector[0] << ", " << subvector[1] << ") "; } std::cout << std::endl; }*/
 }

    //cout<<"Total size: " << totalCount<<endl;
    //adjList totalAdjList(totalCount+1); //Add one to size.
    totalAdjList.resize(totalCount+1);
    totalAdjList[0] ={0,0};
    //int count = 0; mapTuples([&count, &totalAdjList, &totalCount](const std::vector<adjList>& elements) { count++; resultStream<<analyseListOfAdjLists(elements,totalAdjList,totalCount);  },allNumbers); std::cout << "Total tuples: " << count << std::endl;

    //omp_set_num_threads(8);
    //vector<indexedResult> allTuplesRun = mapTuples<decltype(analyseListOfAdjLists), vector<adjList>,analysisResultType>(analyseListOfAdjLists, allNumbers); 

    //Trying to get it work without having to explicitly spell out the template.
    //auto testLambdaFunction =[](const vector<adjList>& element){return analyseListOfAdjLists(element);};
    //vector<string> allTuplesRun = mapTuples<decltype(analyseListOfAdjLists), vector<adjList>(testLambdaFunction, allNumbers); 
    //std::cout << "Total tuples: " << allTuplesRun.size() << std::endl;
    

    map<string,vector<string>> gathered = mapTuplesAndGroup<decltype(analyseListOfAdjLists), vector<adjList>,analysisResultType>(analyseListOfAdjLists, allNumbers); 

    /********************************/
    /*** Now perform output */
    string outputFilename = "/tmp/fastOutput-"+string(argv[2])+"-c"+string(argv[1])+".out";
    ofstream outputFile(outputFilename);
    if (!outputFile) { cout << "Error opening output file: " << outputFilename << endl; return 1; }

    // Concatenate the strings
    std::ostringstream resultStream;
    // Create an ostringstream to store all results for this execution. This is output at the end.

    //map<string,vector<string>> gathered = GroupBy22(allTuplesRun);
    resultStream<< "tensJustNums[{";
    for (const auto& group: gathered) {
            resultStream << "{";
            for (const auto& str: group.second) {
                resultStream << str << "+";
            }
            resultStream.seekp(-1,resultStream.cur);
            resultStream << ",{" << group.first;
            resultStream.seekp(-1,resultStream.cur);
            resultStream << "}},";
        }
        resultStream.seekp(-1,resultStream.cur); //Remove final comma
    resultStream<< "}, obj]";
    outputFile << resultStream.str();
    outputFile.close();

    return 0;
}
