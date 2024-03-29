#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <functional>
#include <stack>
#include <fstream>
#include <sstream>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <omp.h>
//#include <boost/interprocess/file_mapping.hpp>
//#include <boost/interprocess/mapped_region.hpp>

#include "vertexData.h"

using namespace std;

typedef vector<size_t> IndexTuple;
typedef string indexAsKey;
//typedef string transposeRuleAsKey;
typedef vector<int> transposeRuleAsKey;

struct FastVectorHasher { //https://stackoverflow.com/a/53283994
    int operator()(const vector<int> &V) const {
        int hash = V.size();
        for(auto &i : V) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};
struct VectorHasher {
    int operator()(const std::vector<int>& vec) const {
        int hash = 0;
        constexpr int prime = 31;
        for (const auto& num : vec) {
            hash = hash * prime + (num);
        }
        return hash;
    }
};
//typedef unordered_map<indexAsKey, map<int, int>,VectorHasher> umap_indices_to_loopcounter;
typedef unordered_map<indexAsKey, map<int, int>> umap_indices_to_loopcounter;
typedef int vert;
typedef vector<array<int,2>> adjList;
typedef pair<uint, vector<vert>> analysisResultType;
typedef pair<IndexTuple, analysisResultType> indexedResult;

// constexpr indicates a value that's constant and known during compilation.  const indicates a value that's only constant; it's not compulsory to know during compilation.
constexpr int num_colours = 3;

int num_coupling_constants = 0;
vector<array<vert,2>>  contractionRules;
adjList totalAdjList;
vector<vector<int>> ccIntToCoeffIntCanonicaliser;
vector<unordered_map<string,string>> ccIntToMapCoeffCanonicalIntToStr;

pair<vector<vector<vert>>, vector<vector<vert>>> get_cycles_and_sequences_for_fixed_contraction(const adjList& adjacency_dict) {
    
    //So in case your defining constants inside a function and the function is called multiple times, then to reduce the cost of constructing this object every time the function is called you would make it a static variable (Mainly useful when its a constant, Else you may end up changing this value and it soon gets messed up).
    static const int numVerticesTotalPlusOne =adjacency_dict.size();
    // Assumes the vertex IDs are consecutive integers starting from 1.
    vector<vector<vert>> cycles;
    vector<vector<vert>> sequences;
    cycles.reserve(10); //Guess at a reasonable number of cycles.
    sequences.reserve(6*num_coupling_constants);

    //unordered_set<vert> visited(numVerticesTotalPlusOne-1);
    vector<bool> visited(numVerticesTotalPlusOne, false);

    vector<pair<vert, vert>> lengths(numVerticesTotalPlusOne); //Reserve num_indices +1 
    // Iterate from 1, as that is the first vertex. We leave the zero entry empty.
    
    for (int i = 1; i < numVerticesTotalPlusOne; ++i) {
        lengths[i] = make_pair(i, adjacency_dict[i][1]); //This is zero if we are on a valence one index, and nonzero otherwise. Perfect.
    }

    vector<pair<vert,vert>> verticesAndValencies(lengths.begin() + 1, lengths.end()); //Once again, skip the first element, which we have left empty.
    sort(verticesAndValencies.begin(), verticesAndValencies.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    }); // Sort these in place


    // Iterating from smallest valency this way ensures that we always start in the beginning of sequences.
    // This means we never start accidentally in the middle, then having to patch up the two ends afterwards, which would take time.
    for (const auto& vertexAndValency : verticesAndValencies) {
        vert vertex = vertexAndValency.first;
        vert isValencyTwo = vertexAndValency.second; //Valency is 0 if started on vertex of valency one, otherwise nonzero.

        //if (visited.find(vertex) == visited.end()) { //This does indeed check if an object is present.
        if (!visited[vertex]) { //This does indeed check if an object is present.
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
                //visited.insert(current);
                visited[current] = true;

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

/*
 std::vector<int> flattenVector(const std::vector<array<int,2>>& inputVector, const size_t& finalTruleSize) {
     std::vector<int> flattened(finalTruleSize);
     for (size_t i=0; i<finalTruleSize; i+=2) {
         flattened[i] = inputVector[i][0];
         flattened[i+1] = inputVector[i][1];
     }
     return flattened;
 }
*/

std::vector<int> flattenVector(const std::vector<array<int,2>>& inputVector) {
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
    size_t finalTruleSize = out.size();
    std::unordered_map<int, int> mapToNew;
    mapToNew.reserve(finalTruleSize);

    std::vector<int> sortedOut = out;
    std::sort(sortedOut.begin(), sortedOut.end());
    for (size_t i = 0; i < finalTruleSize; i++) {
        mapToNew[sortedOut[i]] = i + 1;
    }
    // Apply transformations and partition the resulting vector. Sort each pair.
    std::vector<array<vert,2>> partitioned(finalTruleSize/2);
    vert firstInPair;
    vert secondInPair;
    for (size_t i = 0; i <finalTruleSize; i += 2) {
        firstInPair = mapToNew[out[i]]; //Note, definitely out here. The point being is that we are mapping out -> (position in sortedOut);
        secondInPair = mapToNew[out[i + 1]];
        if (firstInPair < secondInPair) {
            partitioned[i/2] = {firstInPair,secondInPair};
        } else {
            partitioned[i/2] = {secondInPair,firstInPair};
        }
    }

    // Sort the transformed pairs by modulo operation, i.e. by colour. Break ties by whichever number is smaller.
    std::sort(partitioned.begin(), partitioned.end(),
            [](const array<int,2>& a, const array<int,2>& b) {
            if (a[0] % num_colours == b[0] % num_colours){
                return a[0]<b[0];
            } else{
                return (a[0] - 1) % num_colours < (b[0] - 1) % num_colours; // Subtraction necessary here to ensure that the colour 1 appears before the colour 3
            }

        } // comp: comparison function object (i.e. an object that satisfies the requirements of Compare) which returns true if the first argument is less than (i.e. is ordered before) the second. 
    );

    return flattenVector(partitioned);
} 


// When loading up the full list, do the following
std::vector<adjList> processFullVertex(const vector<vector<int>>& fullVertexTlists, int& totalCount) {
    const int size=fullVertexTlists[0].size();
    vector<adjList> result;
    result.reserve(fullVertexTlists.size());

    vector<array<int, 2>> row(size);

    array<int, 2> subvector;
    subvector[1] = 0;

    for (const auto& oneLine : fullVertexTlists) {

        for (int i =0; i< size; i++) {
            subvector[0] = oneLine[i] + totalCount;
            row[i] = subvector;
        }
        result.push_back(row);
    }

    totalCount += size;
    return result;
}


/************************** NOT IN USE***********************************/
template <typename Function, typename List, typename Result>
vector<pair<IndexTuple,Result>> mapTuples(Function&& function, const std::vector<List>& lists) {
    const size_t numLists = lists.size();
    std::vector<size_t> sizes(numLists);
    IndexTuple indices(numLists, 0);


    size_t totalTuples = 1;
    for (size_t i = 0; i < numLists; i++) {
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
    ostringstream oss;

    std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(oss, ","));
    //for (const auto& value : vec) { oss << value << ","; }
    return oss.str();
}



/*
groupedLeadingIndex = GroupBy[MapIndexed[{#2[[1]], #1[[1]]} &, glPowerOutputs[ho][[1]]], #[[2]] &];
(*Careful of of-by-one errors!*)
#[[2]] & /@ SortBy[List @@ # & /@ Flatten@Values[ Thread[(First /@ #) -> First[#][[1]]] & /@ groupedLeadingIndex], First]

We factor in the off-by-one of C++ counting already here.  */
const array<string,21> variousPowersOfN = {"", "N1", "N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15","N16","N17","N18","N19","N20"};

// Helper function to convert a vector to a string
string makeMMAIndexStringCoeffMap(const IndexTuple& indices) {
    // Vector assumed to be of minimal length 1.
    ostringstream ss;
    ss << (ccIntToCoeffIntCanonicaliser[0][indices[0]]);//Don't add one here, as the map already has the addition of one built in.
    for (size_t i=1; i<indices.size(); ++i) {
        ss <<","<< (ccIntToCoeffIntCanonicaliser[i][indices[i]]); // Likewise
    }
    
    return ss.str();
}

template <typename Function, typename List, typename Result>
unordered_map<transposeRuleAsKey,umap_indices_to_loopcounter,FastVectorHasher> mapTuplesAndGroup(Function&& function, const std::vector<List>& lists) {
    const size_t numLists = lists.size();
    std::vector<size_t> sizes(numLists);

    IndexTuple indices(numLists, 0);
    indexAsKey indicesMapped;

    unordered_map<transposeRuleAsKey,umap_indices_to_loopcounter,FastVectorHasher> groups;

    size_t totalTuples = 1;
    for (size_t i = 0; i < numLists; ++i) {
        sizes[i] = lists[i].size();
        totalTuples *= sizes[i];
    }
    std::vector<typename List::value_type> elements(numLists);

    Result resultNumloopsAndTrules;
    
    for (size_t i = 0; i < totalTuples; ++i) {
        for (size_t j = 0; j < numLists; ++j) {
            elements[j] = lists[j][indices[j]];
        }

        resultNumloopsAndTrules = function(elements); // This returns a number (numLoops) and a vector (the transpose rule that defines the vertex it has been reducd to)
        //indicesMapped =  indices;//makeMMAIndexStringCoeffMap(indices); // Map indices to reduce the coeff complexity.
        indicesMapped = makeMMAIndexStringCoeffMap(indices); // Map indices to reduce the coeff complexity.

        //groups[vectorToString(resultNumloopsAndTrules.second)][indicesMapped][resultNumloopsAndTrules.first]++; //Increment the number of entries for this coefficient at this loop count.
        groups[resultNumloopsAndTrules.second][indicesMapped][resultNumloopsAndTrules.first]++; //Increment the number of entries for this coefficient at this loop count.

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

    // This is the only time contractionRules is used. Adding in contraction rules on top of the extant.
    for (const auto& pair : contractionRules) { 
        int vertex1 = pair[0];
        int vertex2 = pair[1];
        totalAdjList[vertex1][1] = vertex2;
        totalAdjList[vertex2][1] = vertex1; //First entry is always populated.
    }
    //std::cout << "Elements: "; for (const auto& element :totalAdjList) { std::cout << "(" << element[0] << " " << element[1] << ") "; } std::cout << std::endl;

    vector<vector<int>> cycles, sequences;
    tie(cycles, sequences) = get_cycles_and_sequences_for_fixed_contraction(totalAdjList);

    /* cout << "Cycles:" << cycles.size() << endl;
    for (const auto& cycle : cycles) {
        for (int vertex : cycle) { cout << vertex << " "; }
        cout << " So colour is " << ((cycle.front()-1) % num_colours) + 1 << endl; }

    cout << "Sequences:" << endl;
    for (const auto& sequence : sequences) {
        for (int vertex : sequence) { cout << vertex << " ";};
        cout << endl; } */
    //std::array<int,num_colours> vectorOfColourLoopCounts = {}; // Create a vector of integers with length num_colours, initialised to zero by std constructor.

    /*for (const auto& cycle : cycles) {
        vectorOfColourLoopCounts[((cycle.front()-1) % num_colours)+ 1 - 1]++; // Subtract one to get 0-indexed object. The -1 +1 otherwise is related to the mod, to ensure that we get colour=3 rather than colour =0 when colour=num_colours;
        //int whichColourInt =((cycle.front()-1) % num_colours) + 1;
        //std::cout << whichColourInt <<  " " << vectorOfColourLoopCounts[whichColourInt] << std::endl;
    }*/
    
    std::vector<int> transposeFormToBe;
    transposeFormToBe.reserve(sequences.size()* 2);
for (const auto& sequence : sequences) {
            transposeFormToBe.emplace_back(sequence.front());
            transposeFormToBe.emplace_back(sequence.back());
    }
    
    transposeFormToBe = getTransposeForm(transposeFormToBe);
    //Use this to force segfault: transposeFormToBe[10000*currentPos] =10 ;
    return make_pair(cycles.size(), transposeFormToBe);
}


//These functions are only executed O(1) times, so no need to optimise.
vector<vector<int>> mapCCstrToCoeffIntCanonicaliser(const string& input) {
 vector<vector<int>> result;
 for (const auto& ch : input) {
     auto it =ccMapToCoeffIntCanonicaliser.find(ch);
     if (it !=ccMapToCoeffIntCanonicaliser.end()) {
         result.push_back(it->second);
     }
 }
 return result;
}
vector<unordered_map<string,string>> mapCCstrToMapCoeffCanonicalIntToStr(const string& input) {
vector<unordered_map<string,string>> result;
 for (const auto& ch : input) {
     auto it =coeffCanonicalIntToStr.find(ch);
     if (it !=coeffCanonicalIntToStr.end()) {
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

/*
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
    return ss.str();
}

map<string,vector<string>> GroupBy22(const vector<indexedResult>&indexedResults) {
    map<string,vector<string>> groups;
    // Group pairs based on the second element (string)
    for (const auto& pair : indexedResults) {
        string key = vectorToString(pair.second.second);
        groups[key].push_back(makeMMAStringOutput(pair.first, pair.second.first));
    }

    return groups;
}*/


string mapIndexStringToActualVariables(const string& inputString) {
    istringstream iss(inputString);
    string token;
    string result;
    int whichCCOn =0; //Start on the first vertex.
    while (getline(iss, token, ',')) { 
        result +=ccIntToMapCoeffCanonicalIntToStr[whichCCOn].at(token);//https://stackoverflow.com/questions/13354394/reading-object-from-const-unordered-map
        whichCCOn++;
    }
    return result;
}

int main(int argc, char* argv[]) {

    if (argc < 3) {
        cout << "Usage: ./tensorOutputs <contraction 1-4-5-8...> <vertices e.g. gghl> <optional outputfile: otherwise automatic>" << endl;
        return 1;
    }
    contractionRules = parseStringContractionRules(argv[1]);
    ccIntToCoeffIntCanonicaliser = mapCCstrToCoeffIntCanonicaliser(argv[2]);
    ccIntToMapCoeffCanonicalIntToStr = mapCCstrToMapCoeffCanonicalIntToStr(argv[2]);

    // The compile time constant vertToVertAdjList is a map indexed by {'g', 'l', 'h'}, with values consisting of a list of adjacency lists.
    // These are generated by Mathematica. 
    // Each adjacency list would be of the form {1->{2}, 2->{1}, 3->{8}, 8->{3}, ...} Therefore INSTEAD we simply write it as an array arr = {2,1,8, ...}.
    // Thus if I want to know the SINGLE vertex that is attached to vertex vert (an integer) I just query arr[vert-1].
    vector<vector<vector<int>>> couplingConstantAdjLists;
    string whichVerticesString = argv[2];
    num_coupling_constants= 0;
    for (const auto& vertexChar : whichVerticesString) {
        auto it = vertToVertAdjList.find(vertexChar);
        if (it !=  vertToVertAdjList.end()) {
             couplingConstantAdjLists.push_back(it->second);
             num_coupling_constants++;
        } else {
            cout<< "Could not load vertex "<< std::string(vertexChar,1)<< ". Exiting 1"<<endl;
            return 1;
        }
    }

    int totalCount = 0;  // Running total of totalCount
    vector<vector<adjList>> allNumbers;
    
     for (const auto&  couplingCAdjList : couplingConstantAdjLists) {
        std::vector<adjList> output = processFullVertex(couplingCAdjList, totalCount);
         allNumbers.push_back(output);
         /*cout << "Total number of integers in " << filename << ": " << totalCount << endl;
        for (const auto& list : output) { for (const auto& subvector : list) { std::cout << "(" << subvector[0] << ", " << subvector[1] << ") "; } std::cout << std::endl; }*/
     }

    //cout<<"Total size: " << totalCount<<endl;
    //adjList totalAdjList(totalCount+1); //Add one to size.
    totalAdjList.resize(totalCount+1);
    totalAdjList[0] ={0,0};
    //int count = 0; mapTuples([&count, &totalAdjList, &totalCount](const std::vector<adjList>& elements) { count++; rS<<analyseListOfAdjLists(elements,totalAdjList,totalCount);  },allNumbers); std::cout << "Total tuples: " << count << std::endl;
    //omp_set_num_threads(8); //vector<indexedResult> allTuplesRun = mapTuples<decltype(analyseListOfAdjLists), vector<adjList>,analysisResultType>(analyseListOfAdjLists, allNumbers); 

    //Trying to get it work without having to explicitly spell out the template.
    //auto testLambdaFunction =[](const vector<adjList>& element){return analyseListOfAdjLists(element);};
    //vector<string> allTuplesRun = mapTuples<decltype(analyseListOfAdjLists), vector<adjList>(testLambdaFunction, allNumbers); 
    //std::cout << "Total tuples: " << allTuplesRun.size() << std::endl;
    //map<string,vector<string>> gathered = mapTuplesAndGroup<decltype(analyseListOfAdjLists), vector<adjList>,analysisResultType>(analyseListOfAdjLists, allNumbers); 
    
    //unordered_map<transposeRuleAsKey,umap_indices_to_loopcounter> gathered = mapTuplesAndGroup<decltype(analyseListOfAdjLists), vector<adjList>,analysisResultType>(analyseListOfAdjLists, allNumbers);
    unordered_map<transposeRuleAsKey,umap_indices_to_loopcounter,FastVectorHasher> gathered = mapTuplesAndGroup<decltype(analyseListOfAdjLists), vector<adjList>,analysisResultType>(analyseListOfAdjLists, allNumbers);

    /********************************/
    /*** Now perform output */
    string outputFilename;
    if (argc==4){ // i.e. args: exec, contraction, vertices, outputfile
        outputFilename = string(argv[3]);
    } else { // Create own outputfile
        outputFilename = "/tmp/fastOutput-"+string(argv[2])+"-c"+string(argv[1])+".out";
    } 
    ofstream outputFile(outputFilename);
    if (!outputFile) { cout << "Error opening output file: " << outputFilename << endl; return 1; }

    // Create an ostringstream to store all results for this execution. This is output at the end.
    std::ostringstream rS;

    rS<< "tensJustNums[{";
    for (const auto& structAndOverallCoeff: gathered) {
        // structAndOverallCoeff is effectively GroupBy[list, transposeStructure]
        // structAndOverallCoeff.first is the transpose rule; structAndOverallCoeff.second is the map of coefficients to {numLoops->occurence, ..}
        // i.e. now we iterate over{ coeff->{0->5, 1->3, 2->10, 3-> ...}, coeff2->{...}}
            rS << "{";
            for (const auto& coeffAndNFactorsPair: structAndOverallCoeff.second) {
                    // Print c[...] 

                    //rS << "c["<< makeMMAIndexStringCoeffMap(coeffAndNFactorsPair.first) ;rS << "](";
                    rS << mapIndexStringToActualVariables(coeffAndNFactorsPair.first);
                    rS << "(";
                    // Print (a N^0 + b N^3 + c N^4 + ...)+
                    for (const auto& loopCountAndValue :   coeffAndNFactorsPair.second){
                        rS << loopCountAndValue.second << variousPowersOfN[loopCountAndValue.first] << "+";
                    }
                    rS.seekp(-1,rS.cur);
                    rS << ")+";
                }
            rS.seekp(-1,rS.cur);
            //rS << ",{" << structAndOverallCoeff.first;
            rS << ",{" << vectorToString(structAndOverallCoeff.first);
        rS.seekp(-1,rS.cur);
        rS << "}},";
    }
    rS.seekp(-1,rS.cur); //Remove final comma
    rS<< "}, obj]";
    outputFile << rS.str();
    outputFile.close();

    return 0;
}
