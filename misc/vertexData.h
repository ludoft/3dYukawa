// vertexData.h

#ifndef VERTEXDATA_H
#define VERTEXDATA_H

#include <unordered_map>
#include <vector>
using namespace std;

extern const unordered_map<char, vector<vector<int>>> vertToVertAdjList;
extern const unordered_map<char, vector<int>> ccMapToCoeffIntCanonicaliser;
extern const unordered_map<char, unordered_map<string,string>> coeffCanonicalIntToStr;

#endif  // VERTEXDATA_H
