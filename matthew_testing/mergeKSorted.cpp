#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <random>

// Assuming the sorted lists are std::vector<double>

using namespace std;




// Approach 1: Divide and Conquer (Iterative)
// Note: Could do recursion as well
// Time Complexity: O(nlogK)
vector<double> mergeKListsDC(vector<vector<double>> lists) {
    if(lists.empty()) {
        return vector<double>();
    }

    while(lists.size() > 1) { // runs logK times
        vector<vector<double>> mergedLists;

        for(int i=0; i<lists.size(); i+=2) { // merge every two lists together
            vector<double> list1 = lists[i];
            vector<double> list2 = (i+1) < lists.size() ? lists[i+1] : vector<double>();
            
            vector<double> merged;
            merged.resize(list1.size() + list2.size());
            merge(list1.begin(), list1.end(), list2.begin(), list2.end(), merged.begin()); 
            
            mergedLists.emplace_back(merged); // O(n) times
        }
        lists = mergedLists;
    }
    return lists[0];
}


// Approach 2: Heap
vector<double> mergeKListsHeap(vector<vector<double>> lists) {
    if(lists.empty()) {
        return vector<double>{};
    }

    struct HeapStruct { // need to store value, index in overall lists, and index in current list
        double value; // TODO modifications: instead of value, make it pair<int, double>
        // where pair.first (int) stores the column index, heap is sorted based on this value,
        // pair.second (double) stores the actual multiplication calculation
        

        // note if indices are the same, still increment the elementIndex++, except instead of 
        // pushing, alter result.back() by adding the value
        int listIndex;
        int elementIndex; 
        
    };

    vector<double> result;

    auto comp = [](const HeapStruct& a, const HeapStruct& b) {
        return a.value > b.value;
    };

    priority_queue<HeapStruct, vector<HeapStruct>, decltype(comp)> minHeap(comp);

    for(int i=0; i<lists.size(); i++) {     // initially populate with first elements of each list
        if(!lists[i].empty()) {
            minHeap.push({lists[i][0], i, 0});
        }
    }

    while(!minHeap.empty()) {
        HeapStruct min = minHeap.top();
        result.push_back(min.value);
        minHeap.pop();
        int nextIndex = min.elementIndex + 1;
        if(nextIndex < lists[min.listIndex].size()) {
            minHeap.push({lists[min.listIndex][nextIndex], min.listIndex, nextIndex});
        }
    }

    return result;
}


// Adoption for SpGEMM
vector<pair<int,double>> mergeKListsMatrix(vector<vector<pair<int, double>>> lists) {
    if(lists.empty()) {
        return vector<pair<int,double>>{};
    }

    struct HeapStruct {
        pair<int, double> value;
        int listIndex;
        int elementIndex; 
    };

    vector<pair<int, double>> result;

    auto comp = [](const HeapStruct& a, const HeapStruct& b) {
        return a.value.first > b.value.first;
    };

    priority_queue<HeapStruct, vector<HeapStruct>, decltype(comp)> minHeap(comp);

    for(int i=0; i<lists.size(); i++) {     // initially populate with first elements of each list
        if(!lists[i].empty()) {
            minHeap.push({lists[i][0], i, 0});
        }
    }

    while(!minHeap.empty()) {
        HeapStruct min = minHeap.top();
        if(!result.empty() && result.back().first == min.value.first) { // same columnIndex, must add together
            result.back().second += min.value.second;
        }
        else {
            result.push_back(min.value);
        }
        minHeap.pop();
        int nextIndex = min.elementIndex + 1;
        if(nextIndex < lists[min.listIndex].size()) {
            minHeap.push({lists[min.listIndex][nextIndex], min.listIndex, nextIndex});
        }
    }

    return result;
}



int main() {
    // int k = 10; // number of lists
    // int n = 20; // max number of elements per list
    // vector<vector<double>> lists(k);

    // for (int i=0; i<k; i++) {
    //     int size = rand() % n + 1;
    //     vector<double>& list = lists[i];
    //     for (int j = 0; j < size; ++j) {
    //         list.push_back(rand() % 1000 + 1);
    //     }
    //     std::sort(list.begin(), list.end());
    // }

    // vector<double> heapResult = mergeKListsHeap(lists);
    // vector<double> DCResult = mergeKListsDC(lists);

    // if(heapResult == DCResult) {
    //     cout << "Both results match each other!" << endl;
    // }
    // else {
    //     cout << "results do not match" << endl;
    // }


    // we have to combine [ _ 2 __ 4 1 ]
    //                    [ _ _ 10 _ 5 ]
    // --------------------------------------
    // expected result:   [ _ 2 10 4 6 ]  
    // ie [1 2] [2 10] [3 4] [4 6]

    vector<pair<int, double>> row1 = {{1,2}, {3,4}, {4,1}};
    vector<pair<int, double>> row2 = {{2,10}, {4,5}};

    vector<vector<pair<int, double>>> rows = {row1, row2};

    vector<pair<int, double>> result = mergeKListsMatrix(rows);
    cout << "Result:" << endl;
    for(pair<int, double> val: result) {
        cout << "[" << val.first << " " << val.second << "]" << " ";
    }
    cout << endl;

    return 0;
}