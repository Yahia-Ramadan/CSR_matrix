import heapq

def merge_sorted(lists: [[(int, float)]]):
    
    #for c++ implementation, follwing changes should be made:
        # precompute length of res
        # use std::tuple for tuples
        # use either std::make_heap or std::priority_queue
        # avoid pushes to res, use direct access by tracking index
        
    pq = []
    
    for i in range(len(lists)):
        heapq.heappush(pq, (lists[i][0], 0, i))
        
    res = []
    
    while pq:
        pair, list_i, id = heapq.heappop(pq)
        res.append(pair)
        list_i += 1
        if list_i < len(lists[id]):
            heapq.heappush(pq, (lists[id][list_i], list_i, id))
        
    return res



v1 = [4, 1, 5, 7]
c1 = [1, 2, 3, 7]
#i dont need the tuples, i just need an index

v2 = [1, 1, 9, 9]
c2 = [0, 1, 7, 8]

v3 = [4, 5, 9]
c3 = [4, 5, 6]