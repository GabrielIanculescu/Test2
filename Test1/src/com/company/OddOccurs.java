package com.company;

import java.util.HashMap;
import java.util.Map;

public class OddOccurs {
    public int solution(int A[], int N){
        Map<Integer, Integer> myMap = new HashMap<Integer, Integer>();
        for(int i=0; i<N; i++){
            Object o = new Integer(A[i]);
            if(myMap.containsKey(o)){
                Integer val = myMap.get(o);
                if(val == 1) {
                    myMap.remove(o);
                } else {
                    myMap.put(A[i], val+1);
                }

            } else {
                myMap.put(A[i], 1);
            }
        }

        int OddNo =0;
        for ( Map.Entry<Integer, Integer> entry : myMap.entrySet()) {
            Integer key = entry.getKey();
            Integer val = entry.getValue();
            if(val > 0){
                OddNo = key;
            }
        }
        return OddNo;
    }
}
