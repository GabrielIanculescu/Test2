package com.company;

public class TapeEq {
    public int solution(int A[]){
        int N = A.length;
        int Sum = 0;
        for (int i=0; i<N; i++){
            Sum = Sum + A[i];
        }
        int min = Integer.MAX_VALUE;
        int val1 = 0;
        for (int i=0; i<N-1; i++){
            val1 = val1 + A[i];
            int val2 = Sum - val1;
            int diff = Math.abs(val1 -val2);
            if(diff < min){
                min = diff;
            }
        }
        return min;
    }
}
