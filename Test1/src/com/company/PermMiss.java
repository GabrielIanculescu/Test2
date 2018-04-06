package com.company;

public class PermMiss {
    public int solution(int A[]){
        int N = A.length;
        int B[] = new int[N+1];
        for (int i=0; i<N; i++){
            B[A[i]-1] = 1;
        }
        for (int i=0; i<=N; i++){
            if(B[i] == 0) {
                return i+1;
            }
        }
        return 0;
    }
}
