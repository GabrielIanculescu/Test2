package com.company;

public class MissingNo {
    public int solution(int A[]){
        int N = A.length;
        int B[] = new int[N+1];
        for (int i=0; i<N; i++){
            if(A[i]>0 && A[i]<=N)
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
