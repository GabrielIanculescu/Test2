package com.company;

public class PermCount {
    public int solution(int A[]){
        int N = A.length;
        int B[] = new int[N];
        for (int i=0; i<N; i++){
            B[i]=0;
        }
        for (int i=0; i<N; i++){
            if(A[i]<1 || A[i]>N){
                return 0;
            }
            else {
                if(B[A[i]-1] > 0){
                    System.out.println(i + "   " + A[i] + "  " + B[A[i]]);
                    return 0;
                }
                B[A[i]-1] = B[A[i]-1]+1;
            }
        }
        return 1;
    }
}
