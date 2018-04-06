package com.company;

public class Rotation {
    public int[] solution(int A[],int K){

        int N = A.length;
        if(K == 0 || N == 0){
            return A;
        }
        K = K%N;
        int B[] = new int[N];
        for (int i=0; i<K; i++){
            B[i] = A[K+i-1];
            System.out.println(B[i]);
        }
        System.out.println("=====");
        for (int i=K; i<N; i++){
            B[i] = A[i-K];
            System.out.println(B[i]);
        }
        System.out.println("=====");
        for (int i=0; i<N; i++){
           System.out.println(B[i]);
        }
        return A;
    }
}
