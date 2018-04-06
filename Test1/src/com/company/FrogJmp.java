package com.company;

public class FrogJmp {
    int solution(int X, int Y, int D){
        int diff = Y-X;
        int lastJmp = diff % D;
        System.out.println(lastJmp);
        int jmp = (Y-X)/D;
        if(lastJmp > 0) {
            jmp++;
        }
        return jmp;
    }
}
