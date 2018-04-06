package com.company;

public class BinaryGap {

    public int Solution(int N){
        String str = Integer.toBinaryString(N);
        int max =0;
        int gap =0;
        for (int index = 0; index < str.length(); index++) {
            char ch = str.charAt(index);
            if(ch == '0'){
                gap++;
            } else {
                if (max < gap){
                    max = gap;
                }
                gap =0;
            }
        }
        return  max;
    }
}
