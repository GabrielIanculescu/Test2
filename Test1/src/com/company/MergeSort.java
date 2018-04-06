package com.company;

public class MergeSort {
    public void sort(int array[]) {
        mergeSort(array, 0, array.length-1);
    }

    void mergeSort(int[] array, int low, int high){
        if(low < high){
            int middle = (low + high) / 2;
            mergeSort(array, low, middle);
            mergeSort(array, middle+1, high);
            merge(array, low, middle, high);
        }
    }

    void merge(int[] array, int low, int middle, int high){
        int[] tempArr = new int[array.length];
        for (int i = low; i <= high; i++) {
            tempArr[i] = array[i];
        }

        int left = low;
        int right = middle+1;
        int current = low;

        while (left <= middle && right <=high) {
            if(tempArr[left] <= tempArr[right]){
                array[current] = tempArr[left];
                left++;

            }else{
                array[current] = tempArr[right];
                right++;
            }
            current ++;
        }

        int remaining = middle - left;
        for (int i = 0; i <= remaining; i++) {
            array[current+i] = tempArr[left+ i];
        }
    }

}
