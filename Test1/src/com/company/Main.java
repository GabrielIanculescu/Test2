package com.company;
import java.io.PrintWriter;
import java.util.ArrayList;

import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import java.util.*;

import java.awt.*;

public class Main {

    public static class Node {
        String parent;
        String name;
        ArrayList<String> children;
    }

    public static void main(String[] args) {
        int [] array = {4,3,4,4,4,2};
        int [] C = {1, 2};
        int [] D = {5, 4, 4, 3, 3, 4};

       // String  line;
       // Scanner sc = new Scanner(System.in);
       // line = sc.nextLine();
       // int n = Integer.parseInt(line);
        HashMap<String, Node> nodes = new HashMap<String, Node>();
        HashSet<String> nodesFree = new HashSet<String>();

        int n = 5;
                String[] l = new String[5];
        l[0] = "WINDOWS SYSTEM32";
        l[1] = "FOO BAR";
        l[2] = "C: WINDOWS";
        l[3] = "C: FOO";
        l[4] = "BAR XYZZY";

        String[][] s = new String[2][n];
        for (int i=0; i<n; i++){
            //line = sc.nextLine();
            String[] p = l[i].split(" ");
            s[0][i]=p[0];
            s[1][i]=p[1];

            Node nn;
            if(!nodes.containsKey(p[0])) {
                nn = new Node();
                nn.parent = null;
                nn.name = p[0];
                nn.children = new ArrayList<String>();
                nodes.put(p[0], nn);
                if(!p[0].equals("C:")) {
                    nodesFree.add(p[0]);
                }
            } else {
                nn = nodes.get(p[0]);
            }
            Node n1 = new Node();
            n1.parent = nn.name;
            n1.name = p[1];
            n1.children = new ArrayList<String>();
            nn.children.add(p[0]);
            nodes.put(p[1], n1);
        }

        Collection<Node> collectionValues = nodes.values();
        Object[] nnn = nodesFree.toArray();

        do {
            for (int i = 0; i < nnn.length; i++) {
                String ss = (String) nnn[i];
                Node nod = nodes.get(ss);
                if (nod.parent == null) {
                    for (Node nod2 : collectionValues) {
                        if (nod2.children.contains(nod.name)) {
                            nod.parent = nod2.name;
                            nodesFree.remove(ss);
                            break;
                        }
                    }
                } else {
                    nodesFree.remove(ss);
                }
            }
        } while(!nodesFree.isEmpty());

        String maxPath ="C:";
        for(Node nod: collectionValues){
            if(nod.children.isEmpty()){
                String path = nod.name;
                Node nn = nod;
                while(nn.parent!=null){
                    nn = nodes.get(nn.parent);
                    path = nn.name + "\\" + path;
                }
                if(maxPath.length()<path.length()){
                    maxPath = path;
                }
            }
        }

        System.out.println(maxPath);
    }

    public static int CountDiv(int A, int B, int K){
        if(A==0){
            if(B==0) {
                return 1;
            } else {
                return 1+ B/K;
            }
        }

        if(A == B){
            if(A % K ==0){
                return 1;
            } else {
                return 0;
            }
        }

        int n1 = (A-1)/K;
        int n2 = B/K;
        return n2-n1;
    }

    public static int FrogRiverOne(int X, int A[]) {
        int N = A.length;
        int [] array = new int[X];

        int n = 0;
        for (int i = 0; i < N; i++) {
            if (A[i] <= X) {
                if(array[A[i]-1] ==0){
                    array[A[i]-1] =1;
                    n++;
                    if(n == X){
                        return i;
                    }
                }
            }
        }
        return -1;
    }

    public static int PassingCars(int A[]){
        int N = A.length;
        int nEast=0;
        int pairs=0;
        for (int i=0; i<N; i++){
            if(A[i] == 0) {
                nEast++;
            } else {
                pairs = pairs + nEast;
                if(pairs > 1000000000){
                    return -1;
                }
            }
        }
        return pairs;
    }

    public static int Triangle(int A[]) {
        int N = A.length;
        HeapSort ms = new HeapSort();
        ms.sort(A);

        int n = 1;
        for (int i=2; i<N; i++) {
            long a1 = (long)A[i-2];
            long a2 = (long)A[i-1];
            long a3 = (long)A[i];
            if(a1+a2 > a3){
                return 1;
            }
        }
        return 0;
    }

    public static int MaxProductOfThree(int A[]) {
        int N = A.length;
        if(N==3) {
            return A[0]*A[1]*A[2];
        }
        HeapSort ms = new HeapSort();
        ms.sort(A);

        boolean isPositive = false;
        for (int i=0; i<N; i++) {
            if(A[i] >= 0){
                isPositive = true;
                break;
            }
        }

        if(isPositive){
            int max1=0;
            int max2 = 0;
            if(A[0]<0 && A[1]<0) {
                max1 = A[0] * A[1] * A[N-1];
                max2 = A[N-3] * A[N-2] * A[N-1];
                return Math.max(max1, max2);
            }
        }

        return  A[N-3] * A[N-2] * A[N-1];

    }

    public static int Brackets(String s) {
        Stack<Character> stack = new Stack<Character>();
        for(int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case ')':
                    if(stack.isEmpty()){
                        return 0;
                    }
                    if(stack.pop() != '('){
                        return 0;
                    }
                    break;
                case '(':
                    stack.push(c);
                    break;
                case ']':
                    if(stack.isEmpty()){
                        return 0;
                    }
                    if(stack.pop() != '['){
                        return 0;
                    }
                    break;
                case '[':
                    stack.push(c);
                    break;
                case '}':
                    if(stack.isEmpty()){
                        return 0;
                    }
                    if(stack.pop() != '{'){
                        return 0;
                    }
                    break;
                case '{':
                    stack.push(c);
                    break;
            }
        }
        if(stack.isEmpty()){
            return 1;
        }
        return 0;
    }

    public static int Nested(String s) {
        int n=0;
        for(int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if(c == '(') {
                n++;
            } else {
                if(n==0){
                    return 0;
                }
                n--;
            }
        }
        if(n==0){
            return 1;
        }
        return 0;
    }

    public static int EquiLeader(int A[]){
        int N = A.length;
        int n=0;
        Stack<Integer> stack = new Stack<Integer>();
        for (int i=0; i<N; i++) {
            if(stack.isEmpty()){
                stack.push(A[i]);
            } else {
                if(stack.peek() != A[i]){
                    stack.pop();
                } else {
                    stack.push((A[i]));
                }
            }

            if(!stack.isEmpty()){
                int lead = stack.peek();
                int nrLead=0;
                for(int j=0; j<=i; j++) {
                    if(A[j]==lead){
                        nrLead++;
                    }
                }
                if(nrLead > (i+1)/2){
                    nrLead = 0;
                    for(int j=i+1; j<N; j++) {
                        if(A[j]==lead){
                            nrLead++;
                        }
                    }
                    if(nrLead > (N-i-1)/2){
                        n++;
                    }
                }
            }
        }
        return n;
    }

    public static int Dominator(int A[]){
        int N = A.length;
        int n=0;
        Stack<Integer> stack = new Stack<Integer>();
        for (int i=0; i<N; i++) {
            if (stack.isEmpty()) {
                stack.push(A[i]);
            } else {
                if (stack.peek() != A[i]) {
                    stack.pop();
                } else {
                    stack.push((A[i]));
                }
            }
        }

        if(stack.isEmpty()) {
            return -1;
        } else {
            int nrLead=0;
            int idx = 0;
            int lead = stack.peek();
            for (int i=0; i<N; i++) {
                if (A[i] == lead) {
                    nrLead++;
                    idx = i;
                }
            }
            if(nrLead > N/2) {
                return idx;
            }
        }
        return -1;
    }

    public static int MaxSliceSum(int A[]) {
        int N = A.length;
        int maxEnd = A[0];
        int maxSlice = A[0];

        for (int i=1; i<N; i++) {
            maxEnd = Math.max(A[i], maxEnd + A[i]);
            maxSlice = Math.max(maxSlice, maxEnd);
        }
        return maxSlice;
    }

    public static int MaxDoubleSliceSum(int A[]) {
        int N = A.length;
        int[] K1 = new int[N];
        int[] K2 = new int[N];

        for (int i = 1; i < N - 1; i++) {
            K1[i] = Math.max(K1[i - 1] + A[i], 0);
        }
        for (int i = N - 2; i > 0; i--) {
            K2[i] = Math.max(K2[i + 1] + A[i], 0);
        }

        int max = 0;

        for (int i = 1; i < N - 1; i++) {
            max = Math.max(max, K1[i - 1] + K2[i + 1]);
        }

        return max;
    }

    public static int maxProfit(int[] A) {
        if(A.length == 1 || A.length == 0){
            return 0;
        }

        int maxSoFar = 0;
        int maxEndingHere = 0;
        int minPrice = A[0];

        for(int i = 1; i < A.length; i++){
            maxEndingHere = Math.max(0, A[i] - minPrice);
            minPrice = Math.min(minPrice, A[i]);
            maxSoFar = Math.max(maxEndingHere, maxSoFar);
        }

        return maxSoFar;
    }

    public static int[] maxCounters(int N, int[] A) {
        int[] counters = new int[N];
        int max = 0;
        int absMax = 0;

        for (int i=0; i<A.length; i++) {
            if (A[i] == N+1) {
                if ((i < A.length-1 && A[i+1] != N+1) || (i==A.length-1)) {
                    absMax += max;
                    max = 0;
                    counters = new int[N];
                }
            } else {
                counters[A[i]-1]++;
                if (max < counters[A[i]-1])
                    max = counters[A[i]-1];
            }
        }

        for (int i=0; i<counters.length; i++)
            counters[i] += absMax;

        return counters;
    }

    int MinAvgTwoSlice(int A[]) {
        int N = A.length;
        int minVal = Integer.MAX_VALUE;
        int minPos = 0;
        int i;
        int val;
        for (i=0; i<N-1; i++) {
            val = (A[i] + A[i+1])*3;
            if (val < minVal) {
                minVal = val;
                minPos = i;
            }
        }
        for (i=0; i<N-2; i++) {
            val = (A[i] + A[i+1] + A[i+2])*2;
            if (val < minVal) {
                minVal = val;
                minPos = i;
            }
        }
        return minPos;
    }

    public static int[] GenomicRangeQuery(String S, int[] P, int[] Q) {
        int len = S.length();
        int[][] arr = new int[len][4];
        int[] result = new int[P.length];

        for(int i = 0; i < len; i++){
            char c = S.charAt(i);
            if(c == 'A') arr[i][0] = 1;
            if(c == 'C') arr[i][1] = 1;
            if(c == 'G') arr[i][2] = 1;
            if(c == 'T') arr[i][3] = 1;
        }
        // compute prefixes
        for(int i = 1; i < len; i++){
            for(int j = 0; j < 4; j++){
                arr[i][j] += arr[i-1][j];
            }
        }

        for(int i = 0; i < P.length; i++){
            int x = P[i];
            int y = Q[i];

            for(int a = 0; a < 4; a++){
                int sub = 0;
                if(x-1 >= 0) sub = arr[x-1][a];
                if(arr[y][a] - sub > 0){
                    result[i] = a+1;
                    break;
                }
            }

        }
        return result;
    }

    public static int NumberOfDiscIntersection(int[] A){
        int l = A.length;

        long[] A1 = new long[l];
        long[] A2 = new long[l];

        for(int i = 0; i < l; i++){
            A1[i] = (long)A[i] + i;
            A2[i] = -((long)A[i]-i);
        }

        Arrays.sort(A1);
        Arrays.sort(A2);

        long cnt = 0;

        for(int i = A.length - 1; i >= 0; i--){
            int pos = Arrays.binarySearch(A2, A1[i]);
            if(pos >= 0){
                while(pos < A.length && A2[pos] == A1[i]){
                    pos++;
                }
                cnt += pos;
            } else{ // element not there
                int insertionPoint = -(pos+1);
                cnt += insertionPoint;
            }

        }

        long sub = (long)l*((long)l+1)/2;
        cnt = cnt - sub;

        if(cnt > 1e7) return -1;

        return (int)cnt;
    }

    public int StoneWall(int[] H) {
        int len = H.length;
        Stack<Integer> stack = new Stack<Integer>();
        int blockCounter = 0;

        for (int i = 0; i < len; ++i) {
            int element = H[i];
            if (stack.isEmpty()) {
                stack.push(element);
                ++blockCounter;
            } else {
                while (!stack.isEmpty() && stack.peek() > element) {
                    stack.pop();
                }
                if (!stack.isEmpty() && stack.peek() == element) {
                    continue;
                } else {
                    stack.push(element);
                    ++blockCounter;
                }
            }
        }

        return blockCounter;
    }

    public static int fish(int[] A, int[] B) {
        Stack<Integer> s = new Stack<Integer>();

        for(int i = 0; i < A.length; i++){
            int size 	= A[i];
            int dir 	= B[i];
            if(s.empty()){
                s.push(i);
            }
            else{
                while(!s.empty() && dir - B[s.peek()] == -1 && A[s.peek()] < size){
                    s.pop();
                }
                if(!s.empty()){
                    if(dir - B[s.peek()] != -1) s.push(i);
                }
                else{
                    s.push(i);
                }
            }
        }
        return s.size();
    }

    static int equiLeader(int[] a){
        int len = a.length;
        int equi_leaders = 0;

        // first, compute the leader
        int leader = a[0];
        int ctr = 1;

        for(int i = 1; i < a.length; i++){
            if(a[i] == leader) ctr++;
            else ctr--;
            if(ctr == 0){
                ctr = 1;
                leader = a[i];
            }
        }

        // check if it's a leader?
        int total = 0;
        for(int i : a){
            if(i == leader) total++;
        }

        if(total <= (len/2)) return 0; // impossible

        int ldrCount = 0;
        for(int i = 0; i < a.length; i++){
            if(a[i] == leader) ldrCount++;
            int leadersInRightPart = (total - ldrCount);
            if(ldrCount > (i+1)/2   &&   leadersInRightPart > (len-i-1)/2){
                equi_leaders++;
            }
        }

        return equi_leaders;
    }

    public static int MinPerimeterRectangle(int N) {
        int min =2*(N+1);

        for (int i = 1; i <=Math.sqrt(N); ++i) {
            if(N%i==0) {
                int l = N/i;
                int perim = 2*(i+l);
                if(perim < min) {
                    min = perim;
                }
            }
        }
        return min;
    }

    public static int CountFactors(int N) {
        int count =0;

        for (int i = 1; i <Math.sqrt(N); ++i) {
            if(N%i==0) {
               count +=2;
            }
        }
        if(N%Math.sqrt(N)==0){
            count++;
        }
        return count;
    }

    public int Flags(int[] A) {
        if(A.length==1)  return 0;
        ArrayList<Integer> al=new ArrayList<Integer>();
        for(int i=1;i<A.length-1;i++)
        {
            if(A[i]>A[i-1] && A[i]>A[i+1])
            {
                al.add(i);
            }
        }
        int s=al.size();
        if(s==1) return 1;
        if(s==0)return 0;
        s=(int) Math.ceil(Math.sqrt(A.length));
        while(s>=0) {
            int lp=al.get(0);
            int c=1;
            for(int i=1;i<al.size();i++)
            {
                int d=Math.abs(al.get(i)-lp);
                if(d>=s)
                {
                    lp=al.get(i);
                    c++;
                    if(c==s)
                        return c;
                }
            }
            s--;}
        return 0;
    }

    public int[] CountSemiprimes(int N, int[] P, int[] Q){
        boolean[] sieve = new boolean[N+1];
        Arrays.fill(sieve, Boolean.TRUE);
        sieve[0] = false;
        sieve[1] = false;
        for(int i=2; i<(int)(Math.sqrt(N)+1); i++) {
            if(sieve[i]) {
                for(int j=i+i; j<N; j+=i) {
                    sieve[j] = false;
                }
            }
        }
        List<Integer> primes = new ArrayList<Integer>();
        for(int i=0; i<sieve.length; i++) {
            if(sieve[i]) {
                primes.add(i);
            }
        }
        int[] semiprimes = new int[N+1];
        int[] sp = new int[N+1];
        long semiprime = 0;
        for(int i=0; i<primes.size(); i++) {
            for(int j=i; j<primes.size(); j++) {
                semiprime = (long)primes.get(i)*(long)primes.get(j);
                if(semiprime>N) {
                    break;
                }
                semiprimes[(int)semiprime] = 1;
                sp[(int)semiprime] = 1;
            }
        }
        for(int i=1; i<semiprimes.length; i++) {
            semiprimes[i] += semiprimes[i-1];
        }
        int[] count = new int[Q.length];
        for(int i=0; i<Q.length; i++) {
            count[i] = semiprimes[Q[i]] - semiprimes[P[i]] + sp[P[i]];
        }
        return count;
    }
    public int Peaks(int []A){
        int N = A.length;

        // Find all the peaks
        ArrayList<Integer> peaks = new ArrayList<Integer>();
        for(int i = 1; i < N-1; i++){
            if(A[i] > A[i-1] && A[i] > A[i+1]) peaks.add(i);
        }

        for(int size = 1; size <= N; size++){
            if(N % size != 0) continue; // skip if non-divisible
            int find = 0;
            int groups = N/size;
            boolean ok = true;
            // Find whether every group has a peak
            for(int peakIdx : peaks){
                if(peakIdx/size > find){
                    ok = false;
                    break;
                }
                if(peakIdx/size == find) find++;
            }
            if(find != groups) ok = false;
            if(ok) return groups;
        }
        return 0;
    }

    public int[] CountNondivisible(int [] A){
        int nIntCount = A.length;
        int nMaxInt = nIntCount + nIntCount;
        int[] anIntCounts = new int[nMaxInt + 1];
        int[] anDivisorCounts = new int[nMaxInt + 1];
        int[] anNotDivCounts = new int[nIntCount];
        for (int i=0; i<nIntCount;i++){
            anIntCounts[A[i]]++;
        }
        for (int i = 0; i <= nMaxInt; i++)
            if (anIntCounts[i] > 0)
                for (int im = i; im <= nMaxInt; im += i)

                    anDivisorCounts[im] += anIntCounts[i];

        for (int i = 0; i < nIntCount; i++)
            anNotDivCounts[i] = nIntCount - anDivisorCounts[A[i]];

        return anNotDivCounts;
    }

    public int ChocolatesByNumbers(int N, int M) {
        return N / gcd(N, M);
    }

    public int CommonPrimeDivisor(int A[], int B[]) {
        int cnt = 0;
        int Z = A.length;
        int i;
        for (i = 0; i < Z; i++){
            int gcd_ab = gcd(A[i], B[i]);
            if (check(A[i], gcd_ab) && check(B[i], gcd_ab)){
                cnt++;
            }
        }

        return cnt;
    }

    boolean check(int a, int gcd_ab)
    {
        int rest = a / gcd_ab;

        while (gcd_ab % rest != 0){

            int gcd_tmp = gcd(gcd_ab, rest);
            if (gcd_tmp == 1){
                return false;
            }
            rest /= gcd_tmp;
        }
        return true;
    }

    private int gcd(int n, int m) {
        int r = n % m;
        while (r != 0) {
            n = m;
            m = r;
            r = n % m;
        }
        return m;
    }

    public int fibFrog(int [] A){
        ArrayList<Integer> fibonacci = new ArrayList<>();
        fibonacci.add(0); // note: f(0) = 0 (as in the quesion)
        fibonacci.add(1);
        // note: using "while" is better than "for" (avoid errors)
        while(true){
            int temp1 = fibonacci.get( fibonacci.size()-1 );
            int temp2 = fibonacci.get( fibonacci.size()-2 );
            fibonacci.add( temp1 + temp2 );

            // if already bigger than length, then break;
            if(temp1 + temp2 > A.length){
                break;
            }
        }

        // reverse "List": from big to small
        Collections.reverse(fibonacci);

        // use "queue" with "point"
        // point(x,y) = point("position", "number of steps")
        ArrayList<Point> queue = new ArrayList<>();
        queue.add( new Point(-1, 0) ); // position:-1, steps:0

        // index: the current index for queue element
        int index=0;
        while(true){
            // cannot take element from queue anymore
            if(index == queue.size() ){
                return -1;
            }

            // take element from queue
            Point current = queue.get(index);

            // from big to small
            for(Integer n: fibonacci){
                int nextPosition = current.x + n;

                // case 1: "reach the other side"
                if(nextPosition == A.length){
                    // return the number of steps
                    return current.y + 1;
                }

                // case 2: "cannot jump"
                // note: nextPosition < 0 (overflow, be careful)
                else if( (nextPosition > A.length) || (nextPosition < 0)|| (A[nextPosition]==0) ){
                    // note: do nothing
                }

                // case 3: "can jump" (othe cases)
                else{
                    // jump to next position, and step+1
                    Point temp = new Point(nextPosition, current.y + 1);
                    // add to queue
                    queue.add(temp);

                    A[nextPosition] = 0; // key point: for high performance~!!
                }
            }

            index++; // take "next element" from queue
        }
    }

    public int[] Ladders(int[] A, int[] B) {
        int L = A.length;

        int Fib[] = new int[L + 1];
        Fib[0] = 1;
        Fib[1] = 1;
        int p = (1 << 30) - 1;
        for (int i = 2; i < L+1; ++i){
            Fib[i] = (Fib[i-1] + Fib[i-2]) & p;
        }
        int result[] = new int[L];
        for (int i = 0; i < L; ++i){
            p = (1 << B[i]) - 1;
            result[i] = Fib[A[i]] & p;
        }

        return result;
    }

    public int MinMaxDivision(int K, int M, int[] A) {
        int n = A.length;
        int sum = 0, max = 0;
        for (int i = 0; i < n; i++) {
            sum += A[i];
            max = Math.max(max, A[i]);
        }
        int left = max, right = sum;
        while (left <= right) {
            int mid = (left + right) >> 1;
            int intervals = countIntervals(A, mid);
            if (intervals > K) {
                left = mid + 1;
            } else right = mid - 1;
        }
        return left;
    }

    private int countIntervals(int[] A, int target) {
        int sum = 0, count = 0;
        for (int i = 0; i < A.length; i++) {
            sum += A[i];
            if (sum > target) {
                count++;
                sum = A[i];
            }
        }
        return count + (sum > 0 ? 1 : 0);
    }

    int NailingPlanks(int A[], int B[], int C[]) {
        int N = A.length;
        int M = C.length;
        // two dimension array to save the original index of array C
        int[][] sortedNail = new int[M][2];
        for (int i = 0; i < M; i++) {
            sortedNail[i][0] = C[i];
            sortedNail[i][1] = i;
        }
        // use the lambda expression to sort two dimension array
        Arrays.sort(sortedNail, (int x[], int y[]) -> x[0] - y[0]);
        int result = 0;
        for (int i = 0; i < N; i++) {//find the earlist position that can nail each plank, and the max value for all planks is the result
            result = getMinIndex(A[i], B[i], sortedNail, result);
            if (result == -1)
                return -1;
        }
        return result + 1;
    }

    public int getMinIndex(int startPlank, int endPlank, int[][] nail, int preIndex) {
        int min = 0;
        int max = nail.length - 1;
        int minIndex = -1;
        while (min <= max) {
            int mid = (min + max) / 2;
            if (nail[mid][0] < startPlank)
                min = mid + 1;
            else if (nail[mid][0] > endPlank)
                max = mid - 1;
            else {
                max = mid - 1;
                minIndex = mid;
            }
        }
        if (minIndex == -1)
            return -1;
        int minIndexOrigin = nail[minIndex][1];
        //find the smallest original position of nail that can nail the plank
        for (int i = minIndex; i < nail.length; i++) {
            if (nail[i][0] > endPlank)
                break;
            minIndexOrigin = Math.min(minIndexOrigin, nail[i][1]);
            // we need the maximal index of nails to nail every plank, so the
            // smaller index can be omitted
            if (minIndexOrigin <= preIndex)
                return preIndex;
        }
        return minIndexOrigin;
    }

    public int AbsDistinct( int[] A) {
        int n = A.length;
        int sum = 0;
        Set<Integer> set = new HashSet<Integer>();
        for (int i = 0; i < n; i++) {
            int x = Math.abs(A[i]);
            if(!set.contains(x)){
                sum++;
                set.add(x);
            }
        }
        return sum;
    }

    public int CountDistinctSlices(int[] A){
        // array to remember last positions of values
        int vMax = A[0];
        for (int i = 1; i < A.length; ++i) vMax = Math.max(vMax, A[i]);
        int[] vLastPos = new int[vMax + 1];
        for (int i = 0; i < vLastPos.length; ++i) vLastPos[i] = -1;

        // each element adds the same number of slices as is the length of current distinct slice
        int vSlices = 0, vNewStart = -1;
        for (int i = 0; i < A.length; ++i)
        {
            int vVal = A[i];
            int vPrevPos = vLastPos[vVal];
            vSlices += i - Math.max(vPrevPos, vNewStart);
            if (vSlices > 1000000000) return 1000000000;
            if (vPrevPos != -1) vNewStart = Math.max(vNewStart, vPrevPos); // actual start of distinct slice
            vLastPos[vVal] = i;
        }

        return vSlices;
    }

    public int CountTriangles(int[] A) {
        Arrays.sort(A);
        int ans = 0, n = A.length;
        for (int i = 0; i < n - 2; i++) {
            int k = 0;  // k is init here
            for (int j = i + 1; j < n - 1; j++) {
                while (k < n && A[i] + A[j] > A[k]) {
                    k++;
                }
                ans += k - j - 1;
            }
        }
        return ans;
    }

    public int MinAbssumOfTwo(int[] A) {
        int N = A.length;
        Arrays.sort(A);
        int tail = 0;
        int head = N - 1;
        int minAbsSum = Math.abs(A[tail] + A[head]);
        while (tail <= head) {
            int currentSum = A[tail] + A[head];
            minAbsSum = Math.min(minAbsSum, Math.abs(currentSum));
            // If the sum has become
            // positive, we should know that the head can be moved left
            if (currentSum <= 0)
                tail++;
            else
                head--;
        }
        return minAbsSum;
    }

    public int TieRopes(int K, int[] A) {
        // write your code in Java SE 8
        // use the greedy algorithm to solve this problem
        // because we can only tie the adjacent ropes, so for each rope we can
        // count the length of the consecutive ropes and figure out if there is
        // a tied rope that the length is greater than or equal to K.
        int N = A.length;
        int result = 0;
        int length = 0;
        for (int rope : A) {
            length += rope;
            if (length >= K) {
                result++;
                length = 0;
            }
        }
        return result;
    }

    public int MaxNonOverlappingSegmenta(int[] A, int[] B) {
        int ans = 0, pre = -1;
        for (int i = 0; i < A.length; i++) {
            if (A[i] > pre) {
                ans++;
                pre = B[i];
            }
        }
        return ans;
    }

    public int NumberSolitaire(int[] A) {
        int[] store = new int[A.length];
        store[0] = A[0];
        for (int i = 1; i < A.length; i++) {
            store[i] = store[i - 1];
            for (int minus = 2; minus <= 6; minus++) {
                if (i >= minus) {
                    store[i] = Math.max(store[i], store[i - minus]);
                } else {
                    break;
                }
            }
            store[i] += A[i];
        }
        return store[A.length - 1];
    }

    public int MinAbsSum(int[] A) {
        if (A.length == 0) return 0;
        if (A.length == 1) return A[0];
        int sum = 0;
        for (int i=0;i<A.length;i++){
            sum += Math.abs(A[i]);
        }
        int[] indices = new int[A.length];
        indices[0] = 0;
        int half = sum/2;
        int localSum = Math.abs(A[0]);
        int minLocalSum = Integer.MAX_VALUE;
        int placeIndex = 1;
        for (int i=1;i<A.length;i++){
            if (localSum<half){
                if (Math.abs(2*minLocalSum-sum) > Math.abs(2*localSum - sum))
                    minLocalSum = localSum;
                localSum += Math.abs(A[i]);
                indices[placeIndex++] = i;
            }else{
                if (localSum == half)
                    return Math.abs(2*half - sum);

                if (Math.abs(2*minLocalSum-sum) > Math.abs(2*localSum - sum))
                    minLocalSum = localSum;
                if (placeIndex > 1) {
                    localSum -= Math.abs(A[indices[placeIndex--]]);
                    i = indices[placeIndex];
                }
            }
        }
        return (Math.abs(2*minLocalSum - sum));
    }

    public static int LongestPass(String S) {
        // write your code in Java SE 8
        String[] words = S.split(" ");
        int max =-1;
        for(String s : words) {
            int nrLetters=0;
            int nrDigits=0;
            int nr=0;
            for(int i = 0; i < s.length(); i++) {
                char c = s.charAt(i);
                nr++;
                if (Character.isLetter(c)){
                    nrLetters++;
                }
                else if (Character.isDigit(c)) {
                    nrDigits++;
                }
                else {
                    nrDigits =0;
                    break;
                }
            }
            if(nrDigits%2 == 1 && nrLetters%2==0){
                if(nr>max){
                    max = nr;
                }
            }
        }
        return max;
    }

    public int Tenis(int P, int C) {
        // write your code in Java SE 8
        return Math.min(P/2, C);
    }

    public static int Sock(int K, int[] C, int[] D) {
        // write your code in Java SE 8
        int n=0;
        HashSet<Integer> hClean = new HashSet<Integer>();
        HashMap<Integer, Integer> hDirty = new HashMap<Integer, Integer>();

        for (int i=0;i<C.length;i++) {
            if(hClean.contains(C[i])){
                n++;
                hClean.remove(C[i]);
            } else {
                hClean.add(C[i]);
            }
        }
        if(K==0){
            return n;
        }

        int nWash = 0;
        for (int i=0;i<D.length;i++) {
            if(nWash==K){
                break;
            }
            if(hClean.contains(D[i])){
                n++;
                hClean.remove(D[i]);
                nWash++;
            } else {
               int nn=1;
               if(hDirty.containsKey(D[i])) {
                   nn  = hDirty.get(D[i]);
                   nn++;
               }
               hDirty.put(D[i], nn);
            }
        }

        Collection<Integer> collectionValues = hDirty.values();
        for(Integer i: collectionValues){
            if(nWash < K) {
                int x = Math.min(i, K-nWash);
                nWash = nWash+x;
                n = n+ x/2;
            } else {
                break;
            }
        }

        return n;
    }
}
