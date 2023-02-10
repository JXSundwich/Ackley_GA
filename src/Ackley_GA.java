import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

public class Ackley_GA {
    double ackleyN;
    int populationSize;
    int[] xBound;
    int DNASize;
    double mutationRate;
    double recombinationRate;
    List<List<int[]>> DNAs;
    String seedsString;
    long[] seeds;
    int seedsListSize;
    int seedsPointer;

    DecimalFormat df = new DecimalFormat("0.00000000");

    Ackley_GA(int n, int populationSize, int[] xBound, int DNASize, double mutationRate, double recombinationRate) {
        this.ackleyN = n;
        this.populationSize = populationSize;
        this.xBound = xBound;
        this.DNASize = DNASize;
        this.mutationRate = mutationRate;
        this.recombinationRate = recombinationRate;
        //initialize seed's index with my matric number mod 100
        this.seedsPointer = 2201757 % 100-1;
        this.seedsString = "181 586 430 151 940 637 314 954 243 159 962 26 463 884 749 726 882 721 276 135 680 405 488 11 359 199 426 533 594 407 478 673 135 360 961 541 583 143 998 210 540 697 173 240 864 266 116 964 565 684 121 716 569 699 67 963 284 10 471 778 294 89 576 17 911 798 817 800 78 943 673 504 46 493 709 185 779 211 268 793 440 512 871 602 874 91 606 619 476 957 497 939 187 261 1 892 36 570 997 881 97 498 171 613 432 200 247 947 197 796 978 532 852 377 311 894 731 523 654 127 373 537 808 135 726 74 503 12 759 849 119 859 789 284 108 419 331 875 110 938 939 967 852 488 345 214 891 845 696 329 252 154 674 653 333 190 335 874 480 678 333 725 237 671 608 27 834 935 285 834 932 639 728 402 154 200 657 899 743 532 661 880 649 973 921 4 922 791 132 42 7 82 136 529 159 406 516 516 172 169 ";
        //initialize 1st generation
        DNAs = new ArrayList<>();

        //random generate 1st generation
        /*
        for (int i = 0; i < this.populationSize; i++) {
            DNAs.add(new ArrayList<>());
            List<int[]> temp = DNAs.get(i);
            for (int j = 0; j < this.ackleyN; j++) {
               temp.add(randomDNA(DNASize));
            }
        }**/

        for(int i=0;i<populationSize;i++){
            DNAs.add(new ArrayList<>());
        }
        DNAs.get(0).add(new int[]{0, 0, 1, 1, 0, 0, 0, 1, 1, 1});
        DNAs.get(0).add(new int[]{1, 0, 0, 0, 0, 0, 1, 0, 0, 1});
        DNAs.get(0).add(new int[]{1, 0, 1, 0, 1, 1, 1, 1, 0, 1});

        DNAs.get(1).add(new int[]{1, 1, 0, 1, 1, 0, 0, 0, 0, 1});
        DNAs.get(1).add(new int[]{1, 0, 0, 0, 0, 0, 1, 1, 0, 0});
        DNAs.get(1).add(new int[]{1, 0, 0, 0, 1, 1, 1, 1, 0, 0});

        DNAs.get(2).add(new int[]{1, 1, 1, 0, 0, 1, 1, 0, 1, 0});
        DNAs.get(2).add(new int[]{1, 1, 0, 0, 1, 1, 0, 1, 1, 1});
        DNAs.get(2).add(new int[]{1, 0, 0, 0, 1, 1, 0, 0, 1, 0});

        DNAs.get(3).add(new int[]{1, 0, 1, 1, 1, 0, 1, 0, 1, 0});
        DNAs.get(3).add(new int[]{1, 0, 0, 1, 0, 1, 0, 1, 1, 0});
        DNAs.get(3).add(new int[]{1, 0, 0, 1, 0, 1, 0, 0, 1, 0});

        DNAs.get(4).add(new int[]{1, 1, 1, 1, 0, 1, 1, 0, 1, 1});
        DNAs.get(4).add(new int[]{1, 1, 1, 0, 0, 1, 1, 1, 0, 1});
        DNAs.get(4).add(new int[]{0, 0, 1, 0, 1, 1, 1, 1, 0, 0});

        DNAs.get(5).add(new int[]{1, 0, 1, 1, 0, 0, 0, 1, 0, 1});
        DNAs.get(5).add(new int[]{0, 0, 0, 1, 1, 0, 1, 0, 0, 0});
        DNAs.get(5).add(new int[]{1, 1, 0, 1, 1, 1, 1, 0, 1, 0});

        DNAs.get(6).add(new int[]{1, 1, 1, 1, 1, 0, 1, 1, 0, 1});
        DNAs.get(6).add(new int[]{1, 1, 1, 1, 1, 0, 0, 1, 1, 1});
        DNAs.get(6).add(new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

        DNAs.get(7).add(new int[]{1, 0, 1, 0, 1, 0, 0, 0, 1, 1});
        DNAs.get(7).add(new int[]{1, 0, 0, 0, 1, 1, 1, 1, 0, 1});
        DNAs.get(7).add(new int[]{0, 0, 1, 1, 0, 1, 0, 0, 0, 1});

        DNAs.get(8).add(new int[]{1, 0, 0, 0, 0, 1, 1, 1, 1, 0});
        DNAs.get(8).add(new int[]{1, 1, 0, 1, 1, 1, 1, 1, 0, 1});
        DNAs.get(8).add(new int[]{0, 0, 1, 1, 1, 1, 1, 1, 1, 0});

        DNAs.get(9).add(new int[]{0, 0, 0, 0, 0, 0, 1, 1, 0, 0});
        DNAs.get(9).add(new int[]{0, 1, 1, 1, 0, 1, 1, 1, 0, 1});
        DNAs.get(9).add(new int[]{0, 1, 0, 1, 0, 1, 0, 0, 0, 0});


        String[] seedsStrings = seedsString.split(" ");
        this.seedsListSize = seedsStrings.length;
        this.seeds = new long[seedsStrings.length];
        for (int i = 0; i < seedsStrings.length; i++) {
            this.seeds[i] = Long.valueOf(seedsStrings[i]);
        }

        long[] transSeeds=new long[seedsListSize];
        for(int j=0;j<40;j++){
            transSeeds[j]=seeds[j*5];
            transSeeds[j+40]=seeds[j*5+1];
            transSeeds[j+80]=seeds[j*5+2];
            transSeeds[j+120]=seeds[j*5+3];
            transSeeds[j+160]=seeds[j*5+4];
        }
        this.seeds=transSeeds;

    }


    public double calculateAckley(List<int[]> binaryList) {
        List<Double> args = new ArrayList<>();
        for (int i = 0; i < this.ackleyN; i++) {
            int temp = 0;
            int[] tempArray = binaryList.get(i);
            for (int j = 0; j < this.DNASize; j++) {
                temp += Math.pow(2, j) * tempArray[j];
            }
            double div=Double.valueOf(this.df.format(temp / (Math.pow(2, DNASize) - 1)));
            double tempRes = (this.xBound[1] - this.xBound[0]) * div + this.xBound[0];
            args.add(Double.valueOf(this.df.format(tempRes)));
        }
        double res = 0;
        double expValue = 0;
        double cosValue = 0;
        for (int i = 0; i < this.ackleyN; i++) {
            expValue += Math.pow(args.get(i), 2);
            cosValue += Math.cos(2 * Math.PI * args.get(i));
        }
        expValue = -0.2 * Math.pow(1 / this.ackleyN*expValue, 1d / 2);
        cosValue=1 / this.ackleyN * cosValue;
        res = -20 * Math.exp(expValue) - Math.exp(cosValue) + 20;
        BigDecimal bd1=new BigDecimal(res);
        bd1=bd1.add(new BigDecimal(Math.E));
        return Double.valueOf(df.format(bd1.doubleValue()));
    }

    public double[] fitnessToCumulativeProbability(List<List<int[]>> _DNAs) {
        int n = _DNAs.size();
        double sum = 0;
        double[] fitness = new double[n];
        double[] cumulativeProbability = new double[n];
        for (int i = 0; i < n; i++) {
            //The maximum value of ackley's function is less than 30
            //TODO 30
            fitness[i] =Math.round((30 - calculateAckley(_DNAs.get(i))));
            sum += fitness[i];
        }
        //
        double[] selectionProp=new double[n];
        selectionProp[0]=Double.valueOf(df.format(fitness[0] / sum));


        cumulativeProbability[0] = Double.valueOf(df.format(fitness[0] / sum));

        for (int i = 1; i < n; i++) {
            selectionProp[i]=Double.valueOf(df.format(fitness[i] / sum));
            cumulativeProbability[i] = Double.valueOf(df.format(cumulativeProbability[i - 1] + selectionProp[i]));
           // cumulativeProbability[i]=Double.valueOf(df.format(cumulativeProbability[i]));
        }
        cumulativeProbability[n - 1] = 1;
        return cumulativeProbability;
    }

    public void selection(double[] cumulativeProbability) {
        List<List<int[]>> nextGeneration = new ArrayList<>();
        double[] randNums = new double[this.populationSize];

        //use the seeds given by prof to select the next generation
        for (int i = 0; i < this.populationSize; i++) {
            randNums[i] = this.getNextSeed()/1000.0;
        }

        Arrays.sort(randNums);
        int randNumPointer = 0;
        int cupPointer = 0;
        while (randNumPointer < this.populationSize) {
            if (randNums[randNumPointer] <= cumulativeProbability[cupPointer]) {
                nextGeneration.add(this.DNAs.get(cupPointer));
                randNumPointer++;
            } else {
                cupPointer++;
                continue;
            }
        }
        this.DNAs = nextGeneration;
    }


    List<Long> mutationSeed=new ArrayList<>();
    public void mutation() {
        mutationSeed.clear();
        for (int i = 0; i < this.populationSize; i++) {
            for (int j = 0; j < this.ackleyN; j++) {
                for (int x = 0; x < this.DNASize; x++) {
                    long tempSeed=this.getNextSeed();
                    mutationSeed.add(tempSeed);
                    if (tempSeed/1000.0 <= mutationRate) {
                        this.DNAs.get(i).get(j)[x] = 1 - this.DNAs.get(i).get(j)[x];
                    }
                }
            }
        }
        System.out.println("end mutation");
    }

    List<Long> recombinationSeeds=new ArrayList<Long>();
    public void recombination() {
        //flip coins
        recombinationSeeds.clear();
        for (int i = 0; i < populationSize-1; i++) {
            long tempSeed=this.getNextSeed();
            recombinationSeeds.add(tempSeed);
            if (tempSeed/1000.0 <= this.recombinationRate) {
                int position = (i+1)%this.populationSize;
                //get next chromosome
                List<List<int[]>> res = executeRecombination(this.DNAs.get(i), this.DNAs.get(position));
                //update DNAs
                DNAs.set(i, res.get(0));
                DNAs.set(position, res.get(1));
                i++;
            }
        }
        System.out.println("end recombination");
    }

    /**
     * execute recombination:
     * exchange the latter part of each variation in chromosome
     *
     * @param chrom1
     * @param chrom2
     * @return
     */
    public List<List<int[]>> executeRecombination(List<int[]> chrom1, List<int[]> chrom2) {
        List<int[]> res1 = new ArrayList<>();
        List<int[]> res2 = new ArrayList<>();

        for (int i = 0; i < this.ackleyN; i++) {
            int[] array1 = new int[this.DNASize];
            int[] array2 = new int[this.DNASize];
/*

            for (int j = 0; j < this.DNASize / 2; j++) {
                array1[j] = chrom1.get(i)[j];
                array2[j] = chrom2.get(i)[j];
            }
            for (int x = this.DNASize / 2; x < DNASize; x++) {
                array1[x] = chrom2.get(i)[x];
                array2[x] = chrom1.get(i)[x];
            }
            */

            //TODO flip coins for crossover
            for(int j=0;j<this.DNASize;j++){
                long tempSeed=this.getNextSeed();
                recombinationSeeds.add(tempSeed);
                if(tempSeed/1000.0>=0.5){
                    array1[j]=chrom1.get(i)[j];
                    array2[j]=chrom2.get(i)[j];
                }else{
                    array2[j]=chrom1.get(i)[j];
                    array1[j]=chrom2.get(i)[j];
                }
            }

            res1.add(array1);
            res2.add(array2);
        }
        List<List<int[]>> finalRes = new ArrayList<>();
        finalRes.add(res1);
        finalRes.add(res2);
        return finalRes;
    }
/*
    public long getNextSeed() {
        long next = this.seeds[seedsPointer];
        seedsPointer = (seedsPointer + 1) % seedsListSize;
        return next;
    }*/

    public long getNextSeed(){
        long seed=this.seeds[seedsPointer];
        this.seedsPointer=(this.seedsPointer+1)%seedsListSize;
        return seed;
    }

    public int[] randomDNA(int size) {
        //generate 0/1 randomly
        Random rand = new Random();
        int[] res = new int[size];
        for (int i = 0; i < size; i++) {
            //generate 0/1 randomly
             res[i] = rand.nextInt(2);
        }
        return res;
    }

    public double[] getAckleyValuesOfGeneration() {

        //the last value is the minimum value in this generation
        double[] res = new double[this.populationSize+1];
        res[this.populationSize]=30;
        for (int i = 0; i < this.populationSize; i++) {
            res[i] = calculateAckley(DNAs.get(i));
            //update minimum value
            res[this.populationSize]=res[i]<res[this.populationSize]?res[i]:res[this.populationSize];
        }
        return res;
    }

    public void startGA() {
        //parents selection
        selection(fitnessToCumulativeProbability(this.DNAs));
        recombination();
        mutation();
        //offspring selection
       // selection(fitnessToCumulativeProbability(this.DNAs));
    }

    public static void main(String[] args) {
        //The recombination rate is usually between 0.6 and 0.9, we choose 0.7
        //The mutation rate is usually between 1/pop_size(10) and 1/chromosome_size(30) so we choose 0.06
        Ackley_GA ag = new Ackley_GA(3, 10, new int[]{-20, 20}, 10, 0.06, 0.7);

    //    System.out.println("Ackley's function values for first generation: ");
       // double[] values2 = ag.getAckleyValuesOfGeneration();
      //  System.out.println(values2);
       // int[] x=new int[500];
       // double[] y=new double[500];
        for (int i = 0; i < 3; i++) {
            ag.startGA();
            double[] values1 = ag.getAckleyValuesOfGeneration();
//            x[i]=i;
//            y[i]=values1[values1.length-1];
       //     System.out.println(values1);
        }

        double[] values = ag.getAckleyValuesOfGeneration();

        System.out.println(values);
    }
}

