import java.io.*;
import java.lang.reflect.Array;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

//C:\workspaceSchool\gaussian\out\production\gaussian
public class gaussian {

    private int n = 0;
    private double[][] coeff = null;
    private ArrayList<Double> consts = new ArrayList<Double>();
    private ArrayList<Double> sol = new ArrayList<Double>();
    private ArrayList<Integer> ind = new ArrayList<Integer>();

    public static void main(String[] args) throws Exception {
        final long startTime = System.currentTimeMillis();
        File toProcess = null;
        int userChoice = -1;

        if (args.length == 0){
            System.out.println("No arguments were given. " +
                    "\n Please pass in a file name ending in .lin as an argument.");
        } else if (args.length == 1 && args[0].contains("spp")){
            System.out.println("No file name was detected. " +
                    "\n Please pass in a file name ending in .lin as an argument " +
                    "\n to perform special partial pivoting gaussian on.");
        } else if (args.length == 1){
            File check = new File(args[0]);
            if (check.exists()){
                toProcess = check;
                userChoice = 1;

            } else {
                System.out.println("Hmm, the file doesn't seem to exist. " +
                        "\n Please run the program again.");
            }
        } else {
            if (args[0].contains("spp")){
                File check = new File(args[1]);
                System.out.println("This is the file name that is being processed: " + args[1]);
                if (check.exists()){
                    toProcess = check;
                    userChoice = 2;

                } else {
                    System.out.println("Hmm, the file doesn't seem to exist. " +
                            "\n Please run the program again.");
                }
            } else {
                System.out.println("For you to pass 2 arguments for this program, " +
                        "\n the first argument has to be --spp, followed by filename.lin");
            }
        }

        processFile(toProcess, userChoice);
        final long endTime = System.currentTimeMillis();
        if (userChoice == 1) {
            System.out.println("Naive Gaussian Execution Time: " + (endTime - startTime));
        } else {
            System.out.println("SPP Gaussian Execution Time: " + (endTime - startTime));
        }
    }

    public static void processFile(File toProcess, int userChoice) throws IOException {
        Scanner fileScanner = new Scanner(toProcess);
        gaussian gaussian = new gaussian();

        gaussian.n = fileScanner.nextInt();
        gaussian.coeff = new double[gaussian.n][gaussian.n];
        fileScanner.nextLine();

        int rows = gaussian.n;
        int columns = gaussian.n;

        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < columns; ++j) {
                if(fileScanner.hasNextDouble()) {
                    gaussian.coeff[i][j] = fileScanner.nextDouble();
                }
            }
        }

        for (int i = 0; i < rows; ++i){
            if (fileScanner.hasNextDouble()) {
                gaussian.consts.add(i, fileScanner.nextDouble());
            }
        }

        System.out.println(Arrays.deepToString(gaussian.coeff));
        System.out.println(gaussian.consts.toString());


        if (userChoice == 1){
            //Naive Gaussian
            NaiveGaussian(gaussian.coeff, gaussian.consts, gaussian.sol);
        } else if (userChoice == 2){
            //SPP Gaussian
            SPPGaussian(gaussian.coeff, gaussian.consts, gaussian.sol, gaussian.ind);

        }


        System.out.println("*****************************************************");
        System.out.println("*****************************************************");
        System.out.println("*****************************************************");
        System.out.println("*****************************************************");
        System.out.println(gaussian.sol);
        System.out.println("Solution has also been put into sys1.sol");
        System.out.println("*****************************************************");
        System.out.println("*****************************************************");
        System.out.println("*****************************************************");
        System.out.println("*****************************************************");
        fileScanner.close();
        outputToFile(gaussian.sol);
    }

    public static void outputToFile(ArrayList<Double> sol) throws IOException {

        Path currentRelativePath = Paths.get("");
        String s = currentRelativePath.toAbsolutePath().toString();
        File file = new File(s + "//sys1.sol");

        System.out.println("Solution file has been created/updated.");

        FileWriter writer = new FileWriter(file);
        for (int i = 0; i < sol.size(); i++){
            writer.write("X" + (i + 1) + " = " + sol.get(i) + String.format("%n"));

        }
        writer.close();

    }

    /*

    ********************************************************************
    ****NAIVE GAUSSIAN**************************************************
    ********************************************************************

     */
    public static void NaiveGaussian(double[][] coeff, ArrayList<Double> consts, ArrayList<Double> sol) {

        for (int i = 0; i < coeff.length; i++){
            sol.add(i, 0.0);
        }

        FwdElimination(coeff, consts);
        BackSubst(coeff, consts, sol);
    }

    //Forward Elimination
    public static void FwdElimination(double[][] coeff, ArrayList<Double> consts){

        for (int k = 0; k < coeff.length; k++){

            for (int i = k + 1; i < coeff.length; i= i + k + 1){


                double mult = coeff[i][k]/coeff[k][k];
                coeff[i][k] = mult;
                for (int j = k + 1; j < coeff.length; j = j + k + 1){
                    coeff[i][j] = coeff[i][j] - (mult) * coeff[k][j];
                }

                consts.set(i, consts.get(i) - (mult) * consts.get(i));
            }
        }
    }

    //Back Substitution
    public static void BackSubst(double[][] coeff, ArrayList<Double> consts, ArrayList<Double> sol){
//        System.out.println(" sol" + sol.toString());
//        System.out.println(" consts" + consts.toString());
//        System.out.println(" coeff" + Arrays.deepToString(coeff));

        sol.set(coeff.length - 1, consts.get(coeff.length - 1) / coeff[coeff.length - 1][coeff.length - 1]);

        for (int i = coeff.length - 2; i > -1; i--){
            double sum = consts.get(i);
            for (int j = i + 1; j < coeff.length; j = j + i + 1){
                sum = sum - coeff[i][j] * sol.get(j);
            }
            sol.set(i, sum / coeff[i][i]);
        }
    }







    /*

    ********************************************************************
    ****SPP GAUSSIAN****************************************************
    ********************************************************************

     */
    public static void SPPGaussian(double[][] coeff, ArrayList<Double> consts, ArrayList<Double> sol, ArrayList<Integer> ind){

        for (int i = 0; i < coeff.length; i++){
            sol.add(i, 0.0);
        }

        for (int i = 0; i < coeff.length; i++){
            ind.add(i, i);
        }

        SPPFwdElimination(coeff, consts, ind);
        SPPBackSubsts(coeff, consts, sol, ind);
    }

    //Forward SPP Elimination
    public static void SPPFwdElimination(double[][] coeff, ArrayList<Double> consts, ArrayList<Integer> ind){
        ArrayList<Double> scaling = new ArrayList<Double>();

        for (int i = 0; i < coeff.length; i++){
            scaling.add(i, 0.0);
        }

        //Initialize index and scaling vectors
        for (int i = 0; i < coeff.length; i++){
            double smax = 0;
            for (int j = 0; j < coeff.length; j++){
                smax = Math.max(smax, Math.abs(coeff[i][j])); //find coefficient with the greatest absolute value
            }

            scaling.set(i, smax);
        }

        for (int k = 1; k < coeff.length - 2; k++){

            double rmax = Math.abs(coeff[ind.get(k)][k]/scaling.get(ind.get(k)));
            int maxInd = k;

            for (int i = k + 1; i < coeff.length; i++){
                double r = Math.abs(coeff[ind.get(i)][k] / scaling.get(ind.get(i))); //Ratio of coeff to scaling factor

                if (r > rmax){
                    rmax = r;
                    maxInd = i;
                }
            }

            //Swapping indeces in ind list
            int temp = ind.get(maxInd);
            ind.set(maxInd, ind.get(k));
            ind.set(k, temp);

            for (int i = k + 1; i < coeff.length; i++){
                double mult = coeff[ind.get(i)][k] / coeff[ind.get(k)][k];
                coeff[ind.get(i)][k] = mult;

                for (int j = k + 1; j < coeff.length; j++){
                    coeff[ind.get(i)][j] = coeff[ind.get(i)][j] - mult * coeff[ind.get(k)][j];
                }

                consts.set(ind.get(i), consts.get(ind.get(i)) - coeff[ind.get(i)][k] * consts.get(ind.get(k)));
            }
        }
    }


    //Back Substitution
    public static void SPPBackSubsts(double[][] coeff, ArrayList<Double> consts, ArrayList<Double> sol, ArrayList<Integer> ind){
//        System.out.println(" sol" + sol.toString());
//        System.out.println(" consts" + consts.toString());
//        System.out.println(" coeff" + Arrays.deepToString(coeff));
//        System.out.println(" ind" + ind.toString());

        sol.set(coeff.length - 1, consts.get(coeff.length - 1) / coeff[ind.get(coeff.length - 1)][coeff.length - 1]);

        for (int i = coeff.length - 2; i > -1; i--){
            //double sum = consts.get(ind.get(i));

            sol.set(i, consts.get(i));
            for (int j = i + 1; j < coeff.length; j++){
                sol.set(i, sol.get(i) - coeff[ind.get(i)][j] * sol.get(j));
            }
            sol.set(i, sol.get(i) / coeff[ind.get(i)][i]);
        }
    }

}


