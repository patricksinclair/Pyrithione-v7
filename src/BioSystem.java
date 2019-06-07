import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

public class BioSystem {


    Random rand = new Random();

    private int K; //karrying kapacity
    private double alpha, c_max; //steepness of antimicrobial gradient, max concn
    private ArrayList<Microhabitat> microhabitats;
    private double timeElapsed;
    private double tau = 0.01; //timestep used in tau-leaping
    private double immigration_rate =  0.8;
    private double migration_rate = 0.2;
    private double deterioration_rate;
    private double delta_x = 5.;
    private int immigration_index, biofilm_edge_index;
    //private int no_of_detachments = 0;


    public BioSystem(double alpha, double c_max){

        this.K = 120;
        this.alpha = alpha;
        this.c_max = c_max;
        this.microhabitats = new ArrayList<>();
        this.timeElapsed = 0.;
        this.immigration_index = 0;
        this.deterioration_rate = 1.e-9;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, this.c_max, this.alpha, delta_x), migration_rate));
        microhabitats.get(0).setSurface(true);
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }

    public BioSystem(double alpha, double c_max, double detachment_rate){

        this.K = 120;
        this.alpha = alpha;
        this.c_max = 0;
        this.microhabitats = new ArrayList<>();
        this.timeElapsed = 0.;
        this.immigration_index = 0;
        this.deterioration_rate = detachment_rate;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, this.c_max, this.alpha, delta_x), migration_rate));

        microhabitats.get(0).setSurface(true);
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }




    public double getTimeElapsed(){
        return timeElapsed;
    }

    public int getImmigration_index(){return immigration_index;}

    public int getN_i(int index){return microhabitats.get(index).getN();}

    public double getDeterioration_rate(){return deterioration_rate;}

    public int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }

    public int getBiofilmThickness(){
        return microhabitats.size();
    }

    public int getFormedBiofilmThickness(){
        int bf_counter = 0;
        for(Microhabitat m : microhabitats){
            if(m.isBiofilm_region()) bf_counter++;
        }
        return bf_counter;
    }

    public double[] populationDistribution(){
        double[] popsizes = new double[microhabitats.size()];
        for(int i = 0; i < microhabitats.size(); i++){
            popsizes[i] = microhabitats.get(i).getN();
        }
        return popsizes;
    }

    public double[] concnProfile(){
        double[] c_vals = new double[microhabitats.size()];
        for(int i = 0; i < microhabitats.size(); i++){
            c_vals[i] = microhabitats.get(i).getC();
        }
        return c_vals;
    }

    public boolean[] biofilmProfile(){
        boolean[] bf_vals = new boolean[microhabitats.size()];
        for(int i = 0; i < microhabitats.size(); i++){
            bf_vals[i] = microhabitats.get(i).isBiofilm_region();
        }
        return bf_vals;
    }

    public int getBiofilmEdge(){
        int edgeIndex = 0;
        for(int i = 0; i < microhabitats.size(); i++){
            if(microhabitats.get(i).isBiofilm_region()) edgeIndex = i;
        }
        return edgeIndex;
    }




    public void immigrate(int mh_index, int n_immigrants){
        microhabitats.get(mh_index).addARandomBacterium_x_N(n_immigrants);
    }


    public void migrate(ArrayList<Microhabitat> updated_microhabs, int mh_index, int bac_index){

        int biof_edge = immigration_index;
        double migrating_bac = updated_microhabs.get(mh_index).getPopulation().get(bac_index);
        updated_microhabs.get(mh_index).removeABacterium(bac_index);

        if(updated_microhabs.get(mh_index).isSurface()){
            updated_microhabs.get(mh_index+1).addABacterium(migrating_bac);

        }else if(updated_microhabs.get(mh_index).isImmigration_zone()){
            updated_microhabs.get(mh_index-1).addABacterium(migrating_bac);

        }else{
            if(rand.nextBoolean()){
                updated_microhabs.get(mh_index+1).addABacterium(migrating_bac);
            }else{
                updated_microhabs.get(mh_index-1).addABacterium(migrating_bac);
            }
        }
    }



    public void updateBiofilmSize(){
        //once the edge microhabitat is sufficiently populated, this adds another microhabitat onto the system list
        //which is then used as the immigration zone

        if(microhabitats.get(immigration_index).atBiofilmThreshold()){

            microhabitats.get(immigration_index).setBiofilm_region(true);
            microhabitats.get(immigration_index).setImmigration_zone(false);

            int i = microhabitats.size();
            microhabitats.add(new Microhabitat(K, BioSystem.calc_C_i(i, c_max, alpha, delta_x), migration_rate));
            immigration_index = i;
            microhabitats.get(immigration_index).setImmigration_zone(true);
        }
    }


    public void performAction(){
        //make a copy of the microhabitats
        //todo does this actually need to be done?
        /*ArrayList<Microhabitat> updated_microhabs = new ArrayList<>(microhabitats.size());
        for(Microhabitat m : microhabitats){
            //updated_microhabs.add(new Microhabitat(m));
        }*/

        double tau_step = tau;

        int system_size = microhabitats.size(); //this is all the microhabs in the system
        int[][] replication_allocations;
        int[][] death_allocations;
        int[][] migration_allocations;
        int[] detachment_allocations;
        int[] original_popsizes;
        int n_immigrants;

        whileloop:
        while(true){
            replication_allocations = new int[system_size][];
            death_allocations = new int[system_size][];
            migration_allocations = new int[system_size][];
            original_popsizes = new int[system_size];
            detachment_allocations = new int[microhabitats.get(immigration_index).getN()];

            for(int mh_index = 0; mh_index < system_size; mh_index++){
                int mh_pop = microhabitats.get(mh_index).getN();
                int[] n_replications = new int[mh_pop];
                int[] n_deaths = new int[mh_pop];
                int[] n_migrations = new int[mh_pop];

                for(int bac_index = 0; bac_index < mh_pop; bac_index++){

                    ////////// MIGRATIONS //////////////////////
                    n_migrations[bac_index] = new PoissonDistribution(microhabitats.get(mh_index).migrate_rate()*tau_step).sample();

                    if(n_migrations[bac_index] > 1){
                        tau_step /= 2.;
                        continue whileloop;
                    }
                    ////////////////////////////////////////////

                    ///////////// DETACHMENTS /////////////////////////
                    if(mh_index == immigration_index){
                        detachment_allocations[bac_index] = new PoissonDistribution(deterioration_rate*tau_step).sample();
                        //no_of_detachments += detachment_allocations[bac_index];
                        if(detachment_allocations[bac_index] > 1){
                            //no_of_detachments -= detachment_allocations[bac_index];
                            tau_step /= 2.;
                            continue whileloop;
                        }
                        //if a bacteria is detaching then it can't migrate
                        if(detachment_allocations[bac_index] != 0){
                            n_migrations[bac_index] = 0;
                        }
                    }
                    ////////////////////////////////////////////////////////

                    ////////////////// REPLICATIONS AND DEATHS ///////////////////////////
                    double g_or_d_rate = microhabitats.get(mh_index).replicationOrDeathRate(bac_index);

                    if(g_or_d_rate == 0.){

                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = 0;

                    }else if(g_or_d_rate > 0){

                        n_replications[bac_index] = new PoissonDistribution(g_or_d_rate*tau_step).sample();
                        n_deaths[bac_index] = 0;

                    }else{
                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = new PoissonDistribution(Math.abs(g_or_d_rate)*tau_step).sample();

                        if(n_deaths[bac_index] > 1){
                            tau_step /= 2.;
                            continue whileloop;
                        }
                        //if a death is occurring, then that bacteria can't migrate or detach
                        if(n_deaths[bac_index] !=0) {
                            n_migrations[bac_index] = 0;
                            if(mh_index == immigration_index) detachment_allocations[bac_index] = 0;
                        }
                    }
                    /////////////////////////////////////////////////////////////////////////
                }

                replication_allocations[mh_index] = n_replications;
                death_allocations[mh_index] = n_deaths;
                migration_allocations[mh_index] = n_migrations;
                original_popsizes[mh_index] = microhabitats.get(mh_index).getN();
            }

            n_immigrants = new PoissonDistribution(immigration_rate*tau_step).sample();
            break whileloop;
        }


        for(int mh_index = 0; mh_index < system_size; mh_index++){
            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){

                if(death_allocations[mh_index][bac_index]!= 0) microhabitats.get(mh_index).removeABacterium(bac_index);

                else{
                    microhabitats.get(mh_index).replicateABacterium_x_N(bac_index, replication_allocations[mh_index][bac_index]);

                    if(system_size > 1){
                        if(migration_allocations[mh_index][bac_index] != 0) migrate(microhabitats, mh_index, bac_index);
                    }

                    if(mh_index == immigration_index){
                        if(detachment_allocations[bac_index] != 0) microhabitats.get(mh_index).removeABacterium(bac_index);
                    }
                }
            }
        }
        //System.out.println(n_immigrants);
        //if(n_immigrants > 0)System.out.println("n_immigrants: "+n_immigrants);
        //System.out.println("\ndetachments");
        //System.out.println(Arrays.toString(detachment_allocations) + "\n");
        immigrate(immigration_index, n_immigrants);
        updateBiofilmSize();
        timeElapsed += tau_step;


    }

    public static double calc_C_i(int i, double c_max, double alpha, double delta_x){
        return c_max*Math.exp(-alpha*i*delta_x);
    }


    public static void tester(){

        double duration = 240.;

        int K = 120;
        double c_max = 10.;
        double alpha = 1e-4;
        double tau = 0.01;

        BioSystem bs = new BioSystem(alpha, c_max);
        //System.out.println(Arrays.toString(bs.getConcentrationProfile())+"\n");
        System.out.println(bs.deterioration_rate);

        while(bs.getTimeElapsed() <= duration) {

            if(true) {

                System.out.println("-----------------------------------------------------------------------------------------------");

                String output = String.format("time elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tbf_edge pop: %d \tbf_edge fracfull: %.3f \timmig index: %d",
                        bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bs.getBiofilmEdge()), bs.microhabitats.get(bs.getBiofilmEdge()).fractionFull(), bs.getImmigration_index());

                String output2 = String.format("time elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d",
                        bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge());

                System.out.println(output);
                System.out.println("\nconcn profile");
                System.out.println(Arrays.toString(bs.concnProfile()));
                System.out.println("\nN distb");
                System.out.println(Arrays.toString(bs.populationDistribution()) + "\n");
                System.out.println("Biofilm y/n?");
                System.out.println(Arrays.toString(bs.biofilmProfile()) + "\n");
                System.out.println("-----------------------------------------------------------------------------------------------");
            }
            bs.performAction();
        }
        //System.out.println(bs.no_of_detachments);
    }







    public static int getThicknessReachedAfterATime(double duration, int i){
        int K = 120;
        double c_max = 10., alpha = 0.01;

        BioSystem bs = new BioSystem(alpha, c_max);
        System.out.println("detach_rate: "+bs.deterioration_rate);
        int nUpdates = 20;
        double interval = duration/nUpdates;
        boolean alreadyRecorded = false;


        while(bs.timeElapsed <= duration){

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int max_poss_pop = bs.getBiofilmThickness()*K;
                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+bs.getTotalN()+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge());
                alreadyRecorded = true;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;


            bs.performAction();
        }

        return bs.getBiofilmEdge();
    }


    public static void getBiofilmThicknessHistoInParallel(int nReps){
        long startTime = System.currentTimeMillis();

        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nSections = 8; //number of sections the reps will be divided into, to avoid using loadsa resources
        int nRuns = nReps/nSections; //number of runs in each section

        double duration = 1680.; //10 week duration


        int[] mh_index_reached = new int[nReps];
        String index_reached_filename = "pyrithione-bf-thickness_histo-t="+String.valueOf(duration)+"-parallel-drate_miniscule";

        for(int j = 0; j < nSections; j++){
            IntStream.range(j*nRuns, (j+1)*nRuns).parallel().forEach(i -> mh_index_reached[i] = BioSystem.getThicknessReachedAfterATime(duration, i));
        }


        Toolbox.writeHistoArrayToFile(index_reached_filename, mh_index_reached);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }



    public static DataBox getAllData(int i){

        int K = 120;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nMeasurements = 40;
        double duration = 200., interval = duration/nMeasurements;

        BioSystem bs = new BioSystem(alpha, c_max);

        boolean alreadyRecorded = false;
        int timerCounter = 0;

        double[] popSizes = new double[nMeasurements+1];
        double[][] popDistbs = new double[nMeasurements+1][];
        double[] biofilmEdges = new double[nMeasurements+1];
        double[][] avgGenotypeDistbs = new double[nMeasurements+1][];
        double[][] genoStDevs = new double[nMeasurements+1][];

        while(bs.timeElapsed <= duration+0.02*interval){

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int max_poss_pop = bs.getBiofilmThickness()*K;
                int total_N = bs.getTotalN();

                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+total_N+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge());

                /*popSizes[timerCounter] = total_N;
                popDistbs[timerCounter] = bs.getPopulationDistribution();
                biofilmEdges[timerCounter] = bs.getBiofilmEdge();
                avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                genoStDevs[timerCounter] = bs.getStDevOfGenotypeDistribution();*/

                alreadyRecorded = true;
                timerCounter++;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }

        return new DataBox(popSizes, popDistbs, biofilmEdges, avgGenotypeDistbs, genoStDevs);
    }


    public static void getInfoInParallel(){

        long startTime = System.currentTimeMillis();

        int K = 120, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nReps = 16;
        double duration = 200.;

        double[][] allPopSizes = new double[nReps][];
        double[][][] allPopDistbs = new double[nReps][][];
        double[][] allBiofilmEdges = new double[nReps][];
        double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        double[][][] allGenoStDevs = new double[nReps][][];

        String popSizeFilename = "pyrithione-testing-pop_size-t="+String.valueOf(duration)+"-parallel";
        String popDistbFilename = "pyrithione-testing-pop_distb-t="+String.valueOf(duration)+"-parallel";
        String biofilmEdgeFilename = "pyrithione-testing-biofilm_edge-t="+String.valueOf(duration)+"-parallel";
        String avgGenotypeDistbFilename = "pyrithione-testing-avgGenoDistb-t="+String.valueOf(duration)+"-parallel";
        String genoStDevDistbFilename = "pyrithione-testing-genoStDevDistb-t="+String.valueOf(duration)+"-parallel";

        //double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        DataBox[] dataBoxes = new DataBox[nReps];

        IntStream.range(0, nReps).parallel().forEach(i -> dataBoxes[i] = BioSystem.getAllData(i));

        for(int j = 0; j < dataBoxes.length; j++){
            allPopSizes[j] = dataBoxes[j].getPopSizes();
            allPopDistbs[j] = dataBoxes[j].getPopDistbs();
            allBiofilmEdges[j] = dataBoxes[j].getBiofilmEdges();
            allAvgGenotypeDistbs[j] = dataBoxes[j].getAvgGenotypeDistbs();
            allGenoStDevs[j] = dataBoxes[j].getGenoStDevs();
        }

        double[] processedPopSizes = Toolbox.averagedResults(allPopSizes);
        double[][] processedPopDistbs = Toolbox.averagedResults(allPopDistbs);
        double[] processedBiofilmEdges = Toolbox.averagedResults(allBiofilmEdges);
        double[][] processedAvgGenotypeDistbs = Toolbox.averagedResults(allAvgGenotypeDistbs);
        double[][] processedGenoStDevs = Toolbox.averagedResults(allGenoStDevs);

        Toolbox.writeAveragedArrayToFile(popSizeFilename, processedPopSizes);
        Toolbox.writeAveragedDistbsToFile(popDistbFilename, processedPopDistbs);
        Toolbox.writeAveragedArrayToFile(biofilmEdgeFilename, processedBiofilmEdges);
        Toolbox.writeAveragedDistbsToFile(avgGenotypeDistbFilename, processedAvgGenotypeDistbs);
        Toolbox.writeAveragedDistbsToFile(genoStDevDistbFilename, processedGenoStDevs);

        long finishTime = System.currentTimeMillis();

        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);

        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }


    public static int[] optimalDetachmentSubSubroutine(double d_rate, double timelimit, int i){
        //this plays a biosystem to completion of one rep for a specified detachment rate
        //returns the thickness and pop size of this one rep
        double alpha = 0.0, c_max = 0.;
        boolean alreadyRecorded = false;
        int nMeasurements = 10;
        double interval = timelimit/nMeasurements;

        BioSystem bs = new BioSystem(alpha, c_max, d_rate);
        while(bs.timeElapsed <= timelimit){

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int total_N = bs.getTotalN();
                System.out.println("d_rate: "+bs.getDeterioration_rate()+"\trep : "+i+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+total_N+"\tbf_edge: "+bs.getBiofilmEdge());
                alreadyRecorded = true;
            }

            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }
        //System.out.println("TIME ELAPSED: "+bs.timeElapsed);
        return new int[]{bs.getBiofilmThickness(), bs.getTotalN()};
    }


    public static double[] optimalDetachmentSubroutine(double d_rate, double timelimit, int nReps){
        //this runs several reps of a biosystem for a given detachment rate
        //returns the average and stdev of the thickness and pop size of these reps.

        int[][] sub_results = new int[nReps][];
        double[] bf_thicknesses = new double[nReps];
        double[] pop_sizes = new double[nReps];
        double alpha = 0.01, c_max = 0.;

        IntStream.range(0, nReps).parallel().forEach(i -> sub_results[i] = BioSystem.optimalDetachmentSubSubroutine(d_rate, timelimit, i));

        for(int r = 0; r < nReps; r++){
            bf_thicknesses[r] = sub_results[r][0];
            pop_sizes[r] = sub_results[r][1];
        }

        double[] thickness_avg_stDev = Toolbox.averageAndStDevOfArray(bf_thicknesses);
        double[] popsize_avg_stDev = Toolbox.averageAndStDevOfArray(pop_sizes);

        double thickness_avg = thickness_avg_stDev[0], thickness_stDev = thickness_avg_stDev[1];
        double popsize_avg = popsize_avg_stDev[0], popsize_stDev = popsize_avg_stDev[1];

        return new double[]{thickness_avg, thickness_stDev, popsize_avg, popsize_stDev};
    }



    public static void findOptimalDetachmentRate(){

        long startTime = System.currentTimeMillis();

        double min_detachment = 0.0515, max_detachment = 0.0517;
        int n_detachments = 24; //number of detachment rates measured
        double detach_increment = (max_detachment-min_detachment)/(double)n_detachments;
        int nReps = 20;
        double duration = 240.; //10 days



        String filename = String.format("optimal_detach_rates-range=%.5f_%.5f_%.5f-REDO", min_detachment, detach_increment, max_detachment);

        double[] dRateArray = new double[n_detachments+1];
        double[] thickness_array_avg = new double[n_detachments+1];
        double[] thickness_array_stDev = new double[n_detachments+1];
        double[] popsize_array_avg = new double[n_detachments+1];
        double[] popsize_array_stDev = new double[n_detachments+1];

        for(int i = 0; i <= n_detachments; i++){
            dRateArray[i] = min_detachment + i*detach_increment;

            double[] subroutine_findings = BioSystem.optimalDetachmentSubroutine(dRateArray[i], duration, nReps);

            thickness_array_avg[i] = subroutine_findings[0];
            thickness_array_stDev[i] = subroutine_findings[1];
            popsize_array_avg[i] = subroutine_findings[2];
            popsize_array_stDev[i] = subroutine_findings[3];
        }

        double[][] collated_results = new double[][]{dRateArray, thickness_array_avg, thickness_array_stDev, popsize_array_avg, popsize_array_stDev};

        Toolbox.writeMultipleColumnsToFile(filename, collated_results);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }
}
