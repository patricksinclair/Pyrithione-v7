import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class BioSystem {


    Random rand = new Random();

    private int L, K; //length, karrying kapacity
    private double alpha, c_max; //steepness of antimicrobial gradient, max concn
    private ArrayList<Microhabitat> microhabitats;
    private double timeElapsed;
    private double tau; //timestep used in tau-leaping
    private double immigration_rate = 2880.;
    private double migration_rate = 0.2;
    private double attachment_rate = 2000.;
    private double detachment_rate = 1.38;
    private double delta_x = 5.;
    private int immigration_index, biofilm_edge_index;


    public BioSystem(int K, double alpha, double c_max, double tau){

        this.K = K;
        this.alpha = alpha;
        this.c_max = c_max;
        this.tau = tau;
        this.microhabitats = new ArrayList<>();
        this.timeElapsed = 0.;
        this.immigration_index = 0;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, c_max, alpha, delta_x), migration_rate));

        microhabitats.get(0).setSurface(true);
        microhabitats.get(0).addARandomBacterium_x_N(25);
    }


    public double getTimeElapsed(){
        return timeElapsed;
    }

    public int getImmigration_index(){return immigration_index;}

    public int getN_i(int index){return microhabitats.get(index).getN();}

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

    public double[] populationDistribution(){
        double[] popsizes = new double[microhabitats.size()];
        for(int i = 0; i < microhabitats.size(); i++){
            popsizes[i] = microhabitats.get(i).getN();
        }
        return popsizes;
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
        ArrayList<Microhabitat> updated_microhabs = new ArrayList<>(microhabitats.size());
        for(Microhabitat m : microhabitats){
            updated_microhabs.add(new Microhabitat(m));
        }

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
                        detachment_allocations[bac_index] = new PoissonDistribution(detachment_rate*tau_step).sample();
                        if(detachment_allocations[bac_index] > 1){
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
        immigrate(immigration_index, n_immigrants);
        updateBiofilmSize();
        timeElapsed += tau_step;


    }

    public static double calc_C_i(int i, double c_max, double alpha, double delta_x){
        return c_max*Math.exp(-alpha*i*delta_x);
    }


    public static void tester(){

        double duration = 240.;

        int K = 500;
        double c_max = 10.;
        double alpha = 1e-4;
        double tau = 0.01;

        BioSystem bs = new BioSystem(K, alpha, c_max, tau);
        //System.out.println(Arrays.toString(bs.getConcentrationProfile())+"\n");

        while(bs.getTimeElapsed() <= duration) {

            if(true) {

                System.out.println("-----------------------------------------------------------------------------------------------");

                String output = String.format("time elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tbf_edge pop: %d \tbf_edge fracfull: %.3f \timmig index: %d",
                        bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bs.getBiofilmEdge()), bs.microhabitats.get(bs.getBiofilmEdge()).fractionFull(), bs.getImmigration_index());

                String output2 = String.format("time elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d",
                        bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge());

                System.out.println(output);
                System.out.println("\nN distb");
                System.out.println(Arrays.toString(bs.populationDistribution()) + "\n");
                System.out.println("Biofilm y/n?");
                System.out.println(Arrays.toString(bs.biofilmProfile()) + "\n");
                System.out.println("-----------------------------------------------------------------------------------------------");
            }
            bs.performAction();

        }

    }
}
