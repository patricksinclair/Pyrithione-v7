import java.util.ArrayList;

public class PyrithioneMain {
    public static void main(String[] args){
        //javac -cp commons-math3-3.6.1.jar:commons-math3-3.6.1-javadoc.jar:commons-math3-3.6.1-sources.jar:commons-math3-3.6.1-tests.jar:commons-math3-3.6.1-test-sources.jar:commons-math3-3.6.1-tools.jar:. *.java
        //System.out.println("precisest d_rate, second attempt");
        //BioSystem.tester();
        System.out.println("Biofilm thickness histogram - with counters - immig");
        BioSystem.getBiofilmThicknessHistoInParallel(96);
        //System.out.println("Biofilm det rate optimal - with antimicrobial");
        //BioSystem.findOptimalDetachmentRate();
    }
}
