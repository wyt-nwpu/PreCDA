package javasrc;

import java.io.*;
import java.util.*;

import static javasrc.circRNAsim.readCirc2Disease;

public class UniNet {

    final static String CIRCRNA_rootString = "E:/data/github/datasource/";
    final static String CIRCRNA_workplace = CIRCRNA_rootString+"output/";
    /***
     * read circRNA similarity file
     * @param filepath the file path
     * @return HashMap of the circRNA similarity
     */
    public static HashMap<String, Double> readCircRNASimFile(String filepath){
        File file_CircSim = new File(CIRCRNA_workplace+filepath+".txt");
        BufferedReader bufReader_CircSim = null;
        HashMap<String, Double> results = new HashMap<>();
        try {
            bufReader_CircSim = new BufferedReader(new FileReader(file_CircSim));
            String temp = null;
            String[] line_strs;
            while ((temp = bufReader_CircSim.readLine()) != null) {
                line_strs = temp.split("\t");
                results.put(line_strs[0]+"\t"+line_strs[1], Double.parseDouble(line_strs[2]));
            }
        }catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (bufReader_CircSim != null) {
                try {
                    bufReader_CircSim.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        return results;
    }

    /****
     * get the file of the circRNA expression and functional similarity
     * @param exsimfile the circRNA expression similarity file
     * @param cossimfile the circRNA functional similarity file
     * @param output the output file name
     * @param threshold1 expression similarity threshold
     * @param threshold2 functional similarity threshold
     * @param threshold3 circRNA similarity threshold
     */
    public static void getCircRNASimByCosineandSpearman(String exsimfile, String cossimfile, String output, double threshold1,double threshold2,double threshold3){
        HashMap<String, Double> records = new HashMap<String, Double>();

        HashMap<String, Double> CircRNAEXP = readCircRNASimFile(exsimfile);
        HashMap<String, Double> CircRNACos = readCircRNASimFile(cossimfile);

        for (Map.Entry<String, Double> cosentry : CircRNACos.entrySet()){
            String[] crnas = cosentry.getKey().split("\t");
            double cos_val = cosentry.getValue();
            String temp_1 = crnas[0] + "\t" + crnas[1];
            String temp_2 = crnas[1] + "\t" + crnas[0];
            double exp_val = 0.0;
            double uni_cosim = 0.0;
            if(CircRNAEXP.containsKey(temp_1)){
                exp_val = CircRNAEXP.get(temp_1);
            }else if(CircRNAEXP.containsKey(temp_2)){
                exp_val = CircRNAEXP.get(temp_2);
            }
            if (exp_val < threshold1) {
                exp_val = 0.0;
            }
            if(cos_val < threshold2){
                cos_val = 0.0;
            }
            if(exp_val>0.0 && cos_val>0.0){
                //average
                uni_cosim = (exp_val + cos_val)/2;
            }else if(cos_val>0.0){
                uni_cosim = cos_val;
            }
            if(uni_cosim > threshold3) {
                records.put(crnas[0] + "\t" + crnas[1], uni_cosim);
            }
        }

        BufferedWriter bufWriter = null;
        //sort
        ArrayList<Map.Entry<String,Double>> list = new ArrayList<Map.Entry<String,Double>>(records.entrySet());
        Collections.sort(list,new Comparator<Map.Entry<String,Double>>() {
            //down
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return o2.getValue().compareTo(o1.getValue());
            }
        });

        try {
            bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace + output + ".txt"));
            for (int n = 0; n < list.size(); n++) {
                bufWriter.write(list.get(n).getKey() + "\t" + list.get(n).getValue() + "\n");
            }
        }catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (bufWriter != null) {
                try {
                    bufWriter.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /***
     * create the expression similarity file
     * @param exsimfile1 expression similarity from circBase
     * @param exsimfile2 expression similarity from circPedia
     * @param output the output file name
     * @param isall flag: all or join
     * @param method flag: ave; min; max
     * @param threshold1 expression similarity threshold
     */
    public static void combineEXSim(String exsimfile1, String exsimfile2, String output, String isall, String method, double threshold1){
        HashMap<String, Double> records = new HashMap<String, Double>();

        HashMap<String, Double> CircRNAEXP1 = readCircRNASimFile(exsimfile1);
        HashMap<String, Double> CircRNAEXP2 = readCircRNASimFile(exsimfile2);

        Set<String> key1 = CircRNAEXP1.keySet();
        Set<String> key2 = CircRNAEXP2.keySet();
        Set<String> keyset = new HashSet<String>();
        keyset.addAll(key1);
        keyset.addAll(key2);

        for(String circRNAstr : keyset){
            String[] keystrs = circRNAstr.split("\t");
            String circRNAstr2 = keystrs[1] + "\t" +keystrs[0];
            double value1 = 0.0;
            double value2 = 0.0;
            double value3 = 0.0;
            if(key1.contains(circRNAstr)){
                value1 = CircRNAEXP1.get(circRNAstr);
            }
            if(key1.contains(circRNAstr2)){
                value1 = CircRNAEXP1.get(circRNAstr2);
            }
            if(key2.contains(circRNAstr)){
                value2 = CircRNAEXP2.get(circRNAstr);
            }
            if(key2.contains(circRNAstr2)){
                value2 = CircRNAEXP2.get(circRNAstr2);
            }
            if(value1>0.0 && value2>0.0){
                if(method.equals("ave")){
                    value3 = (value1+value2)/2;
                }else if(method.equals("min")){
                    if(value1>value2){
                        value3 = value2;
                    }else{
                        value3 = value1;
                    }
                }else if(method.equals("max")){
                    if(value1<value2){
                        value3 = value2;
                    }else{
                        value3 = value1;
                    }
                }else{
                    System.out.println("ERROR");
                }
            }else if(value1 > threshold1){
                if (isall.equals("all")){
                    records.put(circRNAstr,value1);
                }
            }else if(value2 > threshold1){
                if (isall.equals("all")){
                    records.put(circRNAstr,value2);
                }
            }

            if(value3 > threshold1){
                records.put(circRNAstr,value3);
            }
        }

        BufferedWriter bufWriter = null;
        //sort
        ArrayList<Map.Entry<String,Double>> list = new ArrayList<Map.Entry<String,Double>>(records.entrySet());
        Collections.sort(list,new Comparator<Map.Entry<String,Double>>() {
            //down
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2) {
                return o2.getValue().compareTo(o1.getValue());
            }
        });

        try {
            bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace + output + ".txt"));
            for (int n = 0; n < list.size(); n++) {
                bufWriter.write(list.get(n).getKey() + "\t" + list.get(n).getValue() + "\n");
            }
        }catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (bufWriter != null) {
                try {
                    bufWriter.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public static void main(String args[]) {
        int flag = 1;
        double threshold1 = 0.7;//expression sim threshold
        double threshold2 = 0.0;//functional sim threshold
        double threshold3 = 0.0;//circRNA sim threshold
        //create the expression similarity file
        combineEXSim("CircRNAEXSim_"+flag, "CircRNAEXSimPedia_"+flag, "CircRNAEXSimBasePediaMaxJoin_"+threshold1+"_"+flag, "join", "max", 0.0);

        //get the file of the circRNA expression and functional similarity
        getCircRNASimByCosineandSpearman("CircRNAEXSimBasePediaMaxJoin_"+threshold1+"_"+flag,"CircRNASim_"+flag,"CircRNAUniCosSimBPAJAve_"+threshold1+"_"+threshold2+"_"+threshold3+"_"+flag,threshold1,threshold2,threshold3);
    }
}
