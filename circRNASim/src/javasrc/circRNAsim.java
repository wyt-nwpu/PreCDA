package javasrc;

import java.io.*;
import java.util.*;

/***
 * calculate circRNA functional similarity
 */
public class circRNAsim {

    public static HashMap<String, Double> DoSimlist;
    final static String CIRCRNA_rootString = "PreCDA/datasource/";
    final static String CIRCRNA_workplace = CIRCRNA_rootString+"output/";

    /****
     * read circRNA-disease associations
     * @param Do2Circ_filename the file of circRNA-disease associations
     * @param flag the key of the hashmap
     * @return the hashmap of circRNA-disease associations
     */
    public static HashMap<String, ArrayList<String>> readCirc2Disease(String Do2Circ_filename,String flag) {
        File file_Circ2DO = new File(CIRCRNA_workplace+Do2Circ_filename+".txt");
        BufferedReader bufReader_Circ2DO = null;
        HashMap<String, ArrayList<String>> Circ2DO = new HashMap<String, ArrayList<String>>();
        try {
            bufReader_Circ2DO = new BufferedReader(new FileReader(file_Circ2DO));
            String temp = null;
            String[] line_strs;
            while ((temp = bufReader_Circ2DO.readLine()) != null) {
                line_strs = temp.split("\t");
                if(flag.equals("rna")){
                    if (Circ2DO.containsKey(line_strs[0])) {
                        ArrayList<String> values = Circ2DO.get(line_strs[0]);
                        values.add(line_strs[1]);
                    } else {
                        ArrayList<String> values = new ArrayList<String>();
                        values.add(line_strs[1]);
                        Circ2DO.put(line_strs[0], values);
                    }
                }else{
                    if (Circ2DO.containsKey(line_strs[1])) {
                        ArrayList<String> values = Circ2DO.get(line_strs[1]);
                        values.add(line_strs[0]);
                    } else {
                        ArrayList<String> values = new ArrayList<String>();
                        values.add(line_strs[0]);
                        Circ2DO.put(line_strs[1], values);
                    }
                }
            }
        } catch (Exception e) {
            // TODO: handle exception
            System.out.println("readCirc2Disease Exception e: " + e);
            e.printStackTrace();
        } finally {
            if (bufReader_Circ2DO != null) {
                try {
                    bufReader_Circ2DO.close();
                } catch (IOException e) {
                    e.getStackTrace();
                }
            }
        }
        return Circ2DO;
    }

    /****
     * calculate circRNA functional similarity by cosine similarity
     * @param output file name to be saved
     * @param Do2Circ_filename the file of circRNA-disease associations
     */
    public static void getCircRNASimByCosineSim(String output, String Do2Circ_filename) {
        HashMap<String, ArrayList<String>> Circ2DO = readCirc2Disease(Do2Circ_filename, "rna");
        HashSet<String> circ_list = new HashSet<String>();
        BufferedWriter bufWriter = null;
        HashMap<String, Double> records = new HashMap<String, Double>();

        for (String lnckey: Circ2DO.keySet()) {
            circ_list.add(lnckey);
        }
        ArrayList<String> mapCirc_list = new ArrayList<String>(circ_list);
        System.out.println(mapCirc_list.size());//

        Set<String> DoNameSet = new HashSet<String>();
        for(String dopair : DoSimlist.keySet()){
            DoNameSet.add(dopair.split("vs")[0]);
            DoNameSet.add(dopair.split("vs")[1]);
        }
        System.out.println("DO number is " + DoNameSet.size());//

        try {
            bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace+output+".txt"));

            for (int i = 0; i < mapCirc_list.size(); i++) {
                String lnc_key1 = mapCirc_list.get(i);
                for (int j = i+1; j < mapCirc_list.size(); j++) {
                    String lnc_key2 = mapCirc_list.get(j);
                    Double cross_product = 0.0;
                    Double l1 = 0.0;
                    Double l2 = 0.0;
                    Double sim_sum1 = 0.0;
                    Double sim_sum2 = 0.0;
                    ArrayList<String> DOlist1 = Circ2DO.get(lnc_key1);
                    if (DOlist1 == null) {
                        DOlist1 = new ArrayList<String>();
                    }
                    ArrayList<String> DOlist2 = Circ2DO.get(lnc_key2);
                    if (DOlist2 == null) {
                        DOlist2 = new ArrayList<String>();
                    }
                    for(String doname: DoNameSet){
                        if(DOlist1.contains(doname) && DOlist2.contains(doname)){//a=1;b=1
                            cross_product = cross_product+1.0;
                            l1 = l1 + 1.0;
                            l2 = l2 + 1.0;
                        }else if(DOlist1.contains(doname)){//a=1;b=0
                            Double b = getDoMAXSim(doname, DOlist2);
                            cross_product = cross_product + b;
                            l1 = l1 + 1.0;
                            l2 = l2 + (b * b);//l2 + DoubleUtil.mul(b, b);
                        }else if(DOlist2.contains(doname)){//a=0;b=1
                            Double a = getDoMAXSim(doname, DOlist1);
                            cross_product = cross_product + a;
                            l1 = l1 + (a * a);
                            l2 = l2 + 1.0;
                        }else{
                            Double a = getDoMAXSim(doname, DOlist1);
                            Double b = getDoMAXSim(doname, DOlist2);
                            cross_product = cross_product + (a * b);
                            l1 = l1+(a*a);
                            l2 = l2+(b*b);
                        }
                    }
                    sim_sum1 = Math.sqrt(l1);
                    sim_sum2 = Math.sqrt(l2);
                    Double LNCSim = 0.0;
                    Double down = sim_sum1 * sim_sum2;
                    if(down != 0.0){
                        LNCSim = DoubleUtil.div(cross_product, down, 20);
                    }
                    if(LNCSim > 0){
                        records.put(lnc_key1 + "\t" + lnc_key2,LNCSim);
//                        bufWriter.write(lnc_key1 + "\t" + lnc_key2 + "\t" + LNCSim + "\n");
                        System.out.println(lnc_key1 + "\t" + lnc_key2 + "\t" + LNCSim);
                    }
                }
            }
			ArrayList<Map.Entry<String,Double>> list = new ArrayList<Map.Entry<String,Double>>(records.entrySet());
			Collections.sort(list,new Comparator<Map.Entry<String,Double>>() {
	            public int compare(Map.Entry<String, Double> o1,
	                    Map.Entry<String, Double> o2) {
	                return o2.getValue().compareTo(o1.getValue());
	            }

	        });
			for (int n = 0; n < list.size(); n++) {
				bufWriter.write(list.get(n).getKey() + "\t" + list.get(n).getValue() + "\n");
			}
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (bufWriter != null) {
                try {
                    bufWriter.close();
                } catch (IOException e) {
                    e.getStackTrace();
                }
            }
        }
    }

    /****
     * calculate the association degree between one disease and one circRNA
     * @param Do_key the DOID of disease
     * @param DOlist the list of disease related to circRNA
     * @return
     */
    public static Double getDoMAXSim(String Do_key, ArrayList<String> DOlist) {
        Double DOSim = 0.0;
        int num = 0;
        for (int i = 0; i < DOlist.size(); i++) {
            String key = Do_key + "vs" + DOlist.get(i);
            if(DoSimlist.containsKey(key)){
                if(DOSim < DoSimlist.get(key)){
                    DOSim = DOSim + DoSimlist.get(key);
                    num = num + 1;
                }
            }else if(Do_key.equals(DOlist.get(i))){
                DOSim = 0.999;
            }
        }
        if(num!=0){
            DOSim = DOSim/num;
        }
        return DOSim;
    }

    public static void main(String args[]) {
        //read disease similarity
        DoSimlist = filehandler.getFNSemSim();

        //calculate circRNA functional similarity by cosine similarity
        getCircRNASimByCosineSim("CircRNASim_1", "R2D_1_nodup_mapped");
        getCircRNASimByCosineSim("CircRNASim_2", "R2D_2_nodup_mapped");
        getCircRNASimByCosineSim("CircRNASim_3", "R2D_3_nodup_mapped");
    }
}
