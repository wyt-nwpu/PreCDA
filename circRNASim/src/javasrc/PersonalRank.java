package javasrc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class PersonalRank {

	final static double et = 0.7;//expression sim threshold
	final static double ct = 0.0;//cosin sim threshold
	final static double st = 0.0;//sim threshold
	final static String CIRCRNA_rootString = "E:/data/github/datasource/";
	final static String CIRCRNA_workplace = CIRCRNA_rootString+"output/";
	final static String prefixforCosin = "RNAUniCosSimBPAJAve_"+et+"_"+ct+"_"+st+"_";
	/***
	 * read the circNRA similarity network file
	 * @return hasmap
	 */
	public static HashMap<String, HashMap<String,Double>> readCircRNANet(int flag) {
		
		String NETPath = CIRCRNA_workplace + "Circ" + prefixforCosin + flag + ".txt";
		File hnfile = new File(NETPath);
        BufferedReader bufReader = null;
        HashMap<String, HashMap<String,Double>> hnterms = new HashMap<String, HashMap<String,Double>>();

		String R2DPath = CIRCRNA_workplace + "R2D_"+ flag +"_nodup_mapped.txt";
		File R2Dfile = new File(R2DPath);
		BufferedReader R2DbufReader = null;

        try {
            bufReader = new BufferedReader(new FileReader(hnfile));
            String temp = null;
            String[] line_strs;
            while ((temp = bufReader.readLine()) != null) {
            	line_strs = temp.split("\t");
            	if(!hnterms.containsKey(line_strs[0])){
            		HashMap<String,Double> terms_value = new HashMap<String,Double>();
            		terms_value.put(line_strs[1], Double.parseDouble(line_strs[2]));
            		hnterms.put(line_strs[0], terms_value);
            	}else{
            		HashMap<String,Double> terms_value = hnterms.get(line_strs[0]);
            		terms_value.put(line_strs[1], Double.parseDouble(line_strs[2]));
            	}
            	if(!hnterms.containsKey(line_strs[1])){
            		HashMap<String,Double> terms_value = new HashMap<String,Double>();
            		terms_value.put(line_strs[0], Double.parseDouble(line_strs[2]));
            		hnterms.put(line_strs[1], terms_value);
            	}else{
            		HashMap<String,Double> terms_value = hnterms.get(line_strs[1]);
            		terms_value.put(line_strs[0], Double.parseDouble(line_strs[2]));
            	}
            }

			R2DbufReader = new BufferedReader(new FileReader(R2Dfile));
			while ((temp = R2DbufReader.readLine()) != null) {
				line_strs = temp.split("\t");
				if(hnterms.containsKey(line_strs[0])){
					if(!hnterms.containsKey(line_strs[1])){
						HashMap<String,Double> terms_value = new HashMap<String,Double>();
						terms_value.put(line_strs[0], 1.0);
						hnterms.put(line_strs[1], terms_value);
					}else{
						HashMap<String,Double> terms_value = hnterms.get(line_strs[1]);
						terms_value.put(line_strs[0], 1.0);
					}
				}
			}
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (bufReader != null || R2DbufReader != null) {
                try {
                    bufReader.close();
					R2DbufReader.close();
                } catch (IOException e) {
                	e.printStackTrace();
                }
            }
        }
        return hnterms;
	}

	public static void main(String args[]) {

		int thread = 8;
		int flag = 3;
		int max_step = 200;
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(thread);
		
 		HashMap<String, HashMap<String,Double>> terms = readCircRNANet(flag);

 		try {
			for (Entry<String, HashMap<String, Double>> entry1 : terms.entrySet()) {
				if (entry1.getKey().indexOf("DOID:") != -1) {
					PRank p = new PRank(terms, entry1.getKey(), 0.85, max_step, flag);
					fixedThreadPool.execute(p);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			fixedThreadPool.shutdown();
		}

		fixedThreadPool.shutdown();
		System.out.println("finish the process!");
		
	}
}

class TaskException extends Exception {  
    public TaskException(String message) {  
        super(message);  
    }  
}
