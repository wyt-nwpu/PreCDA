package javasrc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;


public class PRank implements Runnable {
	final static double et = 0.7;//expression sim threshold
	final static double ct = 0.0;//cosin sim threshold
	final static double st = 0.0;//sim threshold
	final static String CIRCRNA_rootString = "E:/data/github/datasource/";
	final static String CIRCRNA_workplace = CIRCRNA_rootString+"output/";
	final static String prefix = "RNAUniCosSimBPAJAve_"+et+"_"+ct+"_"+st+"_";

	private HashMap<String, HashMap<String,Double>> terms;
	private String root;
	private double alpha;
	private int max_step;
	private int flag;
	
	public PRank(HashMap<String, HashMap<String, Double>> terms, String root, double alpha, int max_step, int flag) {
		super();
		this.terms = terms;
		this.root = root;
		this.alpha = alpha;
		this.max_step = max_step;
		this.flag = flag;
	}

	public static class ByValueComparator implements Comparator<String> {
        HashMap<String, Double> base_map;
 
        public ByValueComparator(HashMap<String, Double> base_map) {
            this.base_map = base_map;
        }
 
        public int compare(String arg0, String arg1) {
            if (!base_map.containsKey(arg0) || !base_map.containsKey(arg1)) {
                return 0;
            }
 
            if (base_map.get(arg0) < base_map.get(arg1)) {
                return 1;
            } else if (base_map.get(arg0) == base_map.get(arg1)) {
                return 0;
            } else {
                return -1;
            }
        }
    }

	@Override
	public void run() {
		// TODO Auto-generated method stub
		HashMap<String,Double> rank = new HashMap<String,Double>();
		HashMap<String,Double> tmp = new HashMap<String,Double>();

		BufferedWriter bufWriter = null;
		for(Entry<String, HashMap<String,Double>> entry : terms.entrySet()){
			rank.put(entry.getKey(), 0.0);
		}
		rank.put(root, 1.0);
		try {
			System.out.println("the process is " + root);
			File PRfile = new File(CIRCRNA_workplace + prefix+flag + "/"+max_step+"PR/PR_"+ root.substring(5) + ".txt");
			if(!PRfile.getParentFile().exists()){
				PRfile.getParentFile().mkdirs();
			}

			bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace + prefix+flag + "/"+max_step+"PR/PR_"+ root.substring(5) + ".txt"));
			for(int step = 0; step < max_step; step++){
				for(Entry<String, HashMap<String,Double>> entry : terms.entrySet()){
					tmp.put(entry.getKey(), 0.0);
				}
				for(Entry<String, HashMap<String,Double>> entry : terms.entrySet()){
					HashMap<String,Double> items = entry.getValue();

					for(Entry<String,Double> item : items.entrySet()){
//						System.out.println(item.getKey());
//						if(item.getKey().equals("hsa_circ_0000095")){
//							System.out.println();
//						}
						double tmpvalue = tmp.get(item.getKey()) + (alpha * rank.get(entry.getKey()) * item.getValue() / (1.0 * items.size()));
						tmp.put(item.getKey(), tmpvalue);
					}
				}
				tmp.put(root, tmp.get(root)+1-alpha);
				rank = (HashMap<String,Double>)tmp.clone();
				
			}
			
			//sort
			ArrayList<Map.Entry<String,Double>> list = new ArrayList<Map.Entry<String,Double>>(rank.entrySet());
			Collections.sort(list,new Comparator<Map.Entry<String,Double>>() {
	            //down
	            public int compare(Entry<String, Double> o1,
	                    Entry<String, Double> o2) {
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
                    e.printStackTrace();
                }
            }
        }
	}
}
