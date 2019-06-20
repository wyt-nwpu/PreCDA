package javasrc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import static javasrc.circRNAsim.readCirc2Disease;
import static javasrc.filehandler.createRNAIDMap;

/***
 * calculate circRNA expression similarity
 */
public class ExpressProfile {
	final static String CIRCRNA_rootString = "PreCDA/datasource/";
	final static String CIRCRNA_workplace = CIRCRNA_rootString+"output/";

	/***
	 * unify the name of circRNAs and extract the expression profiles
	 */
	public static void createExpressionMapBycircPedia(){
		HashMap<String,String> results = new HashMap<>();
		File file_circLoc1 = new File(CIRCRNA_rootString+"circRNA/circBase/hsa_hg19_circRNA.txt");
		BufferedReader bufReader_circLoc1 = null;

		File file_circLoc2 = new File(CIRCRNA_rootString+"circRNA/circFunBase/Homo_sapiens_circ.txt");
		BufferedReader bufReader_circLoc2 = null;

		File file_circLoc3 = new File(CIRCRNA_rootString+"circRNA/CircR2Disease_1/CircR2Disease_circRNA-disease_associations.txt");
		BufferedReader bufReader_circLoc3 = null;

		File file_express = new File(CIRCRNA_rootString+"circRNA/circpedia/expression.txt");
		BufferedReader bufReader_express = null;

		BufferedWriter bufWriter = null;
		try {
			bufReader_circLoc1 = new BufferedReader(new FileReader(file_circLoc1));

			String temp;
			String[] line_strs;
			int line = 1;
			while ((temp = bufReader_circLoc1.readLine()) != null) {
				if(line > 1){
					line_strs = temp.split("\t");
					String keyvalue = line_strs[0]+":"+line_strs[1]+"-"+line_strs[2];
					results.put(keyvalue, line_strs[4]);
				}
				line++;
			}

			bufReader_circLoc2 = new BufferedReader(new FileReader(file_circLoc2));
			line = 1;
			while ((temp = bufReader_circLoc2.readLine()) != null) {
				if(line > 1){
					line_strs = temp.split("\t");
					if (!results.containsKey(line_strs[1])){
						results.put(line_strs[1],line_strs[0]);
					}
				}
				line++;
			}

			bufReader_circLoc3 = new BufferedReader(new FileReader(file_circLoc3));
			line = 1;
			while ((temp = bufReader_circLoc3.readLine()) != null) {
				if(line > 1){
					line_strs = temp.split("\t");
					String keyloc = null;

					if(line_strs.length > 8) {
						String species_match = line_strs[8].toLowerCase();
						if(species_match.equals("human")) {
							if ((!line_strs[1].equals("N/A")) && line_strs[1].length() > 7) {
								keyloc = line_strs[1].replaceAll(" ", "").replaceAll("\\|", "-");
							}

							String[] cids = line_strs[0].replaceAll("\\|", "-").split("/");

							if (!results.containsKey(keyloc)) {
								if(cids.length>0){
									results.put(keyloc,cids[0]);
								}
							}
						}
					}
				}
				line++;
			}

			bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_rootString+"circRNA/circpedia/expression_mapped.txt"));

			HashMap<String,String> rnaidmap = createRNAIDMap();
			bufReader_express = new BufferedReader(new FileReader(file_express));
			while ((temp = bufReader_express.readLine()) != null) {
				line_strs = temp.split("\t");
				String loc = line_strs[0];
				Double val = Double.parseDouble(line_strs[1]);
				String sample = line_strs[2];
				if(results.containsKey(loc)){
					if(rnaidmap.containsKey(results.get(loc))){
						bufWriter.write(rnaidmap.get(results.get(loc))+"\t"+val+"\t"+sample+"\n");
					}else{
						bufWriter.write(results.get(loc)+"\t"+val+"\t"+sample+"\n");
					}
				}
			}

		}catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (bufReader_circLoc1 != null || bufReader_circLoc2 != null || bufReader_circLoc3 != null || bufWriter != null || bufReader_express!=null) {
				try {
					bufReader_circLoc1.close();
					bufReader_circLoc2.close();
					bufReader_circLoc3.close();
					bufReader_express.close();
					bufWriter.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/***
	 * read the expression profile of circPedia
	 * @return hashmap of cirdRNA and expression values
	 */
	public static HashMap<String, ArrayList<Double>> readCircRNAExpProfileBycircPedia(){
		File file_CircSAM = new File(CIRCRNA_rootString+"circRNA/circpedia/samples.txt");
		BufferedReader bufReader_CircSAM = null;
		ArrayList<String> samplelist = new ArrayList<String>();

		File file_CircEXP = new File(CIRCRNA_rootString+"circRNA/circpedia/expression_mapped.txt");
		BufferedReader bufReader_CircEXP = null;

		HashMap<String, ArrayList<Double>> RNASetEXPs = new HashMap<String, ArrayList<Double>>();
		try {
			bufReader_CircSAM = new BufferedReader(new FileReader(file_CircSAM));
			String temp = null;
			String[] line_strs;
			while ((temp = bufReader_CircSAM.readLine()) != null) {
				line_strs = temp.split("\t");
				samplelist.add(line_strs[0]);
			}

			bufReader_CircEXP = new BufferedReader(new FileReader(file_CircEXP));
			while ((temp = bufReader_CircEXP.readLine()) != null) {
				line_strs = temp.split("\t");
				String rnaname = line_strs[0];
				Double val = Double.parseDouble(line_strs[1]);
				String sample = line_strs[2];

				if(RNASetEXPs.containsKey(rnaname)){
					ArrayList<Double> rnasamples = RNASetEXPs.get(rnaname);
					for (int i=0; i<samplelist.size(); i++){
						if(samplelist.get(i).equals(sample)){
							rnasamples.set(i,val);
						}
					}
				}else{
					ArrayList<Double> rnasamples = new ArrayList<Double>();
					for (int i=0; i<samplelist.size(); i++){
						if(samplelist.get(i).equals(sample)){
							rnasamples.add(val);
						}else{
							rnasamples.add(0.0);
						}
					}
					RNASetEXPs.put(rnaname, rnasamples);
				}
			}
		} catch (Exception e) {
			// TODO: handle exception
			System.out.println("readLncRNAExpProfile Exception e: " + e);
			e.printStackTrace();
		} finally {
			if (bufReader_CircSAM != null || bufReader_CircEXP != null) {
				try {
					bufReader_CircSAM.close();
					bufReader_CircEXP.close();
				} catch (IOException e) {
					e.getStackTrace();
				}
			}
		}
		return RNASetEXPs;
	}

	/***
	 * read the expression profile of circBase
	 * @return hashmap of cirdRNA and expression values
	 */
	public static HashMap<String, ArrayList<Double>> readCircRNAExpProfileBycircBase() {
		File file_CircSAM = new File(CIRCRNA_rootString+"circRNA/circBase/samples.txt");
		BufferedReader bufReader_CircSAM = null;
		ArrayList<String> samplelist = new ArrayList<String>();
		File file_CircEXP = new File(CIRCRNA_rootString+"circRNA/circBase/expression.txt");
		BufferedReader bufReader_CircEXP = null;
		HashMap<String, ArrayList<Double>> RNASetEXPs = new HashMap<String, ArrayList<Double>>();
		try {
			bufReader_CircSAM = new BufferedReader(new FileReader(file_CircSAM));
			String temp = null;
			String[] line_strs;
			while ((temp = bufReader_CircSAM.readLine()) != null) {
				line_strs = temp.split("\t");
				samplelist.add(line_strs[0]);
			}

			bufReader_CircEXP = new BufferedReader(new FileReader(file_CircEXP));
			while ((temp = bufReader_CircEXP.readLine()) != null) {
				line_strs = temp.split("\t");
				String rnaname = line_strs[0];
				String[] samples = line_strs[1].split(", ");
				ArrayList<Double> rnasamples = new ArrayList<Double>();

				for (int i=0; i<samplelist.size(); i++){
					int iscontain = 0;
					for (int j=0; j<samples.length;j++) {
						if(samplelist.get(i).equals(samples[j])){
							iscontain = 1;
						}
					}
					if(iscontain==1){
						rnasamples.add(1.0);
					}else{
						rnasamples.add(0.0);
					}
				}
				RNASetEXPs.put(rnaname, rnasamples);
			}
		} catch (Exception e) {
			// TODO: handle exception
			System.out.println("readLncRNAExpProfile Exception e: " + e);
			e.printStackTrace();
		} finally {
            if (bufReader_CircSAM != null || bufReader_CircEXP != null) {
                try {
					bufReader_CircSAM.close();
					bufReader_CircEXP.close();
                } catch (IOException e) {
                    e.getStackTrace();
                }
            }
        }
		return RNASetEXPs;
	}
	
	/****
	 * calculate the Spearman's rank correlation coefficient between two circRNAs
	 * @param l1 the expression values of circRNA a
	 * @param l2 the expression values of circRNA b
	 * @return the Spearman's rank correlation coefficient
	 */
	public static Double Spearman(ArrayList<Double> l1, ArrayList<Double> l2) {
		int size1 = l1.size();
		ArrayList<Double> rank1 = new ArrayList<Double>();
		for (int i = 0; i < size1; i++) {
			int up = 0;
			int down = 0;
			for (int j = 0; j < size1; j++) {
				if (i != j) {
					if (l1.get(i) < l1.get(j)) {
						up = up + 1;
					}else if (l1.get(i) > l1.get(j)) {
						down = down + 1;
					}
				}
			}
			int equalnum = size1 - up - down;
			double sum = 0.0;
			for (int k = 1; k < equalnum + 1; k++) {
				sum = sum + up + k;
			}
			rank1.add(sum/equalnum);
		}
		int size2 = l2.size();
		ArrayList<Double> rank2 = new ArrayList<Double>();
		for (int i = 0; i < size2; i++) {
			int up = 0;
			int down = 0;
			for (int j = 0; j < size2; j++) {
				if (i != j) {
					if (l2.get(i) < l2.get(j)) {
						up = up + 1;
					}else if (l2.get(i) > l2.get(j)) {
						down = down + 1;
					}
				}
			}
			int equalnum = size2 - up - down;
			double sum = 0.0;
			for (int k = 1; k < equalnum + 1; k++) {
				sum = sum + up + k;
			}
			rank2.add(sum/equalnum);
		}
		
		if (size1 != size2 || rank1.size() != rank2.size() || rank1.size() != size1) {
			System.out.println("ERROR");
			return 0.0;
		}
		double sq_sum = 0.0;
		for (int i = 0; i < size2; i++) {
			sq_sum = sq_sum + (rank1.get(i) - rank2.get(i))*(rank1.get(i) - rank2.get(i));
		}
		
		double p = 1 - (6*sq_sum)/(size2*size2*size2 - size2);
		return p;
	}

	/****
	 * get the circRNA similarity by Spearman
	 * @param output the output file name
	 * @param Do2Circ_filename the file of the circRNA-disease associations
	 * @param threshold the threshold for the similarity between two circRNAs
	 * @param dbname circRNA database name
	 */
	public static void getCircRNASimBySpearman(String output, String Do2Circ_filename, double threshold, String dbname){
		HashMap<String, ArrayList<Double>> CircRNAEXPs= new HashMap<String, ArrayList<Double>>();
		if(dbname.equals("circPedia")){
			CircRNAEXPs = readCircRNAExpProfileBycircPedia();
		}else{
			CircRNAEXPs = readCircRNAExpProfileBycircBase();
		}
		HashMap<String, ArrayList<String>> Circ2DO = readCirc2Disease(Do2Circ_filename, "rna");

		BufferedWriter bufWriter = null;
		HashMap<String, Double> records = new HashMap<String, Double>();

		ArrayList<String> circ_list = new ArrayList<String>(Circ2DO.keySet());
		System.out.println(circ_list.size());

		for (int i = 0; i < circ_list.size(); i++) {
			String circrna1 = circ_list.get(i);
			for (int j = i + 1; j < circ_list.size(); j++) {
				String circrna2 = circ_list.get(j);
				if (CircRNAEXPs.containsKey(circrna1) && CircRNAEXPs.containsKey(circrna2)) {
					double simbyex = Spearman(CircRNAEXPs.get(circrna1), CircRNAEXPs.get(circrna2));
					if (simbyex > threshold) {
						records.put(circrna1 + "\t" + circrna2, simbyex);
					}
				}
			}
		}
		//sort
		ArrayList<Map.Entry<String,Double>> list = new ArrayList<Map.Entry<String,Double>>(records.entrySet());
		Collections.sort(list,new Comparator<Map.Entry<String,Double>>() {
			//down
			public int compare(Entry<String, Double> o1,
							   Entry<String, Double> o2) {
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
		//extract the expression profiles from circPedia
		createExpressionMapBycircPedia();
		//calculate the circRNA expression similarity from circPedia
		getCircRNASimBySpearman("CircRNAEXSimPedia_1", "R2D_1_nodup_mapped",0.0, "circPedia");
		getCircRNASimBySpearman("CircRNAEXSimPedia_2", "R2D_2_nodup_mapped",0.0, "circPedia");
		getCircRNASimBySpearman("CircRNAEXSimPedia_3", "R2D_3_nodup_mapped",0.0, "circPedia");

		//calculate the circRNA expression similarity from circBase
		getCircRNASimBySpearman("CircRNAEXSim_1", "R2D_1_nodup_mapped",0.0, "circBase");
		getCircRNASimBySpearman("CircRNAEXSim_2", "R2D_2_nodup_mapped",0.0, "circBase");
		getCircRNASimBySpearman("CircRNAEXSim_3", "R2D_3_nodup_mapped",0.0, "circBase");
	}
}
