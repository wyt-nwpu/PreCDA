package javasrc;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

/***
 * Processing source files
 */
public class filehandler {
    final static String CIRCRNA_rootString = "E:/data/github/datasource/";
    final static String CIRCRNA_workplace = CIRCRNA_rootString+"output/";

    /****
     * read disease similarity from FNSemSim
     * @return  disease similarity
     */
    public static HashMap<String, Double> getFNSemSim() {
        File file_FNsim = new File(CIRCRNA_rootString+"Disease/fndo_sim_sort_nomal.txt");
        BufferedReader bufReader_FNsim = null;
        HashMap<String, Double> DOSim = new HashMap<String, Double>();
        try {
            bufReader_FNsim = new BufferedReader(new FileReader(file_FNsim));
            String temp = null;
            String[] line_strs;
            while ((temp = bufReader_FNsim.readLine()) != null) {
                line_strs = temp.split("\t");
                Double values = Double.parseDouble(line_strs[2]);
                String key = line_strs[0] + "vs" + line_strs[1];
                String key_reverse = line_strs[1] + "vs" + line_strs[0];
                DOSim.put(key, values);
                DOSim.put(key_reverse, values);
            }
        } catch (Exception e) {
            // TODO: handle exception
            System.out.println("getFNSemSim Exception e: " + e);
            e.printStackTrace();
        } finally {
            if (bufReader_FNsim != null) {
                try {
                    bufReader_FNsim.close();
                } catch (IOException e) {
                    e.getStackTrace();
                }
            }
        }
        return DOSim;
    }

    /***
     * read the circRNA-disease file and match disease name to DOID
     * @param circ2disesase_filename file path of circRNA-disease association database
     * @return hashmap of circRNA and DOID
     */
    public static HashMap<String, HashSet<String>> getDO2CircRNA(String circ2disesase_filename) {

        File disease_name = new File(CIRCRNA_rootString+"Disease/DO_Medic_Allsynonym.txt");
        BufferedReader bufReader_dname = null;
        HashMap<String, ArrayList<String>> disease2name = new HashMap<String, ArrayList<String>>();
        File file_circDisease = new File(CIRCRNA_rootString+"circRNA/"+circ2disesase_filename);
        BufferedReader bufReader_CircDisease = null;
        HashMap<String, HashSet<String>> circdiseaseMap = new HashMap<String, HashSet<String>>();
        BufferedWriter bufWriter = null;
        BufferedWriter mapbufWriter = null;
        try {
            String temp = null;
            String[] line_strs;
            int line = 1;
            bufReader_dname = new BufferedReader(new FileReader(disease_name));
            while ((temp = bufReader_dname.readLine()) != null) {
                line_strs = temp.split("\t");
                if(!disease2name.containsKey(line_strs[0])){
                    ArrayList<String> values = new ArrayList<String>();
                    values.add(line_strs[1]);
                    disease2name.put(line_strs[0], values);
                }else{
                    ArrayList<String> values = disease2name.get(line_strs[0]);
                    values.add(line_strs[1]);
                }
            }

            bufReader_CircDisease = new BufferedReader(new FileReader(file_circDisease));

            if (circ2disesase_filename.equals("circRNADisease_2/2017-12-25.txt")){
                bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace+"R2D_2.txt"));

                while ((temp = bufReader_CircDisease.readLine()) != null) {
                    line_strs = temp.split("\t");

                    String circid1 = line_strs[5];
                    String circid2 = line_strs[6];
                    String circid3 = line_strs[7];
                    String doid = null;
                    String diseasename_match = line_strs[8].toLowerCase().replaceAll("'s", "").replaceAll("[^a-zA-Z0-9]", "");
                    String species_match = line_strs[10].toLowerCase();

                    if(species_match.equals("human")){
                        Iterator iter = disease2name.entrySet().iterator();
                        ArrayList<String> dnames = null;
                        while (iter.hasNext()) {
                            Map.Entry do_name = (Map.Entry) iter.next();
                            dnames = (ArrayList<String>) do_name.getValue();
                            for (int i = 0; i < dnames.size(); i++) {
                                String matchname = dnames.get(i).toLowerCase().replaceAll("'s", "").replaceAll("[^a-zA-Z0-9]", "");
                                if (diseasename_match.equals(matchname)) {
                                    doid = (String)do_name.getKey();
                                }
                            }
                        }

                        if (doid != null) {
                            if (!circid1.equals("-")) {
                                bufWriter.write(circid1 + "\t" + doid + "\n");
                                if (circdiseaseMap.containsKey(doid)) {
                                    HashSet<String> values = circdiseaseMap.get(doid);
                                    values.add(circid1);
                                }else{
                                    HashSet<String> values = new HashSet<String>();
                                    values.add(circid1);
                                    circdiseaseMap.put(doid, values);
                                }
                            }else if (!circid2.equals("-")) {
                                bufWriter.write(circid2 + "\t" + doid + "\n");
                                if (circdiseaseMap.containsKey(doid)) {
                                    HashSet<String> values = circdiseaseMap.get(doid);
                                    values.add(circid2);
                                }else{
                                    HashSet<String> values = new HashSet<String>();
                                    values.add(circid2);
                                    circdiseaseMap.put(doid, values);
                                }
                            }else if (!circid3.equals("")) {
                                bufWriter.write(circid3 + "\t" + doid + "\n");
                                if (circdiseaseMap.containsKey(doid)) {
                                    HashSet<String> values = circdiseaseMap.get(doid);
                                    values.add(circid3);
                                }else{
                                    HashSet<String> values = new HashSet<String>();
                                    values.add(circid3);
                                    circdiseaseMap.put(doid, values);
                                }
                            }
                        }
                    }

                    System.out.println(line++);
                }
            } else if (circ2disesase_filename.equals("CircR2Disease_1/CircR2Disease_circRNA-disease_associations.txt")){
                bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace+"R2D_1.txt"));

                while ((temp = bufReader_CircDisease.readLine()) != null) {
                    line_strs = temp.split("\t");

                    String circids = line_strs[0];
                    String doid = null;

                    if(line_strs.length > 5){
                        String diseasename_match = line_strs[5].toLowerCase().replaceAll("'s", "").replaceAll("[^a-zA-Z0-9]", "");
                        String species_match = line_strs[8].toLowerCase();

                        if(species_match.equals("human")) {
                            Iterator iter = disease2name.entrySet().iterator();
                            ArrayList<String> dnames = null;
                            while (iter.hasNext()) {
                                Map.Entry do_name = (Map.Entry) iter.next();
                                dnames = (ArrayList<String>) do_name.getValue();
                                for (int i = 0; i < dnames.size(); i++) {
                                    String matchname = dnames.get(i).toLowerCase().replaceAll("'s", "").replaceAll("[^a-zA-Z0-9]", "");
                                    if (diseasename_match.equals(matchname)) {
                                        doid = (String) do_name.getKey();
                                    }
                                }
                            }
                        }
                    }


                    if (doid != null) {
                        String[] ids = circids.replaceAll("\\|", "-").split("/");
                        if (ids.length>0) {
                            for (int k=0; k<ids.length; k++) {
                                bufWriter.write(ids[k] + "\t" + doid + "\n");
                                if (circdiseaseMap.containsKey(doid)) {
                                    HashSet<String> values = circdiseaseMap.get(doid);
                                    values.add(ids[k]);
                                } else {
                                    HashSet<String> values = new HashSet<String>();
                                    values.add(ids[k]);
                                    circdiseaseMap.put(doid, values);
                                }
                            }
                        }
                    }
                    System.out.println(line++);
                }
            }else if(circ2disesase_filename.equals("Circ2Disease_3/Circ2Disease_Association.txt")){
                bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace+"R2D_3.txt"));

                while ((temp = bufReader_CircDisease.readLine()) != null) {
                    line_strs = temp.split("\t");

                    String circid1 = line_strs[0];
                    String[] circids2 = null;
                    if(!line_strs[1].equals("N/A")){
                        circids2=line_strs[1].split(";");
                    }

                    String doid = null;

                    if(line_strs.length > 1){
                        String diseasename_match = line_strs[9].toLowerCase().replaceAll("'s", "").replaceAll("[^a-zA-Z0-9]", "");

                        Iterator iter = disease2name.entrySet().iterator();
                        ArrayList<String> dnames = null;
                        while (iter.hasNext()) {
                            Map.Entry do_name = (Map.Entry) iter.next();
                            dnames = (ArrayList<String>) do_name.getValue();
                            for (int i = 0; i < dnames.size(); i++) {
                                String matchname = dnames.get(i).toLowerCase().replaceAll("'s", "").replaceAll("[^a-zA-Z0-9]", "");
                                if (diseasename_match.equals(matchname)) {
                                    doid = (String)do_name.getKey();
                                }
                            }
                        }
                    }


                    if (doid != null) {
                        bufWriter.write(circid1 + "\t" + doid + "\n");
                        if (circdiseaseMap.containsKey(doid)) {
                            HashSet<String> values = circdiseaseMap.get(doid);
                            values.add(circid1);
                        }else{
                            HashSet<String> values = new HashSet<String>();
                            values.add(circid1);
                            circdiseaseMap.put(doid, values);
                        }

                        if(circids2!=null){
                            for(int i=0; i<circids2.length; i++){
                                bufWriter.write(circids2[0] + "\t" + doid + "\n");
                                if (circdiseaseMap.containsKey(doid)) {
                                    HashSet<String> values = circdiseaseMap.get(doid);
                                    values.add(circids2[0]);
                                }else{
                                    HashSet<String> values = new HashSet<String>();
                                    values.add(circids2[0]);
                                    circdiseaseMap.put(doid, values);
                                }
                            }
                        }

                    }
                    System.out.println(line++);
                }
            }


        } catch (Exception e) {
            // TODO: handle exception
            e.printStackTrace();
        } finally {
            if (bufReader_dname != null || bufReader_CircDisease != null || bufWriter != null) {
                try {
                    bufReader_CircDisease.close();
                    bufReader_dname.close();
                    bufWriter.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        return circdiseaseMap;
    }

    /***
     * create the map of circRNA synonyms
     * @return hashmap of circRNA synonyms
     */
    public static HashMap<String,String> createRNAIDMap(){
        HashMap<String,String> results = new HashMap<>();

        File file_circIDconvert = new File(CIRCRNA_rootString+"circRNA/circBase/hg19_circID_to_name.txt");
        BufferedReader bufReader_circIDconvert = null;

        File file_circLoc1 = new File(CIRCRNA_rootString+"circRNA/circBase/hsa_hg19_circRNA.txt");
        BufferedReader bufReader_circLoc1 = null;
        HashMap<String,String> results1 = new HashMap<>();

        File file_circLoc2 = new File(CIRCRNA_rootString+"circRNA/circFunBase/Homo_sapiens_circ.txt");
        BufferedReader bufReader_circLoc2 = null;

        File file_circIDconvert1 = new File(CIRCRNA_rootString+"circRNA/CircR2Disease_1/CircR2Disease_circRNA-disease_associations.txt");
        BufferedReader bufReader_circIDconvert1 = null;
        BufferedWriter bufWriter = null;
        try {
            bufReader_circIDconvert = new BufferedReader(new FileReader(file_circIDconvert));
            String temp;
            String[] line_strs;
            int line = 1;
            while ((temp = bufReader_circIDconvert.readLine()) != null) {
                if(line > 1){
                    line_strs = temp.split("\t");
                    results.put(line_strs[1],line_strs[0]);
                    results.put(line_strs[0],line_strs[0]);
                }
                line++;
            }
            bufReader_circLoc1 = new BufferedReader(new FileReader(file_circLoc1));
            line = 1;
            while ((temp = bufReader_circLoc1.readLine()) != null) {
                if(line > 1){
                    line_strs = temp.split("\t");
                    String keyvalue = line_strs[0]+":"+line_strs[1]+"-"+line_strs[2];

                    results1.put(keyvalue, line_strs[4]);
                }
                line++;
            }

            bufReader_circLoc2 = new BufferedReader(new FileReader(file_circLoc2));
            line = 1;
            while ((temp = bufReader_circLoc2.readLine()) != null) {
                if(line > 1){
                    line_strs = temp.split("\t");
                    if(results1.containsKey(line_strs[1])){
                        String rnaid = results1.get(line_strs[1]);
                        if(!rnaid.equals(line_strs[0])){
                            results.put(line_strs[0],rnaid);
                        }
                    }
                }
                line++;
            }

            bufReader_circIDconvert1 = new BufferedReader(new FileReader(file_circIDconvert1));
            line = 1;
            while ((temp = bufReader_circIDconvert1.readLine()) != null) {
                if(line > 1){
                    line_strs = temp.split("\t");
                    String keyloc = null;
                    if(line_strs.length > 5) {
                        if ((!line_strs[1].equals("N/A")) && line_strs[1].length() > 7) {
                            keyloc = line_strs[1].replaceAll(" ", "").replaceAll("\\|", "-");
                        }
                        String[] cids = line_strs[0].replaceAll("\\|", "-").split("/");
                        if (results1.containsKey(keyloc)) {
                            String rnaid = results1.get(keyloc);
                            for (int i = 0; i < cids.length; i++) {
                                if (!rnaid.equals(cids[i])) {
                                    results.put(cids[i], rnaid);
                                }
                            }
                        }
                    }
                }
                line++;
            }

            Set<String> keys= results.keySet();
            ArrayList<HashSet<String>> circRNASet = getcircRNASet();

            for(int i=0; i<circRNASet.size();i++){
                HashSet<String> cset = circRNASet.get(i);
                HashSet<String> cl = (HashSet<String>) cset.clone();
                cl.retainAll(keys);
                if(cl.size()>0){
                    String id = null;
                    for(String na : cl){
                        id = results.get(na);
                    }
                    for(String na : cset){
                        results.put(na,id);
                    }
                }else{
                    Iterator it = cset.iterator();
                    String id = (String) it.next();
                    for(String na : cset){
                        results.put(na,id);
                    }
                }
            }

            bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace+"maps.txt"));
            for (String key : results.keySet()){
                bufWriter.write(key+"\t"+results.get(key)+"\n");
            }

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (bufReader_circIDconvert != null || bufReader_circLoc1 != null || bufReader_circLoc2 != null || bufReader_circIDconvert1 != null || bufWriter != null) {
                try {
                    bufReader_circIDconvert.close();
                    bufReader_circLoc1.close();
                    bufReader_circLoc2.close();
                    bufReader_circIDconvert1.close();
                    bufWriter.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        return results;
    }

    /****
     * Unify the nomenclature of circRNAs
     *
     * @param filename circRNA-disease association file
     * @param rnaidmap map of circRNA synonyms
     */
    public static void handleR2Dfile(String filename,HashMap<String,String> rnaidmap){

        File file_R2D = new File(CIRCRNA_workplace+filename+".txt");
        BufferedReader bufReader_R2D = null;
        BufferedWriter bufWriter = null;

        HashSet<String> nodup = new HashSet<>();

        try{
            String temp;
            String[] line_strs;
            int line = 1;
            int num = 0;
            bufWriter = new BufferedWriter(new FileWriter(CIRCRNA_workplace+filename+"_nodup_mapped.txt"));

            bufReader_R2D = new BufferedReader(new FileReader(file_R2D));
            while ((temp = bufReader_R2D.readLine()) != null) {
                line_strs = temp.split("\t");
                String rnaid = line_strs[0];
                if(rnaidmap.containsKey(rnaid)){
                    nodup.add(rnaidmap.get(rnaid)+"\t"+line_strs[1]);
                    num++;
                }else{
                    nodup.add(line_strs[0]+"\t"+line_strs[1]);
                }

                line++;
            }
            System.out.println(filename+" mapped number: "+num);

            for(String str : nodup){
                bufWriter.write(str + "\n");
            }

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (bufReader_R2D != null || bufWriter != null) {
                try {
                    bufReader_R2D.close();
                    bufWriter.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public static ArrayList<HashSet<String>> getcircRNASet() {
        File file_1 = new File(CIRCRNA_rootString+"circRNA/CircR2Disease_1/CircR2Disease_circRNA-disease_associations.txt");
        BufferedReader bufReader_1 = null;

        File file_2 = new File(CIRCRNA_rootString+"circRNA/circRNADisease_2/2017-12-25.txt");
        BufferedReader bufReader_2 = null;

        File file_3 = new File(CIRCRNA_rootString+"circRNA/Circ2Disease_3/Circ2Disease_Association.txt");
        BufferedReader bufReader_3 = null;

        ArrayList<HashSet<String>> circRNASet = new ArrayList<HashSet<String>>();
        try {
            String temp = null;
            String[] line_strs;
            bufReader_1 = new BufferedReader(new FileReader(file_1));
            while ((temp = bufReader_1.readLine()) != null) {
                line_strs = temp.split("\t");
                if(line_strs.length>8){
                    String species_match = line_strs[8].toLowerCase();
                    if(species_match.equals("human")) {
                        String[] idvalues = line_strs[0].replaceAll(" ", "").replaceAll("\\|", "-").split("/");
                        HashSet<String> idset = new HashSet<>();
                        for (int i=0; i<idvalues.length; i++){
                            idset.add(idvalues[i]);
                        }
                        circRNASet.add(idset);
                    }
                }
            }

            bufReader_2 = new BufferedReader(new FileReader(file_2));
            while ((temp = bufReader_2.readLine()) != null) {
                line_strs = temp.split("\t");
                String circid1 = line_strs[5];
                String circid2 = line_strs[6];
                String circid3 = line_strs[7];
                String species_match = line_strs[10].toLowerCase();
                if(species_match.equals("human")){
                    HashSet<String> idset = new HashSet<>();
                    if (!circid1.equals("-")) {
                        idset.add(circid1);
                    }
                    if (!circid2.equals("-")) {
                        idset.add(circid2);
                    }
                    if (!circid3.equals("")) {
                        idset.add(circid3);
                    }

                    ArrayList<Integer> duplist = new ArrayList<>();
                    for (int j=0; j<circRNASet.size();j++){
                        HashSet<String> clj = (HashSet<String>) circRNASet.get(j).clone();
                        HashSet<String> clm = (HashSet<String>) idset.clone();
                        clm.retainAll(clj);
                        if (clm.size()>0){
                            duplist.add(j);
                        }
                    }
                    if (duplist.size()>0){
                        ArrayList<HashSet<String>> circRNASetclone = (ArrayList<HashSet<String>>)circRNASet.clone();
                        circRNASet.clear();
                        for (int n=0; n<circRNASetclone.size(); n++){
                            if(duplist.contains(n)){
                                idset.addAll(circRNASetclone.get(n));
                            }else{
                                circRNASet.add(circRNASetclone.get(n));
                            }
                        }
                    }

                    circRNASet.add(idset);
                }
            }

            bufReader_3 = new BufferedReader(new FileReader(file_3));
            while ((temp = bufReader_3.readLine()) != null) {
                line_strs = temp.split("\t");

                if(line_strs.length > 2) {
                    HashSet<String> idset = new HashSet<>();
                    String cid1=line_strs[0];
                    idset.add(cid1);

                    String[] cids = null;
                    if(!line_strs[1].equals("N/A")){
                        cids = line_strs[1].split(";");
                    }
                    if(cids!=null){
                        for (int i=0; i<cids.length; i++){
                            idset.add(cids[i]);
                        }
                    }

                    ArrayList<Integer> duplist = new ArrayList<>();
                    for (int j=0; j<circRNASet.size();j++){
                        HashSet<String> clj = (HashSet<String>) circRNASet.get(j).clone();
                        HashSet<String> clm = (HashSet<String>) idset.clone();
                        clm.retainAll(clj);
                        if (clm.size()>0){
                            duplist.add(j);
                        }
                    }
                    if (duplist.size()>0){
                        ArrayList<HashSet<String>> circRNASetclone = (ArrayList<HashSet<String>>)circRNASet.clone();
                        circRNASet.clear();
                        for (int n=0; n<circRNASetclone.size(); n++){
                            if(duplist.contains(n)){
                                idset.addAll(circRNASetclone.get(n));
                            }else{
                                circRNASet.add(circRNASetclone.get(n));
                            }
                        }
                    }
                    circRNASet.add(idset);
                }
            }
        } catch (Exception e) {
            // TODO: handle exception
            System.out.println("getFNSemSim Exception e: " + e);
            e.printStackTrace();
        } finally {
            if (bufReader_1 != null || bufReader_2 != null || bufReader_3 != null) {
                try {
                    bufReader_1.close();
                    bufReader_2.close();
                    bufReader_3.close();
                } catch (IOException e) {
                    e.getStackTrace();
                }
            }
        }
        return circRNASet;
    }

    public static void main(String args[]) {
        //read the file of the circRNA-disease association database
        getDO2CircRNA("circRNADisease_2/2017-12-25.txt");
        getDO2CircRNA("CircR2Disease_1/CircR2Disease_circRNA-disease_associations.txt");
        getDO2CircRNA("Circ2Disease_3/Circ2Disease_Association.txt");

        //create the map of circRNA synonyms
        HashMap<String,String> map = createRNAIDMap();

        //create the map between circRNA and DOID
        handleR2Dfile("R2D_1",map);
        handleR2Dfile("R2D_2",map);
        handleR2Dfile("R2D_3",map);
    }
}
