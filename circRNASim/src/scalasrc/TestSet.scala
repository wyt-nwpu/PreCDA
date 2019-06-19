package scalasrc

import java.io.{File, PrintWriter}

import auc.ReadList

import scala.collection.mutable.ArrayBuffer
import scala.io.Source

object TestSet {
  val et = 0.7
  val ct = 0.0
  val st = 0.0
  val prefixforCosin = "RNAUniCosSimBPAJAve_"+et+"_"+ct+"_"+st+"_"

  val CIRCRNA_rootString = "E:/data/github/datasource/output/"

  /***
    * read the circRNA-disease association database
    * @param flag the flag of the circRNA-disease association database
    * @param keytype the key of HashMap
    * @return
    */
  def readDO2RNA(flag: Int, keytype: String): Map[String, ArrayBuffer[String]] = {
    var filepath = ""
    filepath = CIRCRNA_rootString + "R2D_"+flag+"_nodup_mapped.txt"
    val R2Dfile = Source.fromFile(filepath)
    var result : Map[String, ArrayBuffer[String]] = Map()
    var rnalist = new ArrayBuffer[String]()
    for(line <- R2Dfile.getLines()){
      val linestr = line.split("\t")
      if (keytype.equals("do")){
        if(result.contains(linestr(1))){
          val rnas = result(linestr(1))
          rnas += linestr(0)
        }else{
          val rnas = new ArrayBuffer[String]()
          rnas += linestr(0)
          result += (linestr(1) -> rnas)
        }
      }else{
        if(result.contains(linestr(0))){
          val rnas = result(linestr(0))
          rnas += linestr(1)
        }else{
          val rnas = new ArrayBuffer[String]()
          rnas += linestr(1)
          result += (linestr(0) -> rnas)
        }
      }
    }
    R2Dfile.close()
    result
  }

  /****
    * read the file name of the circRNA net
    * @param netpath the file name of the circRNA net
    * @return all the circRNAs
    */
  def getAllNetnodes(netpath : String): ArrayBuffer[String] ={
    val netfile = Source.fromFile(netpath)
    var rnalist = new ArrayBuffer[String]()
    for(line <- netfile.getLines()){
      val linestr = line.split("\t")
      rnalist += linestr(0)
      rnalist += linestr(1)
    }
    netfile.close()
    rnalist.distinct
  }

  def isListEqual(l0: ArrayBuffer[String], l1: ArrayBuffer[String], doname: String, RNA2Do_tar: Map[String,ArrayBuffer[String]], net15: ArrayBuffer[String]): Boolean = {
    var str15 = ""
    var str17 = ""
    var flag1 = 0
    var flag2 = 0
    for (o <- l0) {
      if (!l1.contains(o)) {
        str15 = str15 + ";" + o
        flag1 = 1
      }
    }
    for (o <- l1) {
      if (!l0.contains(o)) {
        if (RNA2Do_tar.contains(o) && net15.contains(o)) {
          str17 = str17 + ";" + o
          flag2 = 1
        }
      }
    }
    if (flag2 == 1) {
      System.out.println("=============================================")
      System.out.println(doname + "-----TAR=" + str15 + "; TEST=" + str17)
      return true
    }
    false
  }

  /***
    * get the circRNAs set to be tested
    * @param tar_flag The number of the circRNA-disease association database being tested
    * @param test_flag The number of the circRNA-disease association database as the reference database
    * @return circRNA set
    */
  def getTestSet(tar_flag:Int, test_flag:Int): Map[String,ArrayBuffer[ArrayBuffer[String]]] ={
//    val tar_flag = 2
//    val test_flag = 3
    val RNA2Do_tar = readDO2RNA(tar_flag, "rna")
    val Do2RNA_tar = readDO2RNA(tar_flag, "do")
    val RNA2Do_test = readDO2RNA(test_flag, "rna")
    val Do2RNA_test = readDO2RNA(test_flag, "do")
    val netnodes = getAllNetnodes(CIRCRNA_rootString+"Circ" + prefixforCosin + tar_flag + ".txt")
    println("target "+tar_flag + ": disease number = " + Do2RNA_tar.size + "; RNA number = " + RNA2Do_tar.size)
    println("test "+test_flag + ": disease number = " + Do2RNA_test.size + "; RNA number = " + RNA2Do_test.size)
    println("RNA netnode number = " + netnodes.size)
    var RNAResultSet : Map[String,ArrayBuffer[ArrayBuffer[String]]] = Map()

    var disease_num = 0
    Do2RNA_tar.foreach(t => {
      val do_tar = t._1
      val rnas_tar = t._2
      var relatednodesArrayList = new ArrayBuffer[ArrayBuffer[String]]()
      if(Do2RNA_test.contains(do_tar)){
        val rnas_test = Do2RNA_test(do_tar)
        var RNATestSet_temp = new ArrayBuffer[String]()

        if(isListEqual(rnas_tar,rnas_test,do_tar,RNA2Do_tar,netnodes)){
          println(do_tar + " TAR RNA number is " + rnas_tar.size + "; TEST RNA number is " + rnas_test.size)
          disease_num = disease_num+1
        }

        for(i <- 0 until rnas_test.size){
          if(RNA2Do_tar.contains(rnas_test(i)) && netnodes.contains(rnas_test(i)) && (!rnas_tar.contains(rnas_test(i)))){
            RNATestSet_temp += rnas_test(i)
          }
        }

        var isTestset = 0
        for (m <- 0 until rnas_tar.size){
          if(netnodes.contains(rnas_tar(m))){
            isTestset = 1
          }
        }

        if(RNATestSet_temp.size > 0 && isTestset==1){
          relatednodesArrayList += RNATestSet_temp
          var positiveList = new ArrayBuffer[String]()
          for(j <- 0 until rnas_test.size){
            if(netnodes.contains(rnas_test(j))){
              positiveList += rnas_test(j)
            }
          }
          relatednodesArrayList += positiveList
          RNAResultSet += (do_tar -> relatednodesArrayList)
          println(do_tar + " tested node number is " + relatednodesArrayList(0).size + "; positive node number is " + relatednodesArrayList(1).size)
        }
      }
    })

    println("=============================================")
    println("disease that can be tested =" + RNAResultSet.size + "; diseases that have novel relation with RNA in TEST=" + disease_num)
    println("net nodes number is " + netnodes.size)
    RNAResultSet
  }

  /***
    * get the list of the circRNAs potentially related to one disease
    * @param tar_flag The number of the circRNA-disease association database being tested
    * @param test_flag The number of the circRNA-disease association database as the reference database
    * @param iternum the number of iterations in PR
    * @param doname the name of the disease
    * @param testnodes the tested circRNA set
    * @param positiveSet the known circRNA set related to this disease
    */
  def getROCFileTargetDisease(tar_flag: Int, test_flag: Int, iternum: Int, doname: String, testnodes: ArrayBuffer[String], positiveSet: ArrayBuffer[String]): Unit = {
    val disease = doname.split(":")(1)
    val PRFilePath = CIRCRNA_rootString + prefixforCosin + tar_flag + "/" + iternum + "PR/PR_" + disease + ".txt"

    val PRFile = Source.fromFile(PRFilePath)
    var result : Map[String, ArrayBuffer[String]] = Map()
    var rnalist = new ArrayBuffer[String]()

    val diseaserocFile : File = new File(CIRCRNA_rootString + prefixforCosin + tar_flag + "/" + tar_flag + "V" + test_flag + "_i" + iternum + "ROC/ROC_" + disease + ".txt")
    if(!diseaserocFile.getParentFile.exists()){
      diseaserocFile.getParentFile.mkdirs()
    }
    val writer_diseaseroc = new PrintWriter(diseaserocFile)

    for(line <- PRFile.getLines()){
      val linestr = line.split("\t")

      if(testnodes.contains(linestr(0))){
        writer_diseaseroc.write(linestr(1)+"\t"+1+"\n")
      }else{
        if(!positiveSet.contains(linestr(0))){
          writer_diseaseroc.write(linestr(1)+"\t"+0+"\n")
        }
      }
    }
    PRFile.close()
    writer_diseaseroc.close()
  }

  def getAUC(tar_flag: Int, test_flag: Int, disease: String, iternum: Int): Double ={
    val localConfusion = ReadList.readFile(CIRCRNA_rootString + prefixforCosin + tar_flag + "/"+ tar_flag + "V" + test_flag + "_i" + iternum + "ROC/ROC_" + disease + ".txt", "list")
    val auc = localConfusion.calculateAUCROC()
    auc
  }

  def main(args: Array[String]): Unit = {
    val tar_flag = 3
    val test_flag = 9
    val iternum = 200


    val CircTestSet = getTestSet(tar_flag,test_flag)
    var sumauc = 0.0;

    val AUCreport : File = new File(CIRCRNA_rootString + prefixforCosin + tar_flag + "/" + tar_flag + "V" + test_flag + "_i" + iternum + "_aucreport.txt")
    if(!AUCreport.getParentFile.exists()){
      AUCreport.getParentFile.mkdirs()
    }
    val writer_rocreport = new PrintWriter(AUCreport)
    writer_rocreport.write(tar_flag + " based on " + test_flag + "\n")

    for (entrykey <- CircTestSet.keySet) {
      val doname = entrykey
      val relatednodesArrayList = CircTestSet(entrykey)
      val testnodes = relatednodesArrayList(0)
      val positiveSet = relatednodesArrayList(1)
      getROCFileTargetDisease(tar_flag, test_flag, iternum, doname, testnodes, positiveSet)

      val auc = getAUC(tar_flag, test_flag, doname.split(":")(1), iternum)
      writer_rocreport.write(doname + " is " + auc + "\n")
      sumauc = sumauc + auc
    }

    writer_rocreport.write("Average AUC is " + (sumauc/CircTestSet.size) + "\n")

    writer_rocreport.close()

  }
}
