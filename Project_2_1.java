package project_2_1;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

public class Project_2_1 {

    public static void main(String args[]) {
    	
    	ArrayList<String> Mutations = new ArrayList<String>();
    	Mutations.add("APC");
    	Mutations.add("TP53");
    	Mutations.add("KRAS");
    	Mutations.add("PIK3CA");
    	Mutations.add("PTEN");
    	Mutations.add("ATM");
    	Mutations.add("MUC4");
    	Mutations.add("SMAD4");
    	Mutations.add("SYNE1");
    	Mutations.add("FBXW7");
        
        HashMap<String, Double> geneMap = new HashMap<String, Double>();
    	
    	ArrayList<String> People = new ArrayList<String>();
    	ArrayList<String> Survivors = new ArrayList<String>();
    	HashMap<String, ArrayList<String>> MutPpl = new HashMap<String, ArrayList<String>>();
    	
    	
    	
        Connection con;
        Statement stmt;
        ResultSet rs;
        ResultSetMetaData rsmd;
        
        /* Database credentials */
        
        String user = "cse4701";
        String password = "datamine";
        String host = "query.engr.uconn.edu";
        String port = "1521";
        String sid = "BIBCI";
        String url = "jdbc:oracle:thin:@" + host + ":" + port + ":" + sid;
        
        
        try {
            DriverManager.registerDriver(new oracle.jdbc.OracleDriver());
            con = DriverManager.getConnection(url, user, password);
            stmt = con.createStatement();

            String sql = "select Patient_id from Mutation ";
            sql = "select patient_id from clinical GROUP BY Patient_id";
            ResultSet ras = stmt.executeQuery(sql);
            ResultSetMetaData rasmd = ras.getMetaData();
            
            
            
            while (ras.next())
            {
            	People.add((String)ras.getObject(1));
            }

            for(String Mut : Mutations)
            {
            	ArrayList<String> Ppl = new ArrayList<String>();
            	sql = "select patient_id from Mutation WHERE Variant_Classification != 'Silent' AND GENE_Symbol = '"+ Mut +"' GROUP BY patient_id";
            	rs = stmt.executeQuery(sql);
            	rsmd = rs.getMetaData();
            	while(rs.next())
            	{
            		Ppl.add((String)rs.getObject(1));
            	}
            	MutPpl.put(Mut, Ppl);
            }
            
            sql = "select * from Clinical";
            
            rs = stmt.executeQuery(sql);
            rsmd = rs.getMetaData();
            
            while (rs.next())
            {
            	if(rs.getObject(6).equals("LIVING"))
            	{
            		Survivors.add((String)rs.getObject(1));
            	}
            }
        } catch (SQLException ex) {
            Logger.getLogger(Project_2_1.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        ArrayList<ArrayList<String>> Final = new ArrayList<ArrayList<String>>();
        for(String Person : People)
        {
        	ArrayList<String> Tuple = new ArrayList<String>();
        	Tuple.add(Person);
        	for(String Mut : Mutations)
        	{
        		if(MutPpl.get(Mut).contains(Person+"-01"))
        		{
        			Tuple.add("1");
        		}
        		else
        		{
        			Tuple.add("0");
        		}
        	}
        	if(Survivors.contains(Person.substring(0, Person.length())))
        	{
        		Tuple.add("1");
        	}
        	else
        	{
        		Tuple.add("0");
        	}
        	Final.add(Tuple);
        }
        
        
        
        try {
        	DriverManager.registerDriver(new com.mysql.jdbc.Driver());
        	Connection Conn = DriverManager.getConnection("jdbc:mysql://localhost:3306/cse4701", "root", "amenhotep4.5");
        	Statement Stmt = Conn.createStatement();
                Statement insertQuery = Conn.createStatement();
        	/* THESE TWO ONLY NEED TO BE RUN ONCE */
        	String sql = "CREATE TABLE IR_READY( Patient_ID VarCHAR(20), APC int, TP53 int, KRAS int, PIK3CA int, PTEN int, ATM int, MUC4 int, SMAD4 int, SYNE1 int, FBXW7 int, Status int, Primary KEY (Patient_ID) )";
        	Stmt.executeUpdate(sql);
        	
        	for(ArrayList<String> tuple : Final)
        	{
        		boolean First = true;
        		String S = "INSERT INTO IR_READY VALUES ('";
        		for(String field : tuple)
        		{
        			S += field;
        			if(First)
        			{
        				S += "'";
        				First = false;
        			}
        			S += ",";
        		}
        		S = S.substring(0, S.length() - 1);
        		S += ")";
        		Stmt.executeUpdate(S);
        	}
        	sql = "CREATE TABLE PROOF (APC int, TP53 int, KRAS int, PIK3CA int, PTEN int, ATM int, MUC4 int, SMAD4 int, SYNE1 int, FBXW7 int, Status int)";
        	Stmt.executeUpdate(sql);
        	sql = "INSERT INTO PROOF VALUES ("+MutPpl.get("APC").size()+","+MutPpl.get("TP53").size()+","+MutPpl.get("KRAS").size()+","+MutPpl.get("PIK3CA").size()+","+MutPpl.get("PTEN").size()+","+MutPpl.get("ATM").size()+","+MutPpl.get("MUC4").size()+","+MutPpl.get("SMAD4").size()+","+MutPpl.get("SYNE1").size()+","+MutPpl.get("FBXW7").size()+","+Survivors.size()+")";
        	Stmt.executeUpdate(sql);
                
                String[] genes = {"APC", "TP53", "KRAS", "PIK3CA", "PTEN", "ATM", "MUC4", "SMAD4", "SYNE1", "FBXW7"};
                
                double k = 627.0;
                double infoD = -(130/k)*(Math.log(130/k)/Math.log(2)) -(497/k)*(Math.log(497/k)/Math.log(2));// I(130,497)
                System.out.println("INFOD: " + infoD);
                int current = 1;
                List<gene> Rankedgene = new ArrayList<gene>();
                
                for (String gene : genes) {
                    int count = MutPpl.get(gene).size();
                    
                    double val = calculateTheIntermediate(gene, Stmt, count, Final, current, infoD);
                    gene thisGene = new gene(gene, val);
                    
                    Rankedgene.add(thisGene);
                    
                    current++;
                }
                
                ArrayList<gene> finalGenes = new ArrayList<gene>();
                Collections.sort(Rankedgene, new Comparator<gene>() {
        @Override
        public int compare(gene p1, gene p2) {
           return Double.valueOf(p1.getIg()).compareTo(p2.getIg());
        }
});
                
             sql = "CREATE TABLE IGTable (Gene_ID VARCHAR(20), IG double)";
               Stmt.executeUpdate(sql);
                
                for (int i = Rankedgene.size()-1; i > 4; i--) {
                    gene thisGene = Rankedgene.get(i);
                    finalGenes.add(thisGene);
                    String insertStuff = "INSERT INTO IGTable VALUES (" + "'"+ thisGene.geneName + "'" + "," + thisGene.ig +  ");  ";
                    
                    insertQuery.executeUpdate(insertStuff);
                    System.out.println("GENE: " + thisGene.geneName + " IG: " + thisGene.ig);
                }

        }catch (SQLException e)
        {
        	throw new IllegalStateException("Cannot connect", e);
        }
        
    }
    
    public static double calculateTheIntermediate(String gene, Statement stmt, int count, ArrayList<ArrayList<String>> Final, int current, double infoD) throws SQLException {
        
        ResultSet rs;
        
        int A = 0, B = 0, C = 0, D = 0;
        for (ArrayList list : Final) {
            String statusOfTheString = (String) list.get(11);
            int status = Integer.parseInt(statusOfTheString);
            String mutStrng = (String) list.get(current);
            int mutation = Integer.parseInt(mutStrng);
            
            
            if (status == 1 && mutation == 1) {
                A++;
            } else if (status == 0 && mutation == 1) {
                C++;
            } else if (status == 1 && mutation == 0) {
                B++;
            } else if (status == 0 && mutation == 0) {
                D++;
            }
            
        }
        
        int subtractTheNumber = 627 - count;
        double theFirstFraction = (count / 627.0);
        double theSecondFraction = (subtractTheNumber / 627.0);
        
        double inf1 = getInfo(A,C);
        double inf2 = getInfo(B,D);
        
        double intermediate = (theFirstFraction * inf1) + (theSecondFraction * inf2);
        
        double finalll = infoD - intermediate;
        
        return finalll;
    }
    
    public static double getInfo(double n, double p) {
        
        if (n == 0) {
            n = 0.0000000000000001;
        }
        if (p == 0) {
            p = 0.0000000000000001;
        }
        double k = n + p;
        
        double part1 = -(n/k)*(Math.log(n/k)/Math.log(2));
        double part2 = (p/k)*(Math.log(p/k)/Math.log(2));
        
        
        double sum = part1 - part2;
        return  sum;
    }
    
    
    
}

