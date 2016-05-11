package QtlSimulator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class QtlSimulator {
	private static CommandLineOptions commandLineOptions = new CommandLineOptions(); 
	public static void main(String[] args) throws Exception {
		/*
		 * Simulate QTL data with 
		 * random number as basic noise + ( B1 value with noise distribution * cell type) + ( value with noise distribution * [-]GT )   lymph
		 * + random number as basic noise + (value with noise distribution * cell type) + ( value with noise distribution * [-]GT ) mono
		 * + random number as basic noise + (value with noise distribution * cell type) + ( value with noise distribution * [-]GT ) neut
		 * + random number as basic noise + (value with noise distribution * cell type) + ( value with noise distribution * [-]GT ) eos
		 *  <- simulated expression 
		 * 
		 * @param args List of command line arguments
		 */
		commandLineOptions.parseCommandLine(args);
		simulateQTLs();
	}

	public static void simulateQTLs() throws IOException{
		PoissonDistribution poisonDistributionNoise= new PoissonDistribution(commandLineOptions.getNoise());
		
		new File(commandLineOptions.getOutfolder()).mkdirs();
		File simulatedExpressionFile = new File(commandLineOptions.getOutfolder()+"/simulatedExpression.csv");
		FileOutputStream simulatedExpressionStream = new FileOutputStream(simulatedExpressionFile);
		BufferedWriter expressionWriter = new BufferedWriter(new OutputStreamWriter(simulatedExpressionStream));
		File simulatedGenotypeFile = new File(commandLineOptions.getOutfolder()+"/simulatedGenotypes.csv");
		FileOutputStream genotypeStream = new FileOutputStream(simulatedGenotypeFile);
		BufferedWriter genotypeWriter = new BufferedWriter(new OutputStreamWriter(genotypeStream));
		File infoFile = new File(commandLineOptions.getOutfolder()+"/info.csv");
		FileOutputStream infoStream = new FileOutputStream(infoFile);
		BufferedWriter infoWriter = new BufferedWriter(new OutputStreamWriter(infoStream));
		infoWriter.write("QTL_name\tcoefficientCelltype\tcoefficientGenotype");
		infoWriter.newLine();
		// write headers with sample names
		for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
			expressionWriter.write("\tsample_"+Integer.toString(s));
			genotypeWriter.write("\tsample_"+Integer.toString(s));
		}
		expressionWriter.newLine();
		genotypeWriter.newLine();
		
		HashMap<String, Double> cellcountsPerSample = writeCellCountFile();
		/*
		 * Simulate expression, genotype and cellcount %. The cellcount % are distributions around the pre-selected means. The expression
		 * is y = error + (B1 * cellcount) + (B2 * genotype)
		 */
		int genotypeCoefficientGroups = 3;
		int qtlCounter = 1;
		//The B2 is varied from -10 to +10 to have 20 different groups of QTL strengths (-20 very strong down, 0 nothing, 20 very strong up
		for(int j = 0; j < genotypeCoefficientGroups; j++){ 
			double coeficientGenotype = j-(genotypeCoefficientGroups/2);
			// Make i number of QTLs
			for (int i = 0; i < commandLineOptions.getNumberOfQtls(); i++){
				// the coefficientCelltype is for simulating some genes being correlated to the celltype, chance that it is correlated is taken from a normal distribution				
				// for now have no correlations between celltype and genes
				// NormalDistribution normalDistribution =  new NormalDistribution(0, 2);
				double coefficientCelltype = 1;
				// Error is in poison distribution with mean from command line option -e
				double noise = poisonDistributionNoise.sample();
				// write the QTL name in first column
				expressionWriter.write("QTL_"+Integer.toString(qtlCounter));
				genotypeWriter.write("QTL_"+Integer.toString(qtlCounter));
				// loop over the number of samples and calculate expression level
				for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
					double simulatedExpression = noise;
					// genotype randomly 0, 1 or 2
					int genotype = new Random().ints(1, 0, 3).findFirst().getAsInt();
					for (int z = 0; z <  commandLineOptions.getCellcountPercentages().length; z++){
						if(z > 0){
							// only have a QTL effect on the first celltype, to simplify matters
							coeficientGenotype = 0;
						}
						simulatedExpression += (coefficientCelltype * cellcountsPerSample.get("sample_"+Integer.toString(s))) + 
											   (coeficientGenotype * genotype);
					}
					simulatedExpression = Math.pow(2, simulatedExpression);
					expressionWriter.write("\t"+Double.toString(simulatedExpression));
					genotypeWriter.write("\t"+Integer.toString(genotype));
				}
				infoWriter.write("QTL_"+Integer.toString(qtlCounter)+"\t"+Double.toString(coefficientCelltype) +
								 "\t"+coeficientGenotype);
				qtlCounter++;
				infoWriter.newLine();
				expressionWriter.newLine();
				genotypeWriter.newLine();
			}
		}
		expressionWriter.close();
		genotypeWriter.close();
		infoWriter.close();
		System.out.printf("Outfiles written to: %s",commandLineOptions.getOutfolder());
	}
	
	public static HashMap<String, Double> writeCellCountFile() throws IOException{
		File simulatedCellcountFile = new File(commandLineOptions.getOutfolder()+"/simulatedCellcounts.csv");
		FileOutputStream cellcountStream = new FileOutputStream(simulatedCellcountFile);
		BufferedWriter cellcountWriter = new BufferedWriter(new OutputStreamWriter(cellcountStream));
		// write headers with celltypes
		for (int c = 0; c < commandLineOptions.getCellcountPercentages().length; c++){
			c ++;
			// only write header once, before first sample
			cellcountWriter.write("\tcelltype_"+Integer.toString(c));
		}
		cellcountWriter.newLine();
		 HashMap<String, Double> cellcountsPerSample = new HashMap<String, Double>();

		for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
			cellcountWriter.write("sample_"+Integer.toString(s)+"\t");
			for (double cellcount : commandLineOptions.getCellcountPercentages()){
				// since this stays the same for every sample, calculate the celltype % per sample now
			
				// TODO: normal distribution with cellcount as mean and SD as cellcount/10, should probably be changed to more appropriate
				NormalDistribution normalDistribution =  new NormalDistribution(cellcount, cellcount/10);
				double normalDistributionCellcount =  normalDistribution.sample();
				// cellcount can not be < 0, if < 0 make it 0 + small random number between 0.01 and 0.91
				if (normalDistributionCellcount < 0){
					normalDistributionCellcount = 0 + Math.random()+0.01;
				}
				cellcountWriter.write(Double.toString(normalDistributionCellcount));
				cellcountWriter.write("\t");
				cellcountsPerSample.put("sample_"+Integer.toString(s), normalDistributionCellcount);
			}
			cellcountWriter.newLine();
		}
		cellcountWriter.close();
		return(cellcountsPerSample);
	}
}
