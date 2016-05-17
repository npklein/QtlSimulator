package QtlSimulator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;

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
		infoWriter.write("QTL_name\tcoefficientCelltype\tcoefficientGenotype\tcoefficientInteraction");
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
		int genotypeCoefficientGroups = 1;
		int qtlCounter = 1;
		//The B2 is varied from -3 to +3 to have 3 different groups of QTL strengths (-1 down, 1 nothing, 1 up
		for(int j = 0; j < genotypeCoefficientGroups; j++){ 
			// make a normal distribution around the given mean distribution so that each gene has a different (altho close to given) mean expression distribution
			// the SD of this distribution is large as genes can have a large difference in expression levels
			NormalDistribution noiseDistribution = new NormalDistribution(commandLineOptions.getNoise(), commandLineOptions.getNoise()/4);
			NormalDistribution coeficientGenotypeDistribution = new NormalDistribution(-3 + (j*3), 0.5);
			double coeficientGenotype = coeficientGenotypeDistribution.sample();
			// Make i number of QTLs
			for (int i = 0; i < commandLineOptions.getNumberOfQtls(); i++){
				// the coefficientCelltype is for simulating celltype specific effects
				// NormalDistribution normalDistribution =  new NormalDistribution(0, 2);
				NormalDistribution coefficientInteractionDistribution = new NormalDistribution(3, 0.5);
				NormalDistribution coefficientCelltypeDistribution = new NormalDistribution(1, 0.5);
				double coefficientCelltype = coefficientCelltypeDistribution.sample();
				double coefficientInteraction = coefficientInteractionDistribution.sample();
				// write the QTL name in first column
				expressionWriter.write("QTL_"+Integer.toString(qtlCounter));
				genotypeWriter.write("QTL_"+Integer.toString(qtlCounter));
				// loop over the number of samples and calculate expression level
				for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
					// Base expression is in poison distribution with mean from command line option -e
					// Get the expression level for current QTL and current Sample. SD is low as between samples should not be high difference of "base" expression level 
					// (variance should be mostly in the genotype and interaction coefficient)
					double noiseSample = noiseDistribution.sample();
					NormalDistribution noise = new NormalDistribution(noiseSample, Math.sqrt(Math.pow(noiseSample/6.0,2)));
					double simulatedExpression = noise.sample();
					// genotype randomly 0, 1 or 2
					int genotype = new Random().ints(1, 0, 3).findFirst().getAsInt();
					for (int z = 0; z <  commandLineOptions.getCellcountPercentages().length; z++){
						//if(z > 0){
							// only have a QTL effect on the first celltype, to simplify matters
						//	coeficientGenotype = 0;
						//}
						double cellcountFactor = (cellcountsPerSample.get("celltype_"+Integer.toString(z))/commandLineOptions.getCellcountPercentages()[z]);
						/*
						System.out.printf("Simulated expression: %f\ncoefficientCelltype: %f\ncoefficientGenotype:%f\ncellcountFactor: %f\ngenotype: %d\ncoefficientInteraction: %f\n"+
										  "genotype: %d\ncellcountFactor: %f\nsimulatedExpression: %f\n", simulatedExpression,coefficientCelltype, coeficientGenotype, cellcountFactor, 
										  genotype, coefficientInteraction, genotype, cellcountFactor, simulatedExpression);
						System.out.printf("(coefficientCelltype * cellcountFactor): %f\n(coeficientGenotype * genotype): %f\n(coefficientInteraction * genotype * cellcountFactor): %f\n", 
								(coefficientCelltype * cellcountFactor),(coeficientGenotype * genotype), (coefficientInteraction * genotype * cellcountFactor), 
								  genotype, coefficientInteraction, genotype, cellcountFactor, simulatedExpression);
						System.out.println("-------------------------------");*/
						simulatedExpression += (coefficientCelltype * cellcountFactor) + 
								   (coeficientGenotype * genotype) + (coefficientInteraction * genotype * cellcountFactor);
					}
					//System.out.printf("Simulated expression: %f\n", simulatedExpression);
					//System.out.println("----------------------------------");
					//System.exit(0);
					// expression evel cant be lower than 0, if negative give it random number between 0 and 1 (reflects real life higher level of 0-1 genes as well
					if (simulatedExpression < 0){
						simulatedExpression = 0 + Math.random();
					}
					//simulatedExpression = Math.pow(2, simulatedExpression);
					expressionWriter.write("\t"+Double.toString(simulatedExpression));
					genotypeWriter.write("\t"+Integer.toString(genotype));
				}
				infoWriter.write("QTL_"+Integer.toString(qtlCounter)+"\t"+Double.toString(coefficientCelltype) +
								 "\t"+Double.toString(coeficientGenotype) +"\t"+Double.toString(coefficientInteraction));
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
			// only write header once, before first sample
			cellcountWriter.write("\tcelltype_"+Integer.toString(c));
		}
		cellcountWriter.newLine();
		 HashMap<String, Double> cellcountsPerSample = new HashMap<String, Double>();

		for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
			cellcountWriter.write("sample_"+Integer.toString(s));
			int c = 0;
			for (double cellcount : commandLineOptions.getCellcountPercentages()){
				// this stays the same for every sample, calculate the celltype % per sample now
				// TODO: normal distribution with cellcount as mean and SD as cellcount/10, should probably be changed to more appropriate
				NormalDistribution normalDistribution =  new NormalDistribution(cellcount, cellcount/10);
				double normalDistributionCellcount =  normalDistribution.sample();
				// cellcount can not be < 0, if < 0 make it 0 + small random number between 0.01 and 0.91
				if (normalDistributionCellcount < 0){
					normalDistributionCellcount = 0 + Math.random()+0.01;
				}
				cellcountWriter.write("\t");
				cellcountWriter.write(Double.toString(normalDistributionCellcount));
				cellcountsPerSample.put("celltype_"+Integer.toString(c), normalDistributionCellcount);
				c++;
			}
			cellcountWriter.newLine();
		}
		cellcountWriter.close();
		return(cellcountsPerSample);
	}
}
