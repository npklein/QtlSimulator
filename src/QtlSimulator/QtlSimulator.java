package QtlSimulator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.LinkedListMultimap;
import com.google.common.collect.Multimap;

public class QtlSimulator {
	private static CommandLineOptions commandLineOptions = new CommandLineOptions(); 
	private static BufferedWriter expressionWriter;
	private static BufferedWriter genotypeWriter;
	private static BufferedWriter infoWriter;
    private static HashMap<String, BufferedWriter> expressionPerCelltypeWriterMap = new HashMap<String, BufferedWriter>();
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

	public static void initializeFiles() throws IOException{
		new File(commandLineOptions.getOutfolder()).mkdirs();
		File simulatedExpressionFile = new File(commandLineOptions.getOutfolder()+"/simulatedExpression.csv");
		FileOutputStream simulatedExpressionStream = new FileOutputStream(simulatedExpressionFile);
		expressionWriter = new BufferedWriter(new OutputStreamWriter(simulatedExpressionStream));
		for(int c = 0; c < commandLineOptions.getCellcountPercentages().length; c++){
			File simulatedExpressionCelltypeFile = new File(commandLineOptions.getOutfolder()+"/simulatedExpression_Celltype_"+Integer.toString(c)+".csv");
			FileOutputStream simulatedExpressionCelltypeStream = new FileOutputStream(simulatedExpressionCelltypeFile);
			BufferedWriter expressionCelltypeWriter = new BufferedWriter(new OutputStreamWriter(simulatedExpressionCelltypeStream));
			expressionPerCelltypeWriterMap.put("celltype_"+Integer.toString(c), expressionCelltypeWriter);
			// write headers with sample names
			for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
				expressionCelltypeWriter.write("\tsample_"+Integer.toString(s));
			}
			expressionCelltypeWriter.newLine();
		}
		File simulatedGenotypeFile = new File(commandLineOptions.getOutfolder()+"/simulatedGenotypes.csv");
		FileOutputStream genotypeStream = new FileOutputStream(simulatedGenotypeFile);
		genotypeWriter = new BufferedWriter(new OutputStreamWriter(genotypeStream));
		File infoFile = new File(commandLineOptions.getOutfolder()+"/info.csv");
		FileOutputStream infoStream = new FileOutputStream(infoFile);
		infoWriter = new BufferedWriter(new OutputStreamWriter(infoStream));
		
		infoWriter.write("QTL_name");
		for(int c = 0; c < commandLineOptions.getCellcountPercentages().length; c++){
			infoWriter.write("\tinteractionCelltype_"+Integer.toString(c));
		}
		infoWriter.write("\tcelltypeCoefficient\tgenotypeCoefficient");
		infoWriter.newLine();
		// write headers with sample names
		for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
			expressionWriter.write("\tsample_"+Integer.toString(s));
			genotypeWriter.write("\tsample_"+Integer.toString(s));
		}

		expressionWriter.newLine();
		genotypeWriter.newLine();

	}
	
	public static void simulateQTLs() throws IOException{		
		initializeFiles();

		/*
		 * Simulate expression, genotype and cellcount %. The cellcount % are distributions around the pre-selected means. The expression
		 * is y = error + (B1 * cellcount) + (B2 * genotype)
		 */
		int genotypeCoefficientGroups = 3;
		int interactionCoefficientGroups = 3;
		int qtlCounter = 1;

		/*
		 * Main part of simulation. Do:
		 * 1. Simulate celltype percentages. Loop over number of samples
		 *   2. Loop over cellcount percentages (from command line)
		 *     3. Make NomralDistribution with mean;SD cellcount;cellcount/10
		 *     4. Sample from the normaldistribution; if sample < 0; sample = random(0,1)
		 * 5. Make a celltypeCoefficientDistribution with mean;SD 1;0.5
		 * 6. Loop over range 0 - <genotypeCoefficientGroups> for the genotypeCoefficient (e.g. -3, 0, 3)             -> j
		 * 	 7. Loop over range 0 - <interactionCoefficientGroups> for the interactionCoefficient (e.g. -5, 0, 5)     -> y
		 *     8. Make NormalDistribution biological noise, genotypeCoefficient, interactionCoefficient, 
		 * 				with mean;SD getNoise();getNoise()/4, j;0.5, y;0.5                                            <- The biological noise distribution is the distribution over the genes (QTL), later is used as mean for distribution over the samples per gene  
		 *     9. Loop over number of QTLs (from cmd)/(# genotype and interaction coefficients)      				  <- so if 100 QTLs and 3 genotype and interaction coefficients, loop for QTLs 33 times 
		 *       10. Make NormalDistribution celltypeCoeeficient with mean+SD 1+0.5
		 *       11. Loop over all the samples
		 *         12. From the previously made distributions, sample celltypeCoefficient, noise, genotypeCoefficient, interactionCoefficient
		 *         13. Make NormalDistributiotn with mean;SD noise;abs(noise/6) (over genes)                           <- Distribution is the samples variance per gene
		 *         14. SimulatedExpression = sample from step 10.
		 *         15. Loop over the number of celltypes
		 *           16. Calculate cellcountFactor <- (cellcount from step 4.)/(cellcount from command line)
		 *           17. SimulatedExpression += (celltypeCoefficient * cellcountFactor) + 
												(genotypeCoefficient * genotype) + 
												(interactionCoefficient * genotype * cellcountFactor);
		 */
		Map<String, HashMap<String, Double>> cellcountsPerSample = writeCellCountFile();
		NormalDistribution celltypeCoefficientDistribution = new NormalDistribution(10, 0.5);
		for(int j = 0; j < genotypeCoefficientGroups; j++){ 
			for (int y = 0; y < interactionCoefficientGroups; y++){
				// make a normal distribution around the given mean distribution so that each gene has a different (altho close to given) mean expression distribution
				// the SD of this distribution is large as genes can have a large difference in expression levels
				NormalDistribution noiseDistribution = new NormalDistribution(commandLineOptions.getNoise(), commandLineOptions.getNoise()/4);
				NormalDistribution genotypeCoeficientDistribution = new NormalDistribution(-3 + (j*3), 0.3);
				NormalDistribution interactioncCefficientDistribution = new NormalDistribution(-3+ (y*3), 0.3);
				// Make i number of QTLs
				for (int i = 0; i < commandLineOptions.getNumberOfQtls()/(genotypeCoefficientGroups*interactionCoefficientGroups); i++){
					// write the QTL name in first column
					expressionWriter.write("QTL_"+Integer.toString(qtlCounter));
					genotypeWriter.write("QTL_"+Integer.toString(qtlCounter));
					// loop over the number of samples and calculate expression level
					double averageCelltypeCoefficient = 0;
					double averageGenotypeCoefficient = 0;
					infoWriter.write("QTL_"+Integer.toString(qtlCounter));
					for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
						double celltypeCoefficient = celltypeCoefficientDistribution.sample();
						averageCelltypeCoefficient += celltypeCoefficient/commandLineOptions.getSampleSize();
						double genotypeCoefficient = genotypeCoeficientDistribution.sample();
						averageGenotypeCoefficient += genotypeCoefficient/commandLineOptions.getSampleSize();
						double interactionCoefficient = interactioncCefficientDistribution.sample();
						// Get the expression level for current QTL and current Sample. SD is low as between samples should not be high difference of "base" expression level 
						// (variance should be mostly in the genotype and interaction coefficient)
						double noiseSample = noiseDistribution.sample();
						NormalDistribution noise = new NormalDistribution(noiseSample, Math.sqrt(Math.pow(noiseSample/6.0,2)));
						double biologicalNoise = noise.sample();
						double simulatedExpression = biologicalNoise;
						// genotype randomly 0, 1 or 2
						int genotype = new Random().ints(1, 0, 3).findFirst().getAsInt();
						for (int c = 0; c <  commandLineOptions.getCellcountPercentages().length; c++){
							String celltypeName = "celltype_"+Integer.toString(c);
							String sampleName = "sample"+Integer.toString(s);
							if(c > 0){
							// only have an interaction effect on the first celltype, to simplify matters
								interactionCoefficient = new NormalDistribution(0, 0.3).sample();
							}
							if(s == 0){
								infoWriter.write("\t"+Double.toString(interactionCoefficient));
								expressionPerCelltypeWriterMap.get(celltypeName).write("QTL_"+Integer.toString(qtlCounter));
							}
							
							double cellcountFactor = new NormalDistribution(scale(cellcountsPerSample.get(celltypeName).get(sampleName),
									cellcountsPerSample.get(celltypeName).get("minCellcount"),	
									cellcountsPerSample.get(celltypeName).get("maxCellcount"), 0, 4), 1).sample();
						/*	
						System.out.printf("Simulated expression: %f\ncelltypeCoefficient: %f\ngenotypeCoefficient:%f\ncellcountFactor: %f\ngenotype: %d\ninteractionCoefficient: %f\n"+
										  "genotype: %d\ncellcountFactor: %f\nsimulatedExpression: %f\n", simulatedExpression,celltypeCoefficient, genotypeCoefficient, cellcountFactor, 
										  genotype, interactionCoefficient, genotype, cellcountFactor, simulatedExpression);
						System.out.printf("(coefficientCelltype * cellcountFactor): %f\n(coeficientGenotype * genotype): %f\n(coefficientInteraction * genotype * cellcountFactor): %f\n", 
								(celltypeCoefficient * cellcountFactor),(genotypeCoefficient * genotype), (interactionCoefficient * genotype * cellcountFactor), 
								  genotype, interactionCoefficient, genotype, cellcountFactor, simulatedExpression);
						System.out.println("-------------------------------");
						 */
							//simulatedExpression += interactionCoefficient * genotype * cellcountFactor;
							// write the simulated expression per celltype
							//expressionPerCelltypeWriterMap.get(celltypeName).write("\t"+Double.toString(biologicalNoise+(celltypeCoefficient * cellcountFactor) + 
							//																												(genotypeCoefficient * genotype) + 
							//																												(interactionCoefficient * genotype * cellcountFactor)));
							double celltypeSpecificExpression = (celltypeCoefficient) +// * cellcountFactor) + 
																(genotypeCoefficient * genotype) + 
																(interactionCoefficient * genotype * cellcountFactor);
							if(celltypeSpecificExpression < 0){
								celltypeSpecificExpression = 0;
							}
							simulatedExpression += celltypeSpecificExpression;
							expressionPerCelltypeWriterMap.get(celltypeName).write("\t"+Double.toString(celltypeSpecificExpression));
						}
						// expression evel cant be lower than 0, if negative give it random number between 0 and 1 (reflects real life higher level of 0-1 genes as well
						if (simulatedExpression < 0){
							simulatedExpression = 0 + Math.random();
						}
						/*
						System.out.printf("Simulated expression: %f\n", simulatedExpression);
						System.out.println("----------------------------------");
						System.exit(0);*/
						//simulatedExpression = Math.pow(2, simulatedExpression);
						expressionWriter.write("\t"+Double.toString(simulatedExpression));
						genotypeWriter.write("\t"+Integer.toString(genotype));
					}
					infoWriter.write("\t"+Double.toString(averageCelltypeCoefficient) +
							"\t"+Double.toString(averageGenotypeCoefficient));
					infoWriter.newLine();
					qtlCounter++;
					expressionWriter.newLine();
					genotypeWriter.newLine();
					for(BufferedWriter w : expressionPerCelltypeWriterMap.values()){
						w.newLine();
					}
				}
			}
		}
		expressionWriter.close();
		genotypeWriter.close();
		infoWriter.close();
		System.out.printf("Outfiles written to: %s",commandLineOptions.getOutfolder());
	}
	public static Map<String, HashMap<String, Double>> writeCellCountFile() throws IOException{
		File simulatedCellcountFile = new File(commandLineOptions.getOutfolder()+"/simulatedCellcounts.csv");
		FileOutputStream cellcountStream = new FileOutputStream(simulatedCellcountFile);
		BufferedWriter cellcountWriter = new BufferedWriter(new OutputStreamWriter(cellcountStream));
		// write headers with celltypes
		for (int c = 0; c < commandLineOptions.getCellcountPercentages().length; c++){
			// only write header once, before first sample
			cellcountWriter.write("\tcelltype_"+Integer.toString(c));
		}
		cellcountWriter.newLine();
		
		Map<String, HashMap<String, Double>> cellcountsPerSample = new HashMap<String, HashMap<String, Double>>();
		for (int s = 0; s < commandLineOptions.getSampleSize(); s++){
			cellcountWriter.write("sample_"+Integer.toString(s));
			int c = 0;
			for (double cellcount : commandLineOptions.getCellcountPercentages()){
				// this stays the same for every sample, calculate the celltype % per sample now
				// TODO: normal distribution with cellcount as mean and SD as cellcount/10, should probably be changed to more appropriate
				NormalDistribution normalDistributionCellcount =  new NormalDistribution(cellcount, cellcount/10);
				double cellcountOfSample =  normalDistributionCellcount.sample();
				// cellcount can not be < 0, if < 0 make it 0 + small random number between 0.01 and 0.91
				if (cellcountOfSample < 0){
					cellcountOfSample = 0 + Math.random()+0.01;
				}
				cellcountWriter.write("\t");
				cellcountWriter.write(Double.toString(cellcountOfSample));
				String celltypeName = "celltype_"+Integer.toString(c);
				String sampleName = "sample"+Integer.toString(s);
				// check if key value already exists. If not, make the new key entry
				if(cellcountsPerSample.containsKey(celltypeName)){
					cellcountsPerSample.get(celltypeName).put(sampleName, cellcountOfSample);
					// find the min and max cellcount values of cellcounts to use for scaling during simulation
					if (cellcountOfSample > cellcountsPerSample.get(celltypeName).get("maxCellcount")){
						cellcountsPerSample.get(celltypeName).put("maxCellcount", cellcountOfSample);
					}
					if  (cellcountOfSample < cellcountsPerSample.get(celltypeName).get("minCellcount")){
						cellcountsPerSample.get(celltypeName).put("minCellcount", cellcountOfSample);
					}
				}
				else{
					HashMap<String, Double> sampleCellCount = new HashMap<String, Double>();
					sampleCellCount.put(sampleName, cellcountOfSample);
					cellcountsPerSample.put(celltypeName, sampleCellCount);
					// initialize min and max, every iteration check if there is a new min/max
					cellcountsPerSample.get(celltypeName).put("minCellcount", cellcountOfSample);
					cellcountsPerSample.get(celltypeName).put("maxCellcount", cellcountOfSample);
				}
				c++;
			}
			cellcountWriter.newLine();
		}
		cellcountWriter.close();
		return(cellcountsPerSample);
	}
	public static final double scale(double value, double min, double max, double limitMin, double limitMax){
		/*
		 *              (b-a)(x - min)
		 * 		f(x) = --------------  + a
         *                max - min
         *                
         *  where min, max comes from cellcount, and a (limitMin), b (limitMax) is what I want to scale it to
		 */
		/*
		System.out.printf("limitMax: %s\tlimitMin: %s\tmin: %s\t max: %s\t value: %s\t result: %s\n", 
				limitMax, limitMin, min, max, value, ((limitMax - limitMin) * (value - min) / (max -min)) + limitMin);
			*/
        return ((limitMax - limitMin) * (value - min) / (max -min)) + limitMin;
	}
}
