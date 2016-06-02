package QtlSimulator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class CommandLineOptions {
	private double[] cellcountPercentages = {70, 20, 7, 3};
	private int numberOfQtls = 1000;
	private String outfolder;
	private int sampleSize = 100;
	private int noise = 10;
	public void parseCommandLine(String[] args) throws ParseException {
		/*
		 * Standard command line parsing.
		 * 
		 * @param args A string vector of all arguments given to the command
		 * line, e.g. for `java -jar Deconvolution.jar -o` args = ["-o"]
		 * 
		 * @return A CommandLine object that includes all the options given to
		 * the command line
		 */
		Options options = new Options();
		Option help = new Option("help", "print this message");
		Option cellcountPercentagesOption = Option.builder("c").required(false).hasArg().longOpt("cellcount_percentages")
				.desc("The cell type percentages around which a distribution of % will be used").build();
		Option noiseOption = Option.builder("e").required(false).hasArg().longOpt("error")
				.desc("The mean of poisson distribution from which error is drawn").build();
		Option numberOfQtlsOption = Option.builder("n").required(false).hasArg().longOpt("number_of_qtls_per_group")
				.desc("The number of QTLs to simulate per genotype group").build();
		Option outfolderOption = Option.builder("o").required(true).hasArg().longOpt("outfolder")
				.desc("Outfolder to write results to").build();
		Option sampleSizeOption = Option.builder("s").required(false).hasArg().longOpt("sample_size")
				.desc("Number of samples to use").build();
		
		options.addOption(noiseOption);
		options.addOption(help);
		options.addOption(cellcountPercentagesOption);
		options.addOption(numberOfQtlsOption);
		options.addOption(outfolderOption);
		options.addOption(sampleSizeOption);
		
		CommandLineParser cmdLineParser = new DefaultParser();
		CommandLine cmdLine = cmdLineParser.parse(options, args);
		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		if (cmdLine.hasOption("help")) {
			formatter.printHelp("deconvolution", options, true);
		}
		parseOptions (cmdLine);
		printArgumentValues(cmdLine);
	}
	
	private void parseOptions(CommandLine cmdLine){
		if(cmdLine.hasOption("cellcount_percentages")){
			String[] items = cmdLine.getOptionValue("cellcount_percentages").split(",");
			cellcountPercentages = new double[items.length];
			double total = 0;
			for (int i = 0; i < items.length; i++) {
				cellcountPercentages[i] = Double.parseDouble(items[i]);
				total += cellcountPercentages[i];
			}
			if(total != 100){
				throw new IllegalArgumentException("Cellcount % have to add up to 100%, added up to: "+Double.toString(total));
			}
		}
		if(cmdLine.hasOption("noise")){
			noise = Integer.parseInt(cmdLine.getOptionValue("noise"));
		}
		if(cmdLine.hasOption("number_of_qtls_per_group")){
			numberOfQtls = Integer.parseInt(cmdLine.getOptionValue("number_of_qtls_per_group"));
		}
		outfolder = cmdLine.getOptionValue("outfolder");
		if(cmdLine.hasOption("sample_size")){
			sampleSize = Integer.parseInt(cmdLine.getOptionValue("sample_size"));
		}
	}
	

	private void printArgumentValues(CommandLine cmdLine){
		System.out.println("======= QtlSimulator paramater settings =======");
		System.out.printf("Cellcount percentages (-c): ");
		for(double cellcount : cellcountPercentages){
			System.out.printf("%s, ", cellcount);
		}
		System.out.println();
		System.out.printf("Number of QTLs to write (-n): %s\n", numberOfQtls);
		System.out.printf("Outfolder: (-o): %s\n", outfolder);
		System.out.printf("Samplesize: (-s): %s\n", sampleSize);
		System.out.printf("Noise: (-e): %s\n", noise);
		System.out.println("=================================================");
	}
	public double[] getCellcountPercentages(){
		return (cellcountPercentages);
	}
	public int getNumberOfQtls(){
		return (numberOfQtls);
	}
	public String getOutfolder(){
		return (outfolder);
	}
	public int getSampleSize(){
		return (sampleSize);
	}
	public int getNoise(){
		return (noise);
	}
}
