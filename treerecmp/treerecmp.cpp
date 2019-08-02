#include <iostream>
#include <sstream>
#include "Poco/Util/Application.h"
#include "Poco/Util/OptionProcessor.h"
#include "Poco/Util/HelpFormatter.h"

using namespace std;
using Poco::Util::Application;
using Poco::Util::Option;
using Poco::Util::OptionSet;
using Poco::Util::OptionCallback;
using Poco::Util::OptionProcessor;
using Poco::Util::HelpFormatter;


class TreeReCmpApp : public Application
{

public:
	TreeReCmpApp() : _helpRequested(false)
	{
	}

	void defineOptions(OptionSet& options)
	{
		setUnixOptions(true);
		string distInfo = "Allow to specify distances (from 1 up to 12) :\n"
			"Metrics for unrooted trees :\n"
			"    ms - the Matching Split metric,\n"
			"    rf - the Robinson - Foulds metric,\n"
			"    pd - the Path Difference metric,\n"
			"    qt - the Quartet metric,\n"
			"    um - the UMAST metric,\n"
			"Metrics for rooted trees :\n"
			"    mc - the Matching Cluster metric,\n"
			"    rc - the Robinson - Foulds metric based on clusters,\n"
			"    ns - the Nodal Splitted metric with L2 norm,\n"
			"    tt - the Triples metric,\n"
			"    mp - the Matching Pair metric,\n"
			"    mt - the RMAST metric,\n"
			"    co - the Cophenetic Metric with L2 norm.\n"
			"Example: -d ms rf";
		options.addOption(
			Option("help", "h", "Display this help information.")
			.required(false)
			.repeatable(false)
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleHelp)));
		options.addOption(
			Option("overlap", "s", "Overlapping pair comparison mode. Every two neighboring trees are compared.")
			.required(false)
			.repeatable(false)
			.group("mode"));
		options.addOption(
			Option("window", "w", "Window comparison mode. Every two trees within a	window are compared.")
			.required(false)
			.repeatable(false)
			.group("mode"));
		options.addOption(
			Option("matrix", "m", "Matrix comparison mode. Every two trees in the input file are compared.")
			.required(false)
			.repeatable(false)
			.group("mode"));
		options.addOption(
			Option("reference", "r", "Referential trees to all input trees mode. Each referential tree is compared to each tree in the input"
				"file.")
			.required(false)
			.repeatable(false)
			.argument(" refTreeFile")
			.group("mode"));
		options.addOption(
			Option("distance", "d", distInfo)
			.required(true)
			.repeatable(false)
			.argument(" distance"));
		options.addOption(
			Option("input", "i", "Specify an input trees file (required).")
			.required(true)
			.repeatable(false)
			.argument(" inputFile"));
		options.addOption(
			Option(" output", "o", "Specify an output file (required).")
			.required(true)
			.repeatable(false)
			.argument(" outputFile"));
		options.addOption(
			Option("NORMALIZE", "N", "Report normalized distances.")
			.required(false)
			.repeatable(false));
		options.addOption(
			Option("INCLUDE", "I", "Include summary section in the output file.")
			.required(false)
			.repeatable(false));
		options.addOption(
			Option("PRUNE", "P", "Prune compared trees if needed (trees can have different leaf sets).")
			.required(false)
			.repeatable(false));
		options.addOption(
			Option("ALIGNMENT", "A", "Generate alignment files(only for MS and MC metrics).Cannot be used with -O option.")
			.required(false)
			.repeatable(false));
		options.addOption(
			Option("OPTIMIZE", "O", "Use MS/MC/MP metrics optimized for similar trees. Cannot be used with -A option.")
			.required(false)
			.repeatable(false));
	}

	void handleHelp(const std::string& name, const std::string& value)
	{
		_helpRequested = true;
		stopOptionsProcessing();
	}

	void displayHelp()
	{
		HelpFormatter helpFormatter(options());
		helpFormatter.setUnixStyle(true);
		helpFormatter.setCommand(commandName());
		helpFormatter.setUsage(
			"-s | -w <size> | -m | -r <refTreeFile> -d <distance(s)> -i <inputFile> -o <outputFile> [-N][-P][-I][-A | -O]");
		helpFormatter.setHeader("Application for comparison and reconstruction of phylogenetic trees in polynomial time.");
		helpFormatter.format(std::cout);
	}

	int main(const std::vector<std::string>& args)
	{
		//if (_helpRequested || args.empty())
		//{
		//	displayHelp();
		//}
		//else
		{
			Option distances = options().getOption("d");
			stringstream msg;
			msg << "Metric:\t\t " << distances.argumentName() << endl;
			//<< "Comparison mode: " << modeName << endl
			//<< "Input file:\t " << inFile << endl
			//<< "Output file:\t " << outFile << endl
			//<< "------------" << endl;
			cout << msg.str();
		}
		return Application::EXIT_OK;
	}

private:
	bool _helpRequested;
};

POCO_APP_MAIN(TreeReCmpApp)