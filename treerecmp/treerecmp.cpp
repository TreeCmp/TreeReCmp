#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Poco/Util/Application.h"
#include "Poco/Util/OptionProcessor.h"
#include "Poco/Util/HelpFormatter.h"
#include "Poco/XML/XMLStreamParser.h"
#include "Poco/FileStream.h"

using namespace std;
using Poco::Util::Application;
using Poco::Util::Option;
using Poco::Util::OptionSet;
using Poco::Util::OptionCallback;
using Poco::Util::OptionProcessor;
using Poco::Util::HelpFormatter;
using namespace Poco::XML;

struct DistDesc {

	DistDesc() : _rooted(false) { }

public:
	string _name;
	string _command_name;
	string _method_name;
	string _description;
	string _unif_data;
	string _yule_data;
	string _aln_file_suffix;
	string _fullname;
	string _complexity;
	string _authors;
	string _publication;
	string _year;
	string _full_description;
	bool _rooted;
};

enum class ConfigElementType
{
	name,
	command_name,
	method_name
};

constexpr unsigned int switchHash(const char *s, int off = 0) {
	return !s[off] ? 5381 : (switchHash(s, off + 1) * 33) ^ s[off];
}

class TreeReCmpApp : public Application
{

public:
	TreeReCmpApp() : _helpRequested(false)
	{
		string configFilePath = "..\\config\\config.xml";

		Poco::FileInputStream istr(configFilePath);		

		XMLStreamParser confParser(istr, "config_parser", XMLStreamParser::RECEIVE_ELEMENTS | XMLStreamParser::RECEIVE_CHARACTERS);
		//WhitespaceFilter wsf(confParser);
		try {
			while (confParser.next() != XMLStreamParser::EV_EOF) {
				if (confParser.peek() == XMLStreamParser::EV_START_ELEMENT && (confParser.localName().compare(0, 6,"metric") == 0)) {
					DistDesc cet;
					string elementXML;
					while (!(confParser.next() == XMLStreamParser::EV_END_ELEMENT && (confParser.localName().compare("metric") == 0))) {
						if (confParser.peek() == XMLStreamParser::EV_CHARACTERS && !all_of(confParser.value().begin(), confParser.value().end(), isspace)) {
							elementXML = confParser.localName();
							switch (switchHash(elementXML.c_str())) {
							case switchHash("name"):
								cet._name = confParser.value();
								break;
							case switchHash("command_name"):
								cet._command_name = confParser.value();
								break;
							case switchHash("method_name"):
								cet._method_name = confParser.value();
								break;
							case switchHash("description"):
								cet._description = confParser.value();
								break;
							case switchHash("unif_data"):
								cet._unif_data = confParser.value();
								break;
							case switchHash("yule_data"):
								cet._yule_data = confParser.value();
								break;
							case switchHash("aln_file_suffix"):
								cet._aln_file_suffix = confParser.value();
								break;
							case switchHash("fullname"):
								cet._fullname = confParser.value();
								break;
							case switchHash("complexity"):
								cet._complexity = confParser.value();
								break;
							case switchHash("authors"):
								cet._authors = confParser.value();
								break;
							case switchHash("publication"):
								cet._publication = confParser.value();
								break;
							case switchHash("year"):
								cet._year = confParser.value();
								break;
							case switchHash("full_description"):
								cet._full_description = confParser.value();
								break;
							case switchHash("rooted"):
								cet._rooted = true;
								break;
							}
						}						
					}
					_distances.push_back(cet);
				}				
			}
			
		}
		catch (XMLStreamParserException e) {
			cout << "Error parsing configuration file: " << e.description() << endl;
			exit(Application::EXIT_CONFIG);
		}
	}

	void defineOptions(OptionSet& options)
	{
		setUnixOptions(true);
		string distInfo = "Allow to specify distances (from 1 up to "+ to_string(_distances.size()) + "):\n";
		distInfo += "Metrics for unrooted trees:\n";
		for (DistDesc dd : _distances) {
			if (!dd._rooted) {
				distInfo += " " + dd._command_name + " - " + dd._fullname + ",\n";
			}
		}
		distInfo += "Metrics for rooted trees:\n";
		for (DistDesc dd : _distances) {
			if (dd._rooted) {
				distInfo += " " + dd._command_name + " - " + dd._fullname + ",\n";
			}
		}
		distInfo += "Example: -d ms rf";
		options.addOption(
			Option("help", "h", "Display this help information.")
			.required(false)
			.repeatable(false)
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleHelp)));
		options.addOption(
			Option("overlap", "s", "Overlapping pair comparison mode. Every two neighboring trees are compared.")
			.required(false)
			.repeatable(false)
			.group("mode")
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleMode)));
		options.addOption(
			Option("window", "w", "Window comparison mode. Every two trees within a	window are compared.")
			.required(false)
			.repeatable(false)
			.group("mode")
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleMode)));
		options.addOption(
			Option("matrix", "m", "Matrix comparison mode. Every two trees in the input file are compared.")
			.required(false)
			.repeatable(false)
			.group("mode")
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleMode)));
		options.addOption(
			Option("reference", "r", "Referential trees to all input trees mode. Each referential tree is compared to each tree in the input"
				"file.")
			.required(false)
			.repeatable(false)
			.argument(" refTreeFile")
			.group("mode")
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleMode)));
		options.addOption(
			Option("distance", "d", distInfo)
			.required(true)
			.repeatable(true)
			.argument(" distance")
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleDistance)));
		options.addOption(
			Option("input", "i", "Specify an input trees file (required).")
			.required(true)
			.repeatable(false)
			.argument(" inputFile")
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleInputFile)));
		options.addOption(
			Option(" output", "o", "Specify an output file (required).")
			.required(true)
			.repeatable(false)
			.argument(" outputFile")
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleOnputFile)));
		options.addOption(
			Option("NORMALIZE", "N", "Report normalized distances.")
			.required(false)
			.repeatable(false)
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleAdditionalOptions)));
		options.addOption(
			Option("INCLUDE", "I", "Include summary section in the output file.")
			.required(false)
			.repeatable(false)
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleAdditionalOptions)));
		options.addOption(
			Option("PRUNE", "P", "Prune compared trees if needed (trees can have different leaf sets).")
			.required(false)
			.repeatable(false)
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleAdditionalOptions)));
		options.addOption(
			Option("ALIGNMENT", "A", "Generate alignment files(only for MS and MC metrics).Cannot be used with -O option.")
			.required(false)
			.repeatable(false)
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleAdditionalOptions)));
		options.addOption(
			Option("OPTIMIZE", "O", "Use MS/MC/MP metrics optimized for similar trees. Cannot be used with -A option.")
			.required(false)
			.repeatable(false)
			.callback(OptionCallback<TreeReCmpApp>(this, &TreeReCmpApp::handleAdditionalOptions)));
	}

	void handleHelp(const std::string& name, const std::string& value)
	{
		_helpRequested = true;
		stopOptionsProcessing();
	}

	void handleDistance(const std::string& name, const std::string& value)
	{
		_chosenDistances.push_back(value);
	}

	void handleMode(const std::string& name, const std::string& value)
	{
		_modeName = name;
	}

	void handleInputFile(const std::string& name, const std::string& value)
	{
		_inFile = value;
	}

	void handleOnputFile(const std::string& name, const std::string& value)
	{
		_outFile = value;
	}

	void handleAdditionalOptions(const std::string& name, const std::string& value)
	{
		_additionalOptions.push_back(name);
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
		if (_helpRequested /*|| args.empty()*/)
		{
			displayHelp();
		}
		else
		{
			string message = "Active options:\n";
			message += "Type of the analysis:\t" + _modeName + "\n";
			message += "Metric:\t\t\t";
			for (string distance : _chosenDistances) {
				message += (distance + " ");
			}
			message += "\n";
			message += "Input file:\t\t" + _inFile + "\n";
			message += "Output file:\t\t" + _outFile + "\n";
			message += "Additional options:\t";
			for (string addOpt : _additionalOptions) {
				message += (addOpt + " ");
			}
			cout << message;

			//Reading trees
			/*cout << "Scanning input file... " << flush;
			vector<TreeTemplate<Node> *> treesIn;
			vector<TreeTemplate<Node> *> trees;
			Newick newickReader(false);
			try {
				newickReader.read(inFile, (vector<Tree *>&) (treesIn));
			}
			catch (exception& e) {
				cout << "Error when reading trees. Application terminated.\n";
				return 0;
			}
			cout << treesIn.size() << " trees found." << endl;*/
			for (string distance : _chosenDistances) {
				
			}
		}
		return Application::EXIT_OK;
	}

private:
	bool _helpRequested;
	vector<DistDesc> _distances;
	vector<string> _chosenDistances;
	vector<string> _additionalOptions;
	string _modeName;
	string _inFile;
	string _outFile;
};

POCO_APP_MAIN(TreeReCmpApp)