#ifndef CSVWRITER_H_
#define CSVWRITER_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "rapidcsv.h"

//#undef PY_MAJOR_VERSION
//#define PY_MAJOR_VERSION 2
#ifdef PLOT
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

class CSVWriter {
public:
	CSVWriter(std::string _header, std::vector <double>& _succs, std::vector <double>& _times, int _start, int _end, int _step, std::string _successFile = "SuccessesTab.csv", std::string _timingFile = "TimingTab.csv", bool _plot = false) {
		bool newSuccFile = false;
		bool newTimeFile = false;

    	if(!fileExists(_successFile)) {
        	std::ofstream ofs(_successFile.c_str());
        	ofs.close();
        	newSuccFile = true;
    	}

    	if(!fileExists(_timingFile)) {
        	std::ofstream ofs(_timingFile.c_str());
        	ofs.close();
        	newTimeFile = true;
    	}

    	rapidcsv::Document docSucc(_successFile, rapidcsv::LabelParams(-1, -1));
    	rapidcsv::Document docTime(_timingFile,  rapidcsv::LabelParams(-1, -1));

    	int colSuccSize = 0;
    	if(newSuccFile)
        	colSuccSize = initTable(&docSucc, _start, _end, _step);
    	else
        	colSuccSize = docSucc.GetColumn<std::string>(0).size()-1;

		int colTimeSize = 0;
    	if(newTimeFile)
        	colTimeSize = initTable(&docTime, _start, _end, _step);
    	else
        	colTimeSize = docTime.GetColumn<std::string>(0).size()-1;


    	writeColumn(&docSucc, _header, _succs, colSuccSize);
    	writeColumn(&docTime, _header, _times, colTimeSize);

    	docSucc.Save();
    	docTime.Save();

#ifdef PLOT
		if(_plot) {
			plotTable(docSucc, _successFile);
			plotTable(docTime, _timingFile, true);
		}
#endif
	}	

#ifdef PLOT
	void plotTable(rapidcsv::Document &_doc, std::string &_file, bool _log = false, bool _show = false) {
		std::vector <std::string> colors = std::vector <std::string> ();

		colors.push_back("black");
		colors.push_back("blue");
		colors.push_back("red");
		colors.push_back("orange");
		colors.push_back("green");
		colors.push_back("gold");
		colors.push_back("lime");
		colors.push_back("magenta");
		colors.push_back("yellow");
		colors.push_back("aqua");
		colors.push_back("purple");
		colors.push_back("olive");
		colors.push_back("saddlebrown");
		colors.push_back("tomato");
		colors.push_back("silver");
		colors.push_back("navy");
		colors.push_back("greenyellow");
		colors.push_back("pink");
		colors.push_back("deepskyblue");

		std::vector <std::string> xCol = _doc.GetColumn<std::string>(0);
		std::vector <int> x = std::vector <int>();
		for(unsigned i = 1; i < xCol.size(); i++) {
			x.push_back(atoi(xCol[i].c_str()));
		}

		plt::figure();

		plt::figure_size(1200, 780);
		plt::title(_file);

		std::vector <std::string> headRow = _doc.GetRow<std::string>(0);
        for(unsigned i = 1; i < headRow.size(); i++) {
            xCol = _doc.GetColumn<std::string>(i);
            std::vector <int> data = std::vector <int>();
            
            for(int j = 1; j < xCol.size(); j++) {
                data.push_back(atoi(xCol[j].c_str()));
            }
            
			if(!_log)
	            plt::plot(x, data, {{"label", headRow[i]}});
			else
				plt::semilogy(x, data, colors[(i-1) % colors.size()], {{"label", headRow[i]}});
        }
		

		/*std::vector <std::string> headRow = _doc.GetRow<std::string>(0);
		for(unsigned i = 1; i < headRow.size(); i++) {
			xCol = _doc.GetColumn<std::string>(i);
			std::vector <int> data = std::vector <int>();

			for(int j = 1; j < xCol.size(); j++) {
				data.push_back(atoi(xCol[j].c_str()));
			}

			plt::named_plot(headRow[i], x, data);
		}*/
		plt::legend();
		std::string filename = (_file.substr(0, _file.find_last_of('.')));
		std::string filenamepdf = filename + ".pdf";
		plt::savefig(filenamepdf);
		std::string filenamepng = filename + ".png";
		plt::savefig(filenamepng);
		if(_show) {
			plt::show();
		}
	}
#endif

	bool fileExists(const std::string& name) {
  		struct stat buffer;   
  		return (stat (name.c_str(), &buffer) == 0);
	}

	unsigned initTable(rapidcsv::Document *_doc, int _start, int _end, int _step, std::string _counter = "#Jobs") {
    	int j = 0;

    	_doc->SetCell<std::string>(0, 0, _counter); 
    	for(int i = _start; i <= _end; i+=_step, j++) {
        	std::stringstream ss;
        	ss << i;
        	_doc->SetCell<std::string>(0, j+1, ss.str());
    	}

    	return j;
	}


	template<typename T >
	void writeColumn(rapidcsv::Document *_doc, std::string _name, std::vector <T> &_data, unsigned _size) {
    	std::vector <std::string> head = _doc->GetRow<std::string>(0);
    	int rowSize = head.size();
    	int i = 0;
    	for(; i < rowSize; i++) {
        	if(_name == head[i])
            	break;
    	}
    	rowSize = i;

    	_doc->SetCell<std::string>(rowSize, 0, _name);
    	for(int i = 0; i < _size; i++) {
        	if(i < _data.size())
            	_doc->SetCell<T>(rowSize, i+1, _data[i]);
        	else
            	_doc->SetCell<T>(rowSize, i+1, -1);
    	}
	}

	template<typename T >
	void writeStringColumn(rapidcsv::Document *_doc, std::string _name, std::vector <std::string> &_data, unsigned _size) {
    	std::vector <std::string> head = _doc->GetRow<std::string>(0);
    	int rowSize = head.size();
    	int i = 0;
    	for(; i < rowSize; i++) {
        	if(_name == head[i])
            	break;
    	}
    	rowSize = i;

    	_doc->SetCell<std::string>(rowSize, 0, _name);
    	for(int i = 0; i < _size; i++) {
        	if(i < _data.size())
            	_doc->SetCell<T>(rowSize, i+1, _data[i]);
        	else
            	_doc->SetCell<T>(rowSize, i+1, "--");
    	}
	}

	void writeMachineLoad(std::string _header, std::vector <std::vector <std::string > >* _machineLoad, std::string _filename = "MachineLoad.csv") {
		bool newFile = false;
    	if(!fileExists(_filename)) {
        	std::ofstream ofs(_filename.c_str());
        	ofs.close();
        	newFile = true;
    	}

    	rapidcsv::Document docML(_filename, rapidcsv::LabelParams(-1, -1));

    	int colSize = 0;
    	if(newFile)
        	colSize = initTable(&docML, 1, _machineLoad->at(0).size(), 1, "machine");
    	else
        	colSize = docML.GetColumn<std::string>(0).size()-1;

		for(int i = 0; i < _machineLoad->size(); i++) {
			std::stringstream header;
		    header << _header;
		    header << "-";
		    header << (i+1);
			
			writeStringColumn<std::string>(&docML, header.str(), _machineLoad->at(i), colSize);
		}
		docML.Save();
	}


};

#endif
