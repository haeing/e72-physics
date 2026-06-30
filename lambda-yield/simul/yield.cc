#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

const int runnumber = 2447;


Long64_t GetTrigDScaler(int runnumber)
{
    string scaler_dir =
        "/gpfs/group/had/sks/E72/JPARC2025Nov/share/scaler";

    string scaler_file =
        Form("%s/scaler_%05d.txt",
             scaler_dir.c_str(), runnumber);

    ifstream fin(scaler_file);

    if (!fin.is_open()) {
        cerr << "Cannot open " << scaler_file << endl;
        return -1;
    }

    string line;

    while (getline(fin, line)) {

        stringstream ss(line);

        string name;
        Long64_t value;

        if (!(ss >> name))
            continue;

        if (name != "TRIG-D")
            continue;

        if (!(ss >> value))
            continue;
        return value;
    }

    cerr << "TRIG-D not found in " << scaler_file << endl;
    return -1;
}
void yield(){
  Long64_t nKbeam = GetTrigDScaler(runnumber);
  
}
