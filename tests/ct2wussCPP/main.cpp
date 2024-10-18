//
//  main.cpp
//  ct2wuss
//
//  Created by Marat Kazanov on 11/07/2023.
//  Copyright © 2023 Marat Kazanov. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "esl_wuss.h"
#include "service.h"

using namespace std;

int LoadCT(string path, int * &ct, vector<string> &nt, string &header)
{
    vector<string> flds;
    string line;
    
    ifstream f(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }

    int cnt = 0;
    getline(f, line);
    while(1)
    {
        getline(f, line);
        if(f.eof())
            break;
        cnt++;
    }

    f.clear();
    f.seekg(0, ios::beg);
    
    ct = new int[cnt+1];
    ct[0] = 0;
    
    getline(f, header);
    int i = 0;
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = split(line);
            ct[i+1] = str2i(flds[4]);
            nt.push_back(flds[1]);
            i++;
        }
    }
    
    return(cnt);
}


int main(int argc, const char * argv[]) {

    int *ct;
    int n;
    char *str;
    int res;
    string pathIn,pathOut1,pathOut2;
    vector<string> nt;
    
    pathIn = argv[1];
    pathOut1 = argv[2];
    pathOut2 = argv[3];
    
    string header;
    n = LoadCT(pathIn, ct, nt, header);
    
    str = new char[n+1];
    
    res = esl_ct2wuss(ct,n,str);
    //res = esl_ct2simplewuss(ct,n,str);
    
    ofstream f;
    f.open(pathOut1.c_str());
    f << header << '\n';
    string seq,ss;
    for(int i=0;i<n;i++)
    {
        f << setfill(' ') << setw(5) << i+1
          << setfill(' ') << setw(2) << nt[i]
          << setfill(' ') << setw(6) << i
          << setfill(' ') << setw(6) << i+2
          << setfill(' ') << setw(6) << ct[i+1]
          << setfill(' ') << setw(6) << i+1
          << setfill(' ') << setw(2) << str[i] <<'\n';
        seq += nt[i];
        ss += str[i];
    }
    f.close();
    
    vector<string> headers;
    string headerName;
    headers = split(header);
    headerName = headers[headers.size()-1];
    f.open(pathOut2.c_str());
    f << '>' << headerName << '\n' << seq << '\n' << ss << '\n';
    f.close();
    
    delete[] ct;
    delete[] str;

    return 0;
}
