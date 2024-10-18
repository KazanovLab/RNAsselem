/*
 *  service.h
 *  bkglibML
 *
 *  Created by Marat Kazanov on 2/3/12.
 *  Copyright 2012 BalalayKrokoG. All rights reserved.
 *
 */

#include <string>
#include <vector>

using namespace std;

vector<string> split(string s);
vector<string> splitd(string s, char delimiter);
int isNumeric(string s);
string Upper(string s);
double str2d(string s);
int str2i(string s);
unsigned long str2ul(string s);
string i2str(int i);
string d2str(double d);
string ul2str(unsigned long ul);
