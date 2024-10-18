/*
 *  service.cpp
 *  bkglibML
 *
 *  Created by Marat Kazanov on 2/3/12.
 *  Copyright 2012 BalalayKrokoG. All rights reserved.
 *
 */
//#include <stdafx.h>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "service.h"

vector<string> split(string s)
{
	vector<string> ret;
	istringstream iss(s);
	copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(ret));
	return(ret);
}

vector<string> splitd(string s, char delimiter)
{
	vector<string> ret;
	string item;
	stringstream ss(s);
	while(getline(ss,item,delimiter))
	{
		ret.push_back(item);
	}
	return(ret);
}

int isNumeric(string s)
{
	float f;
	istringstream ss(s);
	if(!(ss >> f))
		return(0);
	else 
		return(1);
}

string Upper(string s)
{
	//string ss;
	//ss = s;
	transform(s.begin(),s.end(),s.begin(),::toupper);
	return(s);
}

double str2d(string s)
{
	double ret;
	istringstream iss(s);
	iss >> ret;
	return(ret);
}

int str2i(string s)
{
    int ret;
    istringstream iss(s);
    iss >> ret;
    return(ret);
}

unsigned long str2ul(string s)
{
    unsigned long ret;
    istringstream iss(s);
    iss >> ret;
    return(ret);
}

string i2str(int i)
{
	ostringstream oss;
	oss << i;
	return(oss.str());
}

string d2str(double d)
{
	ostringstream oss;
	oss << d;
	return(oss.str());
}

string ul2str(unsigned long ul)
{
    ostringstream oss;
    oss << ul;
    return(oss.str());
}
