/******************************************************************************

                    C++ Code written by Asad Usmani on 25/11/2021
					  A novel DNA encoding and decoding scheme

*******************************************************************************/

#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;


char convertDNATable[4][4] = {  {'A', 'T', 'C', 'G'},
                                {'T', 'C', 'G', 'A'},
                                {'C', 'G', 'A', 'T'},
                                {'G', 'A', 'T', 'C'}
                            };
                            
char rightTable[4][4]  =    {   {'T', 'C', 'T', 'C'},
                                {'C', 'G', 'G', 'A'},
                                {'T', 'G', 'T', 'A'},
                                {'C', 'A', 'A', 'T'}
                            };
                            
char leftTable[4][4]    =   {   {'C', 'G', 'G', 'T'},
                                {'G', 'A', 'A', 'C'},
                                {'G', 'A', 'A', 'T'},
                                {'T', 'C', 'T', 'A'}
                            };
    

int convertNTtoInt(char inNt) {
    
    int outIdx;
    switch (inNt) {
        case 'A' : outIdx = 0; break;
        case 'T' : outIdx = 1; break;
        case 'C' : outIdx = 2; break;
        case 'G' : outIdx = 3; break;
        default: return -1;
    }
    return outIdx;
}

char convertInttoNT(int inIdx) {
   
    char outNt;
    switch (inIdx) {
        case 0 : outNt = 'A'; break;
        case 1 : outNt = 'T'; break;
        case 2 : outNt = 'C'; break;
        case 3 : outNt = 'G'; break;
        default: return outNt;
    }
    return outNt;
}

string convertIntoDNAStr(string str) {
    
    string outDNAStr = "";
    for (int i =0; i < str.length() ;  i++){
        
        switch (str[i]) {
            case '0' : outDNAStr += "A"; break;
            case '1' : outDNAStr += "T"; break;
            case '2' : outDNAStr += "C"; break;
            case '3' : outDNAStr += "G"; break;
            default: return "";
        }
    }
    return outDNAStr;
}

string decodeDNAStr(string inputStr) {
    
    bool evenTurn = true;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, mChar, rChar, eChar;
    string outStr = "", substrTemp = "";
    
    for (int i=0, m=3; m <= inputStr.length();  i+=4, m+=4) {
        
        j = convertNTtoInt(inputStr[m-1]);
        k = convertNTtoInt(inputStr[m+1]);

        substrTemp = inputStr.substr(i, 4);
        
        lChar = substrTemp[0];
        mChar = substrTemp[1];
        rChar = substrTemp[2];
        eChar = substrTemp[3];
        
        g = convertNTtoInt(lChar); 
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {
            
            if (eChar == rightTable[j][k]) {
                substrTemp += inputStr.substr(i, 3);
            }
            else if (leftTable[g][h] == mChar){
                substrTemp += lChar; 
                substrTemp += lChar;
                substrTemp += rChar;
            }
             else{
                substrTemp += lChar; 
                substrTemp += rChar;
                substrTemp += rChar; 
            }
            evenTurn = false;
        }
        else {
            
            if (eChar == leftTable[j][k]) {
                substrTemp += inputStr.substr(i, 3);
            }
            else if (leftTable[g][h] == mChar){
                substrTemp += lChar; 
                substrTemp += lChar;
                substrTemp += rChar;
            }
             else{
                substrTemp += lChar; 
                substrTemp += rChar;
                substrTemp += rChar; 
            }
            evenTurn = true;
        }
        outStr += substrTemp;
        substrTemp = "";
    }
    return outStr;
}

string encodeDNAStr(string inputStr) {
    
    bool evenTurn = true;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, rChar, mChar;
    string outStr = "", substrTemp = "";
    
    for (int i=0, m=2; m <= inputStr.length();  i+=3, m+=3) {
        
        j = convertNTtoInt(inputStr[m]);
        k = convertNTtoInt(inputStr[m+1]);

        substrTemp = inputStr.substr(i, 3);
        
        lChar = substrTemp[0];
        mChar = substrTemp[1];
        rChar = substrTemp[2];
        
        g = convertNTtoInt(lChar); 
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {
            if (lChar == mChar) {
                substrTemp += lChar; 
                substrTemp += leftTable[g][h];
                substrTemp += rChar;
                substrTemp += leftTable [j][k];
            }
            else if (mChar == rChar){
                substrTemp += lChar; 
                substrTemp += rightTable[g][h];
                substrTemp += rChar;
                substrTemp += leftTable [j][k];
            }
             else{
                substrTemp += inputStr.substr(i, 3);
                substrTemp += rightTable[j][k];   
            }
            evenTurn = false;
        }
        else {
            if (lChar == mChar) {
                substrTemp += lChar; 
                substrTemp += leftTable[g][h];
                substrTemp += rChar;
                substrTemp += rightTable [j][k];
            }
            else if (mChar == rChar) {
                substrTemp += lChar; 
                substrTemp += rightTable[g][h];
                substrTemp += rChar;
                substrTemp += rightTable [j][k];
            }
            else {
                substrTemp += inputStr.substr(i, 3);
                substrTemp += leftTable[j][k]; 
            }
            evenTurn = true;
        }
        outStr += substrTemp;
        substrTemp = "";
    }
    return outStr;
}


int main()
{
    int strSize = 60000;
    string inputStr = "";
    int countNT[4] = {0, 0, 0, 0};
    
    srand (time(NULL));
    for (int i =0 ; i < strSize ; i++)
    {
        int rval = rand () % 4; 
        (
            rval == 3 ? inputStr += '3': 
            rval == 2 ? inputStr += '2': 
			rval == 1 ? inputStr += '1': 
			inputStr += '0'
		);
    }
    
    string inputDNAStr = convertIntoDNAStr(inputStr);
    //cout << "\nBefore Encoding DNA String:" << inputDNAStr;
    //cout << "  and length:" << strSize;
    
    string encodeOutStr = encodeDNAStr(inputDNAStr);
    //cout << "\nAfter  Encoding DNA String:" << encodeOutStr;
    
    string decodeOutStr = decodeDNAStr(encodeOutStr);
    //cout << "\nAfter  Decoding DNA String:" << decodeOutStr;
    
    bool homopolymers = false;
    for (int i = 0 ; i < encodeOutStr.length()-4 ; i++)
    {
        if ((encodeOutStr[i] == encodeOutStr[i+1]))
        {
            homopolymers = true;
            cout << "\n Index of duplication!" << i ;
        }
    }
    
    for (int i =0 ; i < encodeOutStr.length() ; i++)
    {
        (
        encodeOutStr[i] == 'A' ? countNT[0] += 1: 
        encodeOutStr[i] == 'T' ? countNT[1] += 1: 
        encodeOutStr[i] == 'C' ? countNT[2] += 1: 
	    countNT[3] += 1
	    );
    }
    
    int sumAT = countNT[0] + countNT[1]; 
    int sumGC = countNT[2] + countNT[3];
   
    if ((inputDNAStr.compare(decodeOutStr)) == 0)
        cout << "\n\n1. Encoding-Decoding Correctness!OK";
    else
        cout << "\n\n1. Something is wrong with Encoding or Decoding!NK";
    
    if (!homopolymers)
        cout << "\n2. Homopolymers                 !OK";
    else
        cout << "\n2. Homopolymers !NK";
        
    if (sumAT > sumGC)
        cout << "\n3. Highly Balanced-GC Contents  !OK";
    else
        cout << "\n3.Less Balanced-GC Contents!NK";
    
    cout << "\n4. Total(A+T+C+G) = " << sumAT + sumGC 
    << ";  A+T=" << sumAT << ", " << (sumAT*100/(sumAT + sumGC)) << "%"
    << ";  C+G=" << sumGC << ", " << (sumGC*100/(sumAT + sumGC)) << "%"
    << ";  Diff=" << (sumAT*100/(sumAT + sumGC)) - (sumGC*100/(sumAT + sumGC)) << "%";
    
    return 0;
}
