/******************************************************************************

                    C++ Code written by Asad Usmani on 25/11/2021
		      A novel DNA encoding and decoding scheme

*******************************************************************************/

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;
int printFlag = 1;

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


string convertIdxtoNTStr(int inIdx) {
   
    string outNt = "";
    int i = inIdx, j = 0, k = 4;
    while (k > 0)
    {
        int p = pow(4,k-1);
        j = i / p;
        i = i % p;
        switch (j) {
            case 0 : outNt += 'A'; break;
            case 1 : outNt += 'T'; break;
            case 2 : outNt += 'C'; break;
            case 3 : outNt += 'G'; break;
            default: return outNt;
        }
        k = k - 1;
    }
    return outNt;
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

string convertIntintoNT(int val) {
    
    string charNt;
        switch (val) {
            case 0 : charNt = "A"; break;
            case 1 : charNt = "T"; break;
            case 2 : charNt = "C"; break;
            case 3 : charNt = "G"; break;
            default: return "";
        }
    return charNt;
}


int ** calcVarianceGC(string encodeOutStr, int sublength)
{
    int prevSubSeq = 0;
    int nextSubSeq = 0;
    int countSubSeq = 0;
    int **matrixSubSeq = new int*[256];
    for (int k = 0 ; k < 256; k++)
    {
        matrixSubSeq[k] = new int[256];
    }
    
    for (int j = 0 ; j < 256; j++)
    {
        for (int k = 0 ; k < 256; k++)
        {
            matrixSubSeq[j][k] = 0;
        }
    }
    int count= 0;
    for (int k = 0; k < sublength; k ++)
    {
        for (int j = k; j < encodeOutStr.length()-sublength+1; j = j+sublength)
        {
            int countNT[4] = {  convertNTtoInt(encodeOutStr[j]), 
                                convertNTtoInt(encodeOutStr[j+1]), 
                                convertNTtoInt(encodeOutStr[j+2]), 
                                convertNTtoInt(encodeOutStr[j+3])
            };
            
            nextSubSeq = countNT[0]*64 + countNT[1]*16 + countNT[2]*4 + countNT[3];
            matrixSubSeq[prevSubSeq][nextSubSeq] = matrixSubSeq[prevSubSeq][nextSubSeq] + 1; 
            prevSubSeq = nextSubSeq;
        }
    }
    return matrixSubSeq;
}


int * calcPercentGC(string encodeOutStr, int parts)
{
    int k = encodeOutStr.length()/parts;
    int *matrixParts = new int[parts];
    
    for (int j = 0; j < parts ; j++)
    {
        int countNT[4] = {0, 0, 0, 0};
        for (int i = j*k; i < (j*k)+k ; i++)
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
        matrixParts[j] = sumGC;
    }
    return matrixParts;
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

bool checkHomopolymers(string encodeDNAStr, int subSeqeunceLength)
{
    bool homopolymers = false;
    for (int i = 0 ; i < encodeDNAStr.length()-subSeqeunceLength ; i++)
    {
        if ((encodeDNAStr[i] == encodeDNAStr[i+1]))
        {
            homopolymers = true;
        }
    } 
    return false;
}

string inputBinaryStrGen(int strSize)
{
    string inputStr = "";
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
    return inputStr;
}

bool checkEncodeDecode(string inputDNAStr, string decodeOutStr)
{
    return inputDNAStr.compare(decodeOutStr);
}

void printMatrixGCContent(int * matrixGCcontent, int parts, int totalLen)
{
    cout << "\n ================= GC Content =====================";
    for (int j = 0; j < parts; j++)
    {
        cout << "\nPart." << j+1 << ".GC Content:" << (matrixGCcontent[j]*100)/(totalLen/parts) << "%";
    }
    cout << "\n ================= GC Content =====================";
}

void printMatrixGCVariance(int ** matrixSubSeq, int subSeqeunceLength)
{
    cout << "\n ================= GC Variation =====================";
    int count = 0, flag = 0;
    int martixLen = pow(4, subSeqeunceLength);
    int *matrixSubSeqCol = new int[martixLen];
    int *matrixSubSeqRow = new int[martixLen];
    
    for (int j = 0 ; j < martixLen; j++)
    {
        flag = 0;
        for (int k = 0 ; k < martixLen; k++)
        {
            if ( matrixSubSeq[j][k] > 0)   
            {
                flag = 1;
            }
        }
        if (flag) 
            matrixSubSeqRow[j] = 1;
        else
            matrixSubSeqRow[j] = 0;
    } 
    
    for (int j = 0 ; j < martixLen; j++)
    {
        flag = 0;
        for (int k = 0 ; k < martixLen; k++)
        {
            if ( matrixSubSeq[k][j] > 0)   
            {
                flag = 1;
            }
        }
        if (flag) 
            matrixSubSeqCol[j] = 1;
        else
            matrixSubSeqCol[j] = 0;
    }
    
    
    cout << "\n";
    for (int j = 0 ; j < martixLen; j++)
    {
        if (printFlag && matrixSubSeqRow[j] && (j%5 == 0)) 
        {
            cout << "," << "\"" <<convertIdxtoNTStr(j) << "\"";
        }
    }
    
    cout << "\n[";
    for (int j = 0 ; j < martixLen; j++)
    {
        if (printFlag && matrixSubSeqRow[j] && j%5 == 0) 
        {
            cout << "\n[";
            int flagTemp = 0;
            for (int k = 0 ; k < martixLen; k++)
            {
                if (k%5 == 0)
                {
                    if ( matrixSubSeq[j][k] > 0)   
                    {
                        if (printFlag) cout << matrixSubSeq[j][k] << ",";
                    }
                    else if (matrixSubSeqCol[k] || k==(martixLen-1))
                    {
                        if (printFlag && k<(martixLen-1)) 
                        {
                            cout << "0,";
                        }
                        else cout << "0";
                    }
                }
            }
            cout << "],";
        }
    }
    cout << "]";
    cout << "\n ================= GC Variation =====================";
}

void mainPrint(bool printFlag , int indNo, string outStr)
{
    switch(indNo)
    {
        case 0: if (printFlag) cout << "\nInput DNA String:" << outStr; 
        break;
        case 1: if (printFlag) cout << "\nEncoded String:" << outStr; 
        break;
        case 2: if (printFlag) cout << "\nDecoded String:" << outStr; 
        break;
        case 3: if (printFlag) cout << "\nEncoding-Decoding Correctness:" << outStr;
        break;
        case 4: if (printFlag) cout << "\nHomopolymers don't exit:" << outStr; 
        break;
        default: if (printFlag) cout << "\nProgram Completed:" << outStr; 
    }
}

int main()
{
    // input binary string length
    int strSize = 15000, contentGCParts=10, subSeqeunceLength=4;
    // enable or disable all the console printing
    bool printAll = false;
    
    //generate a random input binary string of given length
    string inputStr = inputBinaryStrGen(strSize);
    
    //generate DNA seqeunce from input binary string
    string inputDNAStr = convertIntoDNAStr(inputStr);
    
    //Encode the DNA sequence using our proposed method
    string encodeOutStr = encodeDNAStr(inputDNAStr);
    
    //Conversely, Decoding is performed below 
    string decodeOutStr = decodeDNAStr(encodeOutStr);

   
    int * matrixGCcontent = calcPercentGC(encodeOutStr, contentGCParts);
    int ** martixGCVariance = calcVarianceGC(encodeOutStr, subSeqeunceLength);
   
    // printing could be performed using the below lines of code, disable indivitual printing using fasle instead.
    mainPrint(printAll && true , 0 , inputDNAStr);
    mainPrint(printAll && true , 1 , encodeOutStr);
    mainPrint(printAll && true , 2 , decodeOutStr);
    mainPrint(printAll && true , 3 , (checkEncodeDecode(inputDNAStr, decodeOutStr) == 0 ? "OK.": "Not Ok."));
    mainPrint(printAll && true , 4 , (checkHomopolymers(encodeOutStr, subSeqeunceLength) == 0 ? "OK.": "Not Ok."));
    //printMatrixGCContent(matrixGCcontent, contentGCParts, encodeOutStr.length());
    printMatrixGCVariance(martixGCVariance, subSeqeunceLength);
    mainPrint(true && true , 6 , "Done!");
        
    return 0;
}
