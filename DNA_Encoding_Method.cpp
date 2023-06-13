/******************************************************************************

                    C++ Code written by Asad Usmani on 09/06/2023
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
         
char rightTable[4][4]   =   {  {'C', 'C', 'T', 'C'},
                                {'C', 'C', 'G', 'A'},
                                {'T', 'G', 'T', 'A'},
                                {'C', 'A', 'A', 'T'}
                            };
                            
char leftTable[4][4]    =   {  {'G', 'G', 'G', 'T'},
                                {'G', 'G', 'A', 'C'},
                                {'G', 'A', 'A', 'T'},
                                {'T', 'C', 'T', 'A'}
                            };
    
                       
char rightTable1[4][4]  =   {   {'T', 'C', 'T', 'C'},
                                {'C', 'G', 'G', 'A'},
                                {'T', 'G', 'T', 'A'},
                                {'C', 'A', 'A', 'T'}
                            };
                            
char leftTable1[4][4]   =   {   {'C', 'G', 'G', 'T'},
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


string decodeDNAStrWithOption(string inputStr, int nHomopoly) {
    
    bool evenTurn = true;
    int strideCount = 2 * nHomopoly + 1;
    int midIndex = strideCount/2; 
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, mChar, rChar, eChar;
    string outStr = "", substrTemp = "";
    
    for (int i=0, m=strideCount; m <= inputStr.length();  
                i+=strideCount+1, m+=strideCount+1) {
        
        j = convertNTtoInt(inputStr[m-1]);
        k = convertNTtoInt(inputStr[m+1]);

        substrTemp = inputStr.substr(i, strideCount+1);
        
        lChar = substrTemp[midIndex-1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex+1];
        eChar = substrTemp[strideCount];
        
        g = convertNTtoInt(lChar); 
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {
            
            if (eChar == rightTable[j][k]) {
                substrTemp += inputStr.substr(i, strideCount);
            }
            else if (leftTable[g][h] == mChar){
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
             else{
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            evenTurn = false;
        }
        else {
            
            if (eChar == leftTable[j][k]) {
                substrTemp += inputStr.substr(i, strideCount);
            }
            else if (leftTable[g][h] == mChar){
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
             else{
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            evenTurn = true;
        }
        outStr += substrTemp;
        substrTemp = "";
    }
    return outStr;
}

string encodeDNAStrWithOption(string inputStr, int nHomopoly) {
    
    bool evenTurn = true;
    int strideCount = 2 * nHomopoly + 1;
    int midIndex = strideCount/2; 
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, rChar, mChar;
    string outStr = "", substrTemp = "";
    
    for (int i=0, m=strideCount-1; m <= inputStr.length();  
                        i+=strideCount, m+=strideCount) {
        
        j = convertNTtoInt(inputStr[m]);
        k = convertNTtoInt(inputStr[m+1]);

        substrTemp = inputStr.substr(i, strideCount);
        
        lChar = substrTemp[midIndex-1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex+1];
        
        g = convertNTtoInt(lChar); 
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {
            if (lChar == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly); 
                substrTemp += leftTable[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable [j][k];
            }
            else if (mChar == rChar){
                substrTemp += inputStr.substr(i, nHomopoly); 
                substrTemp += rightTable[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable [j][k];
            }
             else{
                substrTemp += inputStr.substr(i, strideCount);
                substrTemp += rightTable[j][k];   
            }
            evenTurn = false;
        }
        else {
            if (lChar == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);  
                substrTemp += leftTable[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += rightTable [j][k];
            }
            else if (mChar == rChar) {
                substrTemp += inputStr.substr(i, nHomopoly);  
                substrTemp += rightTable[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += rightTable [j][k];
            }
            else {
                substrTemp += inputStr.substr(i, strideCount);
                substrTemp += leftTable[j][k]; 
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

bool checkHomopolymers(string encodeDNAStr, int homopoly)
{
    int maxHomopoly = 1, newHomopoly = 1;
    for (int i = 0 ; i < encodeDNAStr.length()-1 ; i++)
    {
        if ((encodeDNAStr[i] == encodeDNAStr[i+1]))
        {
            newHomopoly = newHomopoly + 1;
        }
        else
        {
            if (maxHomopoly < newHomopoly)
            {
                 maxHomopoly = newHomopoly;
            }
            newHomopoly = 1;
        }
    } 
    if (maxHomopoly == homopoly)
    {
        return true;
    }
    else
    {
        return false;
    }
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
    cout << "\n================= GC Content =====================";
    for (int j = 0; j < parts; j++)
    {
        cout << "\nPart." << j+1 << ".GC Content:" << (matrixGCcontent[j]*100)/(totalLen/parts) << "%";
    }
    cout << "\n================= GC Content =====================";
}

void printMatrixGCVariance(int ** matrixSubSeq, int subSeqeunceLength)
{
    cout << "\n================= GC Variation =====================";
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
    cout << "\n================= GC Variation =====================";
}

void mainPrint(bool printFlag , int indNo, string outStr)
{
    switch(indNo)
    {
        case 0: if (printFlag) cout << "\nString Length :" << outStr.length(); 
        break;
        case 1: if (printFlag) cout << "\nInput   String:" << outStr; 
        break;
        case 2: if (printFlag) cout << "\nEncoded String:" << outStr; 
        break;
        case 3: if (printFlag) cout << "\nDecoded String:" << outStr; 
        break;
        case 4: if (printFlag) cout << "\nDecoded Correctly:" << outStr;
        break;
        case 5: if (printFlag) cout << "\nHomopolymers are fine:" << outStr; 
        break;
        default: if (printFlag) cout << "\nProgram Completed:" << outStr; 
    }
}

int main()
{
    // input binary string length
    int strSize = 315, contentGCParts=1, subSeqLen=4, homopoly=1;
    // enable or disable all the console printing
    bool printAll = true;
    
    //generate a random input binary string of given length
    string inputStr = inputBinaryStrGen(strSize);
    
    //generate DNA seqeunce from input binary string
    string inputDNAStr = convertIntoDNAStr(inputStr);
    
    //Encode the DNA sequence using our proposed method
    string encodeOutStr = encodeDNAStrWithOption(inputDNAStr,homopoly);
    
    //Conversely, Decoding is performed below 
    string decodeOutStr = decodeDNAStrWithOption(encodeOutStr,homopoly);

   
    int * matrixGCcontent = calcPercentGC(encodeOutStr, contentGCParts);
    int ** martixGCVariance = calcVarianceGC(encodeOutStr, subSeqLen);
   
    // printing could be performed using the below lines of code, 
    // disable indivitual printing using false instead.
    mainPrint(printAll && true,  0 , inputDNAStr);
    mainPrint(printAll && true,  0 , encodeOutStr);
    mainPrint(printAll && false, 2 , encodeOutStr);
    mainPrint(printAll && false, 1 , inputDNAStr);
    mainPrint(printAll && false, 3 , decodeOutStr);
    mainPrint(printAll && true , 4 , (checkEncodeDecode(inputDNAStr, decodeOutStr) == 0 ? "YES": "OH! NO"));
    mainPrint(printAll && true , 5 , (checkHomopolymers(encodeOutStr, homopoly) == true? "YES": "OH! NO"));
    printMatrixGCContent(matrixGCcontent, contentGCParts, encodeOutStr.length());
    //printMatrixGCVariance(martixGCVariance, subSeqLen);
    mainPrint(true && true , 6 , "Done!");
        
    return 0;
}
