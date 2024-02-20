/******************************************************************************

                    C++ Code written by Asad Usmani on 17/07/2023
                    A novel DNA encoding and decoding scheme

*******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;
int printFlag = 1;

char convertDNATable[4][4] = { 
    {'A', 'T', 'C', 'G'},
    {'T', 'C', 'G', 'A'},
    {'C', 'G', 'A', 'T'},
    {'G', 'A', 'T', 'C'}
};

// 40-60 table entries of GC content
char rightTable2[4][4] = {
    {'T', 'C', 'T', 'C'},
    {'C', 'G', 'G', 'A'},
    {'T', 'G', 'T', 'A'},
    {'C', 'A', 'A', 'T'}
};

char leftTable2[4][4] = {
    {'C', 'G', 'G', 'T'},
    {'G', 'A', 'A', 'C'},
    {'G', 'A', 'A', 'T'},
    {'T', 'C', 'T', 'A'}
};

// 50-50 table entries of GC content
char rightTable[4][4] = {
    {'G', 'C', 'T', 'C'},
    {'C', 'G', 'G', 'A'},
    {'T', 'G', 'T', 'A'},
    {'C', 'A', 'A', 'T'}
};

char leftTable[4][4] = {
    {'C', 'G', 'G', 'T'},
    {'G', 'C', 'A', 'C'},
    {'G', 'A', 'A', 'T'},
    {'T', 'C', 'T', 'A'}
};

// 60-40 table entries of GC content
char rightTable1[4][4] = {
    {'G', 'C', 'T', 'C'},
    {'C', 'G', 'G', 'A'},
    {'T', 'G', 'T', 'A'},
    {'C', 'A', 'A', 'C'}
};
  
char leftTable1[4][4] = { 
    {'C', 'G', 'G', 'T'},
    {'G', 'C', 'A', 'C'},
    {'G', 'A', 'G', 'T'},
    {'T', 'C', 'T', 'A'}
};


int convertNTtoInt(char inNt) {

    int outIdx;
    switch (inNt) {
    case 'A': outIdx = 0; break;
    case 'T': outIdx = 1; break;
    case 'C': outIdx = 2; break;
    case 'G': outIdx = 3; break;
    default: return -1;
    }
    return outIdx;
}


string convertIdxtoNTStr(int inIdx) {

    string outNt = "";
    int i = inIdx, j = 0, k = 4;
    while (k > 0)
    {
        int p = pow(4, k - 1);
        j = i / p;
        i = i % p;
        switch (j) {
        case 0: outNt += 'A'; break;
        case 1: outNt += 'T'; break;
        case 2: outNt += 'C'; break;
        case 3: outNt += 'G'; break;
        default: return outNt;
        }
        k = k - 1;
    }
    return outNt;
}

char convertInttoNT(int inIdx) {

    char outNt = 'A';
    switch (inIdx) {
    case 0: outNt = 'A'; break;
    case 1: outNt = 'T'; break;
    case 2: outNt = 'C'; break;
    case 3: outNt = 'G'; break;
    default: return outNt;
    }
    return outNt;
}

string convertBinarytoDNAStr(string bitStr, int homopoly, int percent)
{
    int j = 0;
    bool flag = false;
    string outDNAStr = "";
    string temp = "";
    int bitLen = bitStr.length();
    if (percent > 50)
    {
        flag = true;
        percent = 100 - percent;
    }
    int perBound = (20 * percent) / 100 + 1;

    for (j = 0; j < bitLen; j += 2)
    {
        temp = bitStr[j];
        temp = temp + bitStr[j + 1];
        if (temp == "00")
            outDNAStr += 'A';
        else if (temp == "01")
            outDNAStr += 'G';
        else if (temp == "10")
            outDNAStr += 'T';
        else if (temp == "11")
            outDNAStr += 'C';
    }
    bitLen = outDNAStr.length();
    if (flag == true)
    {
        for (j = 0; j < bitLen; j += 1)
        {
            if (j % perBound == 0)
                outDNAStr[j] = 'G';
        }
    }
    else
    {
        for (j = 0; j < bitLen; j += 1)
        {
            if (j % perBound == 0)
                outDNAStr[j] = 'A';
        }
    }
    
    while (bitLen % (homopoly*2+1) != 0)
    {
        outDNAStr += 'G';
        bitLen = outDNAStr.length();
    }
    return outDNAStr;
}

string convertIntoDNAStr(string str) {

    string outDNAStr = "";
    for (int i = 0; i < str.length(); i++) {

        switch (str[i]) {
        case '0': outDNAStr += "A"; break;
        case '1': outDNAStr += "T"; break;
        case '2': outDNAStr += "C"; break;
        case '3': outDNAStr += "G"; break;
        default: return "";
        }
    }
    return outDNAStr;
}

string convertIntintoNT(int val) {

    string charNt;
    switch (val) {
    case 0: charNt = "A"; break;
    case 1: charNt = "T"; break;
    case 2: charNt = "C"; break;
    case 3: charNt = "G"; break;
    default: return "";
    }
    return charNt;
}


int** calcGCVariance(string encodeOutStr, int subSeqLen)
{
    int prevSubSeq = 0, nextSubSeq = 0, countSubSeq = 0;
    int noOfCombination = pow(4, subSeqLen);
    int** matrixSubSeq = new int*[noOfCombination];
    for (int k = 0; k < noOfCombination; k++)
    {
        matrixSubSeq[k] = new int[noOfCombination];
    }

    for (int j = 0; j < noOfCombination; j++)
    {
        for (int k = 0; k < noOfCombination; k++)
        {
            matrixSubSeq[j][k] = 0;
        }
    }
    int count = 0;
    for (int k = 0; k < subSeqLen; k++)
    {
        for (int j = k; j < encodeOutStr.length() - subSeqLen + 1; j = j + subSeqLen)
        {
            int countNT[4] = { convertNTtoInt(encodeOutStr[j]),
                                convertNTtoInt(encodeOutStr[j + 1]),
                                convertNTtoInt(encodeOutStr[j + 2]),
                                convertNTtoInt(encodeOutStr[j + 3])
            };
            nextSubSeq = countNT[0] * 64 + countNT[1] * 16 + countNT[2] * 4 + countNT[3];
            matrixSubSeq[prevSubSeq][nextSubSeq] = matrixSubSeq[prevSubSeq][nextSubSeq] + 1;
            prevSubSeq = nextSubSeq;
        }
    }
    return matrixSubSeq;
}


int* calcPercentGC(string encodeOutStr, int parts)
{
    int k = encodeOutStr.length() / parts;
    int* matrixParts = new int[parts];

    for (int j = 0; j < parts; j++)
    {
        int countNT[4] = { 0, 0, 0, 0 };
        for (int i = j * k; i < (j * k) + k; i++)
        {
            (
                encodeOutStr[i] == 'A' ? countNT[0] += 1 :
                encodeOutStr[i] == 'T' ? countNT[1] += 1 :
                encodeOutStr[i] == 'C' ? countNT[2] += 1 :
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

    for (int i = 0, m = 3; m <= inputStr.length(); i += 4, m += 4) {

        j = convertNTtoInt(inputStr[m - 1]);
        k = convertNTtoInt(inputStr[m + 1]);

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
            else if (leftTable[g][h] == mChar) {
                substrTemp += lChar;
                substrTemp += lChar;
                substrTemp += rChar;
            }
            else {
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
            else if (leftTable[g][h] == mChar) {
                substrTemp += lChar;
                substrTemp += lChar;
                substrTemp += rChar;
            }
            else {
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


string decodeDNAStrWithOptionTable2(string inputStr, int nHomopoly) {

    bool evenTurn = true;
    int strideCount = 2 * nHomopoly + 1;
    int midIndex = strideCount / 2;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, mChar, rChar, eChar;
    string outStr = "", substrTemp = "";

    for (int i = 0, m = strideCount; m <= inputStr.length();
        i += strideCount + 1, m += strideCount + 1) {

        j = convertNTtoInt(inputStr[m - 1]);
        k = convertNTtoInt(inputStr[m + 1]);

        substrTemp = inputStr.substr(i, strideCount + 1);

        lChar = substrTemp[midIndex - 1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex + 1];
        eChar = substrTemp[strideCount];

        g = convertNTtoInt(lChar);
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {

            if (eChar == rightTable2[j][k]) {
                substrTemp += inputStr.substr(i, strideCount);
            }
            else if (leftTable2[g][h] == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            else {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            evenTurn = false;
        }
        else {

            if (eChar == leftTable2[j][k]) {
                substrTemp += inputStr.substr(i, strideCount);
            }
            else if (leftTable2[g][h] == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            else {
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

string encodeDNAStrWithOptionTable2(string inputStr, int nHomopoly) {

    bool evenTurn = true;
    int strideCount = 2 * nHomopoly + 1;
    int midIndex = strideCount / 2;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, rChar, mChar;
    string outStr = "", substrTemp = "";

    for (int i = 0, m = strideCount - 1; m <= inputStr.length();
        i += strideCount, m += strideCount) {

        j = convertNTtoInt(inputStr[m]);
        k = convertNTtoInt(inputStr[m + 1]);

        substrTemp = inputStr.substr(i, strideCount);

        lChar = substrTemp[midIndex - 1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex + 1];

        g = convertNTtoInt(lChar);
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {
            if (lChar == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += leftTable2[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable2[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rightTable2[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable2[j][k];
            }
            else {
                substrTemp += inputStr.substr(i, strideCount);
                substrTemp += rightTable2[j][k];
            }
            evenTurn = false;
        }
        else {
            if (lChar == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += leftTable2[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += rightTable2[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rightTable2[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += rightTable2[j][k];
            }
            else {
                substrTemp += inputStr.substr(i, strideCount);
                substrTemp += leftTable2[j][k];
            }
            evenTurn = true;
        }
        outStr += substrTemp;
        substrTemp = "";
    }
    return outStr;
}

string decodeDNAStrWithOptionTable1(string inputStr, int nHomopoly) {

    bool evenTurn = true;
    int strideCount = 2 * nHomopoly + 1;
    int midIndex = strideCount / 2;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, mChar, rChar, eChar;
    string outStr = "", substrTemp = "";

    for (int i = 0, m = strideCount; m <= inputStr.length();
        i += strideCount + 1, m += strideCount + 1) {

        j = convertNTtoInt(inputStr[m - 1]);
        k = convertNTtoInt(inputStr[m + 1]);

        substrTemp = inputStr.substr(i, strideCount + 1);

        lChar = substrTemp[midIndex - 1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex + 1];
        eChar = substrTemp[strideCount];

        g = convertNTtoInt(lChar);
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {

            if (eChar == rightTable1[j][k]) {
                substrTemp += inputStr.substr(i, strideCount);
            }
            else if (leftTable1[g][h] == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            else {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            evenTurn = false;
        }
        else {

            if (eChar == leftTable1[j][k]) {
                substrTemp += inputStr.substr(i, strideCount);
            }
            else if (leftTable1[g][h] == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            else {
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

string encodeDNAStrWithOptionTable1(string inputStr, int nHomopoly) {

    bool evenTurn = true;
    int strideCount = 2 * nHomopoly + 1;
    int midIndex = strideCount / 2;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, rChar, mChar;
    string outStr = "", substrTemp = "";

    for (int i = 0, m = strideCount - 1; m <= inputStr.length();
        i += strideCount, m += strideCount) {

        j = convertNTtoInt(inputStr[m]);
        k = convertNTtoInt(inputStr[m + 1]);

        substrTemp = inputStr.substr(i, strideCount);

        lChar = substrTemp[midIndex - 1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex + 1];

        g = convertNTtoInt(lChar);
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {
            if (lChar == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += leftTable1[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable1[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rightTable1[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable1[j][k];
            }
            else {
                substrTemp += inputStr.substr(i, strideCount);
                substrTemp += rightTable1[j][k];
            }
            evenTurn = false;
        }
        else {
            if (lChar == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += leftTable1[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += rightTable1[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rightTable1[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += rightTable1[j][k];
            }
            else {
                substrTemp += inputStr.substr(i, strideCount);
                substrTemp += leftTable1[j][k];
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
    int midIndex = strideCount / 2;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, mChar, rChar, eChar;
    string outStr = "", substrTemp = "";

    for (int i = 0, m = strideCount; m <= inputStr.length();
        i += strideCount + 1, m += strideCount + 1) {

        j = convertNTtoInt(inputStr[m - 1]);
        k = convertNTtoInt(inputStr[m + 1]);

        substrTemp = inputStr.substr(i, strideCount + 1);

        lChar = substrTemp[midIndex - 1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex + 1];
        eChar = substrTemp[strideCount];

        g = convertNTtoInt(lChar);
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {

            if (eChar == rightTable[j][k]) {
                substrTemp += inputStr.substr(i, strideCount);
            }
            else if (leftTable[g][h] == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            else {
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
            else if (leftTable[g][h] == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += lChar;
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
            }
            else {
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

string reverseStr(string inputStr, bool flag)
{
    char lChar = ' ', rChar = ' ', mChar = ' ';
    string outStr = "", substrTemp = "";
    int i = 0, j = 0, g = 0, h = 0, endLen = inputStr.length();
   
    if (flag == true)
    {
        for (i = 0; i < endLen-1; i += 1)
        {
            if (i % 2 == 1)
            {
                g = convertNTtoInt(inputStr[i]);
                h = convertNTtoInt(inputStr[i + 2]);

                outStr += inputStr.substr(i, 1);
                outStr += rightTable[g][h];
                i += 1;
            }
            else
            {
                outStr += inputStr.substr(i, 1);
            }
        }
    }
    else
    {
        for (i = 0; i < endLen-1; i += 1)
        {
            if (i % 2 == 1)
            {
                g = convertNTtoInt(inputStr[i]);
                h = convertNTtoInt(inputStr[i + 2]);

                outStr += inputStr.substr(i, 1);
                outStr += leftTable[g][h];
                i += 1;
            }
            else
            {
                outStr += inputStr.substr(i, 1);
            }
        }
    }
    return outStr;
}

string encodeDNAStrWithExtendedVersion(string inputStr, int nHairpinSeq) {
    char lChar = ' ', rChar = ' ', mChar = ' ';
    string outStr = "", substrTemp = "";
    int i = 0, g = 0, h = 0, j = nHairpinSeq / 2 - 1, 
        endLen = inputStr.length() - nHairpinSeq / 2 + 1;
    bool flag = true;
    for (i = 0; i <= endLen; i += j)   
    {
        g = convertNTtoInt(inputStr[i + j - 1]);
        h = convertNTtoInt(inputStr[i + j]);

        outStr += inputStr.substr(i, j);
        outStr += leftTable[g][h];
    }

    inputStr = outStr;
    outStr = "";
    j = nHairpinSeq / 2 - 1;
    endLen = inputStr.length() - nHairpinSeq / 2;
    for (i = 0; i <= endLen; i += j)
    {
        if (flag == true)
        {
            outStr += reverseStr(inputStr.substr(i, j + 2), true);
            flag = false;
        }
        else if (flag == false)
        {
            outStr += inputStr.substr(i, j + 1);
            flag = true;
        }
        i += 1;
    }
   
    return outStr;
}

string decodeDNAStrWithExtendedVersion(string inputStr, int nHairpinSeq) {
    char lChar = ' ', rChar = ' ', mChar = ' ';
    string outStr = "", substrTemp = "";
    int i = 0, g = 0, h = 0, j = nHairpinSeq / 2 - 1,
        endLen = inputStr.length() - nHairpinSeq / 2 + 1;
    for (i = 0; i <= endLen; i += j)
    {
        g = convertNTtoInt(inputStr[i + j - 1]);
        h = convertNTtoInt(inputStr[i + j + 1]);

        if (leftTable[g][h] != inputStr[i + j])
        {
            outStr += reverseStr(inputStr.substr(i, j + 2), false);
        }
        else
        {
            outStr += inputStr.substr(i, j + 1);
        }
        i += 1;
        //cout << "\nPart: " << inputStr.substr(i, j) << " " << leftTable[g][h];
    }
    inputStr = outStr;
    outStr = "";
    for (i = 0; i <= endLen; i += j)
    {
        outStr += inputStr.substr(i, j);
        i += 1;
        //cout << "\nPart: " << inputStr.substr(i, j) << " " << leftTable[g][h];
    }
    return outStr;
}

string encodeDNAStrWithOption(string inputStr, int nHomopoly) {

    bool evenTurn = true;
    int strideCount = 2 * nHomopoly + 1;
    int midIndex = strideCount / 2;
    int j = 0, k = 0, g = 0, h = 0;
    char lChar, rChar, mChar;
    string outStr = "", substrTemp = "";

    for (int i = 0, m = strideCount - 1; m <= inputStr.length();
        i += strideCount, m += strideCount) {

        j = convertNTtoInt(inputStr[m]);
        k = convertNTtoInt(inputStr[m + 1]);

        substrTemp = inputStr.substr(i, strideCount);

        lChar = substrTemp[midIndex - 1];
        mChar = substrTemp[midIndex];
        rChar = substrTemp[midIndex + 1];

        g = convertNTtoInt(lChar);
        h = convertNTtoInt(rChar);

        substrTemp = "";
        if (evenTurn == true) {
            if (lChar == mChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += leftTable[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rightTable[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += leftTable[j][k];
            }
            else {
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
                substrTemp += rightTable[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += inputStr.substr(i, nHomopoly);
                substrTemp += rightTable[g][h];
                substrTemp += inputStr.substr(i + nHomopoly + 1, nHomopoly);
                substrTemp += rightTable[j][k];
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

    for (int i = 0, m = 2; m <= inputStr.length(); i += 3, m += 3) {

        j = convertNTtoInt(inputStr[m]);
        k = convertNTtoInt(inputStr[m + 1]);

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
                substrTemp += leftTable[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += lChar;
                substrTemp += rightTable[g][h];
                substrTemp += rChar;
                substrTemp += leftTable[j][k];
            }
            else {
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
                substrTemp += rightTable[j][k];
            }
            else if (mChar == rChar) {
                substrTemp += lChar;
                substrTemp += rightTable[g][h];
                substrTemp += rChar;
                substrTemp += rightTable[j][k];
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

int checkHomopolymers(string encodeDNAStr)
{
    int maxHomopoly = 1, newHomopoly = 1;
    for (int i = 0; i < encodeDNAStr.length() - 1; i++)
    {
        if ((encodeDNAStr[i] == encodeDNAStr[i + 1]))
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
    return maxHomopoly;
}

string inputBinaryStrGen(int strSize)
{
    string inputStr = "";
    srand(time(NULL));
    for (int i = 0; i < strSize; i++)
    {
        int rval = rand() % 4;
        (
            rval == 3 ? inputStr += '3' :
            rval == 2 ? inputStr += '2' :
            rval == 1 ? inputStr += '1' :
            inputStr += '0'
            );
    }
    srand(time(NULL));
    return inputStr;
}

bool checkEncodeDecode(string inputDNAStr, string decodeOutStr)
{
    return inputDNAStr.compare(decodeOutStr);
}


void printMatrixGCVariance(int** matrixSubSeq, int subSeqLen, int matValLen)
{
    int count = 0, flag = 0, valCol = 0, valRow = 0;
    int martixLen = pow(4, subSeqLen);
    int* matrixSubSeqCol = new int[martixLen];
    int* matrixSubSeqRow = new int[martixLen];

    for (int j = 0; j < martixLen; j++)
    {
        flag = 0;
        for (int k = 0; k < martixLen; k++)
        {
            if (matrixSubSeq[j][k] > 0)
            {
                flag = 1;
            }
        }
        if (flag)
            matrixSubSeqRow[j] = 1;
        else
            matrixSubSeqRow[j] = 0;
    }

    for (int j = 0; j < martixLen; j++)
    {
        flag = 0;
        for (int k = 0; k < martixLen; k++)
        {
            if (matrixSubSeq[k][j] > 0)
            {
                flag = 1;
            }
        }
        if (flag)
            matrixSubSeqCol[j] = 1;
        else
            matrixSubSeqCol[j] = 0;
    }

    valCol = 0;

    cout << "\n================= GC Variation =====================\n";
    
    for (int j = 0; j < martixLen && valCol <= matValLen * 2; j += 3)
    {
        if (matrixSubSeqRow[j] == 1)
        {
            for (int k = 0; valCol <= matValLen * 2 && k < martixLen; k += 3)
            {
                if (matrixSubSeqCol[k] == 1)
                {
                    valCol = valCol + 1;
                    cout << convertIdxtoNTStr(k) << ", ";
                    if (valCol % matValLen == 0)
                        cout << "\n";
                }
            }
        }
    }
    
    cout << "\n================= GC Variation =====================\n";

    valCol = 0;
    valRow = 0;
    
    cout << "\n[\n";
    int counter = 0;
    for (int j = 0; valRow <= matValLen * 2 && j < martixLen; j += 3)
    {
        if (matrixSubSeqRow[j])
        {
            cout << "\n" << counter++ << "[";
            valCol = 0;
            valRow = valRow + 1;
            
            int flagTemp = 0;
            for (int k = 0; valCol <= matValLen * 2 && k < martixLen; k += 3)
            {
                if (matrixSubSeqCol[k])
                {
                    if (matrixSubSeq[j][k] > 0)
                    {
                        cout << matrixSubSeq[j][k] << ",";

                    }
                    else
                    {
                        cout << "0" << ",";
                    }
                    valCol = valCol + 1;
                }
            }
            cout << "],";
        }
    }
    cout << "\n]\n";
    cout << "\n================= GC Variation =====================\n";
}

void writeToFile(int percentGC, int encoded, int homopoly)
{
    string filePath = "D:\\PhDGoetheUni\\Works\\PhD_Work\\3.Encoding_Method_Paper\\Results\\DNAEncodingResults\\";
    string fileName = filePath + "filename_GContent_";
    if (homopoly == 0)
        fileName += "O.csv";
    else if (homopoly == 1)
        fileName += "1.csv";
    else if (homopoly == 2)
        fileName += "2.csv";
    else if (homopoly == 3)
        fileName += "3.csv";
    else if (homopoly == 4)
        fileName += "4.csv";
    else
        fileName += "5.csv";
    fstream myFile(fileName, ios_base::app);
    // Create and open a text file
    if (myFile.is_open())
    {
        myFile << endl << percentGC ;
        myFile.close();
    }
    else cerr << "Unable to open file";
}

void printMatrixGCContent(int* matrixGCcontent, int parts, int totalLen, int homopoly, bool encoded)
{
    int max = 0, min = 100;
    if (encoded == true)
    {
        //cout << "\nEncoded String : Homopoly" << homopoly << "\n";
    }
    else
    {
        //cout << "\nOriginal String: Homopoly" << homopoly << "\n";
        return;
    }
    
    for (int j = 0; j < parts; j++)
    {
        int percentGC = (matrixGCcontent[j] * 100) / (totalLen / parts);
        if (min > percentGC)
        {
            min = percentGC;
        }
        if (max < percentGC)
        {
            max = percentGC;
        }
        //cout << j+1 << ":" << percentGC << "%, ";
            if (homopoly <= 1)
                homopoly = 0;
            else
                homopoly = 5;

        //if (encoded and (percentGC > 60 or percentGC < 40))
        //{
        //    cout << "  Out of range!";
        //}
        writeToFile(percentGC, encoded, homopoly);
    }
    cout << max - min << ", ";
}

void mainPrint(bool printFlag, int indNo, string outStr)
{
    switch (indNo)
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


bool checkHairPin(string subSeq, int stemLen, int loopLen)
{
    int totalLen = subSeq.length();
    int countMatch = 0;
    char leftNuct = ' ', rightNuct = ' ';
    int endLen = totalLen - stemLen + 1;
    char compNuct[4] = { 'T', 'A', 'G', 'C' };
    int leftEnd = totalLen / 2 - loopLen / 2 - 1;
    int rightEnd = totalLen / 2 + loopLen / 2;
    for (leftEnd; leftEnd >= 0;)
    {
        leftNuct = compNuct[convertNTtoInt(subSeq[rightEnd])];
        if (subSeq[leftEnd] == leftNuct)
        {
            rightEnd = rightEnd + 1;
            leftEnd = leftEnd - 1;
            countMatch = countMatch + 1;
        }
        else if (subSeq[leftEnd] == subSeq[rightEnd])
        {
            rightEnd = rightEnd + 1;
            leftEnd = leftEnd - 1;
        }
        else
        {
            return false;
        }
    }
    if (countMatch >= stemLen)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int searchForHairPins(string inputStr, int lenSubSeq, int stemLen, int loopLen)
{
    int i = 0, j = 0, countHairPins = 0;
    int endLen = inputStr.length() - lenSubSeq + 1;
    string subSeq;
    for (i = 0; i < endLen; i++)
    {
        subSeq = inputStr.substr(i, lenSubSeq);
        if (checkHairPin(subSeq, stemLen, loopLen) == true)
        {
            countHairPins = countHairPins + 1;
            i += lenSubSeq;
            //cout << "\n" << subSeq;
        }
    }
    return countHairPins;
}

int findHaipinSequence(string inputDNAStr, int strandSize, 
    int subSeqLen, int loopLen, int stemLen)
{
    int countHairPins = 0, countTotalHairPins = 0, hasStrand = 0;

    //Check for the Hairpins using our proposed method
    for (int i = 0; i < inputDNAStr.length(); i = i + strandSize)
    {
        hasStrand = 0;
        string subSeq = inputDNAStr.substr(i, strandSize);
        hasStrand = searchForHairPins(subSeq, subSeqLen, loopLen, stemLen);
        if (hasStrand > 0)
        {
            countHairPins += 1;
            countTotalHairPins += hasStrand;
        }
    }
    return countTotalHairPins;
}

int removeHairpinSequence(string inputDNAStr, int strandSize,
    int subSeqLen, int loopLen, int stemLen)
{
    int countHairPins = 0, countTotalHairPins = 0, hasStrand = 0;

    //Check for the Hairpins using our proposed method
    for (int i = 0; i < inputDNAStr.length(); i = i + strandSize)
    {
        hasStrand = 0;
        string subSeq = inputDNAStr.substr(i, strandSize);
        //hasStrand = searchAndRemoveHairPins(subSeq, subSeqLen, loopLen, stemLen);
        if (hasStrand > 0)
        {
            countHairPins += 1;
            countTotalHairPins += hasStrand;
            cout << "\nRemoved! " << countTotalHairPins;
        }
    }
    return countTotalHairPins;
}


float calculateHammingDistance(string inputStr, int subStrLen, int minHamDist)
{
    float avgDist = 0.0, countDist = 0, sumDist = 0;
    int i = 0, j = 0, k = 0;
    int endIdx = inputStr.length() - subStrLen, endIdy = 0;
    for (j = subStrLen; j < endIdx; j += subStrLen)
    {
        k = j;
        countDist = 0;
        for (i = j - subStrLen; i < j; i += 1)
        {
            if (inputStr[i] != inputStr[k])
            {
                countDist += 1;
            }
            k += 1;
        }
        if (countDist >= minHamDist)
            sumDist += 1;
    }
    endIdx = endIdx / subStrLen;
    avgDist = (sumDist * 100) / endIdx;
    //cout << "\nsumDist: " << sumDist << "  endIdx: " << endIdx;
    cout << "\nSublength: " << subStrLen << "  Hamming Distance: " << avgDist << "%";
    return avgDist;
}


void printResutls(bool printAll,
    string inputDNAStr,
    string encodeOutStr,
    string decodeOutStr,
    int homopoly,
    int subSeqLen,
    int contentGCParts,
    int serialNo,
    int strandSize,
    int gcVarMatLen)
{
    // printing could be performed using the below lines of code, 
    // disable indivitual printing using fasle instead.
    float inputLen = inputDNAStr.length();
    float encodeLen = encodeOutStr.length();
    float infoDensity = (2 * inputLen) / encodeLen;

    //cout << "\n==================================================";
    //cout << "\nString#" << serialNo;
    mainPrint(printAll && false, 0, inputDNAStr);
    mainPrint(printAll && false, 1, inputDNAStr);
    mainPrint(printAll && false, 0, encodeOutStr);
    mainPrint(printAll && false, 2, encodeOutStr);
    mainPrint(printAll && false, 0, decodeOutStr);
    mainPrint(printAll && false, 3, decodeOutStr);
    mainPrint(printAll && false, 4, (checkEncodeDecode(inputDNAStr, decodeOutStr) == 0 ? "YES" : "OH! NO"));
    mainPrint(printAll && false, 5, (checkHomopolymers(encodeOutStr) <= homopoly ? "YES" : "OH! NO"));
    mainPrint(printAll && false, 6, "Done!");
    
    //cout << "\nInformation Density: " << infoDensity << "  and homopoly:" << homopoly;
    // Calculate the GC content parts of the input string
    //if (homopoly == 1)
    //{
    //    int** martixGCVariance = calcGCVariance(encodeOutStr, subSeqLen);
    //    //printMatrixGCVariance(martixGCVariance, subSeqLen, gcVarMatLen);
    //}
    
    contentGCParts = inputDNAStr.length() / strandSize;
    int* matrixGCcontent = calcPercentGC(inputDNAStr, contentGCParts);
    printMatrixGCContent(matrixGCcontent, contentGCParts, inputDNAStr.length(), homopoly, false);

    // Calculate the GC content parts of the DNA endcoded string
    contentGCParts = encodeOutStr.length() / strandSize;
    matrixGCcontent = calcPercentGC(encodeOutStr, contentGCParts);
    printMatrixGCContent(matrixGCcontent, contentGCParts, encodeOutStr.length(), homopoly, true);

    //calculateHammingDistance(encodeOutStr, subSeqLen, subSeqLen - 2);

    //cout << "\n==================================================";

}

void clearAllFiles(string filePath, int strLen)
{
    string fileName = filePath + "filename_GContent_O.csv";
    cout << fileName;
    ofstream myFile1(fileName, ios::trunc);
    if (myFile1.is_open())
    {
        myFile1 << "Sequence Length:" << strLen;
        myFile1.close();
    }
    else cout << "File opening error!";
    fileName = filePath + "filename_GContent_1.csv";
    ofstream myFile2(fileName);
    if (myFile2.is_open())
    {
        myFile2 << "Sequence Length:" << strLen;
        myFile2.close();
    }
    fileName = filePath + "filename_GContent_2.csv";
    ofstream myFile3(fileName);
    if (myFile3.is_open())
    {
        myFile3 << "Sequence Length:" << strLen;
        myFile3.close();
    }
    fileName = filePath + "filename_GContent_3.csv";
    ofstream myFile4(fileName);
    if (myFile4.is_open())
    {
        myFile4 << "Sequence Length:" << strLen;
        myFile4.close();
    }
    fileName = filePath + "filename_GContent_4.csv";
    ofstream myFile5(fileName);
    if (myFile5.is_open())
    {
        myFile5 << "Sequence Length:" << strLen;
        myFile5.close();
    }
    fileName = filePath + "filename_GContent_5.csv";
    ofstream myFile6(fileName);
    if (myFile6.is_open())
    {
        myFile6 << "Sequence Length:" << strLen;
        myFile6.close();
    }

}

int main()
{
    int p = 0;
    setprecision(2);
    string filePath = "D:\\PhDGoetheUni\\Works\\PhD_Work\\3.Encoding_Method_Paper\\Results\\DNAEncodingResults\\";
    string fileName = filePath + "BinaryCodeDataset1_250_32000.txt";
    clearAllFiles(filePath, 16385);
    fstream myFile(fileName);
    cout << "\n";
    if (myFile.is_open())
    {
        // Chose the printing flag status
        bool printAll = true;

        string inputStr = "", tempStr = "01";

        // input binary string length
        int strSize = 0, strandSize = 315, serialNo = 0, gcVarMatLen = 16, noOfInstances = 1;

        // Enter the length of the homopoly
        int contGCParts = 1, subSeqLen = 7, homopoly = 1, range = 80, base = 10;
        int hairpinCount = 0, subSeqLenForHairpin = 26;

        int* randPercent = new int[noOfInstances];

        for (int k = 0; k < noOfInstances; k++)
        {
            randPercent[k] = (rand() % range) + base;
        }
        while (getline(myFile, inputStr) and serialNo < noOfInstances)
        {
            
            while (inputStr.length() % strandSize != 0)
            {   
                inputStr += tempStr[rand() % 2];
            }
            
            homopoly = 1;
            string inputDNAStr = convertBinarytoDNAStr(inputStr, homopoly, randPercent[serialNo++]);
            //calculateHammingDistance(inputDNAStr, subSeqLen, subSeqLen-2);

            cout << "\n============================================:" << serialNo;
            //findHaipinSequence(inputDNAStr, strandSize, true);
           
            while (homopoly < 5)
            {
                if (homopoly == 1)
                {
                    inputDNAStr = "GCTGATCATGGATCAAGG";

                    cout << "\nInput: \n" << inputDNAStr;
                    //Encode the DNA sequence using our proposed method
                    string encodeOutStr = encodeDNAStrWithOption(inputDNAStr, homopoly);
                    cout << "\nEncoded: \n" << encodeOutStr;

                    encodeOutStr = encodeDNAStrWithExtendedVersion(encodeOutStr, subSeqLenForHairpin);
                    cout << "\nExtended: \n" << encodeOutStr;

                    //Conversely, Decoding is performed below 
                    
                    string decodeOutStr = decodeDNAStrWithExtendedVersion(encodeOutStr, subSeqLenForHairpin);
                    cout << "\nDe-Extended: \n" << decodeOutStr;

                    decodeOutStr = decodeDNAStrWithOption(decodeOutStr, homopoly);
                    cout << "\nDecoded: \n" << decodeOutStr;

                    return 0;

                    //cout << "\nLength of the string:" << inputDNAStr.length();
                    //cout << "\nLength of the string:" << encodeOutStr.length();
                    // 
                    hairpinCount = findHaipinSequence(encodeOutStr, strandSize, 20, 6, 6);
                    if (hairpinCount > 0)
                    {
                        cout << "\nHairpin: " << hairpinCount;
                        removeHairpinSequence(encodeOutStr, strandSize, 20, 6, 6);
                    }


                    /*printResutls(
                        printAll,
                        inputDNAStr,
                        encodeOutStr,
                        decodeOutStr,
                        homopoly,
                        subSeqLen,
                        contGCParts,
                        serialNo,
                        strandSize,
                        gcVarMatLen);
                        */

                }
                homopoly += 1;
            }
            
        }
        myFile.close();

    }
    else
        cout << "File is not opening!";
    return 0;
}
