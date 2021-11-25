1. We have written this code in c++ language so we need a c++ compiler for it.
2. We can run this "Code_DNA.cpp" file to check our encoding and decoding scheme.
3. We can change the "strSize" variable in the "main" function for the selecting 
the string's length of our choice.   
4. Currently, please just select the length to a multiple of 3 and 4 both.
e.g, 12, 24, 48, 60, 600 etc.  


--------------------------------------------------------------

Sample Output 1: with strSize=24 and printFlag=true

Before Encoding DNA String:ATGCGGTCCTGTATTGTCTTAAAT  and length:24
After  Encoding DNA String:ATGACAGATGCATGTGACTCGTCATGACAGTC
After  Decoding DNA String:ATGCGGTCCTGTATTGTCTTAAAT

1. Encoding-Decoding Correctness!OK
2. Homopolymers                 !OK
3. Highly Balanced-GC Contents  !OK
4. Total(A+T+C+G) = 32;  A+T=17, 53%;  C+G=15, 46%;  Diff=7%

--------------------------------------------------------------

Sample Output 2: with strSize=24 and printFlag=false

1. Encoding-Decoding Correctness!OK
2. Homopolymers                 !OK
3. Highly Balanced-GC Contents  !OK
4. Total(A+T+C+G) = 32;  A+T=17, 53%;  C+G=15, 46%;  Diff=7%

---------------------------------------------------------------
  
