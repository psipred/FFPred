//protein feature calculations for UNIPROT

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Sequence.cpp"

std::stringstream usage("\n./features -i [fsa file]\n");

int main  (int argc, char* argv[])
{
 
 if(argc < 2)
 {
   std::cerr<<usage.str()<<std::endl;
   exit(1);

 }

 // read in command line

 unsigned int i=1;
 std::string fname  ="-";


 while( i < argc )
 {
   //std::cerr<<argv[i]<<std::endl;
  if( argv[i][0] == '-')
  {
   
    i++;
    switch(argv[i-1][1])
    {
      case 'i': { fname=argv[i]; break;}

      case 'h': { std::cerr<<usage.str()<<std::endl; exit(1);}

      default : { std::cerr<<usage.str()<<std::endl; exit(1);}
    }
  }
  i++;
 }


 std::vector<Sequence> seqs;

 std::ifstream fsa(fname.c_str());

 std::cerr<<"parsing fasta file "<<fname<<std::endl;
 //handle unreadable file
 if( fsa.bad() )
 {
   std::cerr<<usage.str()<<"\t Could not read file : "<< fname <<std::endl;
   exit(1);
 }

 //parse seqs into vector using '>' tokenizer
 std::string fasta;
 
 getline(fsa,fasta,'>');
 
 while( getline(fsa,fasta,'>') )
 {
   std::stringstream line(fasta.c_str());

   std::string tmp,aa;

   Sequence s;
   
   //assign sequence id to first line

   getline(line,tmp,'\n');
   
   s.set_id(tmp);
 
   //assign other lines to be sequence lines
   while(getline(line,tmp)){ aa+=tmp;}
   
   s.set_aa(aa);   

   //set sequence properties
   s.set_properties();
   
   //add sequence to vector
   //seqs.push_back(s);

   s.print();
 }
 
 fsa.close();
 
 return 0;
}
